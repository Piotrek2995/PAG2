import arcpy
from collections import defaultdict
from typing import Dict, Tuple, Optional
import math

# --- Konfiguracja parametrów (pobierane z toolboxa) ---
FC_ROADS = arcpy.GetParameterAsText(0)
gdb_path = arcpy.GetParameterAsText(1)
nodes_name = "nodes_out"
edges_name = "edges_out"

# Pola (dostosuj jeśli w Twojej warstwie są inne nazwy)
FIELD_OID    = "OBJECTID"
FIELD_SHAPE  = "SHAPE@"
FIELD_CLASS  = "KLASA_DROG"
FIELD_DIR_OPT = "kierunkowosc"

# Słowniki wynikowe
vertices: Dict[int, Dict] = {}
edges: Dict[int, Dict] = {}

def _map_klasa_bdot(klasa_txt: Optional[str]) -> str:
    if not klasa_txt: return "G"
    s = klasa_txt.lower()
    if "autostr" in s: return "A"
    if "ekspres" in s: return "S"
    if "ruchu" in s: return "GP"
    if "główn" in s: return "G"
    if "zbior" in s: return "Z"
    if "lokal" in s: return "L"
    if "dojaz" in s: return "D"
    if "inna" in s or "wewn" in s: return "I"
    return "G"

def _snap_key(x: float, y: float, tol: float) -> Tuple[int, int]:
    return (round(x / tol), round(y / tol))

def build_graph_from_fc(fc: str, snap_tol: float = 0.25):
    vertex_ids_by_snap = {}
    next_vid = 1
    next_eid = 1

    fld_names_lower = {f.name.lower() for f in arcpy.ListFields(fc)}
    has_dir = FIELD_DIR_OPT.lower() in fld_names_lower
    fields = [FIELD_OID, FIELD_SHAPE, FIELD_CLASS] + ([FIELD_DIR_OPT] if has_dir else [])

    with arcpy.da.SearchCursor(fc, fields) as cur:
        for row in cur:
            jezdnia_oid = int(row[0])
            geom = row[1]
            klasa_txt = row[2]
            kier = row[3] if has_dir else 0

            if not geom: continue
            first_pt = geom.firstPoint
            last_pt = geom.lastPoint
            if first_pt is None or last_pt is None: continue

            x1, y1 = first_pt.X, first_pt.Y
            x2, y2 = last_pt.X, last_pt.Y
            length_m = float(geom.length)
            klasa_code = _map_klasa_bdot(klasa_txt)
            kier = int(kier) if kier is not None else 0

            k1 = _snap_key(x1, y1, snap_tol)
            k2 = _snap_key(x2, y2, snap_tol)

            if k1 not in vertex_ids_by_snap:
                vertex_ids_by_snap[k1] = next_vid
                vertices[next_vid] = {"x": x1, "y": y1, "edge_out": []}
                next_vid += 1
            if k2 not in vertex_ids_by_snap:
                vertex_ids_by_snap[k2] = next_vid
                vertices[next_vid] = {"x": x2, "y": y2, "edge_out": []}
                next_vid += 1

            u = vertex_ids_by_snap[k1]
            v = vertex_ids_by_snap[k2]

            def add_edge(u_id, v_id):
                nonlocal next_eid
                edges[next_eid] = {
                    "id": next_eid,
                    "id_from": u_id,
                    "id_to": v_id,
                    "edge_length_field": length_m,
                    "kier_auto": kier,
                    "klasa_drogi": klasa_code,
                    "jezdnia_oid": jezdnia_oid
                }
                vertices[u_id]["edge_out"].append(next_eid)
                next_eid += 1

            if kier == 3:
                continue
            elif kier == 0:
                add_edge(u, v); add_edge(v, u)
            elif kier == 1:
                add_edge(u, v)
            elif kier == 2:
                add_edge(v, u)
            else:
                add_edge(u, v); add_edge(v, u)

    return len(vertices), len(edges)

def export_graph_to_gdb(gdb_path: str, nodes_name: str = "nodes_out", edges_name: str = "edges_out"):
    sr = arcpy.Describe(FC_ROADS).spatialReference
    nodes_fc = f"{gdb_path}\\{nodes_name}"
    edges_fc = f"{gdb_path}\\{edges_name}"

    for fc in [nodes_fc, edges_fc]:
        if arcpy.Exists(fc):
            arcpy.management.Delete(fc)

    arcpy.management.CreateFeatureclass(gdb_path, nodes_name, "POINT", spatial_reference=sr)
    arcpy.management.AddField(nodes_fc, "node_id", "LONG")
    with arcpy.da.InsertCursor(nodes_fc, ["SHAPE@XY", "node_id"]) as icur:
        for vid, v in vertices.items():
            icur.insertRow(((v["x"], v["y"]), vid))
    arcpy.AddMessage(f"[EXPORT] Zapisano {len(vertices)} węzłów do {nodes_fc}")

    arcpy.management.CreateFeatureclass(gdb_path, edges_name, "POLYLINE", spatial_reference=sr)
    arcpy.management.AddField(edges_fc, "edge_id", "LONG")
    arcpy.management.AddField(edges_fc, "id_from", "LONG")
    arcpy.management.AddField(edges_fc, "id_to", "LONG")
    arcpy.management.AddField(edges_fc, "length_m", "DOUBLE")
    arcpy.management.AddField(edges_fc, "klasa", "TEXT")
    arcpy.management.AddField(edges_fc, "kier", "SHORT")
    arcpy.management.AddField(edges_fc, "jezdnia_oid", "LONG")

    with arcpy.da.InsertCursor(edges_fc,
        ["SHAPE@", "edge_id", "id_from", "id_to", "length_m", "klasa", "kier", "jezdnia_oid"]) as icur:
        for eid, e in edges.items():
            u = vertices[e["id_from"]]; v = vertices[e["id_to"]]
            arr = arcpy.Array([arcpy.Point(u["x"], u["y"]), arcpy.Point(v["x"], v["y"])])
            poly = arcpy.Polyline(arr, sr)
            icur.insertRow((poly, eid, e["id_from"], e["id_to"],
                            e["edge_length_field"], e["klasa_drogi"], e["kier_auto"], e.get("jezdnia_oid", None)))
    arcpy.AddMessage(f"[EXPORT] Zapisano {len(edges)} krawędzi do {edges_fc}")

    # ustaw outputy toola (indeksy: 2 i 3, zgodnie z parametrami narzędzia)
    arcpy.SetParameterAsText(2, nodes_fc)
    arcpy.SetParameterAsText(3, edges_fc)

    # dodaj warstwy natychmiast do mapy i włącz etykiety
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        m = aprx.activeMap
        lyr_nodes = m.addDataFromPath(nodes_fc)
        if lyr_nodes.supports("SHOWLABELS"):
            lyr_nodes.showLabels = True
            lbl = lyr_nodes.listLabelClasses()[0]
            lbl.showClassLabels = True
            lbl.expression = "$feature.node_id"
        m.addDataFromPath(edges_fc)
        aprx.save()
    except Exception as e:
        arcpy.AddWarning(f"Nie udało się dodać warstw/etykiet do mapy: {e}")

# --- main ---
if __name__ == "__main__":
    nV, nE = build_graph_from_fc(FC_ROADS)
    arcpy.AddMessage(f"Graph built: |V|={nV}, |E|={nE}")
    export_graph_to_gdb(gdb_path, nodes_name, edges_name)
