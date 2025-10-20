import arcpy
from typing import Dict, Tuple, Optional

# Konfiguracja
arcpy.env.workspace = r"C:\Users\piotr\Documents\ArcGIS\Projects\Projekt1_PAG2\Projekt1_PAG2.gdb"
FC_ROADS = "skjzPAG2"
WHERE = None

FIELD_SHAPE = "SHAPE@"
FIELD_CLASS = "KLASA_DROG"
FIELD_DIR_OPT = "kierunkowosc"

# Model grafu
vertices: Dict[int, Dict] = {}
edges: Dict[int, Dict] = {}

# Funkcje pomocnicze
def _map_klasa_bdot(klasa_txt: Optional[str]) -> str:
    if not klasa_txt:
        return "G"
    s = klasa_txt.lower()
    if "autostr" in s:   return "A"
    if "ekspres" in s:   return "S"
    if "ruchu" in s:     return "GP"
    if "główn" in s:     return "G"
    if "zbior" in s:     return "Z"
    if "lokal" in s:     return "L"
    if "dojaz" in s:     return "D"
    if "inna" in s or "wewn" in s: return "I"
    return "G"

# Budowa grafu
def build_graph_from_fc(fc: str, where_clause: Optional[str] = WHERE) -> Tuple[int, int]:
    vertex_ids = {}
    next_vid = 1
    next_eid = 1
    fld_names_lower = {f.name.lower() for f in arcpy.ListFields(fc)}
    has_dir = FIELD_DIR_OPT.lower() in fld_names_lower
    fields = ["OBJECTID", FIELD_SHAPE, FIELD_CLASS] + ([FIELD_DIR_OPT] if has_dir else [])
    with arcpy.da.SearchCursor(fc, fields, where_clause=where_clause) as cur:
        for row in cur:
            geom = row[1]
            klasa_txt = row[2]
            kier = row[3] if has_dir else 0
            if geom is None: continue
            first_pt, last_pt = geom.firstPoint, geom.lastPoint
            if first_pt is None or last_pt is None: continue
            start_xy, end_xy = (first_pt.X, first_pt.Y), (last_pt.X, last_pt.Y)
            length_m = float(geom.length)
            for xy in [start_xy, end_xy]:
                if xy not in vertex_ids:
                    vertex_ids[xy] = next_vid
                    vertices[next_vid] = {"x": xy[0], "y": xy[1], "edge_out": []}
                    next_vid += 1
            u, v = vertex_ids[start_xy], vertex_ids[end_xy]
            klasa_code = _map_klasa_bdot(klasa_txt)
            kier = int(kier) if kier is not None else 0
            def add_edge(u_id, v_id):
                nonlocal next_eid
                edges[next_eid] = {"id": next_eid, "id_from": u_id, "id_to": v_id,
                                   "edge_length_field": length_m, "kier_auto": kier, "klasa_drogi": klasa_code}
                vertices[u_id]["edge_out"].append(next_eid)
                next_eid += 1
            if kier == 3: continue
            elif kier == 0: add_edge(u, v); add_edge(v, u)
            elif kier == 1: add_edge(u, v)
            elif kier == 2: add_edge(v, u)
            else: add_edge(u, v); add_edge(v, u)
    return len(vertices), len(edges)

if __name__ == "__main__":
    nV, nE = build_graph_from_fc(FC_ROADS, where_clause=WHERE)
    print(f"Graph built: |V|={nV}, |E|={nE}")