import arcpy
from heapq import heappush, heappop
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional

# Konfiguracja
arcpy.env.workspace = r"C:\Users\piotr\Documents\ArcGIS\Projects\Projekt1_PAG2\Projekt1_PAG2.gdb" #tu zmien!!!
FC_ROADS = "skjzPAG2" #tu tez zmien!!!
WHERE = None

FIELD_OID   = "OBJECTID"
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

# Funkcja do zaokrąglania współrzędnych wierzchołków
def _snap_key(x: float, y: float, tol: float) -> Tuple[int, int]:
    return (round(x / tol), round(y / tol))

# Budowa grafu
def build_graph_from_fc(
    fc: str,
    where_clause: Optional[str] = WHERE,
    snap_tol: float = 0.25  #  tolerancja 0.25 metra
) -> Tuple[int, int]:

    vertex_ids_by_snap: Dict[Tuple[int, int], int] = {}  # klucz = (_snap_key), wartość = id węzła
    next_vid = 1
    next_eid = 1

    # sprawdź, czy pole kierunkowosc istnieje
    fld_names_lower = {f.name.lower() for f in arcpy.ListFields(fc)}
    has_dir = FIELD_DIR_OPT.lower() in fld_names_lower

    fields = [FIELD_OID, FIELD_SHAPE, FIELD_CLASS] + ([FIELD_DIR_OPT] if has_dir else [])

    with arcpy.da.SearchCursor(fc, fields, where_clause=where_clause) as cur:
        for row in cur:
            jezdnia_oid = int(row[0])   # <-- ID jezdni z warstwy
            geom        = row[1]
            klasa_txt   = row[2]
            kier        = row[3] if has_dir else 0

            if geom is None:
                continue
            first_pt = geom.firstPoint
            last_pt  = geom.lastPoint
            if first_pt is None or last_pt is None:
                continue

            x1, y1 = first_pt.X, first_pt.Y
            x2, y2 = last_pt.X,  last_pt.Y
            length_m = float(geom.length)
            klasa_code = _map_klasa_bdot(klasa_txt)
            kier = int(kier) if kier is not None else 0

            # --- węzły z kwantyzacją (snap do siatki tol)
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

            def add_edge(u_id: int, v_id: int):
                nonlocal next_eid
                edges[next_eid] = {
                    "id": next_eid,
                    "id_from": u_id,
                    "id_to": v_id,
                    "edge_length_field": length_m,
                    "kier_auto": kier,
                    "klasa_drogi": klasa_code,
                    "jezdnia_oid": jezdnia_oid,   # <-- przeniesione ID jezdni do krawędzi
                }
                vertices[u_id]["edge_out"].append(next_eid)
                next_eid += 1

            # kierunki przejezdności
            if kier == 3:
                continue
            elif kier == 0:
                add_edge(u, v)
                add_edge(v, u)
            elif kier == 1:
                add_edge(u, v)
            elif kier == 2:
                add_edge(v, u)
            else:
                add_edge(u, v)
                add_edge(v, u)

    return len(vertices), len(edges)

# Kierunek (jesli bedzie potrzebne)
def czy_dobry_kierunek(direction: int, id_from: int, id_to: int, current_vertex_id: int) -> bool:
    if direction == 3: return False
    if direction == 1: return current_vertex_id == id_from
    if direction == 2: return current_vertex_id == id_to
    return True
#rekonstrukcja sciezki
def reconstruct_path(predecessors, edge_to_vertex, start, goal):
    path_nodes, path_edges, cur = [], [], goal
    while cur is not None:
        path_nodes.append(cur)
        if cur in edge_to_vertex: path_edges.append(edge_to_vertex[cur])
        cur = predecessors[cur]
    path_nodes.reverse(); path_edges.reverse()
    if not path_nodes or path_nodes[0] != start: return [], []
    return path_nodes, path_edges

# Algorytm Dijkstry
def dijkstra(start_vertex_id: int, end_vertex_id: int) -> List[int]:
    INF = float("inf")
    dist, pred = defaultdict(lambda: INF), defaultdict(lambda: None)
    visited: Set[int] = set(); edge_to_vertex = {}; neighbors_checked = 0
    dist[start_vertex_id] = 0.0; pq = [(0.0, start_vertex_id)]
    while pq:
        cur_d, u = heappop(pq)
        if u in visited: continue
        visited.add(u)
        if u == end_vertex_id: break
        for eid in vertices[u]["edge_out"]:
            e = edges[eid]; v = e["id_to"]; neighbors_checked += 1
            if not czy_dobry_kierunek(e["kier_auto"], e["id_from"], e["id_to"], u): continue
            if v in visited: continue
            nd = cur_d + e["edge_length_field"]
            if nd < dist[v]:
                dist[v] = nd; pred[v] = u; edge_to_vertex[v] = e["id"]; heappush(pq, (nd, v))
    if dist[end_vertex_id] == INF:
        print("No path found."); return []
    nodes, eids = reconstruct_path(pred, edge_to_vertex, start_vertex_id, end_vertex_id)
    print("[Dijkstra] nodes:", " -> ".join(map(str, nodes)))
    print("[Dijkstra] length [m]:", dist[end_vertex_id])
    print("[Dijkstra] |S|:", len(visited), "|| neighbors checked:", neighbors_checked)
    return eids


def export_graph_to_gdb(gdb_path: str, nodes_name: str = "nodes_out", edges_name: str = "edges_out"):

    sr = arcpy.Describe(FC_ROADS).spatialReference
    nodes_fc = f"{gdb_path}\\{nodes_name}"
    edges_fc = f"{gdb_path}\\{edges_name}"

    # Usuwa stare
    for fc in [nodes_fc, edges_fc]:
        if arcpy.Exists(fc):
            arcpy.management.Delete(fc)

    #Punkty (węzły)
    arcpy.management.CreateFeatureclass(gdb_path, nodes_name, "POINT", spatial_reference=sr)
    arcpy.management.AddField(nodes_fc, "node_id", "LONG")
    with arcpy.da.InsertCursor(nodes_fc, ["SHAPE@XY", "node_id"]) as icur:
        for vid, v in vertices.items():
            icur.insertRow(((v["x"], v["y"]), vid))
    print(f"[EXPORT] Zapisano {len(vertices)} węzłów do {nodes_fc}")

    #Linie (krawędzie)
    arcpy.management.CreateFeatureclass(gdb_path, edges_name, "POLYLINE", spatial_reference=sr)

    # atrybuty
    arcpy.management.AddField(edges_fc, "edge_id",     "LONG")
    arcpy.management.AddField(edges_fc, "id_from",     "LONG")
    arcpy.management.AddField(edges_fc, "id_to",       "LONG")
    arcpy.management.AddField(edges_fc, "length_m",    "DOUBLE")
    arcpy.management.AddField(edges_fc, "klasa",       "TEXT")
    arcpy.management.AddField(edges_fc, "kier",        "SHORT")
    arcpy.management.AddField(edges_fc, "jezdnia_oid", "LONG")  # <<<<<< DODANE

    with arcpy.da.InsertCursor(
        edges_fc,
        ["SHAPE@", "edge_id", "id_from", "id_to", "length_m", "klasa", "kier", "jezdnia_oid"]  # <<<<<< DODANE
    ) as icur:
        for eid, e in edges.items():
            u = vertices[e["id_from"]]
            v = vertices[e["id_to"]]
            arr = arcpy.Array([arcpy.Point(u["x"], u["y"]), arcpy.Point(v["x"], v["y"])])
            poly = arcpy.Polyline(arr, sr)
            icur.insertRow((
                poly,
                eid,
                e["id_from"],
                e["id_to"],
                e["edge_length_field"],
                e["klasa_drogi"],
                e["kier_auto"],
                e.get("jezdnia_oid", None)  # <<<<<< DODANE
            ))
    print(f"[EXPORT] Zapisano {len(edges)} krawędzi do {edges_fc}")


if __name__ == "__main__":
    nV, nE = build_graph_from_fc(FC_ROADS, where_clause=WHERE, snap_tol=0.25)
    print(f"Graph built: |V|={nV}, |E|={nE}")
    export_graph_to_gdb(arcpy.env.workspace)

    dijkstra(1, nV)