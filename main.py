import arcpy
from heapq import heappush, heappop
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional


# Konfiguracja
arcpy.env.workspace = r"C:\Users\cp24\Documents\temp\MyProject2\MyProject2.gdb" #tu zmien!!!
FC_ROADS = "skjzl" #tu tez zmien!!!
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
    print("[Dijkstra] |S|:", len(visited), "neighbors checked:", neighbors_checked)
    return eids

if __name__ == "__main__":
    nV, nE = build_graph_from_fc(FC_ROADS)
    print(f"Graph built: |V|={nV}, |E|={nE}")
    dijkstra(1, nV)