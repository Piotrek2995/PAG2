import arcpy
from heapq import heappush, heappop
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional
import math

# Konfiguracja
arcpy.env.workspace = r"C:\Users\szymo\OneDrive\Dokumenty\ArcGIS\Projects\PAG1\PAG1.gdb"
FC_ROADS = "skjz_kopia2"
WHERE = None

FIELD_OID    = "OBJECTID"
FIELD_SHAPE  = "SHAPE@"
FIELD_CLASS  = "KLASA_DROG"
FIELD_DIR_OPT = "kierunkowosc"

# Domyślne prędkości
SPEED_KPH = {"A":140, "S":120, "GP":90, "G":50, "Z":50, "L":50, "D":30, "I":10}

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

def _euclid(v1: int, v2: int) -> float:
    dx = vertices[v1]["x"] - vertices[v2]["x"]
    dy = vertices[v1]["y"] - vertices[v2]["y"]
    return math.hypot(dx, dy)

def _snap_key(x: float, y: float, tol: float) -> Tuple[int, int]:
    return (round(x / tol), round(y / tol))

def _mps(kph: float) -> float:
    return kph * 1000.0 / 3600.0

VMAX_MPS = _mps(max(SPEED_KPH.values()))

# Budowa grafu
def build_graph_from_fc(
    fc: str,
    where_clause: Optional[str] = WHERE,
    snap_tol: float = 0.25
) -> Tuple[int, int]:

    vertex_ids_by_snap: Dict[Tuple[int, int], int] = {}
    next_vid = 1
    next_eid = 1

    fld_names_lower = {f.name.lower() for f in arcpy.ListFields(fc)}
    has_dir = FIELD_DIR_OPT.lower() in fld_names_lower

    fields = [FIELD_OID, FIELD_SHAPE, FIELD_CLASS] + ([FIELD_DIR_OPT] if has_dir else [])

    with arcpy.da.SearchCursor(fc, fields, where_clause=where_clause) as cur:
        for row in cur:
            jezdnia_oid = int(row[0])
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

            def add_edge(u_id: int, v_id: int):
                nonlocal next_eid
                edges[next_eid] = {
                    "id": next_eid,
                    "id_from": u_id,
                    "id_to": v_id,
                    "edge_length_field": length_m,
                    "kier_auto": kier,
                    "klasa_drogi": klasa_code,
                    "jezdnia_oid": jezdnia_oid,
                }
                vertices[u_id]["edge_out"].append(next_eid)
                next_eid += 1

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

# Kierunek
def czy_dobry_kierunek(direction: int, id_from: int, id_to: int, current_vertex_id: int) -> bool:
    if direction == 3: return False
    if direction == 1: return current_vertex_id == id_from
    if direction == 2: return current_vertex_id == id_to
    return True

# Rekonstrukcja ścieżki
def reconstruct_path(predecessors, edge_to_vertex, start, goal):
    path_nodes, path_edges, cur = [], [], goal
    while cur is not None:
        path_nodes.append(cur)
        if cur in edge_to_vertex: path_edges.append(edge_to_vertex[cur])
        cur = predecessors[cur]
    path_nodes.reverse(); path_edges.reverse()
    if not path_nodes or path_nodes[0] != start: return [], []
    return path_nodes, path_edges

# Dijkstra
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
    print("[Dijkstra] length [m]:", dist[end_vertex_id])
    # print("[Dijkstra] edges:", eids)
    return eids

# A* (po długości)
def a_star_length(start_vertex_id: int, end_vertex_id: int) -> List[int]:
    INF = float("inf")
    g = defaultdict(lambda: INF)
    pred = defaultdict(lambda: None)
    visited: Set[int] = set()
    neighbors_checked = 0
    edge_to_vertex: Dict[int, int] = {}

    g[start_vertex_id] = 0.0
    pq_len = []
    heappush(pq_len, (_euclid(start_vertex_id, end_vertex_id), start_vertex_id))

    while pq_len:
        _, u = heappop(pq_len)
        if u in visited: continue
        visited.add(u)
        if u == end_vertex_id: break

        gu = g[u]
        for eid in vertices[u]["edge_out"]:
            e = edges[eid]; v = e["id_to"]
            neighbors_checked += 1
            if not czy_dobry_kierunek(e["kier_auto"], e["id_from"], e["id_to"], u): continue
            if v in visited: continue
            tentative = gu + e["edge_length_field"]
            if tentative < g[v]:
                g[v] = tentative; pred[v] = u; edge_to_vertex[v] = e["id"]
                heappush(pq_len, (tentative + _euclid(v, end_vertex_id), v))

    if g[end_vertex_id] == INF:
        print("No path found.")
        return []

    nodes, eids = reconstruct_path(pred, edge_to_vertex, start_vertex_id, end_vertex_id)
    print("[A* length] length [m]:", g[end_vertex_id])
    # print("[A* length] edges:", eids)
    return eids

# A* (po czasie)
def a_star_speed(start_vertex_id: int, end_vertex_id: int) -> List[int]:
    INF = float("inf")
    g_time = defaultdict(lambda: INF)
    pred = defaultdict(lambda: None)
    visited: Set[int] = set()
    neighbors_checked = 0
    edge_to_vertex: Dict[int, int] = {}

    g_time[start_vertex_id] = 0.0
    pq = []
    heappush(pq, (_euclid(start_vertex_id, end_vertex_id) / VMAX_MPS, start_vertex_id))

    while pq:
        _, u = heappop(pq)
        if u in visited: continue
        visited.add(u)
        if u == end_vertex_id: break

        gu = g_time[u]
        for eid in vertices[u]["edge_out"]:
            e = edges[eid]; v = e["id_to"]
            neighbors_checked += 1
            if not czy_dobry_kierunek(e["kier_auto"], e["id_from"], e["id_to"], u): continue
            if v in visited: continue

            cls = e.get("klasa_drogi", "G")
            v_mps = _mps(SPEED_KPH.get(cls, SPEED_KPH["G"]))
            if v_mps <= 0:
                continue
            travel = e["edge_length_field"] / v_mps
            tentative = gu + travel

            if tentative < g_time[v]:
                g_time[v] = tentative
                pred[v] = u
                edge_to_vertex[v] = e["id"]
                f = tentative + _euclid(v, end_vertex_id) / VMAX_MPS
                heappush(pq, (f, v))

    if g_time[end_vertex_id] == INF:
        print("No path found.")
        return []

    nodes, eids = reconstruct_path(pred, edge_to_vertex, start_vertex_id, end_vertex_id)
    total_len = sum(edges[eid]["edge_length_field"] for eid in eids)
    print("[A* speed] time [s]:", g_time[end_vertex_id])
    # print("[A* speed] edges:", eids)
    return eids

# Trasa alternatywna
def alternative_route(start_vertex_id: int, end_vertex_id: int, penalty_factor: float = 1.2):
    """
    Wyznacza trasę alternatywną metodą kary kosztowej:
    1. Znajduje najszybszą trasę A*.
    2. Mnoży koszt każdej krawędzi tej trasy przez penalty_factor.
    3. Oblicza trasę ponownie -> otrzymujemy alternatywę.
    """

    print("\n1. Szukanie trasy podstawowej...")
    primary_edges = a_star_speed(start_vertex_id, end_vertex_id)

    if not primary_edges:
        print("Brak trasy podstawowej.")
        return []

    # koszt podstawowy
    base_cost = sum(edges[eid]["edge_length_field"] /
                    _mps(SPEED_KPH.get(edges[eid]["klasa_drogi"], 50))
                    for eid in primary_edges)

    print(f"Koszt podstawowy (czas): {base_cost:.2f} s")
    target_min_cost = base_cost * penalty_factor
    print(f"2. Koszt alternatywy (czas): {target_min_cost:.2f} s")

    # 1) dodajemy karę czasową na krawędzie trasy pierwotnej
    penalty_backup = {}

    for eid in primary_edges:
        e = edges[eid]
        cls = e["klasa_drogi"]
        v_mps = _mps(SPEED_KPH.get(cls, 50))
        travel = e["edge_length_field"] / v_mps

        penalty_backup[eid] = travel
        penalty_travel = travel * penalty_factor

        # zapisujemy zmodyfikowany czas
        edges[eid]["__penalty_cost"] = penalty_travel

    # 2) wykonujemy A* z wykorzystaniem zmodyfikowanych kosztów
    print("3. Obliczanie trasy alternatywnej...")

    alt_edges = a_star_speed_with_penalty(start_vertex_id, end_vertex_id)

    # 3) przywracamy oryginalne koszty
    for eid in primary_edges:
        if "__penalty_cost" in edges[eid]:
            del edges[eid]["__penalty_cost"]

    print("Trasa alternatywna:", alt_edges)
    return alt_edges

# Zmodyfikowana wersja A* obsługująca kary
def a_star_speed_with_penalty(start_vertex_id: int, end_vertex_id: int):
    INF = float("inf")
    g_time = defaultdict(lambda: INF)
    pred = defaultdict(lambda: None)
    visited: Set[int] = set()
    edge_to_vertex: Dict[int, int] = {}

    g_time[start_vertex_id] = 0.0
    pq = []
    heappush(pq, (_euclid(start_vertex_id, end_vertex_id) / VMAX_MPS, start_vertex_id))

    while pq:
        _, u = heappop(pq)
        if u in visited: continue
        visited.add(u)
        if u == end_vertex_id: break

        gu = g_time[u]
        for eid in vertices[u]["edge_out"]:
            e = edges[eid]
            v = e["id_to"]
            if not czy_dobry_kierunek(e["kier_auto"], e["id_from"], e["id_to"], u): continue
            if v in visited: continue

            # koszt standardowy
            cls = e["klasa_drogi"]
            v_mps = _mps(SPEED_KPH.get(cls, 50))
            base_time = e["edge_length_field"] / v_mps

            # jeśli istnieje kara – użyj jej
            travel = e.get("__penalty_cost", base_time)

            tentative = gu + travel
            if tentative < g_time[v]:
                g_time[v] = tentative
                pred[v] = u
                edge_to_vertex[v] = eid
                heappush(pq, (tentative + _euclid(v, end_vertex_id) / VMAX_MPS, v))

    if g_time[end_vertex_id] == INF:
        print("Nie udało się znaleźć trasy alternatywnej.")
        return []

    _, eids = reconstruct_path(pred, edge_to_vertex, start_vertex_id, end_vertex_id)
    return eids

# Eksport grafu
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

    arcpy.management.CreateFeatureclass(gdb_path, edges_name, "POLYLINE", spatial_reference=sr)
    arcpy.management.AddField(edges_fc, "edge_id",     "LONG")
    arcpy.management.AddField(edges_fc, "id_from",     "LONG")
    arcpy.management.AddField(edges_fc, "id_to",       "LONG")
    arcpy.management.AddField(edges_fc, "length_m",    "DOUBLE")
    arcpy.management.AddField(edges_fc, "klasa",       "TEXT")
    arcpy.management.AddField(edges_fc, "kier",        "SHORT")
    arcpy.management.AddField(edges_fc, "jezdnia_oid", "LONG")

    with arcpy.da.InsertCursor(
        edges_fc,
        ["SHAPE@", "edge_id", "id_from", "id_to", "length_m", "klasa", "kier", "jezdnia_oid"]
    ) as icur:
        for eid, e in edges.items():
            u = vertices[e["id_from"]]
            v = vertices[e["id_to"]]
            arr = arcpy.Array([arcpy.Point(u["x"], u["y"]), arcpy.Point(v["x"], v["y"])])
            poly = arcpy.Polyline(arr, sr)
            icur.insertRow((
                poly, eid, e["id_from"], e["id_to"],
                e["edge_length_field"], e["klasa_drogi"], e["kier_auto"], e["jezdnia_oid"]
            ))

# main
if __name__ == "__main__":
    nV, nE = build_graph_from_fc(FC_ROADS, where_clause=WHERE, snap_tol=0.25)
    print(f"Graph built: |V|={nV}, |E|={nE}")

    export_graph_to_gdb(arcpy.env.workspace)

    start = 1
    goal = nV

    dijkstra(start, goal)
    a_star_length(start, goal)
    a_star_speed(start, goal)

    print("\n--- TRASA ALTERNATYWNA ---")
    alternative_route(start, goal, penalty_factor=1.2)
