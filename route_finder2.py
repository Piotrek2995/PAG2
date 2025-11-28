import arcpy
from heapq import heappush, heappop
from collections import defaultdict
import math
from typing import List, Set, Tuple, Dict

# Parametry
nodes_fc = arcpy.GetParameterAsText(0)
edges_fc = arcpy.GetParameterAsText(1)
# --- Pobieranie START/END z Warstwy Nodes ---
aprx = arcpy.mp.ArcGISProject("CURRENT")
m = aprx.activeMap

# Znajdź layer odpowiadający nodes_fc
nodes_layer = None
for lyr in m.listLayers():
    if lyr.isFeatureLayer and lyr.dataSource.lower() == nodes_fc.lower():
        nodes_layer = lyr
        break

if nodes_layer is None:
    raise arcpy.ExecuteError("Nie znaleziono warstwy węzłów (nodes_fc) w aktywnej mapie.")

sel = nodes_layer.getSelectionSet()

if len(sel) < 2:
    raise arcpy.ExecuteError("Zaznacz na mapie DWIE różne węzły: start i end.")
# Sortujemy OID, żeby wynik był deterministyczny
sel = list(sel)
sel.sort()

start_vid = sel[0]
end_vid   = sel[1]

arcpy.AddMessage(f"[AUTO] Start VID = {start_vid}, End VID = {end_vid}")

algorithm = arcpy.GetParameterAsText(2)  # "Dijkstra", "A* (długość)", "A* (prędkość)"
gdb_path  = arcpy.GetParameterAsText(3)
out_fc    = arcpy.GetParameterAsText(4)

# Prędkości
SPEED_KPH = {"A":140, "S":120, "GP":90, "G":50, "Z":50, "L":50, "D":30, "I":10}
def _mps(kph): return kph * 1000.0 / 3600.0
VMAX_MPS = _mps(max(SPEED_KPH.values()))

#Wczytujemy graf z GDB
vertices: Dict[int, Dict] = {}
edges: Dict[int, Dict] = {}

with arcpy.da.SearchCursor(nodes_fc, ["node_id", "SHAPE@XY"]) as cur:
    for vid, (x, y) in cur:
        vertices[vid] = {"x": x, "y": y, "edge_out": []}

fields = ["edge_id", "id_from", "id_to", "length_m", "klasa", "kier", "jezdnia_oid"]
with arcpy.da.SearchCursor(edges_fc, fields) as cur:
    for (eid, u, v, length, cls, kier, oid) in cur:
        # Struktura jak w main.py
        edges[eid] = {
            "id": eid,
            "id_from": u,
            "id_to": v,
            "edge_length": length,
            "edge_length_field": length,   # to samo co w main
            "klasa": cls,                  # oryginalne pole z GDB
            "klasa_drogi": cls,            # alias jak w main.py
            "kier": kier,
            "kier_auto": kier,             # alias kierunku
            "jezdnia_oid": oid,            # alias jak w main.py
        }
        vertices[u]["edge_out"].append(eid)



def _euclid(a: int, b: int) -> float:
    dx = vertices[a]["x"] - vertices[b]["x"]
    dy = vertices[a]["y"] - vertices[b]["y"]
    return math.hypot(dx, dy)

#Kierunek (bezpiecznik)
def czy_dobry_kierunek(direction: int, id_from: int, id_to: int, current_vertex_id: int) -> bool:
    if direction == 3: return False
    if direction == 1: return current_vertex_id == id_from
    if direction == 2: return current_vertex_id == id_to
    return True

#Rekonstrukcja ścieżki
def reconstruct_path(predecessors, edge_to_vertex, start, goal):
    path_nodes, path_edges, cur = [], [], goal
    while cur is not None:
        path_nodes.append(cur)
        if cur in edge_to_vertex:
            path_edges.append(edge_to_vertex[cur])
        cur = predecessors[cur]
    path_nodes.reverse(); path_edges.reverse()
    if not path_nodes or path_nodes[0] != start:
        return [], []
    return path_nodes, path_edges

#Dijkstra (po długości)
def dijkstra(start_vertex_id: int, end_vertex_id: int) -> List[int]:
    INF = float("inf")
    dist: Dict[int, float] = defaultdict(lambda: INF)
    pred: Dict[int, int] = defaultdict(lambda: None)
    visited: Set[int] = set()
    edge_to_vertex: Dict[int, int] = {}
    neighbors_checked = 0

    dist[start_vertex_id] = 0.0
    pq: List[Tuple[float, int]] = [(0.0, start_vertex_id)]

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
        arcpy.AddError("Brak ścieżki (Dijkstra)")
        return []

    nodes, eids = reconstruct_path(pred, edge_to_vertex, start_vertex_id, end_vertex_id)
    arcpy.AddMessage(f"[Dijkstra] węzły: {' -> '.join(map(str, nodes))}")
    arcpy.AddMessage(f"[Dijkstra] długość [m]: {dist[end_vertex_id]:.2f}")
    arcpy.AddMessage(f"[Dijkstra] |S|: {len(visited)} || sprawdzonych sąsiadów: {neighbors_checked}")
    return eids

#A* (po długości)
def a_star_length(start_vertex_id: int, end_vertex_id: int) -> List[int]:
    INF = float("inf")
    g: Dict[int, float] = defaultdict(lambda: INF)
    pred: Dict[int, int] = defaultdict(lambda: None)
    visited: Set[int] = set()
    neighbors_checked = 0
    edge_to_vertex: Dict[int, int] = {}

    g[start_vertex_id] = 0.0
    pq_len: List[Tuple[float, int]] = []
    heappush(pq_len, (_euclid(start_vertex_id, end_vertex_id), start_vertex_id))

    while pq_len:
        _, u = heappop(pq_len)
        if u in visited: continue
        visited.add(u)
        if u == end_vertex_id: break
        gu = g[u]
        for eid in vertices[u]["edge_out"]:
            e = edges[eid]
            v = e["id_to"]
            neighbors_checked += 1
            if not czy_dobry_kierunek(e["kier_auto"], e["id_from"], e["id_to"], u):
                continue
            if v in visited: continue
            tentative = gu + e["edge_length_field"]
            if tentative < g[v]:
                g[v] = tentative
                pred[v] = u
                edge_to_vertex[v] = e["id"]
                f = tentative + _euclid(v, end_vertex_id)
                heappush(pq_len, (f, v))

    if g[end_vertex_id] == INF:
        arcpy.AddError("Brak ścieżki (A* długość)")
        return []

    nodes, eids = reconstruct_path(pred, edge_to_vertex, start_vertex_id, end_vertex_id)
    arcpy.AddMessage(f"[A* długość] węzły: {' -> '.join(map(str, nodes))}")
    arcpy.AddMessage(f"[A* długość] długość [m]: {g[end_vertex_id]:.2f}")
    arcpy.AddMessage(f"[A* długość] |S|: {len(visited)} || sprawdzonych sąsiadów: {neighbors_checked}")
    return eids

# ─── A* (po czasie / prędkości)
def a_star_speed(start_vertex_id: int, end_vertex_id: int) -> List[int]:
    """
    Minimalizuje czas przejazdu: cost(edge) = length / speed(klasa_drogi).
    Heurystyka: Euclid / VMAX_MPS (admissible i consistent) => zachowuje optymalność.
    """
    INF = float("inf")
    g_time: Dict[int, float] = defaultdict(lambda: INF)   # sekundy
    pred: Dict[int, int] = defaultdict(lambda: None)
    visited: Set[int] = set()
    neighbors_checked = 0
    edge_to_vertex: Dict[int, int] = {}

    g_time[start_vertex_id] = 0.0
    pq: List[Tuple[float, int]] = []
    heappush(pq, (_euclid(start_vertex_id, end_vertex_id) / VMAX_MPS, start_vertex_id))

    while pq:
        _, u = heappop(pq)
        if u in visited:
            continue
        visited.add(u)
        if u == end_vertex_id:
            break

        gu = g_time[u]
        for eid in vertices[u]["edge_out"]:
            e = edges[eid]
            v = e["id_to"]
            neighbors_checked += 1
            if not czy_dobry_kierunek(e["kier_auto"], e["id_from"], e["id_to"], u):
                continue
            if v in visited:
                continue

            # czas na krawędzi
            cls = e.get("klasa_drogi", "G")
            v_mps = _mps(SPEED_KPH.get(cls, SPEED_KPH["G"]))
            if v_mps <= 0:
                continue
            travel = e["edge_length_field"] / v_mps  # sekundy

            tentative = gu + travel
            if tentative < g_time[v]:
                g_time[v] = tentative
                pred[v] = u
                edge_to_vertex[v] = e["id"]
                h = _euclid(v, end_vertex_id) / VMAX_MPS
                f = tentative + h
                heappush(pq, (f, v))

    if g_time[end_vertex_id] == INF:
        arcpy.AddError("Brak ścieżki (A* prędkość)")
        return []

    nodes, eids = reconstruct_path(pred, edge_to_vertex, start_vertex_id, end_vertex_id)
    total_len = sum(edges[eid]["edge_length_field"] for eid in eids)
    arcpy.AddMessage(f"[A* prędkość] węzły: {' -> '.join(map(str, nodes))}")
    arcpy.AddMessage(f"[A* prędkość] czas [s]: {g_time[end_vertex_id]:.2f}")
    arcpy.AddMessage(f"[A* prędkość] długość [m]: {total_len:.2f}")
    arcpy.AddMessage(f"[A* prędkość] |S|: {len(visited)} || sprawdzonych sąsiadów: {neighbors_checked}")
    return eids


#Trasa alternatywna
def alternative_route(start_vertex_id: int, end_vertex_id: int, penalty_factor: float = 1.2) -> List[int]:
    """
    Wyznacza trasę alternatywną metodą kary kosztowej:
    1. Znajduje najszybszą trasę A*.
    2. Mnoży koszt każdej krawędzi tej trasy przez penalty_factor.
    3. Oblicza trasę ponownie -> otrzymujemy alternatywę.
    """
    arcpy.AddMessage("\n[Alternatywa] 1. Szukanie trasy podstawowej (A* prędkość)...")
    primary_edges = a_star_speed(start_vertex_id, end_vertex_id)

    if not primary_edges:
        arcpy.AddError("[Alternatywa] Brak trasy podstawowej – nie mogę wyznaczyć alternatywy.")
        return []

    # koszt podstawowy (czas)
    base_cost = sum(
        edges[eid]["edge_length_field"] / _mps(SPEED_KPH.get(edges[eid]["klasa_drogi"], SPEED_KPH["G"]))
        for eid in primary_edges
    )

    arcpy.AddMessage(f"[Alternatywa] Koszt podstawowy (czas): {base_cost:.2f} s")
    target_min_cost = base_cost * penalty_factor
    arcpy.AddMessage(f"[Alternatywa] Cel (po karze): {target_min_cost:.2f} s (penalty_factor={penalty_factor:.2f})")

    # 1) dodajemy karę czasową na krawędzie trasy pierwotnej
    penalty_backup = {}

    for eid in primary_edges:
        e = edges[eid]
        cls = e["klasa_drogi"]
        v_mps = _mps(SPEED_KPH.get(cls, SPEED_KPH["G"]))
        travel = e["edge_length_field"] / v_mps

        penalty_backup[eid] = travel
        penalty_travel = travel * penalty_factor

        # zapisujemy zmodyfikowany czas
        edges[eid]["__penalty_cost"] = penalty_travel

    # 2) wykonujemy A* z wykorzystaniem zmodyfikowanych kosztów
    arcpy.AddMessage("[Alternatywa] 2. Obliczanie trasy alternatywnej...")
    alt_edges = a_star_speed_with_penalty(start_vertex_id, end_vertex_id)

    # 3) przywracamy oryginalne koszty
    for eid in primary_edges:
        if "__penalty_cost" in edges[eid]:
            del edges[eid]["__penalty_cost"]

    if not alt_edges:
        arcpy.AddWarning("[Alternatywa] Nie udało się znaleźć trasy alternatywnej.")
    else:
        arcpy.AddMessage(f"[Alternatywa] Znaleziono trasę alternatywną. Liczba krawędzi: {len(alt_edges)}")

    return alt_edges

def a_star_speed_with_penalty(start_vertex_id: int, end_vertex_id: int) -> List[int]:
    INF = float("inf")
    g_time: Dict[int, float] = defaultdict(lambda: INF)
    pred: Dict[int, int] = defaultdict(lambda: None)
    visited: Set[int] = set()
    edge_to_vertex: Dict[int, int] = {}
    neighbors_checked = 0

    g_time[start_vertex_id] = 0.0
    pq: List[Tuple[float, int]] = []
    heappush(pq, (_euclid(start_vertex_id, end_vertex_id) / VMAX_MPS, start_vertex_id))

    while pq:
        _, u = heappop(pq)
        if u in visited:
            continue
        visited.add(u)
        if u == end_vertex_id:
            break

        gu = g_time[u]
        for eid in vertices[u]["edge_out"]:
            e = edges[eid]
            v = e["id_to"]
            neighbors_checked += 1
            if not czy_dobry_kierunek(e["kier_auto"], e["id_from"], e["id_to"], u):
                continue
            if v in visited:
                continue

            # koszt standardowy
            cls = e.get("klasa_drogi", "G")
            v_mps = _mps(SPEED_KPH.get(cls, SPEED_KPH["G"]))
            base_time = e["edge_length_field"] / v_mps

            # jeśli istnieje kara – użyj jej
            travel = e.get("__penalty_cost", base_time)

            tentative = gu + travel
            if tentative < g_time[v]:
                g_time[v] = tentative
                pred[v] = u
                edge_to_vertex[v] = eid
                h = _euclid(v, end_vertex_id) / VMAX_MPS
                heappush(pq, (tentative + h, v))

    if g_time[end_vertex_id] == INF:
        arcpy.AddError("Brak ścieżki (A* prędkość z karą)")
        return []

    nodes, eids = reconstruct_path(pred, edge_to_vertex, start_vertex_id, end_vertex_id)
    arcpy.AddMessage(f"[A* prędkość z karą] czas [s]: {g_time[end_vertex_id]:.2f}")
    arcpy.AddMessage(f"[A* prędkość z karą] |S|: {len(visited)} || sprawdzonych sąsiadów: {neighbors_checked}")
    return eids

# Wybór algorytmu
if algorithm == "Dijkstra":
    path_eids = dijkstra(start_vid, end_vid)
elif algorithm == "A* (długość)":
    path_eids = a_star_length(start_vid, end_vid)
elif algorithm == "A* (prędkość)":
    path_eids = a_star_speed(start_vid, end_vid)
elif algorithm == "A* (prędkość – alternatywa)":
    path_eids = alternative_route(start_vid, end_vid, penalty_factor=1.2)
else:
    arcpy.AddError(f"Nieznany algorytm: {algorithm}")
    raise SystemExit


# Tworzymy output polyline
sr = arcpy.Describe(nodes_fc).spatialReference
if arcpy.Exists(out_fc):
    arcpy.management.Delete(out_fc)

arcpy.management.CreateFeatureclass(gdb_path, out_fc.split("\\")[-1], "POLYLINE", spatial_reference=sr)
with arcpy.da.InsertCursor(out_fc, ["SHAPE@"]) as cur:
    arr = arcpy.Array()
    for eid in path_eids:
        u = vertices[edges[eid]["id_from"]]
        arr.add(arcpy.Point(u["x"], u["y"]))
    v_end = vertices[edges[path_eids[-1]]["id_to"]]
    arr.add(arcpy.Point(v_end["x"], v_end["y"]))
    cur.insertRow([arcpy.Polyline(arr, sr)])

arcpy.AddMessage("Wyznaczono trasę i zapisano do: " + out_fc)
