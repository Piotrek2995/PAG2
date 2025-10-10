import math
from urllib.parse import uses_query

import arcpy
import os
import graphlib
from collections import defaultdict


# Najlepsza implementacja grafu skierowanego to slownik list
from collections import defaultdict

def read_graph_dir(filename):
    # Wczytuje graf skierowany z pliku.
    graph = defaultdict(list)
    with open(filename, "r") as f:
        for line in f:
            u, v = line.split()
            graph[u].append(v)
    return graph

def read_graph_undir(filename):
    # Wczytuje graf nieskierowany z pliku.
    graph = defaultdict(list)
    with open(filename, "r") as f:
        for line in f:
            u, v = line.split()
            graph[u].append(v)
            graph[v].append(u)
    return graph

# Użycie

print(read_graph_dir("graph1.txt"))

print(read_graph_undir("graph2.txt"))


# liczba = int(input("Ile razy wyświetlić żart? "))
# for i in range(liczba):
#     print(f"{i+1}. Dlaczego programista nie może znaleźć żony? Bo szuka wśród zer i jedynek!")

# graf wszerz BFS

