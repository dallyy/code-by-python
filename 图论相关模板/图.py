import bisect
import copy
import heapq
import itertools
import math
import sys
from collections import defaultdict, deque
from functools import lru_cache, reduce
from math import gcd as Gcd

read = sys.stdin.read
readline = sys.stdin.readline
readlines = sys.stdin.readlines
sys.setrecursionlimit(10**7)


class Graph:
    def __init__(self, N, edges=False, graph=False, weighted=False):
        self.N = N
        self.weighted = weighted
        if not graph:
            self.edges = edges
            self.graph = [[] for i in range(self.N)]
            if weighted:
                for i, j, d in self.edges:
                    self.graph[i].append((j, d))
            else:
                for i, j in self.edges:
                    self.graph[i].append(j)
        else:
            self.graph = graph
            self.edges = []
            for i in range(self.N):
                if self.weighted:
                    for j, d in self.graph[i]:
                        self.edges.append((i, j, d))
                else:
                    for j in self.graph[i]:
                        self.edges.append((i, j))

    def BFS(self, s):
        seen = [False] * self.N
        seen[s] = True
        queue = deque([s])
        bfs_tour = []
        while queue:
            x = queue.popleft()
            for y in self.graph[x]:
                if self.weighted:
                    y, d = y
                if not seen[y]:
                    seen[y] = True
                    queue.append(y)
                    if self.weighted:
                        bfs_tour.append((x, y, d))
                    else:
                        bfs_tour.append((x, y))
        return bfs_tour

    def DFS(self, s):
        dfs_tour = []
        stack = [(~s, 0), (s, 0)] if self.weighted else [~s, s]
        seen = [False] * self.N
        finished = [False] * self.N
        while stack:
            if self.weighted:
                x, d = stack.pop()
            else:
                x = stack.pop()
            if x >= 0 and not seen[x]:
                seen[x] = True
                if self.weighted:
                    dfs_tour.append((x, d))
                else:
                    dfs_tour.append(x)
                for y in self.graph[x]:
                    if self.weighted:
                        y, d = y
                    if not seen[y]:
                        if self.weighted:
                            stack.append((~y, d))
                            stack.append((y, d))
                        else:
                            stack.append(~y)
                            stack.append(y)
            elif x < 0 and not finished[~x]:
                if self.weighted:
                    dfs_tour.append((x, d))
                else:
                    dfs_tour.append(x)
        return dfs_tour

    def Unweighted_Dist(self, s):
        dist = [float("inf")] * self.N
        dist[s] = 0
        for tpl in self.BFS(s):
            x, y = tpl[:2] if self.weighted else tpl
            dist[y] = dist[x] + 1
        return dist

    def Tree_Diameter(self, weighted=False):
        if not weighted:
            dist = self.Unweighted_Dist(0)
        else:
            dist = [0] * self.N
            for x, y, d in self.BFS(0):
                dist[y] = dist[x] + d
        u = 0
        for i in range(self.N):
            if dist[u] < dist[i]:
                u = i
        if not weighted:
            dist = self.Unweighted_Dist(u)
        else:
            dist = [0] * self.N
            for x, y, d in self.BFS(u):
                dist[y] = dist[x] + d
        v = 0
        for i in range(self.N):
            if dist[v] < dist[i]:
                v = i
        return u, v, dist[v]

    def Parents(self, s):
        parents = [None] * self.N
        parents[s] = s
        for tpl in self.BFS(s):
            x, y = tpl[:2] if self.weighted else tpl
            parents[y] = x
        return parents

    def Bipartite_Graph(self):
        bipartite_graph = [0] * self.N
        seen = [False] * self.N
        for s in range(self.N):
            if seen[s]:
                continue
            for tpl in self.BFS(s):
                x, y = tpl[:2] if self.weighted else tpl
                bipartite_graph[y] = bipartite_graph[x] ^ 1
        for tpl in self.edges:
            i, j = tpl[:2] if self.weighted else tpl
            if not bool(bipartite_graph[i]) ^ bool(bipartite_graph[j]):
                return False
        return bipartite_graph

    def DAG(self):
        outdegree = [0] * self.N
        tp_sort = []
        queue = deque([])
        reverse_graph = [[] for i in range(self.N)]
        for tpl in self.edges:
            i, j = tpl[:2] if self.weighted else tpl
            outdegree[i] += 1
            reverse_graph[j].append(i)
        for i in range(self.N):
            if not outdegree[i]:
                queue.append(i)
                tp_sort.append(i)
        while queue:
            x = queue.popleft()
            for y in reverse_graph[x]:
                outdegree[y] -= 1
                if not outdegree[y]:
                    queue.append(y)
                    tp_sort.append(y)
        if len(tp_sort) == self.N:
            return True, tp_sort[::-1]
        else:
            is_cycle = [bool(i) for i in outdegree]
            return False, is_cycle

    def SCC(self):
        reverse_graph = [[] for i in range(self.N)]
        for tpl in self.edges:
            i, j = tpl[:2] if self.weighted else tpl
            reverse_graph[j].append(i)
        seen = [False] * self.N
        postorder = []
        for s in range(self.N):
            if seen[s]:
                continue
            lst = (
                [~i for i, d in self.DFS(s) if i < 0]
                if self.weighted
                else [~i for i in self.DFS(s) if i < 0]
            )
            for x in lst:
                if not seen[x]:
                    seen[x] = True
                    postorder.append(x)
        seen = [False] * self.N
        scc = []
        for s in postorder[::-1]:
            if seen[s]:
                continue
            queue = deque([s])
            seen[s] = True
            lst = []
            while queue:
                x = queue.popleft()
                lst.append(x)
                for y in reverse_graph[x]:
                    if self.weighted:
                        y = y[0]
                    if not seen[y]:
                        seen[y] = True
                        queue.append(y)
            scc.append(lst)
        return scc

    def Build_LCA(self, s):
        self.euler_tour = [x for x, d in self.DFS(s)] if self.weighted else self.DFS(s)
        self.dfs_in_index = [None] * self.N
        self.dfs_out_index = [None] * self.N
        depth = self.Unweighted_Dist(s)
        self.parents = self.Parents(s)
        for i, x in enumerate(self.euler_tour):
            if x >= 0:
                self.dfs_in_index[x] = i
            else:
                self.dfs_out_index[~x] = i
        self.ST = SegmentTree(2 * self.N, lambda x, y: min(x, y), float("inf"))
        lst = [None] * 2 * self.N
        for i in range(2 * self.N):
            if self.euler_tour[i] >= 0:
                lst[i] = depth[self.euler_tour[i]]
            else:
                lst[i] = depth[self.parents[~self.euler_tour[i]]]
        self.ST.Build(lst)

    def LCA(self, a, b):
        m = min(self.dfs_in_index[a], self.dfs_in_index[b])
        M = max(self.dfs_in_index[a], self.dfs_in_index[b])
        x = self.euler_tour[self.ST.Query_Index(m, M + 1)]
        if x >= 0:
            return x
        else:
            return self.parents[~x]

    def Dijkstra(self, s, route_restoration=False):
        dist = [float("inf")] * self.N
        dist[s] = 0
        hq = [(0, s)]
        if route_restoration:
            prev = [s] * self.N
        while hq:
            dx, x = heapq.heappop(hq)
            if dist[x] < dx:
                continue
            for i, di in self.graph[x]:
                if dist[i] > dx + di:
                    dist[i] = dx + di
                    if route_restoration:
                        prev[i] = x
                    heapq.heappush(hq, (dist[i], i))
        if route_restoration:
            return dist, prev
        else:
            return dist

    def Bellman_Ford(self, s, route_restoration=False):
        dist = [float("inf")] * self.N
        dist[s] = 0
        if route_restoration:
            prev = [s] * self.N
        for _ in range(self.N - 1):
            for i, j, d in self.edges:
                if dist[j] > dist[i] + d:
                    dist[j] = dist[i] + d
                    if route_restoration:
                        prev[j] = i
        negative_cycle = []
        for i, j, d in self.edges:
            if dist[j] > dist[i] + d:
                negative_cycle.append(j)
        if negative_cycle:
            is_negative_cycle = [False] * self.N
            for i in negative_cycle:
                if is_negative_cycle[i]:
                    continue
                else:
                    queue = deque([i])
                    is_negative_cycle[i] = True
                    while queue:
                        x = queue.popleft()
                        for y, d in self.graph[x]:
                            if not is_negative_cycle[y]:
                                queue.append(y)
                                is_negative_cycle[y] = True
                                if route_restoration:
                                    prev[y] = x
            for i in range(self.N):
                if is_negative_cycle[i]:
                    dist[i] = -float("inf")
        if route_restoration:
            return dist, prev
        else:
            return dist

    def Warshall_Floyd(self, route_restoration=False):
        dist = [[float("inf")] * self.N for i in range(self.N)]
        for i in range(self.N):
            dist[i][i] = 0
        if route_restoration:
            prev = [[j for j in range(self.N)] for i in range(self.N)]
        for i, j, d in self.edges:
            if dist[i][j] > d:
                dist[i][j] = d
                if route_restoration:
                    prev[i][j] = i
        for k in range(self.N):
            for i in range(self.N):
                for j in range(self.N):
                    if dist[i][j] > dist[i][k] + dist[k][j]:
                        dist[i][j] = dist[i][k] + dist[k][j]
                        if route_restoration:
                            prev[i][j] = prev[k][j]
        for i in range(self.N):
            if dist[i][i] < 0:
                for j in range(self.N):
                    if dist[i][j] != float("inf"):
                        dist[i][j] = -float("inf")
        if route_restoration:
            return dist, prev
        else:
            return dist

    def Route_Restoration(self, prev, g):
        route = []
        while prev[g] != g:
            route.append(g)
            g = prev[g]
        route = [g] + route[::-1]
        return route

    def Kruskal(self, sorted_edges=False, directed=False):
        UF = UnionFind(self.N)
        if not sorted_edges:
            self.sorted_edges = sorted(self.edges, key=lambda x: x[2])
        minimum_spnning_tree = []
        for i, j, d in self.sorted_edges:
            if not UF.Same(i, j):
                UF.Union(i, j)
                minimum_spnning_tree.append((i, j, d))
        return minimum_spnning_tree

    def Ford_Fulkerson(self, s, t):
        max_flow = 0
        graph = [defaultdict(int) for i in range(self.N)]
        residual_graph = [defaultdict(int) for i in range(self.N)]
        if self.weighted:
            for i, j, d in self.edges:
                graph[i][j] += d
                residual_graph[j][i] += d
        else:
            for i, j in self.edges:
                graph[i][j] += 1
                residual_graph[j][i] += 1

        def Flow(graph):
            flow = 0
            while True:
                prev = [None] * self.N
                prev[s] = s
                seen = [False] * self.N
                seen[s] = True
                queue = deque([s])
                while queue:
                    x = queue.popleft()
                    for y in graph[x].keys():
                        if not seen[y]:
                            seen[y] = True
                            queue.append(y)
                            prev[y] = x
                            if y == t:
                                f = float("inf")
                                t_ = t
                                while prev[t_] != t_:
                                    f = min(f, graph[prev[t_]][t_])
                                    t_ = prev[t_]
                                flow += f
                                t_ = t
                                while prev[t_] != t_:
                                    graph[prev[t_]][t_] -= f
                                    if not graph[prev[t_]][t_]:
                                        graph[prev[t_]].pop(t_)
                                    t_ = prev[t_]
                                break
                    else:
                        continue
                    break
                else:
                    break
            return flow

        max_flow += Flow(graph)
        for i in range(self.N):
            for j in graph[i].keys():
                if residual_graph[j][i] == graph[i][j]:
                    residual_graph[j].pop(i)
                else:
                    residual_graph[j][i] -= graph[i][j]
                residual_graph[i][j] += graph[i][j]
        max_flow += Flow(residual_graph)
        return max_flow


# https://atcoder.jp/contests/abc139/submissions/21523615
# N = int(readline())
# edges = []
# for i in range(N):
#     A = [a - 1 for a in list(map(int, readline().split()))]
#     for j in range(N - 2):
#         i1, i2 = i, i
#         a1, a2 = A[j], A[j + 1]
#         if i1 > a1:
#             i1, a1 = a1, i1
#         if i2 > a2:
#             i2, a2 = a2, i2
#         edges.append((i1 * N + a1, i2 * N + a2))
# N **= 2
# G = Graph(N, edges=edges)
# bl, tp_sort = G.DAG()
# if not bl:
#     ans = -1
# else:
#     dp = [1] * N
#     for i in tp_sort:
#         for j in G.graph[i]:
#             dp[j] = max(dp[j], dp[i] + 1)
#     ans = max(dp)
# print(ans)
