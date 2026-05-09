import typing
from sys import setrecursionlimit, stdin

setrecursionlimit(10**7)

input = lambda: stdin.readline().rstrip()
II = lambda: int(input())
MII = lambda: map(int, input().split())
PY = lambda: print("Yes")
PN = lambda: print("No")


class CSR:
    """压缩稀疏行格式，存储图的边"""

    def __init__(self, n: int, edges: typing.List[typing.Tuple[int, int]]):
        self.start = [0] * (n + 1)
        self.elist = [0] * len(edges)
        for u, _ in edges:
            self.start[u + 1] += 1
        for i in range(1, n + 1):
            self.start[i] += self.start[i - 1]
        counter = self.start[:]
        for u, v in edges:
            self.elist[counter[u]] = v
            counter[u] += 1


class SCCGraph:
    """Tarjan算法求强连通分量"""

    def __init__(self, n: int):
        self.n = n
        self.edges: typing.List[typing.Tuple[int, int]] = []

    def add_edge(self, u: int, v: int) -> None:
        self.edges.append((u, v))

    def scc(self) -> typing.List[typing.List[int]]:
        g = CSR(self.n, self.edges)
        now_ord = group_num = 0
        visited = []
        low = [0] * self.n
        order = [-1] * self.n
        ids = [0] * self.n

        def dfs(v: int):
            nonlocal now_ord, group_num
            low[v] = order[v] = now_ord
            now_ord += 1
            visited.append(v)
            for i in range(g.start[v], g.start[v + 1]):
                to = g.elist[i]
                if order[to] == -1:
                    dfs(to)
                    low[v] = min(low[v], low[to])
                else:
                    low[v] = min(low[v], order[to])
            if low[v] == order[v]:
                while True:
                    u = visited.pop()
                    order[u] = self.n
                    ids[u] = group_num
                    if u == v:
                        break
                group_num += 1

        for i in range(self.n):
            if order[i] == -1:
                dfs(i)

        groups = [[] for _ in range(group_num)]
        for i in range(self.n):
            groups[group_num - 1 - ids[i]].append(i)
        return groups


# def solve():
#     n, m = MII()
#     g = [[i] for i in range(n)]
#     for _ in range(m):
#         u, v = MII()
#         u -= 1
#         v -= 1
#         g[u].append(v)
#         g[v].append(u)

#     w = II()
#     s = [input() for _ in range(n)]
#     scc = SCCGraph(n * w)

#     for u in range(n):
#         for v in g[u]:
#             for i in range(w):
#                 if s[u][i] == s[v][(i + 1) % w] == "o":
#                     f = u * w + i
#                     t = v * w + (i + 1) % w
#                     scc.add_edge(f, t)

#     groups = scc.scc()
#     if max(len(g) for g in groups) > 1:
#         PY()
#     else:
#         PN()


# def main():
#     t = II()
#     for _ in range(t):
#         solve()


# if __name__ == "__main__":
#     main()
