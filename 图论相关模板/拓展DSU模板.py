class UnionFind:
    """
    并查集（Union-Find）数据结构
    用于高效地管理元素的分组和合并操作
    时间复杂度：均摊 O(α(n))，其中 α 是反阿克曼函数，实际上近似 O(1)
    """

    def __init__(self, n):
        """
        初始化 n 个元素的并查集，每个元素初始时独立成组
        参数：
            n: 元素个数
        """
        self.parent = list(range(n))  # 父节点数组，初始时每个节点的父节点是自己
        self.size = [1] * n  # 每个组的元素个数

    def leader(self, x):
        """
        查找 x 所在组的根节点（代表元素）
        使用路径压缩优化，将路径上的节点直接连接到根节点
        参数：
            x: 要查找的元素
        返回：
            x 所在组的根节点
        """
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]  # 路径压缩
            x = self.parent[x]
        return x

    def merge(self, x, y):
        """
        合并 x 所在的组和 y 所在的组
        使用按大小合并优化，将小组合并到大组
        参数：
            x, y: 要合并的两个元素
        返回：
            合并后的根节点，如果 x 和 y 已在同一组则返回 -1
        """
        x, y = self.leader(x), self.leader(y)
        if x == y:
            return -1  # 已在同一组
        if self.size[x] < self.size[y]:  # 保证 x 的元素数更大
            x, y = y, x
        # 将 y 连接到 x
        self.size[x] += self.size[y]
        self.parent[y] = x
        return x

    def issame(self, x, y):
        """
        判断 x 和 y 是否在同一组
        参数：
            x, y: 要判断的两个元素
        返回：
            如果在同一组返回 True，否则返回 False
        """
        return self.leader(x) == self.leader(y)

    def getsize(self, x):
        """
        获取 x 所在组的元素个数
        参数：
            x: 要查询的元素
        返回：
            x 所在组的元素个数
        """
        return self.size[self.leader(x)]

    def all_leaders(self):
        """
        获取所有组的根节点列表
        返回：
            所有根节点的列表
        """
        return [i for i, v in enumerate(self.parent) if i == v]

    def groups(self):
        """
        获取所有连通分量（组）的列表
        返回：
            列表的列表，每个子列表包含一个组的所有元素
        """
        n = len(self.parent)
        gp = [[] for _ in range(n)]
        for v in range(n):
            gp[self.leader(v)].append(v)
        return [lst for lst in gp if lst]


# ===== END: /Suzlib/python/data_structure/unionfind/UnionFind.py =====

# ===== BEGIN: /Suzlib/python/graph/bipartite.py =====
# competitive-verifier: TITLE 二部图判定

# [bundled] from python.data_structure.unionfind.UnionFind import UnionFind


class bipartite(UnionFind):
    """
    二部图判定并查集
    用于动态维护图的二部性，支持添加边和查询二部图着色

    原理：使用 2n 个节点，对于原图的每个节点 v，创建 v 和 v+n 两个节点
          v 和 v+n 代表两种不同的颜色
          如果 u 和 v 有边相连（不同色），则合并 u 和 v+n，以及 u+n 和 v
          如果 u 和 u+n 在同一组，说明存在奇环，不是二部图

    用法示例：
        bg = bipartite(n)  # n 个顶点
        bg.add_edge(u, v)  # 添加边 (u,v)，表示 u 和 v 颜色不同
        bg.add_same(u, v)  # 添加约束，表示 u 和 v 颜色相同
        if bg.is_bipartite:  # 判断是否为二部图
            colors = bg.coloring()  # 获取着色方案
    """

    def __init__(self, n):
        """
        初始化 n 个顶点的二部图并查集
        参数：
            n: 原图的顶点数
        """
        super().__init__(2 * n)  # 创建 2n 个节点
        self.n = n
        self.is_bipartite = True  # 是否为二部图
        self.num_conn_comp = 2 * n  # 并查集上的当前连通分量数
        # 每个根节点包含的原始顶点（0..n-1）的个数
        self.orig_vertex_count = [1] * n + [0] * n

    def merge(self, x, y):
        """
        合并 x 所在的组和 y 所在的组（重写父类方法以维护额外信息）
        参数：
            x, y: 要合并的两个元素
        返回：
            合并后的根节点，如果已在同一组则返回 -1
        """
        x, y = self.leader(x), self.leader(y)
        if x == y:
            return -1
        if self.size[x] < self.size[y]:  # 保证 x 的元素数更大
            x, y = y, x
        # 将 y 连接到 x
        self.size[x] += self.size[y]
        self.parent[y] = x
        self.orig_vertex_count[x] += self.orig_vertex_count[y]
        self.num_conn_comp -= 1
        return x

    def add_edge(self, u, v):
        """
        添加边 (u, v)，表示顶点 u 和 v 颜色不同
        参数：
            u, v: 要连接的两个顶点（0-indexed）
        """
        self.merge(u, v + self.n)  # u 和 v 的另一色合并
        self.merge(u + self.n, v)  # u 的另一色和 v 合并
        if self.issame(u, u + self.n):  # 如果 u 的两个颜色在同一组，说明有矛盾
            self.is_bipartite = False

    def add_same(self, u, v):
        """
        添加约束，表示顶点 u 和 v 颜色相同
        参数：
            u, v: 颜色相同的两个顶点
        """
        self.merge(u, v)  # 同色节点合并
        self.merge(u + self.n, v + self.n)  # 另一色也合并
        if self.issame(u, u + self.n):  # 检查是否产生矛盾
            self.is_bipartite = False

    def coloring(self):
        """
        返回二部图的着色方案（0 或 1）
        返回：
            长度为 n 的列表，每个元素为 0 或 1，表示该顶点的颜色
        """
        col = [0] * self.n
        for v in range(self.n):
            col[v] = 0 if self.leader(v) < self.leader(v + self.n) else 1
        return col

    def connected_component_color_sizes(self, v):
        """
        获取包含顶点 v 的连通分量中，两种颜色的顶点个数
        参数：
            v: 查询的顶点
        返回：
            (颜色0的个数, 颜色1的个数)，颜色规约与 coloring() 相同
        """
        r0 = self.leader(v)
        r1 = self.leader(v + self.n)
        x = self.orig_vertex_count[r0]
        y = self.orig_vertex_count[r1]
        return (x, y) if r0 < r1 else (y, x)

    def all_connected_components(self):
        """
        获取所有连通分量的颜色分布信息
        返回：
            列表，每个元素为 (颜色0的个数, 颜色1的个数)
            颜色规约与 coloring() 相同
        """
        res = []
        used = [0] * (2 * self.n)
        for v in range(self.n):
            r0 = self.leader(v)
            r1 = self.leader(v + self.n)
            rep = r0 if r0 < r1 else r1
            if used[rep]:
                continue
            used[rep] = 1
            x = self.orig_vertex_count[r0]
            y = self.orig_vertex_count[r1]
            res.append((x, y) if r0 < r1 else (y, x))
        return res

    def number_of_connected_component(self):
        """
        获取当前的连通分量数
        返回：
            连通分量的个数
        时间复杂度：O(1)
        """
        return self.num_conn_comp // 2


# ===== END: /Suzlib/python/graph/bipartite.py =====


"""
使用示例：解决二部图动态连边问题

问题描述：
https://atcoder.jp/contests/abc451/tasks/abc451_f

解法：
- 使用 bipartite 类维护二部图性质
- 对于每个连通分量，删除顶点数较少的颜色集合
- 合并连通分量时，更新答案

情况1：u 和 v 已在同一连通分量
    - 只是在连通分量内部加边，不改变连通性
    - 颜色分布不变，ans 不变

情况2：u 和 v 在不同连通分量，需要合并
    - 合并前：
        * u 所在连通分量有 a1 个颜色0，b1 个颜色1，贡献 min(a1, b1)
        * v 所在连通分量有 a2 个颜色0，b2 个颜色1，贡献 min(a2, b2)
        * 当前 ans 包含了这两部分：ans = ... + min(a1,b1) + min(a2,b2)

    - 合并后：
        * 两个连通分量合并成一个
        * 新连通分量有 a 个颜色0，b 个颜色1（其中 a+b = a1+b1+a2+b2）
        * 新的贡献是 min(a, b)

    - 更新答案：
        * 先减去旧的贡献：ans -= min(a1,b1) + min(a2,b2)
        * 再加上新的贡献：ans += min(a, b)
        * 这样 ans 就正确反映了合并后所有连通分量的最小删除数之和


"""

# T = int(readline())
# for _ in range(T):
#     ans = solve()
#     print(ans)

# n = int(readline())
# a = [int(i) for i in readline().split()]
# ab = [[int(i) for i in readline().split()] for _ in range()]
# S = readline().strip()
# b = [readline().strip() for _ in range()]

# readline = sys.stdin.readline
# n, Q = [int(i) for i in readline().split()]
# UF = bipartite(n)
# ans = 0
# for _ in range(Q):
#     u, v = [int(i) for i in readline().split()]
#     u -= 1  # 转换为 0-indexed
#     v -= 1
#     if UF.issame(u, v) or UF.issame(u, v+n):
#         UF.add_edge(u, v)
#     else:
#         # u 和 v 在不同连通分量中，需要合并
#         a1, b1 = UF.connected_component_color_sizes(u)
#         a2, b2 = UF.connected_component_color_sizes(v)
#         # 减去合并前两个连通分量各自的最小删除数（减去旧贡献）
#         ans -= min(a1, b1) + min(a2, b2)
#         UF.add_edge(u, v)
#         # 加上合并后连通分量的最小删除数（加上新贡献）
#         a, b = UF.connected_component_color_sizes(v)
#         ans += min(a, b)
#         #print(a1, b1, a2, b2, a, b)
#     if not UF.is_bipartite:
#         ans = -1
#     print(ans)


# https://atcoder.jp/contests/abc447/tasks/abc447_e


# def s() -> str:
#     """在一行中输入一个字符串"""
#     return input()


# def sl() -> list[str]:
#     """在一行中输入多个字符串（以空格分隔）"""
#     return s().split()


# def ii() -> int:
#     """输入一个整数"""
#     return int(s())


# def il(add_num: int = 0) -> list[int]:
#     """在一行中输入多个整数（常用于输入数组）"""
#     return list(map(lambda i: int(i) + add_num, sl()))


# def li(n: int, func, *args: list[Any]) -> list[list[Any]]:
#     """支持多行输入（读取 n 行数据）"""
#     return [func(*args) for _ in [0] * n]


# N, M = il()
# E = li(M, il, -1)
# UF = UnionFind(N)
# ans = 0
# for i in range(M)[::-1]:
#     u, v = E[i]
#     if UF.issame(u, v):
#         continue
#     u, v = UF.leader(u), UF.leader(v)
#     if UF.getsize(u) + UF.getsize(v) == N:
#         ans += pow(2, i + 1, 998244353)
#         ans %= 998244353
#     else:
#         UF.merge(u, v)

# print(ans)
