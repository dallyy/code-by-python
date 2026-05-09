# 算法模板 (Python)

个人算法竞赛模板库，涵盖图论、数据结构、数学、常见算法等常用模板，使用 Python 编写。

## 📁 项目结构

```
.
├── 图论相关模板/
│   ├── 01bfs(边权为01求最短路).py   # 0-1 BFS 最短路
│   ├── SCG.py                        # 强连通分量 (Tarjan)
│   ├── 二分图最大权匹配.py            # KM 算法 (Kuhn-Munkres)
│   ├── 图.py                         # 图论综合模板 (BFS/DFS/Dijkstra/Bellman-Ford/Floyd/Kruskal/Ford-Fulkerson/LCA)
│   └── 拓展DSU模板.py                 # 并查集 + 二部图判定并查集
├── 常见算法/
│   ├── 倍增.py                        # 倍增法 (Doubling)
│   ├── 决策单调.py                    # 决策单调性优化 DP (分治优化)
│   └── 差分约束.py                    # 差分约束系统 (SPFA/Dijkstra)
├── 数学相关模板/
│   ├── Kitamasa(预测递推系数).py      # Kitamasa 线性递推 (NTT + BM)
│   └── 数学.py                       # 数学工具库 (FFT/NTT/多项式/组合数)
├── 数据结构/
│   └── Sgtree.py                     # 懒标记线段树 (Lazy Segment Tree)
└── 管理.py
```

## 📖 模板速览

### 图论

| 文件 | 内容 | 核心算法 |
|------|------|----------|
| `图.py` | 图论综合模板 | BFS, DFS, Dijkstra, Bellman-Ford, Floyd-Warshall, Kruskal, Ford-Fulkerson, LCA |
| `SCG.py` | 强连通分量 | Tarjan 算法 + CSR 格式存储 |
| `01bfs.py` | 0-1 最短路 | 双端队列 BFS，处理边权为 0 或 1 的最短路 |
| `二分图最大权匹配.py` | 二分图最大权匹配 | KM (Kuhn-Munkres) 算法 |
| `拓展DSU模板.py` | 并查集 + 二部图判定 | Union-Find + 带权并查集判定二部图 |

### 数据结构

| 文件 | 内容 | 核心算法 |
|------|------|----------|
| `Sgtree.py` | 懒标记线段树 | 区间更新、区间查询，支持自定义 monoid 和懒标记 |

### 常见算法

| 文件 | 内容 | 核心算法 |
|------|------|----------|
| `倍增.py` | 倍增法 | 预处理后 O(log N) 跳转，支持二分查找最大步数 |
| `决策单调.py` | 决策单调性优化 | 分治优化 DP (Divide and Conquer DP) |
| `差分约束.py` | 差分约束系统 | SPFA 求最长路/最短路，Dijkstra 优化 |

### 数学

| 文件 | 内容 | 核心算法 |
|------|------|----------|
| `数学.py` | 数学工具库 | FFT/NTT、多项式运算、组合数、阶乘 |
| `Kitamasa.py` | 线性递推 | Berlekamp-Massey + Kitamasa 快速求递推第 N 项 |

## 🚀 使用说明

每个模板文件均为独立可运行的 Python 脚本，包含：
- 完整的算法实现
- 详细的中文注释
- LeetCode / AtCoder 等平台的示例用法（注释形式）

直接导入使用或参考注释中的示例：

```python
# 示例：使用并查集
from code.数据结构.DSU import UnionFind

uf = UnionFind(10)
uf.merge(0, 1)
print(uf.issame(0, 1))  # True
```

## 🏷️ 相关链接

- [AtCoder](https://atcoder.jp/)
- [LeetCode](https://leetcode.cn/)
- [Codeforces](https://codeforces.com/)
