INF = float("inf")


def divide_and_conquer_optimization(k, n, f):
    """
    分治优化DP
    f(i, j, step): 左闭右开区间 [i, j) 的代价 (0 <= i < j <= n)
    """
    # 初始化 dp 数组
    dp = [[INF] * (n + 1) for _ in range(k + 1)]
    dp[0][0] = 0

    for k_ in range(1, k + 1):

        def get_cost(y, x):
            if x >= y or dp[k_ - 1][x] == INF:
                return INF
            return dp[k_ - 1][x] + f(x, y, k_)

        res = monotoneminima(n + 1, n + 1, get_cost)
        for j in range(n + 1):
            dp[k_][j] = res[j][1]

    return dp


def monotoneminima(H, W, f):
    """
    对每个 0 <= i < H 求出 f(i, j) 取得最小值的 (j, f(i, j)) (0 <= j < W)
    """
    dp = [(0, 0)] * H  # 用元组代替 [2]int

    def dfs(top, bottom, left, right):
        if top > bottom:
            return

        mid = (top + bottom) >> 1
        index = -1
        res = 0

        for i in range(left, right + 1):
            tmp = f(mid, i)
            if index == -1 or tmp < res:
                index = i
                res = tmp

        dp[mid] = (index, res)
        dfs(top, mid - 1, left, index)
        dfs(mid + 1, bottom, index, right)

    dfs(0, H - 1, 0, W - 1)
    return dp


# https://leetcode.cn/problems/minimum-partition-score/solutions/3893471/jue-ce-dan-diao-xing-you-hua-dp-de-yi-xi-cfpw/
# class Solution:
#     def minPartitionScore(self, nums: list[int], k: int) -> int:
#         n = len(nums)
#         # 前缀和数组
#         presum = [0] * (n + 1)
#         for i, v in enumerate(nums):
#             presum[i + 1] = presum[i] + v

#         # 代价函数
#         def cost(l: int, r: int, _: int = 0) -> int:
#             s = presum[r] - presum[l]
#             return s * (s + 1) // 2

#         dp = divide_and_conquer_optimization(k, n, cost)
#         return dp[k][n]
