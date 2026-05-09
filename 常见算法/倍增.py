"""
倍增维护`n`个状态的转移编号和值.
每个状态对应一个`编号(0-n-1)`和`值(幺半群)`.
从状态`a`转移一次到达状态`b`，状态`b`对应的值与`边权`进行结合运算.
"""

from typing import Callable, Generic, List, TypeVar

E = TypeVar("E")


class Doubling:
    def __init__(
        self, n: int, max_step: int, e: Callable[[], E], op: Callable[[E, E], E]
    ):
        self._n = n
        self._log = 1
        while 2**self._log <= max_step:
            self._log += 1
        self._e = e
        self._op = op

        size = n * self._log
        self._to = [-1] * size
        self._dp = [e() for _ in range(size)]
        self._is_prepared = False

    def add(self, _from: int, to: int, weight: E) -> None:
        """初始状态(leaves):从 `from` 状态到 `to` 状态，边权为 `weight`. 0 <= from, to < n."""
        if self._is_prepared:
            raise RuntimeError("Doubling is prepared")
        if to < -1 or to >= self._n:
            raise ValueError("to is out of range")
        self._to[_from] = to
        self._dp[_from] = weight

    def build(self) -> None:
        if self._is_prepared:
            return
        self._is_prepared = True
        n = self._n
        for k in range(self._log - 1):
            for v in range(n):
                w = self._to[k * n + v]
                next_idx = (k + 1) * n + v
                if w == -1:
                    self._to[next_idx] = -1
                    self._dp[next_idx] = self._dp[k * n + v]
                    continue
                self._to[next_idx] = self._to[k * n + w]
                self._dp[next_idx] = self._op(self._dp[k * n + v], self._dp[k * n + w])

    def jump(self, _from: int, step: int):
        """从 `from` 状态开始，执行 `step` 次操作，返回最终状态的编号和值. 0 <= from < n. 如果最终状态不存在，返回 [-1, e()]."""
        if not self._is_prepared:
            self.build()
        if step >= 2**self._log:
            raise ValueError("step is over max step")
        value = self._e()
        to = _from
        for k in range(self._log):
            if to == -1:
                break
            div = 2**k
            if (step // div) & 1:
                pos = k * self._n + to
                value = self._op(value, self._dp[pos])
                to = self._to[pos]
        return to, value

    def max_step(self, _from: int, check: Callable[[E], bool]) -> int:
        """求从 `from` 状态开始转移 `step` 次，满足 `check` 为 `true` 的最大的 `step`. 0 <= from < n."""
        if not self._is_prepared:
            self.build()
        res = self._e()
        step = 0
        for k in range(self._log - 1, -1, -1):
            pos = k * self._n + _from
            to = self._to[pos]
            next_val = self._op(res, self._dp[pos])
            if check(next_val):
                step += 2**k
                _from = to
                res = next_val
        return step


# https://leetcode.cn/problems/maximize-value-of-function-in-a-ball-passing-game/
class Solution:
    def getMaxFunctionValue(self, receiver: list[int], k: int) -> int:
        n = len(receiver)
        db = Doubling(n, k + 10, lambda: 0, lambda a, b: a + b)
        for i in range(n):
            db.add(i, receiver[i], i)
        db.build()

        res = 0
        for i in range(n):
            _, value = db.jump(i, k + 1)
            res = max(res, value)
        return res
