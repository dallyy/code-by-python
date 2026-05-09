import bisect
import heapq
import math
import operator
import os
import string
import sys
from array import array
from collections import Counter, defaultdict, deque
from copy import copy, deepcopy
from functools import cache, lru_cache
from heapq import heappop, heappush
from io import BytesIO, IOBase
from itertools import combinations, permutations
from typing import *

INF = float("inf")
BUFSIZE = 8192


class FastIO(IOBase):
    newlines = 0

    def __init__(self, file):
        self._file = file
        self._fd = file.fileno()
        self.buffer = BytesIO()
        self.writable = "x" in file.mode or "r" not in file.mode
        self.write = self.buffer.write if self.writable else None

    def read(self):
        while True:
            b = os.read(self._fd, max(os.fstat(self._fd).st_size, BUFSIZE))
            if not b:
                break
            ptr = self.buffer.tell()
            self.buffer.seek(0, 2), self.buffer.write(b), self.buffer.seek(ptr)
        self.newlines = 0
        return self.buffer.read()

    def readline(self):
        while self.newlines == 0:
            b = os.read(self._fd, max(os.fstat(self._fd).st_size, BUFSIZE))
            self.newlines = b.count(b"\n") + (not b)
            ptr = self.buffer.tell()
            self.buffer.seek(0, 2), self.buffer.write(b), self.buffer.seek(ptr)
        self.newlines -= 1
        return self.buffer.readline()

    def flush(self):
        if self.writable:
            os.write(self._fd, self.buffer.getvalue())
            self.buffer.truncate(0), self.buffer.seek(0)


class IOWrapper(IOBase):
    def __init__(self, file):
        self.buffer = FastIO(file)
        self.flush = self.buffer.flush
        self.writable = self.buffer.writable
        self.write = lambda s: self.buffer.write(s.encode("ascii"))
        self.read = lambda: self.buffer.read().decode("ascii")
        self.readline = lambda: self.buffer.readline().decode("ascii")


sys.stdin, sys.stdout = IOWrapper(sys.stdin), IOWrapper(sys.stdout)
input = sys.stdin.buffer.readline

ask = lambda *x: print("?", *x, flush=True)
reply = lambda *x: print("!", *x, flush=True)

RI = lambda: int(sys.stdin.readline())
RF = lambda: float(sys.stdin.readline())
RS = lambda: sys.stdin.readline().strip()
RFF = lambda: map(float, sys.stdin.readline().split())
RII = lambda: map(int, sys.stdin.readline().split())
RSS = lambda: map(str, sys.stdin.readline().strip().split())
RIL = lambda: list(RII())
RFL = lambda: list(RFF())
RSL = lambda: list(RSS())

from types import GeneratorType


def bootstrap(f, stack=[]):
    def wrappedfunc(*args, **kwargs):
        if stack:
            return f(*args, **kwargs)
        else:
            to = f(*args, **kwargs)
            while True:
                if type(to) is GeneratorType:
                    stack.append(to)
                    to = next(to)
                else:
                    stack.pop()
                    if not stack:
                        break
                    to = stack[-1].send(to)
            return to

    return wrappedfunc


class LazySegTree:
    def __init__(self, data, op, e, mapping, composition, id):

        self.op = op
        self.e = e
        self.mapping = mapping
        self.composition = composition
        self.id = id

        self._n = len(data)
        v = data

        self.size = 1
        while self.size < self._n:
            self.size <<= 1
        self.log = self.size.bit_length() - 1

        self.d = [self.e for _ in range(self.size << 1)]
        self.lz = [self.id for _ in range(self.size)]

        for i in range(self._n):
            self.d[self.size + i] = v[i]
        for i in range(self.size - 1, 0, -1):
            self._update(i)

    def _update(self, k):
        self.d[k] = self.op(self.d[k << 1], self.d[k << 1 | 1])

    def _apply(self, k, f):
        """对节点k应用标记f（不下推）"""
        self.d[k] = self.mapping(f, self.d[k])
        if k < self.size:
            self.lz[k] = self.composition(f, self.lz[k])

    def _push(self, k):
        """下推节点k的懒标记到子节点"""
        self._apply(k << 1, self.lz[k])
        self._apply(k << 1 | 1, self.lz[k])
        self.lz[k] = self.id

    def set(self, p, x):
        """将位置p的值设置为x"""
        assert 0 <= p < self._n
        p += self.size
        for i in range(self.log, 0, -1):
            self._push(p >> i)
        self.d[p] = x
        for i in range(1, self.log + 1):
            self._update(p >> i)

    def get(self, p):
        """获取位置p的值"""
        assert 0 <= p < self._n
        p += self.size
        for i in range(self.log, 0, -1):
            self._push(p >> i)
        return self.d[p]

    def prod(self, l, r):
        """返回区间[l, r]的聚合值（包含两端）"""
        assert 0 <= l <= r < self._n
        l += self.size
        r += self.size + 1

        for i in range(self.log, 0, -1):
            if ((l >> i) << i) != l:
                self._push(l >> i)
            if ((r >> i) << i) != r:
                self._push((r - 1) >> i)

        sml = self.e
        smr = self.e
        while l < r:
            if l & 1:
                sml = self.op(sml, self.d[l])
                l += 1
            if r & 1:
                r -= 1
                smr = self.op(self.d[r], smr)
            l >>= 1
            r >>= 1
        return self.op(sml, smr)

    def all_prod(self):
        """返回整个数组的聚合值"""
        return self.d[1]

    def apply_point(self, p, f):
        """对位置p应用标记f"""
        assert 0 <= p < self._n
        p += self.size
        for i in range(self.log, 0, -1):
            self._push(p >> i)
        self.d[p] = self.mapping(f, self.d[p])
        for i in range(1, self.log + 1):
            self._update(p >> i)

    def apply_range(self, l, r, f):
        """对区间[l, r]应用标记f"""
        assert 0 <= l <= r < self._n
        if l > r:
            return
        l += self.size
        r += self.size + 1

        for i in range(self.log, 0, -1):
            if ((l >> i) << i) != l:
                self._push(l >> i)
            if ((r >> i) << i) != r:
                self._push((r - 1) >> i)

        l2, r2 = l, r
        while l < r:
            if l & 1:
                self._apply(l, f)
                l += 1
            if r & 1:
                r -= 1
                self._apply(r, f)
            l >>= 1
            r >>= 1
        l, r = l2, r2

        for i in range(1, self.log + 1):
            if ((l >> i) << i) != l:
                self._update(l >> i)
            if ((r >> i) << i) != r:
                self._update((r - 1) >> i)

    def __repr__(self):
        return "LazySegTree({0})".format([self.get(p) for p in range(self._n)])


"""
data : 初始数组
op   : 合并函数，op(a:S, b:S) -> S
e    : 单位元 ,op(a,e)==a
mapping : 应用标记到值，mapping(f:F, a:S) -> S
composition : 合并标记，composition(f:F, g:F) -> F   (先g后f)
id   : 标记单位元，id -> F , comp(f,id)==f, op(id,a)==a
"""


def main():
    n, q = RII()

    MOD = 998244353
    e = (0, 0, 0)

    def op(a, b):
        return ((a[0] + b[0]) % MOD, (a[1] + b[1]) % MOD, a[2] + b[2])

    def map(f, a):
        return (
            (a[0] + a[2] * f) % MOD,
            (a[1] + a[2] * f * f + 2 * f * a[0]) % MOD,
            a[2],
        )

    def comp(f, g):
        return (f + g) % MOD

    id = 0

    d = LazySegTree([(0, 0, 1)] * n, op, e, map, comp, id)

    inv2 = pow(2, -1, MOD)
    for _ in range(q):
        l, r, x = RII()
        l -= 1
        r -= 1
        d.apply_range(l, r, x)
        res = d.prod(l, r)
        ans = res[0] * res[0] - res[1]
        ans *= inv2
        ans %= MOD
        print(ans)


if __name__ == "__main__":
    main()
