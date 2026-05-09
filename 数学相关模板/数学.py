# 常用字母表和方向向量
alp_low = "abcdefghijklmnopqrstuvwxyz"  # 小写字母
alp_up = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"  # 大写字母
dij = [[0, 1], [1, 0], [0, -1], [-1, 0]]  # 四个方向：右、下、左、上

import itertools
import platform
import sys

sys.setrecursionlimit(10**9)


def factorial(n):
    """
    计算n的阶乘，结果对mod取模
    """
    ans = 1
    for i in range(1, n + 1):
        ans *= i
        ans %= mod
    return ans


def binom(n, r):
    """
    计算二项式系数 C(n,r) = n! / (r! * (n-r)!)，结果对mod取模
    使用费马小定理计算模逆元
    """
    if r >= mod:
        raise ValueError("r is too big")
    if n < 0:
        return 0
    if r > n:
        return 0
    if r < 0:
        return 0
    ans = factorial(n) * pow(factorial(r), -1, mod) * pow(factorial(n - r), -1, mod)
    return ans % mod


# ↓以下库由 Chat GPT 生成
def hilbert_order(x: int, y: int) -> int:
    """
    计算希尔伯特曲线上点(x,y)的顺序值
    希尔伯特曲线是一种空间填充曲线，可以将2D坐标映射到1D顺序
    常用于优化查询顺序，提高缓存局部性

    参数:
        x, y: 二维坐标
    返回:
        该点在希尔伯特曲线上的位置索引
    """
    # 根据x, y的最大值确定所需的2^k
    maxc = max(x, y)
    k = maxc.bit_length()
    maxn = 1 << k

    d = 0
    s = maxn >> 1
    while s:
        rx = 1 if (x & s) else 0
        ry = 1 if (y & s) else 0

        d += s * s * ((rx * 3) ^ ry)

        if ry == 0:
            if rx == 1:
                x = maxn - 1 - x
                y = maxn - 1 - y
            x, y = y, x

        s >>= 1

    return d


def multiply_naive(a, b):
    """
    朴素多项式乘法（O(n*m)复杂度）
    计算两个多项式a和b的乘积，结果对mod取模
    """
    n = len(a)
    m = len(b)
    out_len = max(n, m)
    res = [0] * out_len
    for i in range(n):
        ai = a[i] % mod
        if ai == 0:
            continue
        for j in range(m):
            k = i + j
            if k >= out_len:
                continue
            res[k] = (res[k] + ai * b[j]) % mod
    return res


def multiply_sparse(a, b):
    """
    稀疏多项式乘法优化版本
    仅对非零项进行计算，适用于稀疏多项式
    """
    n = len(a)
    m = len(b)
    out_len = max(n, m)

    nz_a = [(i, a[i] % mod) for i in range(n) if a[i] % mod != 0]
    nz_b = [(j, b[j] % mod) for j in range(m) if b[j] % mod != 0]

    res = [0] * out_len
    for i, ai in nz_a:
        for j, bj in nz_b:
            k = i + j
            if k >= out_len:
                continue
            res[k] = (res[k] + ai * bj) % mod
    return res


def fps_div_naive(a, b, deg=-1):
    """
    朴素形式幂级数除法
    计算多项式a除以b的商，精度为deg项
    """
    if deg == -1:
        deg = max(len(a), len(b))
    if not b or b[0] % mod == 0:
        raise ZeroDivisionError("b[0] == 0 mod mod")

    inv_b0 = pow(b[0] % mod, mod - 2, mod)
    q = [0] * deg

    for n in range(deg):
        s = a[n] % mod if n < len(a) else 0
        upper = min(n, len(b) - 1)
        for i in range(1, upper + 1):
            s = (s - b[i] * q[n - i]) % mod
        q[n] = s * inv_b0 % mod

    return q


def fps_div_sparse(a, b, deg=-1):
    """
    稀疏形式幂级数除法优化版本
    仅对b的非零项进行计算
    """
    if deg == -1:
        deg = max(len(a), len(b))
    if not b or b[0] % mod == 0:
        raise ZeroDivisionError("b[0] == 0 mod mod")

    inv_b0 = pow(b[0] % mod, mod - 2, mod)

    b_nz = {i: bi % mod for i, bi in enumerate(b) if bi % mod != 0}

    q = [0] * deg
    for n in range(deg):
        s = a[n] % mod if n < len(a) else 0
        for i, bi in b_nz.items():
            if i == 0 or i > n:
                continue
            s = (s - bi * q[n - i]) % mod
        q[n] = s * inv_b0 % mod

    return q


# ↑以上库由 Chat GPT 生成

# NTT（数论变换）实现
# 参考: https://judge.yosupo.jp/submission/140558
MOD = 998244353  # 常用的NTT模数，满足 MOD = c * 2^k + 1
# NTT预计算的常数
_IMAG = 911660635  # 虚数单位i在模MOD下的表示
_IIMAG = 86583718  # 虚数单位-i在模MOD下的表示
# 预计算的旋转因子，用于加速NTT计算
_rate2 = (
    0,
    911660635,
    509520358,
    369330050,
    332049552,
    983190778,
    123842337,
    238493703,
    975955924,
    603855026,
    856644456,
    131300601,
    842657263,
    730768835,
    942482514,
    806263778,
    151565301,
    510815449,
    503497456,
    743006876,
    741047443,
    56250497,
    867605899,
    0,
)
_irate2 = (
    0,
    86583718,
    372528824,
    373294451,
    645684063,
    112220581,
    692852209,
    155456985,
    797128860,
    90816748,
    860285882,
    927414960,
    354738543,
    109331171,
    293255632,
    535113200,
    308540755,
    121186627,
    608385704,
    438932459,
    359477183,
    824071951,
    103369235,
    0,
)
_rate3 = (
    0,
    372528824,
    337190230,
    454590761,
    816400692,
    578227951,
    180142363,
    83780245,
    6597683,
    70046822,
    623238099,
    183021267,
    402682409,
    631680428,
    344509872,
    689220186,
    365017329,
    774342554,
    729444058,
    102986190,
    128751033,
    395565204,
    0,
)
_irate3 = (
    0,
    509520358,
    929031873,
    170256584,
    839780419,
    282974284,
    395914482,
    444904435,
    72135471,
    638914820,
    66769500,
    771127074,
    985925487,
    262319669,
    262341272,
    625870173,
    768022760,
    859816005,
    914661783,
    430819711,
    272774365,
    530924681,
    0,
)


def _fft(a):
    """NTT的前向变换（内部实现）"""
    n = len(a)
    h = (n - 1).bit_length()
    le = 0
    for le in range(0, h - 1, 2):
        p = 1 << (h - le - 2)
        rot = 1
        for s in range(1 << le):
            rot2 = rot * rot % MOD
            rot3 = rot2 * rot % MOD
            offset = s << (h - le)
            for i in range(p):
                a0 = a[i + offset]
                a1 = a[i + offset + p] * rot
                a2 = a[i + offset + p * 2] * rot2
                a3 = a[i + offset + p * 3] * rot3
                a1na3imag = (a1 - a3) % MOD * _IMAG
                a[i + offset] = (a0 + a2 + a1 + a3) % MOD
                a[i + offset + p] = (a0 + a2 - a1 - a3) % MOD
                a[i + offset + p * 2] = (a0 - a2 + a1na3imag) % MOD
                a[i + offset + p * 3] = (a0 - a2 - a1na3imag) % MOD
            rot = rot * _rate3[(~s & -~s).bit_length()] % MOD
    if h - le & 1:
        rot = 1
        for s in range(1 << (h - 1)):
            offset = s << 1
            l = a[offset]
            r = a[offset + 1] * rot
            a[offset] = (l + r) % MOD
            a[offset + 1] = (l - r) % MOD
            rot = rot * _rate2[(~s & -~s).bit_length()] % MOD


def _ifft(a):
    """NTT的逆向变换（内部实现）"""
    n = len(a)
    h = (n - 1).bit_length()
    le = h
    for le in range(h, 1, -2):
        p = 1 << (h - le)
        irot = 1
        for s in range(1 << (le - 2)):
            irot2 = irot * irot % MOD
            irot3 = irot2 * irot % MOD
            offset = s << (h - le + 2)
            for i in range(p):
                a0 = a[i + offset]
                a1 = a[i + offset + p]
                a2 = a[i + offset + p * 2]
                a3 = a[i + offset + p * 3]
                a2na3iimag = (a2 - a3) * _IIMAG % MOD
                a[i + offset] = (a0 + a1 + a2 + a3) % MOD
                a[i + offset + p] = (a0 - a1 + a2na3iimag) * irot % MOD
                a[i + offset + p * 2] = (a0 + a1 - a2 - a3) * irot2 % MOD
                a[i + offset + p * 3] = (a0 - a1 - a2na3iimag) * irot3 % MOD
            irot = irot * _irate3[(~s & -~s).bit_length()] % MOD
    if le & 1:
        p = 1 << (h - 1)
        for i in range(p):
            l = a[i]
            r = a[i + p]
            a[i] = l + r if l + r < MOD else l + r - MOD
            a[i + p] = l - r if l - r >= 0 else l - r + MOD


def ntt(a) -> None:
    """数论变换（Number Theoretic Transform）- FFT的模运算版本"""
    if len(a) <= 1:
        return
    _fft(a)


def intt(a) -> None:
    """逆数论变换"""
    if len(a) <= 1:
        return
    _ifft(a)
    iv = pow(len(a), MOD - 2, MOD)
    for i, x in enumerate(a):
        a[i] = x * iv % MOD


def multiply(s: list, t: list) -> list:
    """
    使用NTT进行快速多项式乘法
    对于小规模输入使用朴素算法，大规模使用NTT优化
    """
    n, m = len(s), len(t)
    l = n + m - 1
    if min(n, m) <= 60:
        a = [0] * l
        for i, x in enumerate(s):
            for j, y in enumerate(t):
                a[i + j] += x * y
        return [x % MOD for x in a]
    z = 1 << (l - 1).bit_length()
    a = s + [0] * (z - n)
    b = t + [0] * (z - m)
    _fft(a)
    _fft(b)
    for i, x in enumerate(b):
        a[i] = a[i] * x % MOD
    _ifft(a)
    a[l:] = []
    iz = pow(z, MOD - 2, MOD)
    return [x * iz % MOD for x in a]


def pow2(s: list) -> list:
    """
    计算多项式的平方（自身与自身相乘）
    优化版本，避免重复计算
    """
    n = len(s)
    l = (n << 1) - 1
    if n <= 60:
        a = [0] * l
        for i, x in enumerate(s):
            for j, y in enumerate(s):
                a[i + j] += x * y
        return [x % MOD for x in a]
    z = 1 << (l - 1).bit_length()
    a = s + [0] * (z - n)
    _fft(a)
    for i, x in enumerate(a):
        a[i] = x * x % MOD
    _ifft(a)
    a[l:] = []
    iz = pow(z, MOD - 2, MOD)
    return [x * iz % MOD for x in a]


def ntt_doubling(a: list) -> None:
    """
    NTT倍增技巧
    将NTT结果的长度加倍，用于某些高级算法
    """
    M = len(a)
    b = a[:]
    intt(b)
    r = 1
    zeta = pow(3, (MOD - 1) // (M << 1), MOD)
    for i, x in enumerate(b):
        b[i] = x * r % MOD
        r = r * zeta % MOD
    ntt(b)
    a += b


# 形式幂级数（Formal Power Series）库
# 参考: https://nyaannyaan.github.io/library/fps/formal-power-series.hpp
def shrink(a: list) -> None:
    """移除多项式尾部的零项"""
    while a and not a[-1]:
        a.pop()


def fps_add(a: list, b: list) -> list:
    """形式幂级数加法"""
    if len(a) < len(b):
        res = b[::]
        for i, x in enumerate(a):
            res[i] += x
    else:
        res = a[::]
        for i, x in enumerate(b):
            res[i] += x
    return [x % MOD for x in res]


def fps_add_scalar(a: list, k: int) -> list:
    """形式幂级数加上标量（常数项加k）"""
    res = a[:]
    res[0] = (res[0] + k) % MOD
    return res


def fps_sub(a: list, b: list) -> list:
    """形式幂级数减法"""
    if len(a) < len(b):
        res = b[::]
        for i, x in enumerate(a):
            res[i] -= x
        res = fps_neg(res)
    else:
        res = a[::]
        for i, x in enumerate(b):
            res[i] -= x
    return [x % MOD for x in res]


def fps_sub_scalar(a: list, k: int) -> list:
    """形式幂级数减去标量"""
    return fps_add_scalar(a, -k)


def fps_neg(a: list) -> list:
    """形式幂级数取负"""
    return [MOD - x if x else 0 for x in a]


def fps_mul_scalar(a: list, k: int) -> list:
    """形式幂级数乘以标量"""
    return [x * k % MOD for x in a]


def fps_matmul(a: list, b: list) -> list:
    """形式幂级数逐项相乘（未验证）"""
    return [x * b[i] % MOD for i, x in enumerate(a)]


def fps_div(a: list, b: list) -> list:
    """形式幂级数除法"""
    if len(a) < len(b):
        return []
    n = len(a) - len(b) + 1
    cnt = 0
    if len(b) > 64:
        return multiply(a[::-1][:n], fps_inv(b[::-1], n))[:n][::-1]
    f, g = a[::], b[::]
    while g and not g[-1]:
        g.pop()
        cnt += 1
    coef = pow(g[-1], MOD - 2, MOD)
    g = fps_mul_scalar(g, coef)
    deg = len(f) - len(g) + 1
    gs = len(g)
    quo = [0] * deg
    for i in range(deg)[::-1]:
        quo[i] = x = f[i + gs - 1] % MOD
        for j, y in enumerate(g):
            f[i + j] -= x * y
    return fps_mul_scalar(quo, coef) + [0] * cnt


def fps_mod(a: list, b: list) -> list:
    """形式幂级数取模运算"""
    res = fps_sub(a, multiply(fps_div(a, b), b))
    while res and not res[-1]:
        res.pop()
    return res


def fps_divmod(a: list, b: list):
    """形式幂级数除法和取模，同时返回商和余数"""
    q = fps_div(a, b)
    r = fps_sub(a, multiply(q, b))
    while r and not r[-1]:
        r.pop()
    return q, r


def fps_eval(a: list, x: int) -> int:
    """计算形式幂级数在x处的值（霍纳法则）"""
    r = 0
    w = 1
    for v in a:
        r += w * v % MOD
        w = w * x % MOD
    return r % MOD


def fps_inv(a: list, deg: int = -1) -> list:
    """
    形式幂级数求逆
    计算1/a(x) mod x^deg
    要求a[0] != 0
    """
    if deg == -1:
        deg = len(a)
    res = [0] * deg
    res[0] = pow(a[0], MOD - 2, MOD)
    d = 1
    while d < deg:
        f = [0] * (d << 1)
        tmp = min(len(a), d << 1)
        f[:tmp] = a[:tmp]
        g = [0] * (d << 1)
        g[:d] = res[:d]
        ntt(f)
        ntt(g)
        for i, x in enumerate(g):
            f[i] = f[i] * x % MOD
        intt(f)
        f[:d] = [0] * d
        ntt(f)
        for i, x in enumerate(g):
            f[i] = f[i] * x % MOD
        intt(f)
        for j in range(d, min(d << 1, deg)):
            if f[j]:
                res[j] = MOD - f[j]
            else:
                res[j] = 0
        d <<= 1
    return res


def fps_pow(a: list, k: int, deg=-1) -> list:
    """
    形式幂级数的k次幂
    计算a(x)^k mod x^deg
    """
    n = len(a)
    if deg == -1:
        deg = n
    if k == 0:
        if not deg:
            return []
        ret = [0] * deg
        ret[0] = 1
        return ret
    for i, x in enumerate(a):
        if x:
            rev = pow(x, MOD - 2, MOD)
            ret = fps_mul_scalar(
                fps_exp(
                    fps_mul_scalar(fps_log(fps_mul_scalar(a, rev)[i:], deg), k), deg
                ),
                pow(x, k, MOD),
            )
            ret[:0] = [0] * (i * k)
            if len(ret) < deg:
                ret[len(ret) :] = [0] * (deg - len(ret))
                return ret
            return ret[:deg]
        if (i + 1) * k >= deg:
            break
    return [0] * deg


def fps_exp(a: list, deg=-1) -> list:
    """
    形式幂级数的指数函数
    计算exp(a(x)) mod x^deg
    要求a[0] == 0
    """
    if deg == -1:
        deg = len(a)
    inv = [0, 1]

    def inplace_integral(F: list) -> list:
        n = len(F)
        while len(inv) <= n:
            j, k = divmod(MOD, len(inv))
            inv.append((-inv[k] * j) % MOD)
        return [0] + [x * inv[i + 1] % MOD for i, x in enumerate(F)]

    def inplace_diff(F: list) -> list:
        return [x * i % MOD for i, x in enumerate(F) if i]

    b = [1, (a[1] if 1 < len(a) else 0)]
    c = [1]
    z1 = []
    z2 = [1, 1]
    m = 2
    while m < deg:
        y = b + [0] * m
        ntt(y)
        z1 = z2
        z = [y[i] * p % MOD for i, p in enumerate(z1)]
        intt(z)
        z[: m >> 1] = [0] * (m >> 1)
        ntt(z)
        for i, p in enumerate(z1):
            z[i] = z[i] * (-p) % MOD
        intt(z)
        c[m >> 1 :] = z[m >> 1 :]
        z2 = c + [0] * m
        ntt(z2)
        tmp = min(len(a), m)
        x = a[:tmp] + [0] * (m - tmp)
        x = inplace_diff(x)
        x.append(0)
        ntt(x)
        for i, p in enumerate(x):
            x[i] = y[i] * p % MOD
        intt(x)
        for i, p in enumerate(b):
            if not i:
                continue
            x[i - 1] -= p * i % MOD
        x += [0] * m
        for i in range(m - 1):
            x[m + i], x[i] = x[i], 0
        ntt(x)
        for i, p in enumerate(z2):
            x[i] = x[i] * p % MOD
        intt(x)
        x.pop()
        x = inplace_integral(x)
        x[:m] = [0] * m
        for i in range(m, min(len(a), m << 1)):
            x[i] += a[i]
        ntt(x)
        for i, p in enumerate(y):
            x[i] = x[i] * p % MOD
        intt(x)
        b[m:] = x[m:]
        m <<= 1
    return b[:deg]


def fps_log(a: list, deg=-1) -> list:
    """
    形式幂级数的对数函数
    计算log(a(x)) mod x^deg
    要求a[0] == 1
    """
    if deg == -1:
        deg = len(a)
    return fps_integral(multiply(fps_diff(a), fps_inv(a, deg))[: deg - 1])


def fps_integral(a: list) -> list:
    """形式幂级数的积分"""
    n = len(a)
    res = [0] * (n + 1)
    if n:
        res[1] = 1
    for i in range(2, n + 1):
        j, k = divmod(MOD, i)
        res[i] = (-res[k] * j) % MOD
    for i, x in enumerate(a):
        res[i + 1] = res[i + 1] * x % MOD
    return res


def fps_diff(a: list) -> list:
    """形式幂级数的微分"""
    return [i * x % MOD for i, x in enumerate(a) if i]


# --------------------------------------------
# 试验性库（如果不喜欢可以删除）
class Cumulative_Sum:
    """
    累加和（前缀和）数据结构
    支持O(1)时间查询区间和
    """

    def __init__(self, a: list[int]):
        self.a = a
        self.rui = [0]
        for i in a:
            self.rui.append(i + self.rui[-1])

    def sum(self, l, r):
        """返回区间[l,r)的和"""
        return self.rui[r] - self.rui[l]


class Interval:
    """
    半开整数区间 [l, r)

    用法:
        Interval(r)     -> [0, r)
        Interval(l, r)  -> [l, r)

    支持的操作:
        &: 交集
        |: 并集（仅当重叠或相邻时）
        迭代
    """

    # 代码由 ChatGPT 生成
    __slots__ = ("l", "r")

    def __init__(self, l, r=None):
        if r is None:
            l, r = 0, l
        if r < l:
            l = r = l  # 空区間に正規化
        self.l = l
        self.r = r

    def __repr__(self):
        return f"[{self.l},{self.r})"

    # 迭代器支持
    def __iter__(self):
        return iter(range(self.l, self.r))

    # 交集
    def __and__(self, other):
        l = max(self.l, other.l)
        r = min(self.r, other.r)
        return Interval(l, r)

    # 并集（仅当重叠或相邻时）
    def __or__(self, other):
        if self.r < other.l or other.r < self.l:
            raise ValueError("disjoint intervals")
        return Interval(min(self.l, other.l), max(self.r, other.r))

    # 空区间判定
    def empty(self):
        return self.l >= self.r

    # 长度
    def __len__(self):
        return self.r - self.l

    # in 运算符
    def __contains__(self, x):
        return self.l <= x < self.r


# --------------------------------------------
mod = 998244353


def nin():
    """读取一行输入并转换为整数列表"""
    return list(map(int, input().split()))


def deq(x):
    """将列表中所有元素减1（用于将1-indexed转换为0-indexed）"""
    return [i - 1 for i in x]


# def main():
#     (n,) = nin()  # 读取点的数量
#     p = [(nin(), i) for i in range(n)]  # 读取n个点，格式: [([x,y], 原始索引), ...]
#     p.sort(key=lambda i: hilbert_order(i[0][0], i[0][1]))  # 按希尔伯特曲线顺序排序
#     ans = [i[1] + 1 for i in p]  # 提取排序后的索引并转换为1-indexed
#     ans.remove(1)  # 移除索引1
#     print(*([1] + ans))  # 将1放在最前面并输出


# if __name__ == "__main__":
#     main()
