import sys

# 设置递归深度以防万一
sys.setrecursionlimit(200000)

MOD = 998244353
G = 3

def qpow(a, b):
    return pow(a, b, MOD)

def get_rev(cnt):
    rev = [0] * cnt
    limit = cnt.bit_length() - 1
    for i in range(cnt):
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (limit - 1))
    return rev

def ntt(a, inv):
    n = len(a)
    rev = get_rev(n)
    for i in range(n):
        if i < rev[i]:
            a[i], a[rev[i]] = a[rev[i]], a[i]
    
    length = 2
    while length <= n:
        mid = length // 2
        # 原根步长
        w_n = qpow(G if inv == 1 else qpow(G, MOD - 2), (MOD - 1) // length)
        for i in range(0, n, length):
            w = 1
            for j in range(mid):
                u = a[i + j]
                v = (w * a[i + j + mid]) % MOD
                a[i + j] = (u + v) % MOD
                a[i + j + mid] = (u - v + MOD) % MOD
                w = (w * w_n) % MOD
        length <<= 1
    
    if inv == -1:
        inv_n = qpow(n, MOD - 2)
        for i in range(n):
            a[i] = (a[i] * inv_n) % MOD

def poly_inv(f, n):
    """多项式求逆"""
    if n == 1:
        return [qpow(f[0], MOD - 2)]
    
    res = poly_inv(f, (n + 1) // 2)
    
    # 扩展长度到 2 的幂
    length = 1
    while length < (n << 1):
        length <<= 1
    
    tmp_f = f[:n] + [0] * (length - len(f[:n]))
    res += [0] * (length - len(res))
    
    ntt(tmp_f, 1)
    ntt(res, 1)
    
    for i in range(length):
        res[i] = res[i] * (2 - tmp_f[i] * res[i] % MOD + MOD) % MOD
    
    ntt(res, -1)
    return res[:n]

def poly_mul(f, g, length=None):
    """多项式乘法"""
    n = len(f)
    m = len(g)
    if length is None:
        length = 1
        while length < n + m:
            length <<= 1
    
    f_cp = f + [0] * (length - n)
    g_cp = g + [0] * (length - m)
    
    ntt(f_cp, 1)
    ntt(g_cp, 1)
    for i in range(length):
        f_cp[i] = f_cp[i] * g_cp[i] % MOD
    ntt(f_cp, -1)
    return f_cp

def poly_mod(f, g, g_inv_rev, k):
    """多项式取模 f % g, 其中 k 为 g 的度"""
    n = len(f) - 1
    if n < k:
        return f[:k]
    
    # Q_rev = F_rev * G_inv_rev
    f_rev = f[::-1]
    q_rev = poly_mul(f_rev[:n - k + 1], g_inv_rev[:n - k + 1])
    q = q_rev[:n - k + 1][::-1]
    
    # R = F - Q * G
    qg = poly_mul(q, g)
    res = [(f[i] - qg[i] + MOD) % MOD for i in range(k)]
    return res

def berlekamp_massey(arr):
    """BM算法求递推系数"""
    n = len(arr)
    best_poly = []
    last_poly = []
    best_fail_idx = -1
    delta_last = 0
    
    cur_poly = []
    for i in range(n):
        delta = arr[i]
        for j in range(len(cur_poly)):
            delta = (delta - cur_poly[j] * arr[i - j - 1]) % MOD
        
        if delta == 0:
            continue
        
        if not cur_poly:
            cur_poly = [0] * (i + 1)
            last_poly = [0]
            best_fail_idx = i
            delta_last = delta
            continue
        
        # mul = delta / delta_last
        mul = delta * qpow(delta_last, MOD - 2) % MOD
        
        # 构造新多项式
        tmp_poly = [0] * (i - best_fail_idx - 1) + [mul]
        for x in last_poly:
            tmp_poly.append((MOD - mul * x % MOD) % MOD)
            
        new_poly = cur_poly[:]
        if len(new_poly) < len(tmp_poly):
            new_poly += [0] * (len(tmp_poly) - len(new_poly))
        
        for j in range(len(tmp_poly)):
            new_poly[j] = (new_poly[j] + tmp_poly[j]) % MOD
            
        if len(cur_poly) - i < len(last_poly) - best_fail_idx:
            last_poly = cur_poly
            best_fail_idx = i
            delta_last = delta
            
        cur_poly = new_poly
        
    return cur_poly

def solve():
    # 模拟 C++ 的输入
    # 第一行: k (已知项数), n (目标项索引从0开始)
    # 第二行: 前 k 项的值
    input_data = sys.stdin.read().split()
    if not input_data: return
    
    k_in = int(input_data[0])
    n_target = int(input_data[1])
    initial_terms = [int(x) % MOD for x in input_data[2:2+k_in]]
    
    # 1. 使用 BM 算法求递推系数
    coeffs = berlekamp_massey(initial_terms)
    k = len(coeffs)
    
    # 2. 构造特征多项式 G(x) = x^k - c1*x^{k-1} - ... - ck
    # 在本算法中，我们计算 x^n mod G(x)
    g = [0] * (k + 1)
    g[k] = 1
    for i in range(k):
        g[k - 1 - i] = (MOD - coeffs[i]) % MOD
    
    # 预计算 G 的逆，用于多项式取模加速
    g_rev = g[::-1]
    g_inv_rev = poly_inv(g_rev, k + 1)
    
    # 3. 快速幂求 x^n mod G(x)
    res = [1] # 初始为 1 (多项式)
    base = [0, 1] # 初始为 x
    
    p = n_target
    while p > 0:
        if p & 1:
            res = poly_mul(res, base)
            res = poly_mod(res, g, g_inv_rev, k)
        base = poly_mul(base, base)
        base = poly_mod(base, g, g_inv_rev, k)
        p >>= 1
    
    # 4. 最终答案 = sum(res[i] * initial_terms[i])
    ans = 0
    for i in range(min(len(res), len(initial_terms))):
        ans = (ans + res[i] * initial_terms[i]) % MOD
        
    print(ans)

if __name__ == "__main__":
    solve()