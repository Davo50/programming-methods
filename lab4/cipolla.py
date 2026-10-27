import random

def sign(power):
    return -1 if power % 2 == 1 else 1


def jacobi(a, n):
    if n <= 0 or n % 2 == 0:
        raise ValueError("n должно быть нечётным натуральным числом")

    result = 1
    while n > 1:
        b = a % n
        if b > n // 2:
            b -= n
        if b == 0:
            return 0

        t = 1 if b < 0 else 0
        c = abs(b)
        k = 0
        while c % 2 == 0:
            c //= 2
            k += 1

        if t == 1:
            result *= sign((n - 1) // 2)
        if k % 2 == 1:
            result *= sign((n * n - 1) // 8)

        result *= sign((c - 1) // 2 * (n - 1) // 2)
        a, n = n, c

    return result


def is_prime(n):
    if n < 2:
        return False
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in bases:
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    for a in bases:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def legendre(a, p):
    return jacobi(a, p)


def check_input(a, p):
    if p < 3 or p % 2 == 0 or not is_prime(p):
        raise ValueError("p должно быть нечётным простым числом")
    if a % p == 0:
        raise ValueError("a должно быть взаимно просто с p")
    if legendre(a, p) != 1:
        raise ValueError(f"({a}/{p}) = -1, уравнение x^2 = {a} в Z_{p} неразрешимо")


def decompose(p):
    q = p - 1
    m = 0
    while q % 2 == 0:
        q //= 2
        m += 1
    return m, q

def poly_mul(u, v, b, a, p):
    u1, u0 = u
    v1, v0 = v
    high = u1 * v1 % p
    w1 = (u1 * v0 + u0 * v1 + high * b) % p
    w0 = (u0 * v0 - high * a) % p
    return w1, w0


def poly_pow(base, k, b, a, p):
    result = (0, 1)
    for bit in bin(k)[2:]:
        result = poly_mul(result, result, b, a, p)
        if bit == "1":
            result = poly_mul(result, base, b, a, p)
    return result


def cipolla(a, p, verbose=False):
    check_input(a, p)
    a %= p

    attempts = 0
    while True:
        b = random.randrange(p)
        d = (b * b - 4 * a) % p
        attempts += 1
        if legendre(d, p) == -1:
            break
        if verbose:
            print(f"  b = {b}: b^2 - 4a = {d}, ({d}/{p}) = {legendre(d, p)} — не подходит")

    if verbose:
        print(f"  b = {b}: b^2 - 4a = {d}, ({d}/{p}) = -1 — подходит (попыток: {attempts})")
        print(f"  Шаг 2: f(y) = y^2 - {b}y + {a}")

    x1, x0 = poly_pow((1, 0), (p + 1) // 2, b, a, p)

    if verbose:
        print(f"  Шаг 3: y^(({p}+1)/2) mod f(y) = {x1}y + {x0}")

    assert x1 == 0, "остаток не является константой"
    return x0


def cipolla_attempts(a, p):
    a %= p
    attempts = 0
    while True:
        b = random.randrange(p)
        attempts += 1
        if legendre((b * b - 4 * a) % p, p) == -1:
            return attempts


def probabilistic_sqrt(a, p, verbose=False):
    check_input(a, p)
    a %= p
    m, q = decompose(p)

    if verbose:
        print(f"  p - 1 = 2^{m} * {q}")

    while True:
        b = random.randrange(1, p)
        if legendre(b, p) == -1:
            break

    if verbose:
        print(f"  Шаг 1: b = {b}, ({b}/{p}) = -1")

    a_seq = [a]
    k_seq = []
    while True:
        a_i = a_seq[-1]
        k = 0
        t = pow(a_i, q, p)
        while t != 1:
            t = t * t % p
            k += 1
        k_seq.append(k)
        if verbose:
            print(f"  Шаг 2: a_{len(a_seq)} = {a_i}, k_{len(k_seq)} = {k}")
        if k == 0:
            break
        a_seq.append(a_i * pow(b, 1 << (m - k), p) % p)

    n = len(a_seq)
    r = pow(a_seq[n - 1], (q + 1) // 2, p)
    if verbose:
        print(f"  Шаг 3: r_{n} = {r}")
    for i in range(n - 1, 0, -1):
        factor = pow(b, 1 << (m - k_seq[i - 1] - 1), p)
        r = r * pow(factor, p - 2, p) % p
        if verbose:
            print(f"  Шаг 3: r_{i} = {r}")

    return r


def probabilistic_iterations(a, p):
    a %= p
    m, q = decompose(p)
    while True:
        b = random.randrange(1, p)
        if legendre(b, p) == -1:
            break

    a_i = a
    n = 0
    while True:
        n += 1
        k = 0
        t = pow(a_i, q, p)
        while t != 1:
            t = t * t % p
            k += 1
        if k == 0:
            return n
        a_i = a_i * pow(b, 1 << (m - k), p) % p


def sqrt_mod_p_simple(a, p):
    if p % 4 != 3:
        raise ValueError("формула применима только при p = 3 (mod 4)")
    return pow(a, (p + 1) // 4, p)


def roots(x0, p):
    return sorted({x0 % p, (-x0) % p})
