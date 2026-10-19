import random
import sys
import time
from math import gcd

from montgomery import Montgomery, mon_pro_base2


def random_modulus(bits, b):
    while True:
        N = random.getrandbits(bits) | (1 << (bits - 1)) | 1
        if gcd(N, b) == 1:
            return N


def demo():
    M = Montgomery(79, 10)
    print(f"N = 79, b = 10")
    print(f"n = {M.n}, R = {M.R}, R' = {M.R_inv}, N' = {M.N_prime}, R^2 mod N = {M.R2}")

    x, y = 17, 26
    print(f"\nАлгоритм 2 для x = {x}, y = {y}:")
    z = M.mon_pro(x, y, verbose=True)
    print(f"xyR^-1 mod N = {z}")
    print(f"xy mod N     = {M.mul_mod(x, y)}")
    print(f"xyR mod N    = {M.mont_mul(x, y)}")
    print(f"x^13 mod N   = {M.pow_binary(x, 13)}")
    print()


def test(rounds=500):
    for _ in range(rounds):
        b = random.choice([2, 10, 2**16, 2**32, 2**64])
        N = random_modulus(random.choice([8, 64, 256, 1024]), b)
        M = Montgomery(N, b)
        x = random.randrange(N)
        y = random.randrange(N)
        k = random.randrange(1000)

        assert M.reduce(x * y) == x * y * M.R_inv % N
        assert M.mon_pro(x, y) == x * y * M.R_inv % N
        assert M.mul_mod(x, y) == x * y % N
        assert M.mont_mul(x, y) == x * y * M.R % N
        assert M.pow_simple(x, k % 50) == pow(x, k % 50, N)
        assert M.pow_binary(x, k) == pow(x, k, N)
        if b == 2:
            assert mon_pro_base2(x, y, N) == M.mon_pro(x, y)

    print(f"Все {rounds} тестов пройдены\n")


def bench():
    b = 2**32
    N = random_modulus(1024, b)
    M = Montgomery(N, b)
    x = random.randrange(N)

    print("Возведение в степень, N — 1024 бита")
    print(f"{'k':>6} {'pow_simple, мс':>15} {'pow_binary, мс':>15}")
    for k in [10, 100, 1000, 10000]:
        start = time.perf_counter()
        M.pow_simple(x, k)
        middle = time.perf_counter()
        M.pow_binary(x, k)
        end = time.perf_counter()
        print(f"{k:>6} {(middle - start) * 1000:>15.2f} {(end - middle) * 1000:>15.2f}")


if __name__ == "__main__":
    demo()
    test()
    if "--no-bench" not in sys.argv:
        bench()
