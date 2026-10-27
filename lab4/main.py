import random
import sys
import time

from cipolla import (
    cipolla,
    cipolla_attempts,
    decompose,
    is_prime,
    legendre,
    probabilistic_iterations,
    probabilistic_sqrt,
    roots,
)


def random_prime(bits):
    while True:
        p = random.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(p):
            return p


def random_prime_with_m(bits, m):
    """Простое p, для которого p - 1 = 2^m * q, q нечётное."""
    while True:
        q = random.getrandbits(bits - m) | (1 << (bits - m - 1)) | 1
        p = (q << m) + 1
        if is_prime(p) and decompose(p)[0] == m:
            return p


def random_residue(p):
    x = random.randrange(1, p)
    return x * x % p


def demo():
    for a, p in [(10, 13), (8, 41), (12346, 65537)]:
        print(f"Уравнение x^2 = {a} в Z_{p}, ({a}/{p}) = {legendre(a, p)}")

        print(" Алгоритм Чипполы:")
        x0 = cipolla(a, p, verbose=True)
        print(f"  x0 = {x0}, x0^2 mod p = {x0 * x0 % p}, решения: {roots(x0, p)}")

        print(" Вероятностный алгоритм:")
        y0 = probabilistic_sqrt(a, p, verbose=True)
        print(f"  x0 = {y0}, x0^2 mod p = {y0 * y0 % p}, решения: {roots(y0, p)}")
        print()

    p, a = 13, 5
    print(f"Для невычета алгоритм неприменим: ({a}/{p}) = {legendre(a, p)}, "
          f"квадраты по модулю {p}: {sorted({x * x % p for x in range(1, p)})}\n")


def test():
    for p in range(3, 400, 2):
        if not is_prime(p):
            continue
        for a in range(1, p):
            if legendre(a, p) != 1:
                continue
            for algorithm in (cipolla, probabilistic_sqrt):
                x0 = algorithm(a, p)
                assert x0 * x0 % p == a % p
                assert roots(x0, p) == sorted({x for x in range(p) if x * x % p == a % p})

    for bits in [64, 256, 512]:
        for _ in range(5):
            p = random_prime(bits)
            for _ in range(5):
                a = random_residue(p)
                assert cipolla(a, p) ** 2 % p == a
                assert probabilistic_sqrt(a, p) ** 2 % p == a

    for m in [8, 16]:
        p = random_prime_with_m(128, m)
        for _ in range(5):
            a = random_residue(p)
            assert cipolla(a, p) ** 2 % p == a
            assert probabilistic_sqrt(a, p) ** 2 % p == a

    print("Все тесты пройдены\n")


def bench_random_primes():
    print("Случайные простые, среднее время на одно извлечение корня, мс")
    print(f"{'бит':>6} {'Чиппола':>10} {'вероятн.':>10} {'m':>4}")
    for bits in [64, 128, 256, 512, 1024]:
        p = random_prime(bits)
        values = [random_residue(p) for _ in range(20)]

        start = time.perf_counter()
        for a in values:
            cipolla(a, p)
        middle = time.perf_counter()
        for a in values:
            probabilistic_sqrt(a, p)
        end = time.perf_counter()

        t1 = (middle - start) / len(values) * 1000
        t2 = (end - middle) / len(values) * 1000
        print(f"{bits:>6} {t1:>10.2f} {t2:>10.2f} {decompose(p)[0]:>4}")
    print()


def bench_by_m():
    print("Простые вида p = 2^m*q + 1 (512 бит), среднее время, мс")
    print(f"{'m':>4} {'Чиппола':>10} {'вероятн.':>10} {'итераций':>10}")
    for m in [1, 2, 8, 32, 64, 128]:
        p = random_prime_with_m(512, m)
        values = [random_residue(p) for _ in range(20)]

        start = time.perf_counter()
        for a in values:
            cipolla(a, p)
        middle = time.perf_counter()
        for a in values:
            probabilistic_sqrt(a, p)
        end = time.perf_counter()

        iterations = sum(probabilistic_iterations(a, p) for a in values) / len(values)
        t1 = (middle - start) / len(values) * 1000
        t2 = (end - middle) / len(values) * 1000
        print(f"{m:>4} {t1:>10.2f} {t2:>10.2f} {iterations:>10.1f}")
    print()


def bench_attempts():
    print("Среднее число попыток на шаге 1 алгоритма Чипполы (ожидается ~2)")
    print(f"{'бит':>6} {'попыток':>10}")
    for bits in [16, 64, 256]:
        p = random_prime(bits)
        values = [random_residue(p) for _ in range(200)]
        average = sum(cipolla_attempts(a, p) for a in values) / len(values)
        print(f"{bits:>6} {average:>10.2f}")


if __name__ == "__main__":
    demo()
    test()
    if "--no-bench" not in sys.argv:
        bench_random_primes()
        bench_by_m()
        bench_attempts()
