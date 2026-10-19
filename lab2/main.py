import random
import sys
import time

from jacobi import is_prime, jacobi, legendre


def euler_criterion(a, p):
    r = pow(a, (p - 1) // 2, p)
    if r == 0:
        return 0
    return 1 if r == 1 else -1


def factorize(n):
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def jacobi_by_definition(a, n):
    result = 1
    for p in factorize(n):
        result *= euler_criterion(a, p)
    return result


def random_prime(bits):
    while True:
        p = random.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(p):
            return p


def demo():
    for a, p in [(219, 383), (1001, 9907)]:
        print(f"({a}/{p}):")
        value = legendre(a, p, verbose=True)
        print(f"Результат: {value}, критерий Эйлера: {euler_criterion(a, p)}\n")

    print("Символ Якоби равен 1, но сравнение x^2 = 2 (mod 15) решений не имеет:")
    squares = sorted({x * x % 15 for x in range(15)})
    print(f"(2/15) = {jacobi(2, 15)}, квадраты по модулю 15: {squares}\n")


def test():
    for p in range(3, 500, 2):
        if not is_prime(p):
            continue
        squares = {x * x % p for x in range(1, p)}
        for a in range(1, p):
            assert legendre(a, p) == euler_criterion(a, p)
            assert (legendre(a, p) == 1) == (a in squares)
            assert legendre(-a, p) == euler_criterion(-a, p)

    for n in range(3, 300, 2):
        for a in range(-n, n):
            assert jacobi(a, n) == jacobi_by_definition(a, n)

    for bits in [64, 256, 1024]:
        for _ in range(5):
            p = random_prime(bits)
            for _ in range(50):
                a = random.randrange(-p, p)
                assert jacobi(a, p) == euler_criterion(a, p)

    print("Все тесты пройдены\n")


def bench():
    print("Время на одно вычисление, мкс")
    print(f"{'бит':>6} {'Якоби':>10} {'Эйлер':>10}")
    for bits in [64, 256, 512, 1024, 2048]:
        p = random_prime(bits)
        numbers = [random.randrange(1, p) for _ in range(200)]

        start = time.perf_counter()
        for a in numbers:
            jacobi(a, p)
        middle = time.perf_counter()
        for a in numbers:
            euler_criterion(a, p)
        end = time.perf_counter()

        jacobi_time = (middle - start) / len(numbers) * 1e6
        euler_time = (end - middle) / len(numbers) * 1e6
        print(f"{bits:>6} {jacobi_time:>10.1f} {euler_time:>10.1f}")


if __name__ == "__main__":
    demo()
    test()
    if "--no-bench" not in sys.argv:
        bench()
