import random
from math import gcd

# Ответы тестов простоты (в терминологии лекции)
COMPOSITE = "N — составное"
UNKNOWN = "неизвестно"


# --- Символ Якоби (из лабораторной работы 2) ---


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


# --- Вспомогательные функции ---


def check_odd(N):
    if N < 3 or N % 2 == 0:
        raise ValueError("N должно быть нечётным числом, N >= 3")


def decompose(n):
    """n = 2^t * u, u нечётное."""
    u = n
    t = 0
    while u % 2 == 0:
        u //= 2
        t += 1
    return t, u


def factorize(n):
    """Каноническое разложение перебором делителей (эталон для проверки)."""
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def trial_division(n):
    """Детерминированная проверка простоты перебором делителей."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def euler_phi(n):
    result = n
    for p in factorize(n):
        result = result // p * (p - 1)
    return result


# --- Тест простоты на основе малой теоремы Ферма ---


def fermat_round(N, a, verbose=False):
    """Один раунд теста Ферма с основанием a."""
    check_odd(N)

    # Шаг 1. Вычислить (a, N) = d. Если d > 1 — N составное.
    d = gcd(a, N)
    if verbose:
        print(f"  Шаг 1: a = {a}, d = ({a}, {N}) = {d}")
    if d > 1:
        return COMPOSITE

    # Шаг 2. Проверить сравнение a^(N-1) = 1 (mod N).
    r = pow(a, N - 1, N)
    if verbose:
        print(f"  Шаг 2: a^(N-1) mod N = {r}")
    return UNKNOWN if r == 1 else COMPOSITE


# --- Тест Соловея–Штрассена ---


def solovay_strassen_round(N, a, verbose=False):
    """Один раунд теста Соловея-Штрассена с основанием a."""
    check_odd(N)

    # Шаг 1.
    d = gcd(a, N)
    if verbose:
        print(f"  Шаг 1: a = {a}, d = ({a}, {N}) = {d}")
    if d > 1:
        return COMPOSITE

    # Шаг 2. Проверить сравнение a^((N-1)/2) = (a/N) (mod N).
    r = pow(a, (N - 1) // 2, N)
    symbol = jacobi(a, N)
    if verbose:
        print(f"  Шаг 2: a^((N-1)/2) mod N = {r}, ({a}/{N}) = {symbol} = {symbol % N} (mod N)")
    return UNKNOWN if r == symbol % N else COMPOSITE


# --- Тест Миллера–Рабина ---


def miller_rabin_round(N, a, verbose=False):
    """Один раунд теста Миллера-Рабина с основанием a."""
    check_odd(N)

    # Шаг 1.
    d = gcd(a, N)
    if verbose:
        print(f"  Шаг 1: a = {a}, d = ({a}, {N}) = {d}")
    if d > 1:
        return COMPOSITE

    # Шаг 2. N - 1 = 2^t * u, вычислить r_k = a^(2^k * u) mod N, k = 0, ..., t-1.
    t, u = decompose(N - 1)
    r = pow(a, u, N)
    if verbose:
        print(f"  Шаг 2: N - 1 = 2^{t} * {u}, r_0 = {r}")

    if r == 1:
        return UNKNOWN
    for k in range(t):
        if r == N - 1:
            if verbose:
                print(f"         r_{k} = -1 (mod N) — ответ «{UNKNOWN}»")
            return UNKNOWN
        r = r * r % N
        if verbose and k + 1 < t:
            print(f"         r_{k + 1} = {r}")
    return COMPOSITE


ROUNDS = {
    "Ферма": fermat_round,
    "Соловея-Штрассена": solovay_strassen_round,
    "Миллера-Рабина": miller_rabin_round,
}


# --- Применение теста для s различных оснований ---


def probable_prime(round_test, N, s=10, bases=None):
    """Тест для s случайных (или заданных) оснований.

    Возвращает False, если найдено основание, доказывающее составность N,
    и True, если все проверенные условия выполнены (N простое с некоторой
    вероятностью).
    """
    if N < 2:
        return False
    if N % 2 == 0:
        return N == 2
    if N == 3:
        return True

    if bases is None:
        if N - 1 <= s:
            bases = range(1, N)
        else:
            chosen = set()
            while len(chosen) < s:
                chosen.add(random.randrange(1, N))
            bases = chosen
    for a in bases:
        if round_test(N, a) == COMPOSITE:
            return False
    return True


def random_prime(bits, round_test=miller_rabin_round, s=10):
    while True:
        N = random.getrandbits(bits) | (1 << (bits - 1)) | 1
        if probable_prime(round_test, N, s):
            return N


# --- Свидетели лжи и числа Кармайкла ---


def liars(round_test, N):
    """Основания a, при которых тест отвечает «неизвестно» (свидетели лжи)."""
    return [a for a in range(1, N) if round_test(N, a) == UNKNOWN]


def liar_fraction(round_test, N):
    """Доля оснований из {1, ..., N-1}, не доказывающих составность N."""
    return len(liars(round_test, N)) / (N - 1)


def is_carmichael(N):
    """Критерий Корсельта: N составное, свободно от квадратов и (p-1) | (N-1)."""
    if N < 2 or N % 2 == 0 or trial_division(N):
        return False
    factors = factorize(N)
    if len(factors) < 2:
        return False
    return all(power == 1 and (N - 1) % (p - 1) == 0 for p, power in factors.items())


def carmichael_numbers(limit):
    return [N for N in range(3, limit, 2) if is_carmichael(N)]
