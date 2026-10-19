def sign(power):
    return -1 if power % 2 == 1 else 1


def symbol_minus_one(n):
    return sign((n - 1) // 2)


def symbol_two(n):
    return sign((n * n - 1) // 8)


def reciprocity_sign(a, n):
    return sign((a - 1) // 2 * (n - 1) // 2)


def jacobi(a, n, verbose=False):
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
            result *= symbol_minus_one(n)
        if k % 2 == 1:
            result *= symbol_two(n)

        if verbose:
            print(f"  ({a}/{n}): b = {b} = (-1)^{t} * 2^{k} * {c}")

        result *= reciprocity_sign(c, n)
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


def legendre(a, p, verbose=False):
    if p < 3 or not is_prime(p):
        raise ValueError("p должно быть нечётным простым числом")
    return jacobi(a, p, verbose)
