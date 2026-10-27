import random
import sys
import time

from primality import (
    COMPOSITE,
    ROUNDS,
    UNKNOWN,
    carmichael_numbers,
    euler_phi,
    factorize,
    fermat_round,
    is_carmichael,
    liar_fraction,
    liars,
    miller_rabin_round,
    probable_prime,
    random_prime,
    solovay_strassen_round,
    trial_division,
)


def product(factors):
    return " * ".join(f"{p}^{k}" if k > 1 else str(p) for p, k in sorted(factors.items()))


def demo():
    print("Пошаговая работа тестов\n")
    for N, a in [(97, 35), (91, 17), (561, 2)]:
        kind = "простое" if trial_division(N) else f"составное = {product(factorize(N))}"
        print(f"N = {N} ({kind}), основание a = {a}")
        for name, round_test in ROUNDS.items():
            print(f" Тест {name}:")
            answer = round_test(N, a, verbose=True)
            print(f"  Ответ: «{answer}»")
        print()

    print("Числа Кармайкла до 10000:", carmichael_numbers(10000), "\n")

    print("Доля оснований a из {1, ..., N-1}, НЕ доказывающих составность N")
    print(f"{'N':>6} {'разложение':>16} {'Ферма':>10} {'Соловей':>10} {'Миллер':>10}")
    for N in [25, 91, 2047, 561, 1105, 1729, 2465]:
        fractions = [liar_fraction(round_test, N) for round_test in ROUNDS.values()]
        mark = " <- Кармайкл" if is_carmichael(N) else ""
        print(f"{N:>6} {product(factorize(N)):>16}"
              + "".join(f"{value:>10.3f}" for value in fractions) + mark)
    print("Теория: для Соловея-Штрассена доля <= 1/2, для Миллера-Рабина <= 1/4;")
    print("для чисел Кармайкла тест Ферма бессилен: лжецами являются все a из Z*_N.\n")

    N = 561
    print(f"Для N = {N}: phi(N) = {euler_phi(N)}, лжецов Ферма = {len(liars(fermat_round, N))},")
    print(f"то есть P0 = 1 - {euler_phi(N)}/{N - 1} = {1 - euler_phi(N) / (N - 1):.3f} "
          f"(ситуация 3 из лекции)\n")

    N, a = 2047, 2
    print(f"Сильно псевдопростое число: N = {N} = {product(factorize(N))}, a = {a}")
    print(f"  Миллер-Рабин отвечает «{miller_rabin_round(N, a)}», хотя N составное\n")


def test():
    # Корректность: тест не должен объявлять составным простое число
    for N in range(3, 3000, 2):
        prime = trial_division(N)
        for round_test in ROUNDS.values():
            for a in range(1, min(N, 60)):
                answer = round_test(N, a)
                if prime:
                    assert answer == UNKNOWN, (N, a)
                if answer == COMPOSITE:
                    assert not prime, (N, a)

    # Вложенность множеств лжецов: Миллер-Рабин <= Соловей-Штрассен <= Ферма
    for N in range(9, 1500, 2):
        if trial_division(N):
            continue
        mr = set(liars(miller_rabin_round, N))
        ss = set(liars(solovay_strassen_round, N))
        fm = set(liars(fermat_round, N))
        assert mr <= ss <= fm, N

        # Теоретические границы на долю лжецов
        assert len(ss) <= (N - 1) / 2, N
        assert len(mr) <= (N - 1) / 4, N

    # Числа Кармайкла: лжецами Ферма является вся группа Z*_N
    for N in carmichael_numbers(3000):
        assert len(liars(fermat_round, N)) == euler_phi(N), N

    # Серия раундов на случайных числах
    for _ in range(300):
        N = random.randrange(3, 20000, 2)
        expected = trial_division(N)
        assert probable_prime(miller_rabin_round, N, s=20) == expected, N
        assert probable_prime(solovay_strassen_round, N, s=30) == expected or not expected

    # Большие числа: все три теста согласованы на заведомо простых
    for bits in [64, 256, 512]:
        p = random_prime(bits)
        for round_test in ROUNDS.values():
            assert probable_prime(round_test, p, s=10)
        assert not probable_prime(miller_rabin_round, p * p, s=5)

    print("Все тесты пройдены\n")


def bench_speed():
    print("Время одного раунда теста, мс")
    print(f"{'бит':>6} {'Ферма':>10} {'Соловей':>10} {'Миллер':>10}")
    for bits in [64, 256, 512, 1024, 2048]:
        p = random_prime(bits)
        bases = [random.randrange(1, p) for _ in range(20)]
        times = []
        for round_test in ROUNDS.values():
            start = time.perf_counter()
            for a in bases:
                round_test(p, a)
            times.append((time.perf_counter() - start) / len(bases) * 1000)
        print(f"{bits:>6}" + "".join(f"{value:>10.3f}" for value in times))
    print()


def bench_error():
    print("Эмпирическая вероятность ошибки на числе Кармайкла N = 1729 (10000 запусков)")
    print(f"{'s':>4} {'Ферма':>10} {'Соловей':>10} {'Миллер':>10} {'1/2^s':>10} {'1/4^s':>10}")
    N = 1729
    attempts = 10000
    for s in [1, 2, 3, 5]:
        shares = []
        for round_test in ROUNDS.values():
            errors = sum(probable_prime(round_test, N, s) for _ in range(attempts))
            shares.append(errors / attempts)
        print(f"{s:>4}" + "".join(f"{value:>10.4f}" for value in shares)
              + f"{0.5 ** s:>10.4f}{0.25 ** s:>10.4f}")
    print("Границы 1/2^s и 1/4^s выполняются для Соловея-Штрассена и Миллера-Рабина;")
    print("тест Ферма им не подчиняется — для чисел Кармайкла его ошибка не убывает так быстро.\n")

    print("Доля оснований, не доказывающих составность, для чисел Кармайкла")
    print(f"{'N':>6} {'Ферма':>10} {'Соловей':>10} {'Миллер':>10}")
    for N in carmichael_numbers(3000):
        fractions = [liar_fraction(round_test, N) for round_test in ROUNDS.values()]
        print(f"{N:>6}" + "".join(f"{value:>10.3f}" for value in fractions))
    print()


def bench_generation():
    print("Поиск случайного простого числа (10 раундов на кандидата)")
    print(f"{'бит':>6} {'Ферма, мс':>12} {'Соловей, мс':>12} {'Миллер, мс':>12}")
    for bits in [128, 256, 512]:
        times = []
        for round_test in ROUNDS.values():
            # одинаковая последовательность кандидатов для всех трёх тестов
            random.seed(bits)
            start = time.perf_counter()
            for _ in range(5):
                random_prime(bits, round_test, s=10)
            times.append((time.perf_counter() - start) / 5 * 1000)
        print(f"{bits:>6}" + "".join(f"{value:>12.1f}" for value in times))
    random.seed()


if __name__ == "__main__":
    demo()
    test()
    if "--no-bench" not in sys.argv:
        bench_speed()
        bench_error()
        bench_generation()
