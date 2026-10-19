def ext_gcd(a, b):
    x0, x1 = 1, 0
    y0, y1 = 0, 1
    while b != 0:
        q = a // b
        a, b = b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


def digits_count(x, base):
    count = 0
    while x > 0:
        x //= base
        count += 1
    return count


class Montgomery:
    def __init__(self, N, b=2**32):
        g, _, _ = ext_gcd(N, b)
        if N < 2 or g != 1:
            raise ValueError("N должно быть >= 2 и взаимно просто с b")

        self.N = N
        self.b = b
        self.n = digits_count(N, b)
        self.R = b ** self.n

        _, u, _ = ext_gcd(self.R, N)
        self.R_inv = u % N
        self.N_prime = (self.R * self.R_inv - 1) // N
        self.R2 = self.R * self.R % N

    def reduce(self, x):
        m = x * self.N_prime % self.R
        t = (x + m * self.N) // self.R
        if t >= self.N:
            t -= self.N
        return t

    def mon_pro(self, x, y, verbose=False):
        z = 0
        for i in range(self.n):
            x_i = x % self.b
            x //= self.b
            u = (z + x_i * y) % self.b
            v = u * self.N_prime % self.b
            z = (z + x_i * y + v * self.N) // self.b
            if verbose:
                print(f"  i = {i}: x_i = {x_i}, u = {u}, v = {v}, z = {z}")
        if z >= self.N:
            z -= self.N
        return z

    def phi(self, x):
        return self.mon_pro(x, 1)

    def phi_inv(self, x):
        return self.mon_pro(x, self.R2)

    def mul_mod(self, x, y):
        z = self.mon_pro(x, y)
        return self.phi_inv(z)

    def mont_mul(self, x, y):
        a = self.phi_inv(x)
        b = self.phi_inv(y)
        return self.mon_pro(a, b)

    def pow_simple(self, x, k):
        if k == 0:
            return 1
        y1 = self.phi_inv(x)
        y = y1
        for _ in range(k - 1):
            y = self.mon_pro(y, y1)
        return self.phi(y)

    def pow_binary(self, x, k):
        x_r = self.phi_inv(x)
        result = self.R % self.N
        for bit in bin(k)[2:]:
            result = self.mon_pro(result, result)
            if bit == "1":
                result = self.mon_pro(result, x_r)
        return self.phi(result)


def mon_pro_base2(x, y, N):
    z = 0
    for i in range(N.bit_length()):
        if (x >> i) & 1:
            z += y
        if z % 2 == 1:
            z += N
        z //= 2
    if z >= N:
        z -= N
    return z
