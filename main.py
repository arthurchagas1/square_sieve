import math   # Usado para funções matemáticas como sqrt e log.
import numpy as np  # Necessário para manipulação de arrays e solução de sistemas mod 2.


# Função para gerar os primeiros primos de maneira eficiente
def generate_primes(limit):
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for start in range(2, int(limit**0.5) + 1):
        if sieve[start]:
            for i in range(start * start, limit + 1, start):
                sieve[i] = False
    return [num for num, is_prime in enumerate(sieve) if is_prime]

def legendre_symbol(a, p):
    return pow(a, (p - 1) // 2, p)

# Algoritmo de Tonelli-Shank para resolver congruencias x2=y2 mod p
def tonelli_shanks(n, p):
    assert legendre_symbol(n, p) == 1, "n não é um resíduo quadrático módulo p"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    
    z = 2
    while legendre_symbol(z, p) != -1:
        z += 1
    
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)
    
    while t != 0 and t != 1:
        t2i = t
        i = 0
        for i in range(1, m):
            t2i = pow(t2i, 2, p)
            if t2i == 1:
                break
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = pow(b, 2, p)
        t = (t * c) % p
        r = (r * b) % p
    return r

# Crivo quadratico para fatorar primos gigantes
def quadratic_sieve(N):
    B = int(math.exp(0.5 * math.log(N) * math.log(math.log(N))))
    primes = generate_primes(B)

    # Calcula os valores de f(j) e os armazena
    def f(j):
        return (j + int(math.sqrt(N)))**2 - N

    smooth_numbers = []
    for j in range(1, B + 1):
        val = f(j)
        if all(val % p == 0 for p in primes):
            smooth_numbers.append((j, val))

    # Encontra x e y
    x = smooth_numbers[0][1]
    y = smooth_numbers[1][1]
    return x, y, math.gcd(x - y, N), math.gcd(x + y, N)

# Solução de sistemas lineares em Z2
def solve_mod_2(A, b):
    A = np.array(A) % 2
    b = np.array(b) % 2
    solution = np.linalg.solve(A, b)
    return solution % 2

def main(N1, N2):
    for N in [N1, N2]:
        print(f"\nFatorando N = {N}")
        x, y, gcd1, gcd2 = quadratic_sieve(N)
        print(f"x: {x}, y: {y}")
        print(f"mdc(x-y, N): {gcd1}")
        print(f"mdc(x+y, N): {gcd2}")

if __name__ == "__main__":
    N1 = int(input("Digite o primeiro número N a ser fatorado: "))
    N2 = int(input("Digite o segundo número N a ser fatorado: "))
    main(N1, N2)
