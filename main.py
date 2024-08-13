import math
import numpy as np
from sympy.ntheory import primerange
from scipy.linalg import null_space

def generate_primes(limit):
    """Gera uma lista de primos até o limite fornecido"""
    return list(primerange(2, limit))

def is_perfect_square(N):
    """Verifica se N é um quadrado perfeito"""
    root = int(math.isqrt(N))
    return root * root == N

def quadratic_sieve(N):
    """Implementa o Crivo Quadrático para fatorar o número N"""
    if is_perfect_square(N):
        root = int(math.isqrt(N))
        return root, root, 1, N

    B = 2000  # Fixar o valor de B em 2000 para garantir primos suficientes
    primes = generate_primes(B)
    sqrtN = int(math.isqrt(N))
    smooth_numbers = []
    x_values = []
    j = 1

    while len(smooth_numbers) < len(primes) and j < 10 * B:  # Limitar o valor de j para evitar loop infinito
        x = sqrtN + j
        x2 = x * x - N
        factors = []
        for prime in primes:
            while x2 % prime == 0:
                x2 //= prime
                factors.append(prime)
        if x2 == 1:  # É um número smooth
            smooth_numbers.append(factors)
            x_values.append(x)
        j += 1

    if len(smooth_numbers) < 2:
        print(f"Não foram encontradas relações suaves suficientes para N = {N}.")
        return None, None, None, None

    # Construção da matriz A (equivalente à forma escalonada sobre Z2)
    M = len(primes)
    A = np.zeros((len(smooth_numbers), M), dtype=int)

    for i, factors in enumerate(smooth_numbers):
        for factor in factors:
            A[i, primes.index(factor)] += 1
        A[i] %= 2  # Forma escalonada sobre Z2

    # Busca do núcleo da matriz A
    ns = null_space(A)

    if ns.size == 0:
        print(f"Nenhuma solução encontrada para N = {N}.")
        return None, None, None, None

    # Extrair fatores a partir do núcleo
    for solution in ns.T:
        x_indices = np.where(solution != 0)[0]
        if len(x_indices) == 0:
            continue
        
        x = np.prod([x_values[i] for i in x_indices]) % N

        # Para y2, devemos aplainar a lista de fatores antes de multiplicar
        y2_factors = [factor for i in x_indices for factor in smooth_numbers[i]]
        y2 = np.prod(y2_factors) % N
        y = int(math.isqrt(y2))

        # Verificar se y^2 é igual a y2 (para evitar erros de raiz quadrada)
        if y * y != y2:
            continue
        
        factor1 = math.gcd(x - y, N)
        factor2 = math.gcd(x + y, N)
        if 1 < factor1 < N:
            return x, y, factor1, factor2

    return None, None, None, None

def main(N1, N2):
    """Função principal para fatorar dois números fornecidos"""
    for N in [N2, N1]:
        print(f"\nFatorando N = {N}")
        x, y, factor1, factor2 = quadratic_sieve(N)
        if x and y:
            print(f"x = {x}, y = {y}")
            print(f"mdc(x - y, N) = {factor1}")
            print(f"mdc(x + y, N) = {factor2}")
        else:
            print(f"Não foi possível fatorar N = {N} com os parâmetros atuais.")

if __name__ == "__main__":
    N1 = int(input("Digite o primeiro número N a ser fatorado: "))
    N2 = int(input("Digite o segundo número N a ser fatorado: "))
    main(N1, N2)
