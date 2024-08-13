import math

def TS(p, n):
    """Implementação do algoritmo de Tonelli-Shanks para resolver y^2 ≡ n (mod p)"""
    if int(math.pow(n, (p - 1) // 2)) % p != 1:
        return "No solutions"

    # Encontrar o maior valor de s tal que 2^s divide p-1
    s = 0
    while (p - 1) % math.pow(2, s) == 0:
        s += 1
    s -= 1
    q = int((p - 1) / math.pow(2, s))  # p-1 = q * 2^s

    # Selecionar um z tal que z seja um resíduo quadrático não-resíduo módulo p
    z = 1
    while int(math.pow(z, (p - 1) // 2)) % p != p - 1:
        z += 1

    c = int(math.pow(z, q)) % p
    r = int(math.pow(n, (q + 1) // 2)) % p
    t = int(math.pow(n, q)) % p
    m = s

    while t % p != 1:
        i = 0
        t_temp = t
        while t_temp != 1:
            t_temp = int(math.pow(t_temp, 2)) % p
            i += 1
        
        b = int(math.pow(c, int(math.pow(2, m - i - 1)))) % p
        r = (r * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        m = i

    return r

def gcd(a, b):
    """Calcula o máximo divisor comum (MDC) usando o algoritmo de Euclides"""
    while b:
        a, b = b, a % b
    return a

def find_x_y(N):
    """Calcula x e resolve a congruência para y"""
    x = math.isqrt(N) + 1  # Começar com x = ⌈√N⌉
    while True:
        n = (x * x - N) % N  # Calcular n = x^2 - N (mod N)
        y = TS(N, n)
        if isinstance(y, str):
            x += 1  # Se não encontrar solução, tentar próximo x
        else:
            return x, y

def main():
    """Função principal para calcular x, y e os MDCs para dois valores de N"""
    N1 = int(input("Digite o primeiro valor de N: "))
    N2 = int(input("Digite o segundo valor de N: "))
    
    for N in [N1, N2]:
        print(f"\nFatorando N = {N}")
        x, y = find_x_y(N)
        factor1 = gcd(x - y, N)
        factor2 = gcd(x + y, N)
        print(f"x = {x}, y = {y}")
        print(f"mdc(x - y, N) = {factor1}")
        print(f"mdc(x + y, N) = {factor2}")

if __name__ == "__main__":
    main()

