import numpy as np
import matplotlib.pyplot as plt
from math import gcd, log
import os

# Create plots folder if needed
os.makedirs('plots', exist_ok=True)

def get_primes(n):
    if n < 2:
        return np.array([])
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(np.sqrt(n)) + 1):
        if sieve[i]:
            for j in range(i * i, n + 1, i):
                sieve[j] = False
    return np.array([i for i in range(2, n + 1) if sieve[i]])

max_n = 10**6  # Increase to 10**7 for better accuracy (slower)
primes = get_primes(max_n)
log_primes = np.log(primes.astype(float))

def chebyshev_psi(x):
    if x < 2:
        return 0.0
    psi = 0.0
    for i, p in enumerate(primes):
        if p > x:
            break
        lp = log_primes[i]
        max_k = int(log(x) / log(p))
        psi += lp * max_k
    return psi

# Experiment A: Light-cone scaling
t_values = np.linspace(2, np.log(max_n), 200)
F, G, H = [], [], []
for t in t_values:
    x = np.exp(t)
    psiv = chebyshev_psi(x)
    f = psiv - x
    g = np.exp(-t / 2) * f
    h = g / (t ** 2)
    F.append(f)
    G.append(g)
    H.append(h)

print('G stats:', np.min(G), np.max(G), np.mean(G), np.std(G))
print('H stats:', np.min(H), np.max(H), np.mean(H), np.std(H))
variances = [np.var(H[i:i+20]) for i in range(len(H)-20)]
print('Mean sliding variance of H:', np.mean(variances))

plt.figure(figsize=(10, 8))
plt.subplot(3, 1, 1)
plt.plot(t_values, F, label='F(t)')
plt.title('Prime Fluctuations')
plt.ylabel('F(t)')
plt.subplot(3, 1, 2)
plt.plot(t_values, G, label='G(t)')
plt.ylabel('G(t)')
plt.subplot(3, 1, 3)
plt.plot(t_values, H, label='H(t)')
plt.ylabel('H(t)')
plt.xlabel('t = log x')
plt.tight_layout()
plt.savefig('plots/light_cone.png')
plt.close()

# Experiment B: Base-collapse
bases = [3, 8, 10, 12, 30, 60, 101]
base_data = {}
for b in bases:
    k_max = int(np.log(max_n) / np.log(b))
    data = []
    coprimes = [r for r in range(b) if gcd(r, b) == 1]
    phi_b = len(coprimes)
    uniform = 1.0 / phi_b
    prime_factors_b = {p for p in primes if p <= b and b % p == 0}
    for k in range(1, k_max + 1):
        N = b ** k
        primes_upto_N = primes[primes <= N]
        relevant_primes = [p for p in primes_upto_N if p not in prime_factors_b]
        if len(relevant_primes) < 10:
            continue
        residues = [p % b for p in relevant_primes]
        count_arr = np.bincount(residues, minlength=b)
        total = len(relevant_primes)
        tv = 0.5 * sum(abs(count_arr[r] / total - uniform) for r in coprimes)
        t = np.log(N)
        data.append((t, tv))
    base_data[b] = data

# Compute collapse score
all_points = [p for data in base_data.values() for p in data]
all_points.sort(key=lambda x: x[0])
bin_edges = np.linspace(min(p[0] for p in all_points), max(p[0] for p in all_points), 11)
bin_means, bin_counts = [], []
for i in range(len(bin_edges) - 1):
    vals = [p[1] for p in all_points if bin_edges[i] <= p[0] < bin_edges[i+1]]
    bin_means.append(np.mean(vals) if vals else 0)
    bin_counts.append(len(vals))
deviations = []
for p in all_points:
    for i in range(len(bin_edges) - 1):
        if bin_edges[i] <= p[0] < bin_edges[i+1] and bin_counts[i] > 1:
            deviations.append((p[1] - bin_means[i]) ** 2)
            break
rms = np.sqrt(np.mean(deviations)) if deviations else 0
print('Collapse RMS score:', rms)

# Plot
plt.figure(figsize=(10, 6))
for b, data in base_data.items():
    if data:
        t, d = zip(*data)
        plt.plot(t, d, label=f'Base {b}')
plt.legend()
plt.xlabel('t = log N')
plt.ylabel('D_b (TV Distance)')
plt.title('Base-Collapse: Residue Deviation vs. Multiplicative Time')
plt.savefig('plots/base_collapse.png')
plt.close()

# Experiment C: Residue-class horizon
moduli = [3, 5, 7, 11, 10, 12, 30, 101]
epsilon = 0.1
for q in moduli:
    coprimes = [r for r in range(q) if gcd(r, q) == 1]
    phi_q = len(coprimes)
    uniform = 1.0 / phi_q
    prime_factors_q = {p for p in primes if p <= q and q % p == 0}
    x_grid = np.logspace(1, np.log10(max_n), 50)
    d_values = []
    first_below = None
    for X in x_grid:
        primes_upto_X = primes[primes <= X]
        relevant_primes = [p for p in primes_upto_X if p not in prime_factors_q]
        if len(relevant_primes) < 10:
            continue
        residues = [p % q for p in relevant_primes]
        count_arr = np.bincount(residues, minlength=q)
        total = len(relevant_primes)
        tv = 0.5 * sum(abs(count_arr[r] / total - uniform) for r in coprimes)
        t = np.log(X)
        d_values.append((t, tv))
        if tv < epsilon and first_below is None:
            first_below = t
    print(f'Modulus {q}: first t below {epsilon}: {first_below if first_below else "None"}, 2 log q: {2 * np.log(q)}')

    # Plot per modulus
    if d_values:
        t, d = zip(*d_values)
        plt.figure(figsize=(8, 5))
        plt.plot(t, d)
        plt.axvline(2 * np.log(q), color='r', linestyle='--', label='Horizon t ≈ 2 log q')
        plt.xlabel('t = log X')
        plt.ylabel('D_q (TV Distance)')
        plt.title(f'Residue Equidistribution for Modulus {q}')
        plt.legend()
        plt.savefig(f'plots/horizon_mod_{q}.png')
        plt.close()

print('All experiments complete. Plots saved to plots/ folder.')
