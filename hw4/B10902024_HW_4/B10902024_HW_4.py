import math
import argparse
import ipdb
import numpy as np

def binomial_tree(S, X, r, s, dt, H, L, n, option_type):
    """
    Calculate the price of a double barrier option using a binomial tree model.
    """
    # Calculate parameters
    u = math.exp(s * math.sqrt(dt))
    d = 1 / u
    R = math.exp(r * dt)
    p = (R - d) / (u - d)

    # Calculate the number of up and down steps needed to reach H and L
    h = round(math.log(H / S) / math.log(u))  # Steps to reach H
    l = round(math.log(L / S) / math.log(u))  # Steps to reach L

    # Initialize option values at maturity
    option_values = [0] * (n + 1)

    for i in range(n + 1):
        if option_type == "call":
            option_values[i] = max(0, S * (u ** (n - i)) * (d ** i) - X)
        elif option_type == "put":
            option_values[i] = max(0, X - S * (u ** (n - i)) * (d ** i))

    if (n - h) % 2 == 0 and 0 <= (n - h) / 2 <= n:
        option_values[int((n - h) / 2)] = 0

    if (n - l) % 2 == 0 and 0 <= (n - l) / 2 <= n:
        option_values[int((n - l) / 2)] = 0

    for j in range(n - 1, -1, -1):
        for i in range(j + 1):
            option_values[i] = (p * option_values[i] + (1 - p) * option_values[i + 1]) / R

        if (j - h) % 2 == 0 and 0 <= (j - h) / 2 <= j:
            option_values[int((j - h) / 2)] = 0

        if (j - l) % 2 == 0 and 0 <= (j - l) / 2 <= j:
            option_values[int((j - l) / 2)] = 0

    return option_values[0]

def double_barrier_option(S, X, r, s, T, H, L, k, option_type):
    T = T / 365.0
    dt = (math.log(H / L) / (2 * s * k))**2
    periods = math.floor(T / dt) - 1
    dt_ = T - periods * dt
    u = math.exp(s * math.sqrt(dt))
    d = 1 / u
    p = (math.exp(r * dt) - d) / (u - d)

    # Let node B be the node whose logarithmic return is closest to mu amonbg all nodes at time t+dt_
    mu = (r - 0.5 * s**2) * dt_
    var = s**2 * dt_
    delta_V = 2 * s * math.sqrt(dt)
    mu_hat = round((mu - math.log(L / S)) / delta_V) * delta_V + math.log(L / S)
    logS = np.log(S)
    log_prices = np.array([
        logS + mu_hat + delta_V,
        logS + mu_hat,
        logS + mu_hat - delta_V
    ])
    S_nodes = np.exp(log_prices)

    # Solve for three branching probabilities
    beta = mu_hat - mu
    alpha = beta + delta_V
    gamma = beta - delta_V
    A = np.array([
        [alpha, beta, gamma],
        [alpha**2, beta**2, gamma**2],
        [1, 1, 1]
    ])
    b = np.array([0, var, 1])

    probs = np.linalg.solve(A, b)

    # Calculate option prices at the final nodes
    values = []
    for S0 in S_nodes:
        val = binomial_tree(S0, X, r, s, dt, H, L, periods, option_type)
        values.append(val)    

    return np.exp(-r * dt_) * np.dot(probs, values)


def main():
    parser = argparse.ArgumentParser(description="Calculate double barrier option prices and deltas.")
    parser.add_argument("S", nargs="?", type=float, default=95, help="Current stock price")
    parser.add_argument("X", nargs="?", type=float, default=100, help="Strike price")
    parser.add_argument("r", nargs="?", type=float, default=0.10, help="Risk-free interest rate")
    parser.add_argument("s", nargs="?", type=float, default=0.25, help="Volatility")
    parser.add_argument("T", nargs="?", type=int, default=365, help="Time to maturity in days")
    parser.add_argument("H", nargs="?", type=float, default=140, help="Upper barrier")
    parser.add_argument("L", nargs="?", type=float, default=90, help="Lower barrier")
    parser.add_argument("k", nargs="?", type=int, default=50, help="Number of time steps")

    args = parser.parse_args()

    S = args.S
    X = args.X
    r = args.r
    s = args.s
    T = args.T
    H = args.H
    L = args.L
    k = args.k

    call = double_barrier_option(S, X, r, s, T, H, L, k, "call")
    call_up = double_barrier_option(S * 1.01, X, r, s, T, H, L, k, "call")
    call_down = double_barrier_option(S * 0.99, X, r, s, T, H, L, k, "call")
    call_delta = (call_up - call_down) / (2 * S * 0.01)
    put = double_barrier_option(S, X, r, s, T, H, L, k, "put")
    put_up = double_barrier_option(S * 1.01, X, r, s, T, H, L, k, "put")
    put_down = double_barrier_option(S * 0.99, X, r, s, T, H, L, k, "put")
    put_delta = (put_up - put_down) / (2 * S * 0.01)

    print(f"{call:.6f}, {call_delta:.6f}, {put:.6f}, {put_delta:.6f}")


if __name__ == "__main__":
    main()