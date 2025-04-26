import argparse
import math
import ipdb

def binomial_tree_barrier_options(S, X, r, s, T, H, n):
    # Convert T to years
    T = T / 365.0

    # Calculate parameters
    dt = T / n
    u = math.exp(s * math.sqrt(dt))
    d = 1 / u
    R = math.exp(r * dt)
    p = (R - d) / (u - d)

    # Set h so that H = S * (u ** h)
    h = math.ceil(math.log(H / S) / math.log(u))

    up_out_call = [0] * (n + 1)
    up_out_put = [0] * (n + 1)
    call = [0] * (n + 1)
    put = [0] * (n + 1)

    # Initialize stock prices at maturity
    for i in range(n + 1):
        up_out_call[i] = max(0, S * (u ** (n - i)) * (d ** i) - X)
        up_out_put[i] = max(0, X - S * (u ** (n - i)) * (d ** i))
        call[i] = max(0, S * (u ** (n - i)) * (d ** i) - X)
        put[i] = max(0, X - S * (u ** (n - i)) * (d ** i))

    if (n - h) % 2 == 0 and 0 <= (n - h) / 2 <= n:
        up_out_call[int((n - h) / 2)] = 0
        up_out_put[int((n - h) / 2)] = 0

    for j in range(n - 1, -1, -1):
        for i in range(j + 1):
            up_out_call[i] = (p * up_out_call[i] + (1 - p) * up_out_call[i + 1]) / R
            up_out_put[i] = (p * up_out_put[i] + (1 - p) * up_out_put[i + 1]) / R
            call[i] = (p * call[i] + (1 - p) * call[i + 1]) / R
            put[i] = (p * put[i] + (1 - p) * put[i + 1]) / R

        if (j - h) % 2 == 0 and 0 <= (j - h) / 2 <= j:
            up_out_call[int((j - h) / 2)] = 0
            up_out_put[int((j - h) / 2)] = 0

    # Return option prices
    return up_out_call[0], up_out_put[0], call[0] - up_out_call[0], put[0] - up_out_put[0]

def main():
    parser = argparse.ArgumentParser(description="Calculate the price of a Bermudan option using a binomial tree model.")
    parser.add_argument("S", nargs="?", type=float, default=100, help="Stock price")
    parser.add_argument("X", nargs="?", type=float, default=110, help="Strike price")
    parser.add_argument("r", nargs="?", type=float, default=0.03, help="Continuously compounded annual interest rate")
    parser.add_argument("s", nargs="?", type=float, default=0.3, help="Volatility of the stock")
    parser.add_argument("T", nargs="?", type=int, default=60, help="Time to maturity")
    parser.add_argument("H", nargs="?", type=float, default=120, help="up-and-out barrier")
    parser.add_argument("n", nargs="?", type=int, default=100, help="Number of time steps")
    
    args = parser.parse_args() 

    up_out_call_price, up_out_put_price, up_in_call_price, up_in_put_price = binomial_tree_barrier_options(args.S, args.X, args.r, args.s, args.T, args.H, args.n)
    print(f"{up_out_call_price:.6f}, {up_out_put_price:.6f}, {up_in_call_price:.6f}, {up_in_put_price:.6f}")

if __name__ == "__main__":
    main()