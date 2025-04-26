import argparse
import math
import ipdb

def bermudan_option_pricing(S, X, r, s, T, m, E):
    """
    Calculate the price of a Bermudan option using a binomial tree model.

    :param S: Initial stock price
    :param X: Strike price
    :param r: Continuously compounded annual interest rate
    :param s: Volatility of the stock
    :param T: Time to maturity (in days)
    :param m: Number of periods per day for the tree
    :param E: Exercise dates from now, a list of integers
    :return: Tuple containing the put and call option prices
    """
    
    # Calculate parameters for the binomial tree
    u = math.exp(s * math.sqrt(1 / (m * 365))) # Up factor
    d = math.exp(-s * math.sqrt(1 / (m * 365)))  # Down factor
    R = math.exp(r / (m * 365)) 
    p = (R - d) / (u - d)  # Risk-neutral probability
    n = T * m  # Total number of periods

    # Initialize stock prices at maturity
    stock_prices = [S * (u ** (n - i)) * (d ** i) for i in range(n + 1)]
    
    # Initialize option values at maturity
    call_values = [max(0, price - X) for price in stock_prices]
    put_values = [max(0, X - price) for price in stock_prices]

    # Backward induction to calculate option values at earlier nodes
    for j in range(n - 1, -1, -1):
        for i in range(j + 1):
            call_values[i] = (p * call_values[i] + (1 - p) * call_values[i + 1]) / R
            put_values[i] = (p * put_values[i] + (1 - p) * put_values[i + 1]) / R

            day = j // m
            if day in E:
                call_values[i] = max(call_values[i], S * (u ** (j - i)) * (d ** i) - X)
                put_values[i] = max(put_values[i], X - S * (u ** (j - i)) * (d ** i))

    return put_values[0], call_values[0]

def main():
    parser = argparse.ArgumentParser(description="Calculate the price of a Bermudan option using a binomial tree model.")
    parser.add_argument("S", nargs="?", type=float, default=100, help="Initial stock price")
    parser.add_argument("X", nargs="?", type=float, default=110, help="Strike price")
    parser.add_argument("r", nargs="?", type=float, default=0.03, help="Continuously compounded annual interest rate")
    parser.add_argument("s", nargs="?", type=float, default=0.3, help="Volatility of the stock")
    parser.add_argument("T", nargs="?", type=int, default=60, help="Time to maturity")
    parser.add_argument("m", nargs="?", type=int, default=5, help="Number of periods per day for the tree")
    parser.add_argument("E", nargs="*", type=int, default=[10, 20, 30, 40, 50], help="Exercise dates from now, a list of intergers")
    
    args = parser.parse_args() 

    put, call = bermudan_option_pricing(args.S, args.X, args.r, args.s, args.T, args.m, args.E)
    print(f"{put:.6f}, {call:.6f}")

if __name__ == "__main__":
    main()
