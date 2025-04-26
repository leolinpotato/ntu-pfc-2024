import argparse
from scipy.optimize import fsolve, brentq

def calculate_loan_details(cash, n, m, r1, r2):
    """
    Calculate the loan details and bond investment returns.
    
    :param cash: The initial cash available for investment.
    :param n: The number of years for the loan.
    :param m: The number of years for the bond investment.
    :param r1: The annual interest rate for the loan (as a decimal).
    :param r2: The annual interest rate for the bond (as a decimal).
    :return: A tuple containing the maximum loan amount, total interest paid on the loan,
             total interest received from the bond, and the annualized internal rate of return.
    """
    v1 = int(cash*r2/(r1/(1-(1+r1/m)**(-n*m))-r2))
    c = v1*r1/m/(1-(1+r1/m)**(-n*m))
    total_interest_paid = round(c*n*m-v1, 6)
    total_interest_received = round((cash+v1)*(r2/m)*n*m, 6)
    
    def irr_equation(irr):
        return sum(((cash+v1)*(r2/m)-c)/(1+irr)**(t+1) for t in range(n*m))+(cash+v1)/(1+irr)**(n*m)-cash
    

    '''
    # Use fsolve to find the IRR
    irr_guess = 0.04  # Initial guess for IRR
    irr_period = fsolve(irr_equation, irr_guess)[0]
    irr = round(irr_period*m, 6)
    '''

    # Use Brent's method, searching for IRR in a reasonable range (0% to 100%)
    try:
        irr_period = brentq(irr_equation, 0, 1)
        irr = round(irr_period * m, 6)  # Convert to annualized IRR
    except ValueError:
        irr = float('nan')  # Handle failure gracefully

    return v1, total_interest_paid, total_interest_received, irr 

def main():
    parser = argparse.ArgumentParser(description="Analyze bond investment financing through a loan.")
    parser.add_argument("cash", type=float, help="The initial cash available for investment.")
    parser.add_argument("n", type=int, help="The number of years for the loan.")
    parser.add_argument("m", type=int, help="The number of years for the bond investment.")
    parser.add_argument("r1", type=float, help="The annual interest rate for the loan (as a decimal).")
    parser.add_argument("r2", type=float, help="The annual interest rate for the bond (as a decimal).")
    
    args = parser.parse_args()
    
    v1, total_interest_paid, total_interest_received, irr = calculate_loan_details(
        args.cash, args.n, args.m, args.r1, args.r2)
    
    print(f'{v1}, {total_interest_paid:.6f}, {total_interest_received:.6f}, {irr:.6f}')

if __name__ == "__main__":
    main()
