# DEF CONSTANTS:

# t0 <- 12/19/2022

# t1 <- 12/27/2022

# t2 <- 08/19/2024

# t3 <- 08/21/2024

# problem 1

from math import log, sqrt, pi, exp

from scipy.stats import norm

from datetime import datetime, date



# d1 = lambda S, K, r, div, sigma, t0, t2: math.log(S_t0/K) + (r - div + 0.5 * sigma ** 2) * ((t2 - t0).days / 365) / (sigma * math.sqrt(((t2 - t0).days / 365)))

def d1(S,K,T,r,sigma, d):
    return(log(S/K)+(r - d + sigma**2/2)*T)/(sigma*sqrt(T))

def d2(S,K,T,r,sigma, d):
    return d1(S,K,T,r,sigma, d)-sigma*sqrt(T)

def call_value(S_t0, K, r_t1_t3, r_t0_t2, div, sigma, t0, t2, t1, t3):

    # div is div(t0, t2), by default

    T_02 = (t2 - t0).days / 365

    T_13 = (t3 - t1).days / 365

    d_1 = d1(S = S_t0,K = K,T = T_02,r = r_t0_t2,sigma = sigma, d = div)

    d_2 = d_1 - sigma * sqrt(T_02)

    risk_free_discount = exp(-r_t1_t3 * T_13)

    expected_payoff = S_t0 * exp((r_t0_t2 - div) * T_02) * norm.cdf(d_1) - K * norm.cdf(d_2)

    return risk_free_discount * expected_payoff

def put_value(S_t0, K, r_t1_t3, r_t0_t2, div, sigma, t0, t2, t1, t3):
    # div is div(t0, t2), by default
    T_02 = (t2 - t0).days / 365

    T_13 = (t3 - t1).days / 365

    d_1 = d1(S = S_t0,K = K,T = T_02,r = r_t0_t2,sigma = sigma, d = div)

    d_2 = d_1 - sigma * sqrt(T_02)

    risk_free_discount = exp(-r_t1_t3 * T_13)

    expected_payoff = -S_t0 * exp((r_t0_t2 - div) * T_02) * norm.cdf(-d_1) + K * norm.cdf(-d_2)

    return risk_free_discount * expected_payoff

def bond_price(b, discount):
    return b * discount


S_t0 = 3817.66 # stock price

K = -1 # strike price

div = -1 

div_t0_t2 = (1.867 / 100)

# div_t1_t3 = -1/100

sigma = -1

discount_t1_t3 = 0.925283

discount_t0_t2 = 0.927082

t0 = date(2022, 12, 19)

t1 = date(2022, 12, 27)

t2 = date(2024, 8, 19)

t3 = date(2024, 8, 21)

T_02 = (t2 - t0).days / 365

T_13 = (t3 - t1).days / 365

r_t1_t3 = log(discount_t1_t3) / (-T_13)

r_t0_t2 = log(discount_t0_t2) / (-T_02)

# print(r_t1_t3, r_t0_t2)

sigma_short_put = 28.445 / 100 # K = 0.8 S

sigma_long_call = 24.453 / 100 # 

sigma_short_call = 21.424 / 100 # k = 1.127 S

long_call = -call_value(S_t0 = S_t0, K = 1 * S_t0, r_t1_t3 = r_t1_t3, r_t0_t2 = r_t0_t2, div = div_t0_t2, sigma = sigma_long_call, t0 = t0, t2 = t2, t1 = t1, t3 = t3)

short_call = call_value(S_t0 = S_t0, K = 1.127 * S_t0, r_t1_t3 = r_t1_t3, r_t0_t2 = r_t0_t2, div = div_t0_t2, sigma = sigma_short_call, t0 = t0, t2 = t2, t1 = t1, t3 = t3)

short_put = put_value(S_t0 = S_t0, K = 0.8 * S_t0, r_t1_t3 = r_t1_t3, r_t0_t2 = r_t0_t2, div = div_t0_t2, sigma = sigma_short_put, t0 = t0, t2 = t2, t1 = t1, t3 = t3)
# one zero coupon bond, N_1 short put, N_2 long call, N_3 short call
N_1, N_2, N_3 = 0.032743 * 10, 0.052388 * 10, 0.052388 * 10

bond = -bond_price(b = 1000, discount = discount_t1_t3)
# discount = exp(-rT) so r = log(discount) / (-T)

print(f"The value of this portfolio is: ${N_1 * short_put + N_2 * long_call + N_3 * short_call + bond}")

# Problem 2: Explain

# Problem 3:

S_t2 = 4008.543

ratio_0_2 = S_t2/S_t0 # 1.05

Pay = (10 / 100) / T_02

if Pay > r_t0_t2:
    print(Pay)

else:
    print("Bad")








