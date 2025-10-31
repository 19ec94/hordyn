import math

c_a = 0.25
c_gamma = 0.15
sigma = math.sqrt(0.002)

factor = 1/(2 * math.pi * pow(sigma,2))


a = 0.25
gamma = 0.125 * 6 + 0.125/2 
a_contribution = pow(a-c_a,2) / pow(sigma,2)
gamma_contribution = pow(gamma-c_gamma,2) / pow(sigma,2)

# print(factor * math.exp(-0.5 * (a_contribution + gamma_contribution)))

# print(0.00341622 * 2 + 2)

u = 0.00341622
u_bar = 0.02
tau_h = 0.7
gamma = 0.125 * 5 + 0.125/2
c_1 = 0.68
c_2 = 0.08
print(tau_h * (-pow(gamma,2) + (c_1 * gamma + c_2) * (1- math.exp(-u/u_bar)) ))


gamma_s = 0.75
Lambda_bar = 0.1
gamma_bar = 0.01
U = 0.999748
print(Lambda_bar * math.exp(-pow(gamma - gamma_s,2)/gamma_bar) * (1- U))
