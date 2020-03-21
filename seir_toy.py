import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy import integrate
plt.style.use('fivethirtyeight')


def dSdt(t, S, I, beta, N, L, mu):
    return L - mu * S - beta * (I/N) * S


def dEdt(t, S, E, I, N, gamma, mu, a):
    return beta * (I/N) * S - (mu + a) * E


def dIdt(t, E, I, N, gamma, mu, a):
    return a * E - (gamma + mu) * I


def dRdt(t, I, R, gamma, mu):
    return gamma * I - mu * R


def dydt(y, t, N, beta, gamma, mu, a, L):

    S, E, I, R = y
    rhs = [dSdt(t, S, I, beta, N, L, mu),
           dEdt(t, S, E, I, N, gamma, mu, a),
           dIdt(t, E, I, N, gamma, mu, a),
           dRdt(t, I, R, gamma, mu)]

    return rhs


beta = 1  # infection rate
gamma = 0.0001  # recovery rate
N = 1000  # total populus
L = 1  # birth rate
mu = 0.001  # death rate
a = 0.1  # inverse incubation period


I0 = 100
R0 = 0
E0 = 0

y0 = [N-I0-E0, E0, I0, R0]
t = np.linspace(0, 2000, 101)

sol = integrate.odeint(dydt, y0, t, args=(N, beta, gamma, mu, a, L))

plt.style.use('fivethirtyeight')
plt.figure(figsize=(12, 4))
plt.ylim([0, N])
plt.plot(t, sol[:, 0], label='Susceptible')
plt.plot(t, sol[:, 1], label='Infected')
plt.plot(t, sol[:, 2], label='Exposed')
plt.plot(t, sol[:, 3], label='Recovered')
plt.legend()
plt.show()
