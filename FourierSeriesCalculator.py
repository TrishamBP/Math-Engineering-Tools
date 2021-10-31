#Fourier Series Calculator
import sympy as sym

import numpy as np

import matplotlib.pyplot as plt

sym.init_session()
sym.init_printing()


xi, n, L ,m = sym.symbols('xi , n, L , m')

sym.printing.pprint(xi)
sym.printing.pprint(n)

def f(x):
    return 7*x**2+19*x

sym.Q.integer(n)

a_0=sym.integrate((2/L)*f(xi),(xi,0,L))

sym.printing.pprint(a_0)

a_n=sym.integrate((2/L)*f(xi)*sym.cos(2*sym.pi*n*xi/L),(xi,0,L))

sym.printing.pprint(a_n)


b_n=sym.integrate((2/L)*f(xi)*sym.sin(2*sym.pi*n*xi/L),(xi,0,L))

sym.printing.pprint(b_n)


a0=sym.lambdify(L,a_0,modules='numpy')
an=sym.lambdify((n,L),a_n,modules='numpy')
bn=sym.lambdify((n,L),b_n,modules='numpy')

f_l=sym.lambdify((xi,L),f(xi),modules='numpy')

N=100;


P=11

n_list=np.arange(N)+1

a_list=np.array([an(i,P) for i in n_list])
b_list=np.array([bn(i,P) for i in n_list])


xx=np.linspace(-2*P,2*P,10000)

def g(x):
    uu=np.array([a_list[i]*np.cos(2*np.pi*n_list[i]*xx/P)+b_list[i]*np.sin(2*np.pi*n_list[i]*xx/P) for i in np.arange(N)])
    return a0(P)/2+np.sum(uu,axis=0)

y_FS=g(xx)


fig1 = plt.figure('Fourier Series')

plt.clf()

ax1 = fig1.gca()


# ax1.set_xticks(np.linspace(-2*P, 2*P, 10))
# ax1.set_yticks(np.linspace(-2*P, 2*P, 10))

plt.grid()


ax1.plot(xx,f_l(xx,P))

ax1.plot(xx,g(xx))





ax1.set_xlabel('x')
# ax1.set_xlim(-P, P)
ax1.set_ylabel('y')
# ax1.set_ylim(-P, P)






plt.show()
