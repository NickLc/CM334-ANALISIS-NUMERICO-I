import numpy as np
import matplotlib.pyplot as plt


coef = np.array([1, -8, 28, -56, 70, -56, 28, -8, 1])

def f(x) :
    r = 0.0
    for i in np.arange(9) :
        r += coef[i]*pow(x, i)
    return (r)

def g(x) :
    r = 1.0
    for i in np.arange(1, 9) :
        r *= x
        r += coef[i]
    return (r)

def h(x) :
    return (pow(x-1.0, 8))

points = 101
startPoint = 0.99
endPoint = 1.01

t = np.linspace(startPoint, endPoint, points)
a = np.linspace(startPoint, endPoint, points)
b = np.linspace(startPoint, endPoint, points)
c = np.linspace(startPoint, endPoint, points)

for i in np.arange(points) :
    a[i] = f(a[i])
    b[i] = g(b[i])
    c[i] = h(c[i])

plt.plot(t, a, 'r', label="f(x)")
plt.plot(t, b, 'g', label="g(x)")
plt.plot(t, c, 'b', label="h(x)")

plt.show()