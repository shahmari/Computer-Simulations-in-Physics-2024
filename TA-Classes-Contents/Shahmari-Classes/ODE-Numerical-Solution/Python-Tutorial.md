### Tutorial on Numerical Methods for Solving Differential Equations in Python

Numerical methods are essential for solving differential equations that cannot be solved analytically. This tutorial will introduce you to some of the most popular numerical methods for solving differential equations and demonstrate how to implement them in Python using basic programming constructs and the `scipy` library.

#### Prerequisites
- Basic understanding of differential equations.
- Familiarity with Python programming.
- Installed `numpy` and `scipy` libraries (`pip install numpy scipy`).

### Popular Numerical Methods

1. **Euler's Method**
2. **Improved Euler's Method (Heun's Method)**
3. **Runge-Kutta Methods**
4. **Using `scipy.integrate.odeint`**

#### 1. Euler's Method

Euler's method is the simplest numerical method for solving ordinary differential equations (ODEs).

**Algorithm:**
1. Start at the initial condition \((t_0, y_0)\).
2. Use the formula \( y_{n+1} = y_n + h \cdot f(t_n, y_n) \) to find \( y_{n+1} \).
3. Increment \( t_n \) by the step size \( h \): \( t_{n+1} = t_n + h \).
4. Repeat steps 2 and 3 for the desired number of steps or until the end of the interval.

**Python Implementation:**

```python
import numpy as np
import matplotlib.pyplot as plt

def euler_method(f, t0, y0, h, n_steps):
    t = t0
    y = y0
    t_values = [t]
    y_values = [y]
    
    for _ in range(n_steps):
        y = y + h * f(t, y)
        t = t + h
        t_values.append(t)
        y_values.append(y)
    
    return np.array(t_values), np.array(y_values)

# Example differential equation dy/dt = -2y
def f(t, y):
    return -2 * y

# Initial conditions
t0 = 0
y0 = 1
h = 0.1
n_steps = 50

t_values, y_values = euler_method(f, t0, y0, h, n_steps)

plt.plot(t_values, y_values, label="Euler's Method")
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()
```

#### 2. Improved Euler's Method (Heun's Method)

This method improves the accuracy of Euler's method by averaging the slopes at the beginning and the end of the interval.

**Python Implementation:**

```python
def improved_euler_method(f, t0, y0, h, n_steps):
    t = t0
    y = y0
    t_values = [t]
    y_values = [y]
    
    for _ in range(n_steps):
        k1 = f(t, y)
        y_predict = y + h * k1
        k2 = f(t + h, y_predict)
        y = y + (h / 2) * (k1 + k2)
        t = t + h
        t_values.append(t)
        y_values.append(y)
    
    return np.array(t_values), np.array(y_values)

t_values, y_values = improved_euler_method(f, t0, y0, h, n_steps)

plt.plot(t_values, y_values, label="Improved Euler's Method")
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()
```

#### 3. Runge-Kutta Methods

The fourth-order Runge-Kutta method (RK4) is one of the most widely used methods due to its balance between accuracy and computational efficiency.

**Python Implementation:**

```python
def runge_kutta_4(f, t0, y0, h, n_steps):
    t = t0
    y = y0
    t_values = [t]
    y_values = [y]
    
    for _ in range(n_steps):
        k1 = h * f(t, y)
        k2 = h * f(t + h / 2, y + k1 / 2)
        k3 = h * f(t + h / 2, y + k2 / 2)
        k4 = h * f(t + h, y + k3)
        y = y + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        t = t + h
        t_values.append(t)
        y_values.append(y)
    
    return np.array(t_values), np.array(y_values)

t_values, y_values = runge_kutta_4(f, t0, y0, h, n_steps)

plt.plot(t_values, y_values, label="Runge-Kutta Method")
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()
```

#### 4. Using `scipy.integrate.odeint`

The `scipy.integrate.odeint` function provides a robust and efficient way to solve ODEs.

**Python Implementation:**

```python
from scipy.integrate import odeint

# Example differential equation dy/dt = -2y
def f(y, t):
    return -2 * y

# Initial conditions
t0 = 0
y0 = 1
t = np.linspace(t0, 5, 100)  # t values from 0 to 5

# Solve ODE
y_values = odeint(f, y0, t)

plt.plot(t, y_values, label="scipy.integrate.odeint")
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()
```