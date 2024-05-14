### Tutorial on Numerical Methods for Solving Differential Equations in Julia

Julia is a high-performance programming language that's particularly well-suited for numerical and scientific computing. In this tutorial, we will explore some of the most popular numerical methods for solving differential equations and demonstrate how to implement them in Julia. We will also use the `DifferentialEquations.jl` package, which provides a powerful suite of tools for solving differential equations.

#### Prerequisites
- Basic understanding of differential equations.
- Familiarity with Julia programming.
- Installed `DifferentialEquations.jl` package (`using Pkg; Pkg.add("DifferentialEquations")`).

### Popular Numerical Methods

1. **Euler's Method**
2. **Improved Euler's Method (Heun's Method)**
3. **Runge-Kutta Methods**
4. **Using `DifferentialEquations.jl`**

#### 1. Euler's Method

Euler's method is the simplest numerical method for solving ordinary differential equations (ODEs).

**Algorithm:**
1. Start at the initial condition \((t_0, y_0)\).
2. Use the formula \( y_{n+1} = y_n + h \cdot f(t_n, y_n) \) to find \( y_{n+1} \).
3. Increment \( t_n \) by the step size \( h \): \( t_{n+1} = t_n + h \).
4. Repeat steps 2 and 3 for the desired number of steps or until the end of the interval.

**Julia Implementation:**

```julia
using Plots

function euler_method(f, t0, y0, h, n_steps)
    t_values = [t0]
    y_values = [y0]
    t = t0
    y = y0
    
    for _ in 1:n_steps
        y = y + h * f(t, y)
        t = t + h
        push!(t_values, t)
        push!(y_values, y)
    end
    
    return t_values, y_values
end

# Example differential equation dy/dt = -2y
f(t, y) = -2 * y

# Initial conditions
t0 = 0.0
y0 = 1.0
h = 0.1
n_steps = 50

t_values, y_values = euler_method(f, t0, y0, h, n_steps)

plot(t_values, y_values, label="Euler's Method", xlabel="t", ylabel="y")
```

#### 2. Improved Euler's Method (Heun's Method)

This method improves the accuracy of Euler's method by averaging the slopes at the beginning and the end of the interval.

**Julia Implementation:**

```julia
function improved_euler_method(f, t0, y0, h, n_steps)
    t_values = [t0]
    y_values = [y0]
    t = t0
    y = y0
    
    for _ in 1:n_steps
        k1 = f(t, y)
        y_predict = y + h * k1
        k2 = f(t + h, y_predict)
        y = y + (h / 2) * (k1 + k2)
        t = t + h
        push!(t_values, t)
        push!(y_values, y)
    end
    
    return t_values, y_values
end

t_values, y_values = improved_euler_method(f, t0, y0, h, n_steps)

plot!(t_values, y_values, label="Improved Euler's Method")
```

#### 3. Runge-Kutta Methods

The fourth-order Runge-Kutta method (RK4) is one of the most widely used methods due to its balance between accuracy and computational efficiency.

**Julia Implementation:**

```julia
function runge_kutta_4(f, t0, y0, h, n_steps)
    t_values = [t0]
    y_values = [y0]
    t = t0
    y = y0
    
    for _ in 1:n_steps
        k1 = h * f(t, y)
        k2 = h * f(t + h / 2, y + k1 / 2)
        k3 = h * f(t + h / 2, y + k2 / 2)
        k4 = h * f(t + h, y + k3)
        y = y + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        t = t + h
        push!(t_values, t)
        push!(y_values, y)
    end
    
    return t_values, y_values
end

t_values, y_values = runge_kutta_4(f, t0, y0, h, n_steps)

plot!(t_values, y_values, label="Runge-Kutta Method")
```

#### 4. Using `DifferentialEquations.jl`

The `DifferentialEquations.jl` package provides a robust and efficient way to solve ODEs.

**Julia Implementation:**

```julia
using DifferentialEquations

# Example differential equation dy/dt = -2y
function f!(dy, y, p, t)
    dy .= -2 .* y
end

# Initial conditions
t0 = 0.0
y0 = 1.0
tspan = (t0, 5.0)
u0 = [y0]

# Solve ODE
prob = ODEProblem(f!, u0, tspan)
sol = solve(prob, Tsit5(), saveat=0.1)

plot(sol, label="DifferentialEquations.jl", xlabel="t", ylabel="y")
```