# Problem Set 4
## 1-
Solve the following differential equation numerically (by discretizing the time axis) exactly or approximately (repetition method). The purpose of this question is your attempt to apply creative methods (according to your previous learning or new creativity).
$$md^2x/dt^2 + kx = ax^3/6$$

## 2-
Calculate the energy and wave function of the ground and excited states of a particle in a two-dimensional potential well (a square with dimensions L).

## 3-
In the following example, I have used the Metropolis-Hasting algorithm to generate x-values based on the normal distribution function and accordingly compute $\langle x \rangle$ and $\langle x^2 \rangle$. As you see, the histogram of generated x values is quite close to the normal distribution function.

Repeat the same procedure for the following integrals

$$\langle f(x) \rangle = \int f(x)w(x) dx\int w(x) dx = \int f(x) p(x)$$
where $p(x) = w(x)/\int w(x) dx$,  via the Metropolis-Hasting algorithm susing the $r_{12} = W(x_2)/W(x_1)$ as the transition ration and for:

(a) $f(x) = x,  w(x) = exp(-|x|)$

(b) $f(x) = x^2, w(x) = 1/(1+x^2)$

(c) $f(x) = x^2, w(x) = cos(x)$

(d) $f(x) = x^2, w(x) = cos(x)exp(-x^2/2)$

(e) $f(x) = cos(10x) , w(x) = exp(-x^2/2)$


Use direct integration methods (analytically or using the Riemann integral evaluation numerically) and answer these questions besides the Monte-Carlo method.

Do we have sign problem in the above problems? In which? How does these expectation values change by increasing the number of samplings? Using this information, study the convergence rate of the the above four cases (qualitatively) and determine which one converges slowly and which one converges rapidly.

```matlab
%=========================================================================================================================================%
%%%%% In the following sample code we examine the Metropolis-Hasting %%% Monte-Carlo algorithm to generate n_sampling x values according to
%%%%% the normal distribution with mean 0 and variance 1 and compute several
%%%%% expectation values accordingly.

clear;

tic;
% defining an inline function (a.k.a. anonymous function) for the Gaussian distribution function
p=@(x) exp(-x.^2/2)/sqrt(2*pi);

dx = 10^(-1);
n_sampling = 10^7;
n_warmup = round(n_sampling/5); % we will ignore n_warmup first generated x values

x = zeros(1,n_sampling);


x(1) = randn(1);
for i=2:n_sampling
    x0 = x(i-1);
    del_x = dx*randn(1);
    x1 = x0 + del_x;
    r10 = p(x1)/p(x0);
    if rand(1) < r10
       x(i) = x0+del_x;
    else
       x(i) = x0;
    end
end

x = x(n_warmup+1:end);
n_bins = 100;
[h,c] = hist(x,n_bins);
dc = c(2) - c(1);

% computing <x>
x_avg = mean(x)

% computing <x^2>
x2_avg = mean(x.^2)

% computing <cos(x)>
cosx_avg = mean(cos(x))

% computing <cosh(x)>
coshx_avg = mean(cosh(x))

%%%% to plot the Gaussian normal distribution
x0 = -5:dc:5;
y0 = p(x0);


r_pdf = sum(h)/sum(y0);
h = h/r_pdf; % maxing sure h(x) describes a probability distribution function (pdf)

figure(1);
hold on;
plot(x0,y0,'.');
plot(c,h);

toc
```

## 4-

### a.
Write a code that uses the Lanczos method and finds the lowest eigenvalue and eigenvector of $M = A + A^\dag$ for a random nxn matrix A. Test whether or not your code generates the correct result for $n=2000$ by comparing it with matlab's (or python's or Julia's) diagoanlization method. Please report and compare their CPU times as well

Now consider $M = 2A_1 \otimes B_1 + 3A_2 \otimes B_2 + h.c$. where $\otimes$ denotes the Kronecker tesnsor product (in LATEX) and $A_1$ and $A_2$ are $n\times n$ matrices ($n=50$) and $B_1$ and $B_2$ are $m\times m$ ($m=100$).  

### b.
Construct M explicitly and use part a to obtains its lowest eigenvalues and eigenvectors.

### c.
Write a code based on the Lanczos method to obtain the lowest eigenvalue and eigenvectors without explicitly constructing matrix M (as explained in the class). Compare the results of part b and c (including their CPU times)