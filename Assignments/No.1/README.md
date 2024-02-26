# Assignment 1:

Convert the following code from MATLAB to Python:

Code 1:

```matlab
clc;
close all;
clear;


dt = 1/10; % tme step
t = 0:dt:15; % t vector (a row vector)
t = reshape(t,[],1); % transforming t into a column vector. alternatively we could use t = reshape(t,[],1) or t = t'; 

N = length(t);
D1 = (-diag(ones(1,N),-1)+diag(ones(1,N+1),0))/dt;
D1(1,1) = 1; D1(1,2) = 0;



% solving df/dt = g=cos(t) equation assuming f(0)= 1 
f0 = 0;
g = cos(1*t(1:end));

g_tilde = [f0; g];

f = inv(D1)*g_tilde; % alternatively we could use f = D1\A;


figure(1);
hold on;
plot(t(1:end),g,'-m','LineWidth',4);
plot(t,f(2:end),'-c','LineWidth',3);




% second order derivative
D2 = (-2*diag(ones(1,N+1),0)+diag(ones(1,N),1)+diag(ones(1,N),-1))/dt^2;
```

Code 2:

```MATLAB
%%%%%%%=======================================================%%%%%%%
%%%%%%%        Solving Ordinary Differential Equations        %%%%%%%
%%%%%%%=======================================================%%%%%%%

clc;
close all;
clear;

% discretization of x axis
dx = 0.05;
x = (0:dx:30)';
n_x = length(x);

% representation of the 1st derivative w.r.t. x axis in a finite space
D1_x = (diag(ones(1,n_x-1),1)-diag(ones(1,n_x),0))/dx;
D1_x(n_x,1) = 1/dx;



% representation of the 2nd derivative w.r.t. x axis in a finite space
D2_x = -(D1_x'*D1_x+D1_x*D1_x')/2;


% modifying 1st & 2nd derivatives such that their determinants are nonzero,
% hence both become invertible
D1_x_tilde = D1_x; 
D1_x_tilde(:,1) = 0; D1_x_tilde(1,:) = 0; 
D1_x_tilde(1,1) = 1;

D2_x_tilde = D2_x;
D2_x_tilde(1,:) = 0; D2_x_tilde(:,1) = 0;
D2_x_tilde(1,1) = 1;
D2_x_tilde(2,2) = 1/dx;
D2_x_tilde(2,3) = -1/dx;


% Examples:


% example 1: solving df(x)/dx = -2x equation 
g_x_tilde = -2*x;
g_x_tilde(1) = 0;
g_x_tilde(2) = 0;

% D1_x_tilde * f = g_x_tilde ==> f = inv(D1_x_tilde) * g_x_tilde;
f_x_tilde = inv(D1_x_tilde) * g_x_tilde;
figure(1);
plot(x(2:end),f_x_tilde(2:end), 'Linewidth', 2);
xlabel('x');
ylabel('f(x)');
title('plot 1: solution of df/dx = -2x equation');



% example 2: solving df(x)/dx = -f(x) equation
% trick: we have to modify it as df/dx - f(x) = a * delta(x) --> (d/dx -
% id) f = a * delta(x) where a is arbitrary and can be fixed using f(x=0)
% value. Here, we use a = 1

delta_x = zeros(size(x));
delta_x(1) = 1/dx;
% (D1_x - eye(n_x)/4)f(x) = -delta(x)
f_x_tilde = inv(D1_x - 1/4*eye(n_x)) * (-delta_x);
figure(2);
% plot(x(2:end),log(abs(f_x_tilde(2:end))), 'Linewidth', 2);
plot(x(2:end),f_x_tilde(2:end), 'Linewidth', 2);
xlabel('x');
ylabel('f(x)');
title('plot 2: solution of df/dx = f/4 equation');



% example 3: solving d^2f/dx^2 + f(x) = delta(x) equation
% trick: we modify it as: d^2f/dx^2 + f(x) = delta(x) ==> (d^2/dx^2 + eye(n_x)) f = delta(x) ==> 
% f = inv(d^2/dx^2 + eye(n_x)) * delta(x)
f_x_tilde = inv(D2_x + eye(n_x)) * delta_x;
figure(3);
plot(x(2:n_x),f_x_tilde(2:n_x), 'Linewidth', 2);
xlabel('x');
ylabel('f(x)');
title('plot 3: simple harmonic oscilator');


% example 4: solving d^2f/dx^2 + gamma * df/dx + f(x) = 0 equation (damped harmonic oscillator)
% trick: we modify it as d^2f/dx^2 + gamma * df(x)/dx + f(x) = delta(x) ==> (d^2/dx^2 + gamma * dx/dx + eye(n_x)) f = delta(x) ==> 
% f = inv(d^2/dx^2 + gamma * dx/dx + eye(n_x)) * delta(x)
gamma = 0.4 ;
f_x_tilde = inv(D2_x + gamma * D1_x + eye(n_x)) * delta_x;
figure(4);
plot(x(2:n_x),f_x_tilde(2:n_x), 'Linewidth', 2);
xlabel('x');
ylabel('f(x)');
title('plot 4: damped harmonic oscilator');



% example 5: solving d^2f/dx^2 + gamma * df/dx + f(x) = g(x)) equation (driven damped harmonic oscillator)
% f = inv(d^2/dx^2 + gamma * dx/dx + eye(n_x)) * g(x)
g_x = sin(1*x);

gamma = 0.5 ;
f_x_tilde = inv(D2_x + gamma * D1_x + eye(n_x)) * g_x;
figure(5);
plot(x(2:n_x),f_x_tilde(2:n_x), 'Linewidth', 2);
xlabel('x');
ylabel('f(x)');
title('plot 5: driven harmonic oscilator');
```

Code 3:

```MATLAB
%%%%%%%=======================================================%%%%%%%
%%%%%%%   Eigenvalue problem & Solvig Schrodinger Equation    %%%%%%%
%%%%%%%=======================================================%%%%%%%

clc;
close all;
clear;


dx = 0.01;
x = (-10:dx:10)';
omega = 1; % omega = sqrt(k/m) = frequency of oscillations
m = 1; % mass
hbar = 1; % Planck constant

n_x = length(x);
D1_x = (diag(ones(1,n_x-1),1)-diag(ones(1,n_x),0))/dx;
D1_x(n_x,1) = 1;
P1_x = (hbar/1i) * D1_x; % momentum operator = hbar/i d/dx

V_x_diag = m * omega^2 * x.^2/2; % potential function for the simple harmonics oscilator 

% potential function for a quantum well stretched from x=-3 to x=+4
% V_x_diag = 1e6*ones(1,n_x); % 1e6 = 10^6
% V_x_diag(abs(x)<3) = 0;

% V_x_diag = (1-cos(x/10)); % cosine potential function


% V_x_diag = abs(x); % linear distance potential function



figure(1);
plot(V_x_diag,'LineWidth',4)
xlabel('x');
ylabel('V(x)');
title('Potential functioon');

V_x = diag(V_x_diag);
H = (P1_x'*P1_x + P1_x * P1_x')/2/(2*m) + V_x; % Hamiltonian operator

[psi, Hd] = eig(full(H)); % diagonalizing Hamiltonian operator

E = diag(Hd); % energy eigenvalues
E1 = E(1) % lowest energy (ground-state energy)
psi_1 = psi(:,1); % ground-state energy wavefunction
figure(2);
hold on;
plot(x,abs(psi_1).^2, 'Linewidth', 2) % plotting probability distribution function (pdf)

E2 = E(2) % 2nd energy (1st excited state energy)
psi_2 = psi(:,2); % 1st excited state wave function
plot(x,abs(psi_2).^2, 'Linewidth', 2) % plotting probability distribution function (pdf)
xlabel('x');
ylabel('|\psi(x)|^2');
legend('ground state p.d.f.','1st excited-state p.d.f.')
title('probability distribution function');


% plotting energy eigenvalues (theory: E(n) = (n+1/2) \hbar \omega),
% n=0,1,2,...
figure(3);
nE = 8;
a = zeros(nE,1);
b = E(1:nE);
plot(a,b, '-s','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6]) % plotting probability distribution function (pdf)
set(gca,'FontSize',16)
ylabel('E(n)');
title('Lowest 6 energy eigenvalues');



%%%% some sanity checks in linear algebra (and quatum mechanics)

% randn(n,m) --> generate an nxm matrix with a normal distribution with mean 0, st. dev. 1
% randn(n,m) --> generate an nxm matrix with a unifrom distribution with mean 0.5, st. dev. 1
M0 = 10 * (randn (3,3) + 1i * rand(3,3)) % M0 not Hermitian (sulf-conjugate)

M1 = (M0 + M0')/2 % M1 is Hermitian (sulf-conjugate)
% M1 = M1 * M1'
% M1 = M1' * M1


[U,Md] = eig(full(M1))


prod_U_Md_Udag = U * Md * U'

v1 = U(:,1) % |v_1>

v2 = U(:,2) % |v_2>

v3 = U(:,3) % |v_2>



% dot(A,B) = A' * B
overlap_11 = dot(v1,v1) % <v_1|v_1>
overlap_12 = dot(v1,v2) % <v_1|v_2>
overlap_13 = dot(v2,v3) % <v_2|v_1>
overlap_21 = dot(v2,v1) % <v_2|v_1>
overlap_22 = dot(v2,v2) % <v_1|v_1>
overlap_23 = dot(v2,v3) % <v_1|v_2>
overlap_31 = dot(v3,v1) % <v_2|v_1>
overlap_32 = dot(v3,v2) % <v_2|v_1>
overlap_33 = dot(v3,v3) % <v_1|v_1>

prod_U_Udag = U * U' % U * U^dag

prod_Udag_U = U' * U % U * U^dag



id_rersolutio = v1 * v1' +v2 * v2' + v3 * v3' % \sum_{j=1}^{n} |v_j> < v_j|
```

Code 4:

```MATLAB
%%%%%%%=======================================================%%%%%%%
%%%%%%%  Eigenvalue problem & Solvig Schrodinger Equation  %%%%%%%
%%%%%%%=======================================================%%%%%%%

clc;
close all;
clear;


dx = 0.2;
x = sparse(-4:dx:4);
y = x;
omega = 1; % omega = sqrt(k/m) = frequency of oscillations
m = 1; % mass
hbar = 1; % Planck constant

n_x = length(x);
n_y = n_x;
D1 = sparse((diag(ones(1,n_x-1),1)-diag(ones(1,n_x),0))/dx);
D1(n_x,1) = 1;
p1_x = (hbar/1i) * D1; % momentum operator = hbar/i d/dx
p1_y = p1_x;

p2_x = (p1_x'*p1_x+p1_x*p1_x')/2;
p2_y = p2_x;

id_x = speye(n_x);
id_y = speye(n_y);

X = kron(diag(x),id_y);
Y = kron(id_x, diag(y));

P1_x = kron(p1_x, id_y);
P1_y = kron(id_x, p1_y);

P2_x = kron(p2_x, id_y); % (P1_x'*P1_x+P1_x*P1_x')/2;
P2_y = kron(id_x, p2_y); % (P1_y'*P1_y+P1_y*P1_y')/2;


H = (P1_x'*P1_x+P1_x*P1_x')/(4*m) + (P1_y'*P1_y+P1_y*P1_y')/(4*m) + m*omega^2*(X^2+Y^2)/2;
Vxy = reshape(diag(m*omega^2*(X^2+Y^2)/2),n_y,n_x);

figure(1);
surface(x,y,Vxy,'LineWidth',4)
xlabel('x');
ylabel('y');
zlabel('V(x,y)')
title('Potential function');


[psi, Hd] = eigs(H,8,'sa'); % diagonalizing Hamiltonian operator

E = diag(Hd); % energy eigenvalues
E1 = E(1) % lowest energy (ground-state energy)
psi_1 = psi(:,1); % ground-state energy wavefunction
psi_1 = reshape(psi_1,n_y,[]);
figure(2);
surf(x,y,abs(psi_1).^2) % plotting probability distribution function (pdf)legend('ground state p.d.f.')
title('probability distribution function of the ground state');

E2 = E(2) % 1st excited energy
psi_2 = psi(:,2); % 1st excited energy
psi_2 = reshape(psi_2,n_y,[]);
figure(3);
surf(x,y,abs(psi_2).^2) % plotting probability distribution function (pdf)legend('1st excited-state p.d.f.')
title('probability distribution function of the 1st excited state');



% plotting energy eigenvalues (theory: E(n) = (n+1/2) \hbar \omega),
% n=0,1,2,...
figure(4);
nE = 8;
a = zeros(nE,1);
b = E(1:nE);
% b = b+0.1*rand(size(b));
plot(a,b, '-s','MarkerSize',10,...
      'MarkerEdgeColor','red',...
      'MarkerFaceColor',[1 .6 .6]) % plotting probability distribution function (pdf)
set(gca,'FontSize',16)
ylabel('E(n)');
title('Lowest 6 energy eigenvalues');
```

Solution is provided in the `ex1-python.ipynb` notebook.