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