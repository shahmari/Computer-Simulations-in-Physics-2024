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