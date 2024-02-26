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