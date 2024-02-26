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