# Assignment 2:

Convert the following code from MATLAB to Python:

Code 1:

```matlab
clc;
clear;
close all;

%%%%% Below, we solve the fullowing Hamiltonian for Np fermions
%%%%% H = \sum_{j1=1:Np} h_1p^{j1} + V1 * \sum_{i1}\sum_{j1<j2}
%%%%% n_1p{i1}^{j1}(n_1p{i1+1}^{j2}+n_1p{i1-1}^{j2})


Np = 3; % number of particles
L = 10; % (1D) system size 
V1 = 5; % interaction strength (coupling constant)
id = speye(L,L); % sigle particle identity operator
bc = 1; % boundary condition
T = sparse(diag(ones(1,L-1),1) + bc*diag(1,-(L-1))); % translation operator

% h_1p = 1-particle kietic energy. In general h_1p = -hbar^2*(T+T'-2*id)/(2*M*dx^2)
% here, we have assumed hbar^2/(M*dx^2) = 1; and have removed 2*id part
h_1p = -(T+T');  

n_1p = cell(1,L); % saving 1-particle density operator into a cell struncture
empty =  sparse(L,L); % empty matrix
for i1=1:L
    n_1p{i1} = empty; % 1-particle density operator at site i1
    n_1p{i1}(i1,i1) = 1;
end



%%% warmup: 2 fermion problem

H = kron(h_1p,id) + kron(id,h_1p);
for i1 = 1:L
    i1_p = mod(i1,L)+1; %i1_p = i1 + 1 mod L       
    i1_m = mod(i1-2,L)+1; %i1_p = i1 - 1 mod L
    H = H + V1*kron(n_1p{i1},n_1p{i1_p}+n_1p{i1_m});
end

% making permutation operator
index_12 = reshape(1:L^2,L,L);
index_21 = reshape(permute(index_12,[2,1]),[],1);
A = speye(L^2,L^2);
P12 = A(:,index_21);

%%%%% diagonalization + imposing fermionic statistics
HP = H + 10^4*P12;
HP = (HP + HP')/2; % ensuring HP is hermitian (the non-hermiticity might be caused by the computer's round-off error)
nE = min(nchoosek(L,2), size(HP,1));  % number of desired eigenvalues
[UH,DH] = eigs(HP, nE, 'sr','Tolerance',1e-16);
DH = real(diag(DH));
% sorting eigenavlues and the eigenvectors accordingly
[DH, S] = sort(DH,'ascend');
UH = sparse(UH(:,S));
E = diag(real(UH' * H * UH)); % energy eigenvalues

V1

L


E_2particle = full(E(1:min(10,length(E))))  % energy eigenvalues (for fermions)

dE_2particle = E_2particle(2:end)-E_2particle(1)  % excitation energies


%%%%% general problem: Np > 2 particles
%%%%% Below, we generate and update the Hamiltonian and other operators iteratively

% we start with 1 particle Hamiltonian
H = h_1p;
Id = id;
Id_prev = 1; % will be used to generate P_{i,i+1}
N_tot = cell(L,1); % N_tot{i1} = sum_{j} n_1p{i1}^{j} the sum of denisty for particles 1:j
for i1 = 1:L
    N_tot{i1} = n_1p{i1}; 
end

Perm = cell(Np,1); % A cell containg adjacent electrons permutation operator
for np = 2:Np
    
    
    %%%%%%%
    H = kron(H,id)+ kron(Id,h_1p);
    for i1 = 1:L
        i1_p = mod(i1,L)+1; %i1_p = i1 + 1 mod L       
        i1_m = mod(i1-2,L)+1; %i1_p = i1 - 1 mod L
        H = H + V1*kron(N_tot{i1},n_1p{i1_p}+n_1p{i1_m});
    end
    
    for i1 = 1:L
        N_tot{i1} = kron(N_tot{i1},id) + kron(Id, n_1p{i1});
    end
    
    Perm{np} = kron(Id_prev, P12); % Perm{p} = P_{p-1,p} 
    for ip = 2:np-1
        Perm{ip} = kron(Perm{ip},id);
    end
    
    Id_prev = Id; % identity operator in the current step
    Id = kron(Id, id);
    
    %%%%% diagonalization + imposing fermionic statistics
    
    lambda = 1e4; % lamda = 10^4
    
    P_tot = sparse(0); % P_tot = sum_{i=1}^{np} P_{i-1,i}
    for ip = 2:np
        P_tot = P_tot + Perm{ip};
    end
    HP = H + lambda*P_tot;
    HP = (HP + HP')/2; % ensuring HP is hermitian (the non-hermiticity might be caused by the computer's round-off error)
    nE = min(nchoosek(L,np), size(HP,1)); % number of desired eigenvalues
    [UH,DH] = eigs(HP, nE, 'sr','Tolerance',1e-16);
    DH = real(diag(DH));
    % sorting eigenavlues and the eigenvectors accordingly
    [DH, S] = sort(DH,'ascend'); 
    UH = sparse(UH(:,S));
    E_Np_particle = full(diag(real(UH' * H * UH)));  % energy eigenvalues
    
end

Np


E_Np_particle = E_Np_particle(1:min(10,length(E_Np_particle))) % energy eigenvalues (for fermions)

dE_Np_particle = E_Np_particle(2:end)- E_Np_particle(1) % excitation energies
```

Code 2:

```MATLAB
clear; % clearig RAM
clc; % clearing screen (command widow)
close all; % closing all open plots/figures

% Pauli matrics

L = 16;
J1 = 1; % Nearest neighbor Heisenberg coupling
J2 = 0.0; % Next nearest neighbor Heisenberg coupling
bnd = 1; % boundary condition : 1--> PBC  0 --> OBC  -1 --> APBC
hbar = 1;

tic; % start time 

% pauli matrices (defined as sparse matrices)
sigma_0 = speye(2);
sigma_z = sparse([1,0; 0,-1]);
sigma_x = sparse([0, 1; 1,0]);
sigma_y = sparse([0,-1i; 1i,0]);
sigma_p = (sigma_x + 1i * sigma_y)/2;
sigma_m = (sigma_x - 1i * sigma_y)/2;
n = (sigma_p * sigma_m); % number operator



% spin operators for two sites problem
sx_1 = kron(hbar/2*sigma_x, sigma_0);
sx_2 = kron(sigma_0, hbar/2*sigma_x);

sy_1 = kron(hbar/2*sigma_y, sigma_0);
sy_2 = kron(sigma_0, hbar/2*sigma_y);

sz_1 = kron(hbar/2*sigma_z, sigma_0);
sz_2 = kron(sigma_0, hbar/2*sigma_z);

% Heisenberg Hamiltonian for two sites
H = J1 * (sx_1*sx_2 + sy_1*sy_2 + sz_1*sz_2);

[U, D] = eig(full(H)); % digonalizing Hamiltonian



% for three sites
sx_1 = kron(hbar/2*sigma_x, kron(sigma_0, sigma_0));
sx_2 = kron(sigma_0, kron(hbar/2*sigma_x, sigma_0));
sx_3 = kron(sigma_0, kron(sigma_0, hbar/2*sigma_x));

sy_1 = kron(hbar/2*sigma_y, kron(sigma_0, sigma_0));
sy_2 = kron(sigma_0, kron(hbar/2*sigma_y, sigma_0));
sy_3 = kron(sigma_0, kron(sigma_0, hbar/2*sigma_y));

sz_1 = kron(hbar/2*sigma_z, kron(sigma_0, sigma_0));
sz_2 = kron(sigma_0, kron(hbar/2*sigma_z, sigma_0));
sz_3 = kron(sigma_0, kron(sigma_0, hbar/2*sigma_z));


H = J1 * ((sx_1*sx_2 + sy_1*sy_2 + sz_1*sz_2) + (sx_2*sx_3 + sy_2*sy_3 + sz_2*sz_3));

[U, D] = eig(full(H));

% spin operators for three sites
sx_1 = kron(hbar/2*sigma_x, kron(sigma_0, kron(sigma_0,sigma_0)));
sx_2 = kron(sigma_0, kron(hbar/2*sigma_x, kron(sigma_0,sigma_0)));
sx_3 = kron(sigma_0, kron(sigma_0, kron(hbar/2*sigma_x,sigma_0)));
sx_4 = kron(sigma_0, kron(sigma_0, kron(sigma_0,hbar/2*sigma_x)));



% general problem: spin operators for n sites
sx{1} = hbar/2 * sigma_x;
sy{1} = hbar/2 * sigma_y;
sz{1} = hbar/2 * sigma_z;
id = sigma_0;
for num_sites = 2: L
    for pos = 1:num_sites-1
        sx{pos} = kron(sx{pos}, sigma_0);
        sy{pos} = kron(sy{pos}, sigma_0);
        sz{pos} = kron(sz{pos}, sigma_0);
    end
    sx{num_sites} = kron(id, hbar/2*sigma_x);
    sy{num_sites} = kron(id, hbar/2*sigma_y);
    sz{num_sites} = kron(id, hbar/2*sigma_z);
    id = kron(id, sigma_0);
end

H = sparse(0);
% nearest neighbor interaction
for pos = 1:L-1
    H = H + J1 * (sx{pos}*sx{pos+1} + sy{pos}*sy{pos+1} + sz{pos}*sz{pos+1});
end
if bnd == 1
    H = H + bnd * J1 * (sx{1}*sx{L} + sy{1}*sy{L} + sz{1}*sz{L});
end

% next nearest neighbor interaction
for pos = 1:L-2
    H = H + J2 * (sx{pos}*sx{pos+2} + sy{pos}*sy{pos+2} + sz{pos}*sz{pos+2});
end
if bnd == 1
    H = H + bnd * J2 * (sx{1}*sx{L-1} + sy{1}*sy{L-1} + sz{1}*sz{L-1} + sx{2}*sx{L} + sy{2}*sy{L} + sz{2}*sz{L});
end

% all eigenvalues ( indeed for > 2-3% of eigenvalues are needed)
% [U, D] = eig(full(H));

% Finding a number of lowest eigenvalues and eigenvectors only
n_eigs = min(10, size(H,1));
[U, D] = eigs(H, n_eigs, 'sa'); % 'sa', 'sr', 'la', 'lr'
E = diag(D)

psi = U(:,1); % |psi> ground-state wavefunction
% computing spin-spin correlationn functions
SS_corr = zeros(L);
for i1 = 1:L
    for i2 = i1+0:L 
        O = sx{i1}*sx{i2}+sy{i1}*sy{i2}+sz{i1}*sz{i2};
        SS_corr(i1,i2) = psi' * O * psi; % <psi| s_2.s_8 |psi>
        SS_corr(i2,i1) = SS_corr(i1,i2);
    end
end
SS_corr

% plotting correlation functions
x = (1:L)-round(L/2);
y = SS_corr(round(L/2),:);
plot(x,y,'-d');
xlabel('distance')
ylabel('Spin-Spin correlation')
elapsed_time = toc % finish time
```
