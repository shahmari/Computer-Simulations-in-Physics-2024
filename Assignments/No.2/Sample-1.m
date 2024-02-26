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