# Assignment 2:

Convert the following code from MATLAB to Python:

## Code 1:

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

An example is available in the file 'Example-1.py'. The code is also available below:

```python
# Credits to: Ali Ashtari Jafari

#%matplotlib qt
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import time # I wanted to evaluate the run time of my program
t0=time.time()
# Solving following hamiltonian for fermions
# H = \sum_{j1=1:Np} h_1p^{j1} + V1 * \sum_{i1}\sum_{j1<j2}
Np=3 # number of particles
L=10 # size of system
V1= 5 # interaction strength
V2= 0 # for part 2
V3= -5 # for part 2
Id=sp.eye(L, L, k=0, format='csr') # single particle identity operator
bc=1 # boundary condition
T=sp.csr_matrix(np.diag(np.ones(L-1), 1) + np.diag([bc], -(L-1))) # translation operator
h_1p=-(T+T.T) # in fact, it should be T+T.T.conjugate() accorfing to the main code, but it's real, so its adjoint is itselt
n_1p=[]
empty=sp.csr_matrix(np.zeros((L, L)))
for i in range(L):
    n_1p.append(empty.copy())
    n_1p[i][i, i]=1

# warmup: 2 fermion problem
H=sp.kron(h_1p, Id) + sp.kron(Id, h_1p) # forming the hamiltonian
for i in range(L):
    # imposing periodic boundary condition for next point
    if i==L-1:
        i1_p=0
    else:
        i1_p=i%L + 1
    # imposing periodic boundary condition for previous point
    if i==0:
        i1_m= L-1
    elif i==1:
        i1_m=0
    else:
        i1_m=(i-2)%L + 1
    H = H + V1*sp.kron(n_1p[i], n_1p[i1_p] + n_1p[i1_m])
# making permutation operator
index_12= np.reshape(np.arange(L**2), (L, L))
index_21= np.reshape(np.transpose(index_12, (1, 0)), -1)
A=np.eye(L**2, L**2)
P12=sp.csr_matrix(A[:, index_21])
# diagonalization and imposing fermionic statistic
HP= H + ((10**4) * P12)
HP=(HP + HP.T.conjugate())/2 # ensuring HP is hermitian
nE= min(np.math.comb(L, 2), HP.shape[0])
eigenvalues, eigenvectors= sp.linalg.eigs(HP, k=nE, which='SR', tol=1e-16) # this is sorted with ascending order
indices=eigenvalues.argsort()
eigenvalues=eigenvalues[indices]
eigenvectors=sp.csr_matrix(eigenvectors[:, indices])
UH=eigenvectors
E=np.real(np.diag(((UH.T.conjugate()).dot(H.dot(UH))).toarray())) # energy eigenvalues
E_2particle=E[: min(10, len(E))] # 10 lowest energy eigenvalues
dE_2particle=E_2particle[1:] - E_2particle[0] # 10 lowest excitation energies
# plotting 10 lowest excitation energies
dE_labels=[]
for i in range(len(dE_2particle)):
    dE_labels.append(fr'$\Delta E_{i+1} = {dE_2particle[i]:.4f}$')
plt.figure()
plt.bar(dE_labels, dE_2particle)
plt.xlabel('eigenenergies')
plt.ylabel('values')
plt.show()

### General problem: Np>2 particles
H=h_1p
H2=h_1p # for part 2
H3=h_1p # for part 2
ID=Id
Id_prev=1
N_tot=[] # The structure for saving density operators
for i1 in range(L):
    N_tot.append(n_1p[i1])
Perm={} # The dictionary structure for saving fermionic permutation operators
for n_p in range(2, Np+1):
    H=sp.kron(H, Id) + sp.kron(ID, h_1p)
    H2=sp.kron(H2, Id) + sp.kron(ID, h_1p)
    H3=sp.kron(H3, Id) + sp.kron(ID, h_1p)
    for i1 in range(L):
        if i1==L-1:
            i1_p=0
        else:
            i1_p=i1%L + 1
        if i1==0:
            i1_m=L-1
        elif i1==1:
            i1_m=0
        else:
            i1_m=(i1-2)%L + 1
        H = H + V1*sp.kron(N_tot[i1], n_1p[i1_p] + n_1p[i1_m])
        H2 = H2 + V2*sp.kron(N_tot[i1], n_1p[i1_p] + n_1p[i1_m]) # for part 2
        H3 = H3 + V3*sp.kron(N_tot[i1], n_1p[i1_p] + n_1p[i1_m]) # for part 2
    for i1 in range(L):
        N_tot[i1]=sp.kron(N_tot[i1], Id) + sp.kron(ID, n_1p[i1])
    Perm[n_p]= sp.kron(Id_prev, P12)
    for ip in range(2, n_p):
        Perm[ip]= sp.kron(Perm[ip], Id)
    Id_prev=ID
    ID=sp.kron(ID, Id)
    # diagonalizing the hamiltonian and imposing fermionic statistics
    Lambda=1e4
    P_tot=sp.csr_matrix(Perm[2].shape)
    for ip in range(2, n_p+1):
        P_tot= P_tot+Perm[ip]
    HP= H + Lambda*P_tot
    HP=(HP.T.conjugate() + HP)/2 # ensuring HP is hermitian
    HP2= H2 + Lambda*P_tot
    HP2=(HP2.T.conjugate() + HP2)/2 # for part 2
    HP3= H3 + Lambda*P_tot
    HP3=(HP3.T.conjugate() + HP3)/2 # for part 2
    nE=min(np.math.comb(L, n_p), HP.shape[0])
    eigenvalues1, eigenvectors1=sp.linalg.eigs(HP, k=nE, which='SR')
    # evaluating and sorting eigenvalues1 ascending
    indices1=eigenvalues1.argsort()
    eigenvalues1=eigenvalues1[indices1]
    eigenvectors1=sp.csr_matrix(eigenvectors1[:, indices1])
    # evaluating and sorting eigenvalues2 ascending
    eigenvalues2, eigenvectors2=sp.linalg.eigs(HP2, k=nE, which='SR') # for part 2
    indices2=eigenvalues2.argsort()
    eigenvalues2=eigenvalues2[indices2]
    eigenvectors2=sp.csr_matrix(eigenvectors2[:, indices2])
    # evaluating and sorting eigenvalues3 ascending
    eigenvalues3, eigenvectors3=sp.linalg.eigs(HP3, k=nE, which='SR') # for part 2
    indices3=eigenvalues3.argsort()
    eigenvalues3=eigenvalues3[indices3]
    eigenvectors3=sp.csr_matrix(eigenvectors3[:, indices3])
    #
    UH1=eigenvectors1
    E1=np.real(np.diag(((UH1.T.conjugate()).dot(H.dot(UH1))).toarray()))
    E_np_particles=E1[:min(10, len(E1))]
    dE_np_particles=E_np_particles[1:] - E_np_particles[0]
    print(dE_np_particles)
    dE_np_labels=[]
    for i in range(len(dE_np_particles)):
        dE_np_labels.append(fr'$\Delta E_{i+1} = {dE_np_particles[i]:.4f}$')
    plt.figure()
    plt.bar(dE_np_labels, dE_np_particles)
    plt.xlabel('eigenenergies')
    plt.ylabel('values')
    plt.title(f'Excitation energies for Np={n_p} particles')
    plt.show()
    
### Computing average number of particles and correlation function of density_density

# Because in this part we are dealing with Mp=3 particles, and the last iteration in the previous part gives us N_tot for Np=3 particle, I use N_tot from the previous part.
# Using ground state wave functions to calculating the desired expectation values
psi_1=eigenvectors1[:, 0]
psi_2=eigenvectors2[:, 0]
psi_3=eigenvectors3[:, 0]
# Average number of particles
N1_av, N2_av, N3_av= [], [], [] 
for i in range(L):
    N1_av.append(np.real((((psi_1.T.conjugate()).dot(N_tot[i].dot(psi_1))).toarray())[0, 0]))
    N2_av.append(np.real((((psi_2.T.conjugate()).dot(N_tot[i].dot(psi_2))).toarray())[0, 0]))
    N3_av.append(np.real((((psi_3.T.conjugate()).dot(N_tot[i].dot(psi_3))).toarray())[0, 0]))
print(N1_av)
print(N2_av)
print(N3_av)
X=[i for i in range(L)]
plt.figure()
plt.plot(X, N1_av, 'o', label='average number of particles for V1= 5')
plt.plot(X, N2_av, 'o',  label='average number of particles for V2= 0')
plt.plot(X, N3_av, 'o',  label='average number of particles for V3= -5')
plt.xlabel('X')
plt.ylabel(r'$\langle N_i \rangle$')
plt.grid()
plt.legend()

# Correlation function
def corr(A, B, v):
    return np.real(((v.T.conjugate()).dot((A.dot(B)).dot(v)) - ((v.T.conjugate()).dot(A.dot(v)))*((v.T.conjugate()).dot(B.dot(v)))).toarray()[0, 0])
Corr_1, Corr_2, Corr_3=[], [], []
for i in range(L):
    Corr_1.append(corr(N_tot[0], N_tot[i], psi_1))
    Corr_2.append(corr(N_tot[0], N_tot[i], psi_2))
    Corr_3.append(corr(N_tot[0], N_tot[i], psi_3))
print(Corr_1)
print(Corr_2)
print(Corr_3)
plt.figure()
plt.plot(X, Corr_1, 'o', label='correlation function for V1= 5')
plt.plot(X, Corr_2, 'o',  label='correlation function for V2= 0')
plt.plot(X, Corr_3, 'o',  label='correlation function for V3= -5')
plt.xlabel('X')
plt.ylabel(r'$Corr (N_i, N_j)$')
plt.grid()
plt.legend()

t1=time.time()
print(f'Run time of the program is: {t1 - t0} s')
```

## Code 2:

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
