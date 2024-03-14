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