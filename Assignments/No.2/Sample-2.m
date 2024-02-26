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