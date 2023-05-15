function A = CreateGluedMatrix(m, p, s, r, t)
% Example 2 matrix from [Smoktunowicz et. al., 2006]. Generates a glued
% matrix of size m x ps, with parameters r and t specifying the powers of
% the largest condition number of the first stage and second stages of the
% matrix, respectively:
%
% Stage 1: A = U * Sigma * V'
%
% Stage 2: A = A * kron(I, Sigma_block) * kron(I, V_block)

%%
n = p * s;
U = orth(randn(m,n));
V = orth(randn(n,m));
Sigma = diag(logspace(0, r, n));
A = U * Sigma * V';

Sigma_block = diag(logspace(0, t, s));
V_block = orth(randn(s,s));

ind = 1:s;
for i = 1:p
    A(:,ind) = A(:,ind) * Sigma_block * V_block';
    ind = ind + s;
end
end