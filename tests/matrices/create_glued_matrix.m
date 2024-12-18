function A = create_glued_matrix(m, p, s, r, t)
% A = CREATE_GLUED_MATRIX(m, p, s, r, t) generates matrices according to
% Example 2 from [Smoktunowicz et. al., 2006]. The resutling matrix is of
% size m x ps, with parameters r and t specifying the powers of the largest
% condition number of the first stage and second stages of the matrix,
% respectively:
%
% Stage 1: A = U * Sigma * V'
%
% Stage 2: A = A * kron(I, Sigma_block) * kron(I, V_block)
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
n = p * s;
rng(1); U = orth(randn(m,n));
rng(2); V = orth(randn(n,m));
Sigma = diag(logspace(0, r, n));
A = U * Sigma * V';

Sigma_block = diag(logspace(0, t, s));
rng(3); V_block = orth(randn(s,s));

ind = 1:s;
for i = 1:p
    A(:,ind) = A(:,ind) * Sigma_block * V_block';
    ind = ind + s;
end
end