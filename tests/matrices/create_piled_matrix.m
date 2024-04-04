function A = create_piled_matrix(m, p, s, c, tol)
% A = CREATE_PILED_MATRIX(m, p, s, t, tol) generated matrices designed to
% challenged BCGS-type algorithms. The resutling matrix is of size m x ps,
% structured as follows:
%
%       A = [A_1, A_2, ..., A_p],
%       with A_{i+1} = A_i + tol*U*D*V',
%
% where U and V are unitary, and D is a diagonal matrix with condtion
% number 10^c. The default value of tol is 1.

%%
if nargin == 4
    tol = 1;
end

rng(1); U = orth(randn(m, s));
rng(100); V = qr(randn(s, s));
D = diag(logspace(0, 4, s));

n = p*s;
A = zeros(m, n);
ind = 1:s;
A(:, ind) = U * D * V';
ind = ind + s;
for i = 2:p
    rng(i); U = orth(randn(m, s));
    rng(100*i); V = orth(randn(s, s));
    D = diag(logspace(0, c, s));
    A(:, ind) = A(:, ind-s) + tol * U * D * V';
    ind = ind + s;
end
end