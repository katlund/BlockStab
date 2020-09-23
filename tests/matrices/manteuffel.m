% Generate Manteuffel matrix
function A = manteuffel(n, beta)
if nargin == 0
    n = 25;
    beta = 1;
elseif nargin == 1
    beta = 1;
end
h = 1/(n+1);
M = gallery('poisson', n);  % matrix will be n^2
e = ones(n,1);
z = zeros(n,1);
N0 = spdiags([-e z e], -1:1, n, n);
N1 = spdiags(e, 0, n, n);
T = spdiags(e, 1, n, n);
N = kron(eye(n), N0) + kron(T, N1) + kron(T', -N1);

A = M / h^2 + beta * N / (2*h);
end