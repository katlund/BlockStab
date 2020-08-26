function [Q, R] = cgs_sror(X, rpltol)
% [Q, R] = CGS_SROR(X, rpltol) performs Classical Gram-Schmidt with
% Selective ReOrthogonalization and Replacement on the m x s matrix X as
% described in [Stewart 2008]. The core part of the routine is contained in
% CGS_STEP_SROR.

%%
% Pre-allocate memory for Q and R
[m, s] = size(X);
Q = zeros(m, s);
R = zeros(s);

% Default for rpltol
if isempty(rpltol)
    rpltol = 1;
end

% Classical Gram-Schmidt with CGS_STEP_SROR
[Q(:,1), ~, R(1,1)] = cgs_step_sror(zeros(m,0), X(:,1), rpltol);
for k = 1:s-1
    [Q(:,k+1), R(1:k,k+1), R(k+1,k+1)] =...
        cgs_step_sror(Q(:,1:k), X(:,k+1), rpltol);
end
end
