function [y, r, rho, northog] = cgs_step_sror(Q, x, nu, rpltol)
% [y, r, rho, northog] = CGS_STEP_SROR(Q, x, nu, rpltol) orthogonalizes x
% against the columns of Q using the the Classical Gram-Schmidt method.
% Reorthogonalization is performed as necessary to ensure orthogonality.
% Specifically, given an orthonormal matrix Q and a vector x, the function
% returns a vectors y and r, and a scalar rho satisfying
% 
% (*)   x = Q*r + rho*y
% 
% where y is a normalized vector orthogonal to the column space of Q.
% The matrix Q may have zero columns.
% 
% The optional argument nu is the norm of the original value of x, to be
% used when x has been subject to previous orthogonalizaton steps.  If
% nu is absent or is less than or equal to norm(x), is is set to norm(x).
% 
% The optional argument rlptol controls when the current vector y is
% replaced by a random vector.  Its default value is 1.  If it is set
% to a value greater than one, the relation (*) will be compromised
% somewhat, but the number of orthogonalizations may be decreased.
% 
% The optional output argument northog, if present, contains the number
% of orthogonalizations.
%
% Originally written by G. W. Stewart, 2008. Modified by Kathryn Lund,
% 2020.

%%
% Initialize
[n, nq] = size(Q);
r = zeros(nq, 1);
nux = norm(x);
if nargout >= 4
   northog = 0;
end

% If Q has no columns, return normalized x if x~=0, otherwise return a
% normalized random vector.
if nq == 0
    if nux == 0
        y = rand(n,1) - 0.5;
        y = y / norm(y);
        rho = 0;
    else
        y = x / nux;
        rho = nux;
    end
    return
end

%  If required, intitalize the optional arguments nu and rpltol.
if nargin < 3
  nu = nux;
else
  if nu < nux
     nu = nux;
  end
end

if nargin < 4
    rpltol = 1;
end

% If norm(x)==0, set it to a nonzero vector and proceed, noting the fact in
% the flag zeronorm.
if nux ~= 0
    zeronorm = false;
    y = x / nux;
    nu = nu / nux;
else
    zeronorm = true;
    y = rand(n,1) - 0.5;
    y = y / norm(y);
    nu = 1;
end


% Main orthogonalization loop.
nu1 = nu;
while true
    if nargout == 4
        northog = northog + 1;
    end
    s = Q' * y;
    r = r + s;
    y = y - Q * s;
    nu2 = norm(y);

    % Return if y is orthogonal.
    if nu2 > 0.5*nu1
        break
    end

    % Continue orthogonalizing if nu2 is not too small.
    if nu2 > rpltol * nu * eps
        nu1 = nu2;
    else % Replace y by a random vector.
        nu = nu * eps;
        nu1 = nu;
        y = rand(n,1) - 0.5;
        y = nu * (y / norm(y));
    end
end

% Calculate rho and normalize y.
if ~zeronorm
    rho = norm(y);
    y = y / rho;
    rho = rho * nux;
    r = r * nux;
else
    y = y / norm(y);
    r = zeros(nq, 1);
    rho = 0;
end