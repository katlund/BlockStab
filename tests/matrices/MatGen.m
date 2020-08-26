function [XX, XXstr, XXprops] = MatGen(matstr, XXdim)
% [XX, XXstr, XXprops] = MATGEN(matstr, XXdim) generates the type of matrix
% specified by matstr with dimensions m x ps and block partitions of size s
% given by XXdim = [m, p, s].  It also saves the output in a .mat file in
% the format
%
%   sprintf('%s_m%d_p%d_s%d', matstr, m, p, s)
%
% and generates a string (XXstr) and a struct (XXprops) for use in RUNTEST.
%
% Possible options for matstr include
%   'rand_uniform' - random entries drawn from uniform distribution
%   'rand_normal' - random entries drawn from normal distribution
%   'rank_def' - like rand_uniform but with a block vector set to 100 times
%       another, in order to force rank deficiency
%   'laeuchli' - the classic L�uchli matrix
%   'monomial' - each block vector spans a Krylov subspace of size s
%   'stewart' - a matrix with a geometric sequence of singular values,
%       ranging from 1 to 10^-t, t = 20; 1st column set to 25 column, and
%       35th column set to 0
%   'stewart_extreme' - a matrix with a geometric sequence of singular
%       values, ranging from 1 to 10^-t, t = 10; second half of singular
%       values set to identically zero
%   'hilbert' - the Hilbert matrix, generated by the built-in function
%      HILB(m,p*s)
%   's-step' and 'newton' - similar to 'monomial', but both generate
%       matrices that span a Krylov subspace of size sp
%   'glued' - 

%%
% Extract dimensions
m = XXdim(1); p = XXdim(2); s = XXdim(3);
n = p*s;
switch matstr
    case {'rand_uniform'}
        XX = rand(m,n);
        XXstr = sprintf('Random matrix from uniform distribution');
        
    case {'rand_normal'}
        XX = randn(m,n);
        XXstr = sprintf('Random matrix from normal distribution');
        
    case {'rank_def'}
        XX = randn(m,n);
        XX(:,1:s) = 100*XX(:,end-s+1:end);
        XXstr = sprintf('Rank-deficient matrix from normal distribution');
        
    case {'laeuchli'}
        XX = zeros(m,n);
        XX(1,:) = 1;
        for i = 2:n
            XX(i,i-1) = eps + (sqrt(eps)-eps)*rand(1);
        end
        XXstr = sprintf('Laeuchli matrix');
        
    case {'monomial'}
        A = spdiags(linspace(.1,10,m)',0,m,m);                              % mimics an operator
        pp = 1:p;
        XXhat = rand(m,p);
        XXhat(:,pp) = XXhat(:,pp)/norm(XXhat(:,pp));
        for i = 2:s
            pp = pp + p;
            XXhat(:,pp) = A*XXhat(:,pp - p);
        end
        % Reshape XX
        XX = zeros(m,n);
        ind = 1:s:n;
        jj = 1:p;
        for j = 1:s
            XX(:,ind) = XXhat(:,jj);
            ind = ind + 1;
            jj = jj + p;
        end
        XXstr = sprintf('Monomial matrix (s-step-like block vectors)');
        
    case {'stewart'}        
        [U,~] = qr(randn(m,n),0);
        [V,~] = qr(randn(n,m),0);
        t = 20;                                                             % higher t leads to worse conditioning
        S = diag(10.^(linspace(0,-t,n))');
        XX = U*S*V';
        XX(:,25) = XX(:,1);
        XX(:,35) = 0;
        XXstr = sprintf('Stewart matrix');
        
    case {'stewart_extreme'}        
        [U,~] = qr(randn(m,n),0);
        [V,~] = qr(randn(n,m),0);
        t = 10;                                                             % higher t leads to worse conditioning
        S = diag([10.^(linspace(0,-t,n/2)) zeros(1,n/2)]');
        XX = U*S*V';
        XXstr = sprintf('Stewart extreme matrix');
        
    case {'hilbert'}
        XX = hilb(m); XX = XX(:,1:n);
        XXstr = sprintf('Hilbert matrix');
        
    case {'s-step'}
        A = spdiags(linspace(.1,10,m)',0,m,m);
        x = rand(m,1); x = x/norm(x);
        XX = zeros(m,n);
        XX(:,1) = x;
        ind = 1:s;
        for j = 1:p
            for i = 2:s
                x = A*XX(:,ind(i-1));
                XX(:,ind(i)) = x/norm(x);
            end
            if j < p
                X = XX(:,ind);
%                 [X, ~] = qr(X,0);       % hard scaling
                x = A*X(:,end);
                ind = ind + s;
                XX(:,ind(1)) = x/norm(x);
            end
        end
        XXstr = sprintf('s-step matrix');
        
    case {'newton'}
        evec = linspace(.1,10,m);
        A = spdiags(evec',0,m,m);
        x = rand(m,1);
        alp = lejaorder(evec,p*s+1);

        XX(:,1) = (A-alp(1).*eye(m))*x;
        XX(:,1) = XX(:,1)./norm(XX(:,1));

        for i = 0:p-1
            for j = 1:s
                XX(:,i*s+j+1) = (A-alp(i*s+j+1)*eye(m))*XX(:,i*s+j);
                XX(:,i*s+j+1) = XX(:,i*s+j+1)./norm(XX(:,i*s+j+1));
            end
        end
        XX = XX(:,1:p*s); 
        XXstr = sprintf('Newton matrix');        
end
XXprops.cond = cond(XX);
XXprops.sv = svd(XX);
XXprops.normF = norm(XX,'fro');

savestr = sprintf('%s_m%d_p%d_s%d', matstr, m, p, s);
save(savestr, 'XX', 'XXstr', 'XXprops');
end

%% Auxiliary functions
function z = lejaorder(x,n)
     % Leja order of the first n components of the vector x
     m = length(x);
     n = min(n,m);
     z = zeros(size(x));
     [v index] = max(abs(x));
     z(1) = x(index);
     temp = abs(x-z(1));
     [v index] = max(temp);
     z(2) = x(index);
     for k = 2:n-1
         for i = 1:m
             temp(i) = temp(i)*abs(x(i)-z(k));
         end
         [v index] = max(temp);
         z(k+1) = x(index);
     end
     z = z(1:n);
end