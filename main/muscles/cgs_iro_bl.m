function [Q, R] = cgs_iro_bl(X, verbose)
% Taken from https://github.com/dbielich/DCGS2 and modified with our
% verbose statements

%%
[m,n] = size(X);
Q = X;
R = zeros(n,n);

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
end

for k = 2:n
    [Q(1:m,k-1:k), R(1:k-1,k-1:k)] = orth_dcgs2_qr(Q(1:m,1:k), R(1:k-2,k-1));

    if verbose
        fprintf('%3.0d:', k-1);
        fprintf('  %2.4e  |',...
            norm( eye(k-1) - Q(:,1:k-1)' * Q(:,1:k-1) ) );
        fprintf('  %2.4e\n',...
            norm( X(:,1:k-1) - Q(:,1:k-1) * R(1:k-1,1:k-1) ) / norm(X(:,1:k-1)) );
    end
end
[Q(1:m,n), R(1:n,n)] = orth_dcgs2_qr_cleanup(Q(1:m,1:n), R(1:n-1,n));
if verbose
    fprintf('%3.0d:', k);
    fprintf('  %2.4e  |', norm( eye(n) - Q' * Q ) );
    fprintf('  %2.4e\n', norm( X - Q * R ) / norm(X) );
end
end

%% Subroutines that add one vector at a time
function [q, r] = orth_dcgs2_qr(Q, RR)
m = size(Q,1);
k = size(Q,2);

r = zeros(k-1,2);
q = zeros(m,2);
if k == 2
    r(1,1) = Q(1:m,1)' * Q(1:m,1);
    r(1,2) = Q(1:m,1)' * Q(1:m,2); 
    r(1,1) = sqrt(r(1,1));
    r(1,2) = r(1,2) / r(1,1);
    q(1:m,1) = Q(1:m,1) / r(1,1);
    q(1:m,2) = Q(1:m,2) - q(1:m,1) * r(1,2); 
end

if k >= 3
    tmp1(1:k-1,1:2) = Q(1:m,1:k-1)' * Q(1:m,k-1:k);
    tmp2(1:k-2,1) = tmp1(1:k-2,1);

    r(1:k-1,2) = tmp1(1:k-1,2);
    r(k-1,1) = tmp1(k-1,1);
    
    r(k-1,2) = r(k-1,2) - tmp2(1:k-2,1)' * r(1:k-2,2);
    r(k-1,1) = r(k-1,1) - tmp2(1:k-2,1)' * tmp2(1:k-2,1); 
    r(1:k-2,1) = RR(1:k-2,1) + tmp2(1:k-2,1); 
    r(k-1,1) = sqrt( r(k-1,1) );
    r(k-1,2) = r(k-1,2) / r(k-1,1);
    
    Q(1:m,k-1) = Q(1:m,k-1) - Q(1:m,1:k-2) * tmp2(1:k-2,1); 
    Q(1:m,k) = Q(1:m,k) - Q(1:m,1:k-2) * r(1:k-2,2);
    
    q(1:m,1) = Q(1:m,k-1) / r(k-1,1);
    q(1:m,2) = Q(1:m,k) - q(1:m,1) * r(k-1,2);
end
end

function [q, r] = orth_dcgs2_qr_cleanup(Q, R)
m = size(Q,1);
n = size(Q,2); 

r = zeros(n,1);
q = zeros(m,1);

r(1:n,1) = Q(1:m,1:n)' * Q(1:m,n); 
r(n,1) = r(n,1) - r(1:n-1,1)' * r(1:n-1,1);
r(n,1) = sqrt( r(n,1) );  
      
Q(1:m,n) = Q(1:m,n) - Q(1:m,1:n-1) * r(1:n-1,1);  
r(1:n-1,1) = R(1:n-1,1) + r(1:n-1,1);  

q(1:m,1) = Q(1:m,n) / r(n,1); 
end 