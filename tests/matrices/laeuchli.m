function XX = laeuchli(m, n, eps_val)
% XX = LAEUCHLI(m, n, eps_val) generates a Laeuchli matrix with the
% specified dimensions and epsilon value
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
XX = zeros(m,n);
XX(1,:) = 1;
for i = 2:n
    XX(i,i-1) = eps_val;
end
end