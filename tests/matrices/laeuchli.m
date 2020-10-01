function XX = laeuchli(m, n, eps_val)
% XX = LAEUCHLI(m, n, eps_val) generates a Laeuchli matrix with the
% specified dimensions and epsilon value

%%
XX = zeros(m,n);
XX(1,:) = 1;
for i = 2:n
    XX(i,i-1) = eps_val;
end
end