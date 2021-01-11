function [XX,RR] = breducegs(XX, s, local_m, LMstr)
% [XX,RR] = breducegs(XX, s, LGSstr, ReduceGSstr, verbose) performs block reduce 
% Gram-Schmidt algorithm
% TODO: Document parameters

addpath(genpath('../'))
addpath(genpath('../main/'))

X00 = XX;

function [X1, r1] = bmgs_step(Q1, X1)
  X0 = X1;
  [m1, s1] = size(X1);
  [m1, n1] = size(Q1);
  p1 = n1/s1;
  r1 = zeros(n1+s1, s1);
  for k1 = 1:p1
    kk1 = (1:s1) + (k1-1)*s1;
    r1(kk1, :) = (Q1(:,kk1)'*X1);
    X1 = X1 - Q1(:,kk1)*r1(kk1,:);
  end
  %disp(norm(Q1'*X1))
  kk1 = (1:s1) + n1;
  [X1, r1(kk1,:)] = IntraOrtho(X1, LMstr);
  %disp(norm(eye(s1,s1)-X1'*X1))
  %disp(norm(X0 - Q1*r1(1:n1,:) - X1*r1(n1+1:end,:)))
end


[m, n] = size(XX);
assert(mod(n,s) == 0)
assert(mod(m,local_m) == 0)
p = n/s;
procs = m/local_m;

RR = zeros(n,n);
reduced_R = zeros(n * procs, n);

for k = 1:p
  kk = (1:s) + s*(k-1);
  % local orthogonalization
  for proc = 1:procs
    rows = (proc-1)*local_m+1 : local_m*proc;
    reduced_rows = (proc-1)*n+1:((proc-1)*n+s*k);
    [XX(rows, kk),reduced_R(reduced_rows, kk)] = bmgs_step(XX(rows, 1:s*(k-1)), XX(rows, kk));
  end
  % global orthogonalization of reduced R
  [reduced_R(:,kk), RR(1:s*(k),kk)] = bmgs_step(reduced_R(:,1:s*(k-1)),reduced_R(:,kk));
end

  % application of orthogonalized reduced R
for proc = 1:procs
  rows = (proc-1)*local_m+1 : local_m*proc;
  reduced_rows = (proc-1)*n+1:((proc-1)*n+s*k);
  XX(rows, :) = XX(rows, :)*(reduced_R(reduced_rows, :));
end
%disp(norm(eye(s,s) - XX(:,kk)'*XX(:,kk)))
%disp(norm(X00 - XX*RR))


end
