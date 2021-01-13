function [XX,RR] = breducegs(XX, s, LMstr)
% [XX,RR] = breducegs(XX, s, LGSstr, ReduceGSstr, verbose) performs block reduce 
% Gram-Schmidt algorithm
% TODO: Document parameters

addpath(genpath('../'))
addpath(genpath('../main/'))
addpath(genpath('../skeletons/'))

[m, n] = size(XX);
ortho = ReduceBlockGS(2*n);

for k = 0:s:n-1
  kk = (1:s) + k;
  [XX(:,kk), RR(1:(k+s),kk)] = ortho.orthogonalize(XX(:,kk));
end

end