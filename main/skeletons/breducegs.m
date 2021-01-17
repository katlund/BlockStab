function [XX,RR] = breducegs(XX, s, IOstr)
%%% [XX,RR] = breducegs(XX, s) performs the recursive Gram-Schmidt algorithm on
%%% the m x n matrix XX with p=n/s block partitions of size s with inner
%%% orthogonalization precudure determined by IOstr. Both orthogonalizations,
%%% the local and reduction orthogonalization use reiterated CGS with 2
%%% iterations. The lower bound for the recursion end is determined by the
%%% number of columns n (maximal recursion depth).

  addpath(genpath('../'))
  addpath(genpath('../main/'))
  addpath(genpath('../skeletons/'))

  [m, n] = size(XX)
  Q = [];
  para.local_gs_it = 2;
  para.local_mus = IOstr;
  para.reduction_gs_it = 2;
  para.reduction_mus = IOstr;
  para.recursionFactor = 2;
  para.recursionThreshold = n;
  RR = zeros(n,n);

  for k = 0:s:n-1
    ii = (1:s) + k;
    [Q, RR(1:k+s,ii)] = recursivegs_step(Q, XX(:,ii), para);
  end

  XX = recursive2global(Q);

  function [X, r] = cgs_step(Q, X, it, mus)
    %% [X, r] = cgs_step Adds the block vector X to the orthonormal basis Q by
    %% appling CGS using the inner orthogonalization procedure determined by mus.
    if nargin < 3
      it = 2;
    end
    if nargin < 4
      mus = 'mgs'
    end
    k = size(Q,2);
    [m,s] = size(X);
    r = zeros(k, s);
    rd = eye(s,s);
    for i = 1:it
      if k > 0
        r1 = InnerProd(Q,X,mus);
        X = X - Q*r1;
        r = r + r1*rd;
      end
      [X, rd1] = IntraOrtho(X, mus);
      rd = rd1*rd;
    end
    r = vertcat(r,rd);
    X = horzcat(Q,X);
  end

  function [Q, r] = recursivegs_step(Q, X, para)
    %% [Q,r,data] = recursivegs_step Q is a recursive basis representation,
    %% i.e. it constists of local orthonormal bases (maybe recursive again) and
    %% factors to combine this local basis to a global orthonormal basis. The
    %% struct para determines the behavior of the recursion and
    %% orthogonalization procedures.
    [m,s] = size(X);
    rf = para.recursionFactor;
    if m/rf <= para.recursionThreshold
      [Q,r] = cgs_step(Q,X, para.local_gs_it, para.local_mus);
    else
      if ~isfield(Q,"Rbasis")
        Q.Rbasis = zeros(rf*s,0);
        Q.Qbasis  = cell(rf, 1);
      end
      k = size(Q.Rbasis,2);
      mr = floor(m/rf);
      local_r = {};
      for i = 1:rf-1
        [Q.Qbasis{i}, local_r{i}] =
        recursivegs_step(Q.Qbasis{i}, X(1+(i-1)*mr:i*mr,:), para);
      end
      [Q.Qbasis{rf}, local_r{rf}] =
      recursivegs_step(Q.Qbasis{rf}, X(1+(rf-1)*mr:end,:), para);
      Rstack = cell2mat(local_r');
      if k > 0
        Rbasis_stretched = zeros(0,k);
        for i = 1:rf
          Rbasis_stretched =
          vertcat(Rbasis_stretched,Q.Rbasis(1+(i-1)*(k):i*k,:),
                  zeros(s, k));
        end
        Q.Rbasis = Rbasis_stretched;
      end
      [Q.Rbasis, r] = cgs_step(Q.Rbasis, Rstack, para.reduction_gs_it,
                               para.reduction_mus);
    end
  end

  function [Qg] = recursive2global(Q)
    %% [Qg] = recursive2global(Q) Convertes the recursive basis Q to a globally
    %% orthonormal basis Qg which is a m x n matrix.
    if ~isstruct(Q)
      Qg = Q;
    else
      k = length(Q.Qbasis);
      m = size(Q.Rbasis,2);
      Qgcell = cell();
      for i = 1:k
        Ql = recursive2global(Q.Qbasis{i});
        Qgcell{end+1} = Ql*Q.Rbasis(1+(i-1)*m:i*m,:);
      end
      Qg = cell2mat(Qgcell');
    end
  end

end
