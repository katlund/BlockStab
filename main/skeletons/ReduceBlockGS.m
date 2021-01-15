classdef ReduceBlockGS  < handle

  properties
    QX = []
    QR = {}
    GS1
    GS2
    recursionEnd = 200
    gs_red_it = 2 % iterations for reduction gram-schmidt method
    gs_loc_it = 2% iterations for local gram-schmidt method
  end

  methods
    function self = ReduceBlockGS(recursionEnd, gs_red_it, gs_loc_it)
      self.recursionEnd = recursionEnd;
      if nargin >= 2
        self.gs_red_it = gs_red_it
      end
      if nargin >= 2
        self.gs_loc_it = gs_loc_it
      end
    end

    function [X, R] = orthogonalize(self,X)
      [n, s] = size(X);
      k = size(self.QX, 2)/s;
      if n < self.recursionEnd % recursion end (BMGS)
        R = vertcat(zeros(k*s, s),eye(s,s));
        for j = 1:self.gs_loc_it
          for i = 1:k
            ii = (1:s) + (i-1)*s;
            Rtmp = self.QX(:,ii)'*X;
            R(ii,:) = R(ii,:) + Rtmp*R(end-s+1:end,end-s+1:end);
            X =  X - self.QX(:,ii)*Rtmp;
          end
          [X, Rtmp] = IntraOrtho(X, "mgs");
          R(end-s+1:end,:) = Rtmp*R(end-s+1:end,:);
        end
        self.QX = horzcat(self.QX,X);
      else
        % initialization
        if isempty(self.GS1)
          self.GS1 = ReduceBlockGS(self.recursionEnd);
          self.GS2 = ReduceBlockGS(self.recursionEnd);
        end

        % divide
        m = floor(n/2);
        [X(1:m,:), R1] = self.GS1.orthogonalize(X(1:m,:));
        [X(m+1:end, :), R2] = self.GS2.orthogonalize(X(m+1:end,:));
        self.QX = horzcat(self.QX,X);

        % and conquer
        R = vertcat(zeros(k*s,s), eye(s,s));
        for j = 1:self.gs_red_it
          for i = 1:k
            ii = (1:s) + (i-1)*s;
            Rtmp = (self.QR{i}{1}'*R1(1:(i*s),:) + self.QR{i}{2}'*R2(1:(i*s),:));
            R(ii,:) = R(ii,:) + Rtmp*R(end-s+1:end,end-s+1:end);
            R1(1:i*s,:) = R1(1:i*s,:) - self.QR{i}{1}*Rtmp;
            R2(1:i*s,:) = R2(1:i*s,:) - self.QR{i}{2}*Rtmp;
          end
          R12 = vertcat(R1,R2);
          [R12, Rtmp] = IntraOrtho(R12, "mgs");
          R(end-s+1:end,:) = Rtmp*R(end-s+1:end,:);
          R1 = R12(1:(k+1)*s,:);
          R2 = R12((k+1)*s+1:end,:);
        end

        X(1:m,:) = self.QX(1:m,:)*R1;
        X(m+1:end,:) = self.QX(m+1:end,:)*R2;
        self.QR{end+1} = {R1,R2};
      end
    end
  end
end
