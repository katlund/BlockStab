classdef ReduceBlockGS  < handle
  
  properties
    QX = []
    QR = {}
    GS1
    GS2
    recursionEnd = 200
  end

  methods
    function self = ReduceBlockGS(recursionEnd)
      self.recursionEnd = recursionEnd;
    end

    function [X, R] = orthogonalize(self,X)
      [n, s] = size(X); 
      k = size(self.QX, 2)/s;
      if n < self.recursionEnd % recursion end (BMGS)
        R = zeros(k*s+s, s);
        for i = 1:k
          ii = (1:s) + (i-1)*s;
          R(ii,:) = self.QX(:,ii)'*X;
          X =  X - self.QX(:,ii)*R(ii,:);
        end
        [X, R(end-s+1:end,:)] = IntraOrtho(X, "mgs");
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
        R = zeros((k+1)*s,s);
        for i = 1:k
          ii = (1:s) + (i-1)*s;
          R(ii,:) = self.QR{i}{1}'*R1(1:(i*s),:) + self.QR{i}{2}'*R2(1:(i*s),:);
          R1(1:i*s,:) = R1(1:i*s,:) - self.QR{i}{1}*R(ii,:);
          R2(1:i*s,:) = R2(1:i*s,:) - self.QR{i}{2}*R(ii,:);
        end        
        R12 = vertcat(R1,R2);
        [R12, R(end-s+1:end,:)] = IntraOrtho(R12, "mgs");
        R1 = R12(1:(k+1)*s,:);
        R2 = R12((k+1)*s+1:end,:);

        X(1:m,:) = self.QX(1:m,:)*R1;
        X(m+1:end,:) = self.QX(m+1:end,:)*R2;
        self.QR{end+1} = {R1,R2};
      end
    end
  end
end
