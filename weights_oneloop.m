function [W] = weights_oneloop(X,n,r)

%% making the r-NN
% or really, only the non-neighbors
W = sparse(n,n);
C = transpose(X);
rnot = sparse(n,n-r);
tic();
for i = 1:n
   dii = pdist2(X(i,:),X,'euclidean');
   % outputs sorted indices for ii point, closest to furthest
   [~,ridx] = sort(dii);
   % ii is not a neighbor of itself
   rnot(i,1) = i; 
   % others (rneigh = 2:r+1)
   rnot(i,2:end) = ridx(r+2:end);
%% using lsqlin to make W
    % solves min(0.5||Cx - d||_2^2
    % C = transpose of X, d = X_i
    % x = W_i, ith column of W
    d = transpose(X(i,:));
    
    % constraint Aeq*x = beq
    % Aeq is (N-r+1) x N)
    % first N-r rows are indicator vectors of not neigbors
    % last row is ones vector
    Aeq = sparse(n-r+1,n);
    % first rows ensure weight of not-NN is zero
    Aeq((1:n-r)+(n-r+1)*(rnot(i,:)-1))=1;
    % last ensures sum of W_i to one
    Aeq(n-r+1,:) = 1;
    beq = sparse(n-r+1,1);
    beq(n-r+1) = 1;
    
    x = lsqlin(C,d,[],[],Aeq,beq,[],[]);

    % adding column
    W(i,:) = sparse(x);
end
toc();