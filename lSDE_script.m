% you can run individual sections
% or the whole code
% helps if we want to track run-times for later
% or if we want to skip steps sine they take too long

% REMEMBER TO RUN cvx_setup prior to this
%% loading previous data if we wnat
% since W and Q calcs take a while for large n
% load 
%% initial parameters
clear all;
n = 200; % # data points, already works for 200
m = 40; % # landmarks
r = 12; % # NN for W
k = 4; %4; % # NN for L
noise = false;

%% generating  Swiss data
[X,c] = swiss(n,m,noise);
%[X,c] = swiss1k(1000); % for n = 100 0 run
% save('X.mat','X')
% save('c.mat','c')
%% generating weight matrix W
%[W] = weights_twoloops(X,n,r);
[W] = weights_oneloop(X,n,r);
% save('W.mat','W')
clear r

%% generating Q matrix
[Q] = lintrans_Q(X,W,n,m,c);
% save('Q.mat,'Q')

%% solving for L matrix
[L] = l_sdp(X,Q,n,m,k);
% pause(10);
[L] = l_sdp_reduced(X,Q,n,m,k);

%% plotting low-dimesnional embedding
[evec,lambda] = eig(L);

% check PSD
check = diag(lambda) > 0;
sum(check)
[R,p] = chol(L);
p;
% ans = m, p = 0 --> PSD! we're good

% low-dim landmarks
lam1 = lambda(m,m);
lam2 = lambda(m-1,m-1);
lam3 = lambda(m-2,m-2); % to check third not being significant
evec1 = evec(:,m);
evec2 = evec(:,m-1);
lowmarks = [sqrt(lam1)*evec1 sqrt(lam2)*evec2];

% low-dim output
Y = Q*lowmarks;

% plotting
figure(3);
scatter(Y(:,1),Y(:,2),[],c,'filled');
tit = ['lSDE Swiss Roll 2-D Embedding: n=', num2str(n),', m=' num2str(m), ', k=', num2str(k)];
title(tit)
