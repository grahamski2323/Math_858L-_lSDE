function [Q] = lintrans_Q(X,W,n,m,c)
% making Phi
Phi = transpose(eye(n) - W)*(eye(n) - W);

% making Q, X parittioned -> so is W -> so is Phi
Q = zeros(n,m);
Q(1:m,:) = eye(m);
% -Phi_uu^(-1)*Phi_ul
Q(m+1:end,:) = -Phi(m+1:end,m+1:end)\Phi(m+1:end,1:m);

% making reconstruction
Xcon = Q*X(1:m,:);

%% plotting
tit = ['Reconstructed Swiss Roll: n=', num2str(n),', m=', num2str(m)];
figure(2);
scatter3(Xcon(:,1),Xcon(:,2),Xcon(:,3),[],c,'fill','MarkerEdgeColor','k');
title(tit);
view(-20,5);
% tit = ['Reconstructed Swiss Roll w/noise dimension, n=', num2str(n),', m=', num2str(m)];
% figure(12);
% scatter3(Xcon(:,1),Xcon(:,2),Xcon(:,4),[],c,'fill','MarkerEdgeColor','k');
% title(tit);
% view(-20,5);

% printing norm of X - Xcon
sprintf(['L2 norm of (Original - Reconstruction) is ', num2str(norm(X-Xcon)), ' for n=', num2str(n), ' and m=', num2str(m)])