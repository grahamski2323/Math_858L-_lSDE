function [L] = l_sdp(X,Q,n,m,k) %[L,time] = l_sdp(X,Q,n,m,k)

%% getting k-NN
mu = sparse(n,n); % the k-NN indicator matrix
for ii = 1:n % all points
    dii = pdist2(X(ii,:),X,'euclidean');
   [~,sidx] = sort(dii);
    % taking k closest, aside from the point itslef
    mu(ii,sidx(2:k+1)) = 1; % k- nearest, put into row
end
%nnz(mu)
clear dii sidx

for ii = 1:n-1
    ii_nn = find(mu(ii,:) == 1);
    for jj = ii+1:n
        jj_nn = find(mu(jj,:) == 1);
        % finding all neighbors which are shared
        shared = intersect(ii_nn,jj_nn);
        % making mu(ii,jj) = 1 if there are any shared as well
        if ~isempty(shared)
            mu(ii,jj) = 1; % since they share at least one NN
        end
    end
end
%nnz(mu)

clear shared jj_nn ii_nn ii jj
vecidx = find(mu == 1);
% converting to matrix index pairs
[i,j] = ind2sub(size(mu),vecidx);
nn_con = length(i);
% should be 1346 for the non-random marks total n=200
clear mu vecidx
%% CVX SPECS - make sure to have run cvx_setup prior to this
cvx_solver sdpt3 %for this, not semidefinite fully, one negative eginval
% cvx_solver sedumi % check both, oe works sometimes for ertain data wile the other doesnt
% N = 200 only works for sdpt3, sedumi does not have semidefinite sdp and
% fails isometry 88 times
cvx_solver
cvx_precision high
clear k
%% CVX loop, non-sdp (the one we weill probably use)
cvx_begin
    %tic();
    variable L(m,m) symmetric;
    expression K(n,n)
    expression sumK % for elemntwise sum
    sumK=0;
    % QLQ^T, not sure if I need this as expression but its good practice
    K = Q*L*transpose(Q);% QLQ^T
    % objective
    maximize(trace(K)); 
    subject to
        % 1)
        L == semidefinite(m);
        % sum constraint 2)
        
        % double-sum is less memory sum sum(K(:))
        %sum(sum(K)) == 0;
        % element-wise is slower, but it'll work memory-wise
        % for large N
        for l = 1:n
            for h = 1:n
                sumK = sumK + K(l,h);
            end
        end
        sumK == 0;
        
        % all the 3) NN constraints for isometry
        for ii = 1:nn_con
            % A little more memory, not an issue at all
            K(i(ii),i(ii)) -2*K(i(ii),j(ii)) + K(j(ii),j(ii)) <= pow_abs(vecnorm(X(i(ii),:) - X(j(ii),:)),2);
            % re-writing as element-wise multiplication
            %(K(i(ii),i(ii)) -2*K(i(ii),j(ii)) + K(j(ii),j(ii))) <= ((X(i(ii),1) - X(j(ii),1))*(X(i(ii),1) - X(j(ii),1)) +(X(i(ii),2) - X(j(ii),2))*(X(i(ii),2) - X(j(ii),2)) + (X(i(ii),3) - X(j(ii),3))*(X(i(ii),3) - X(j(ii),3)));
        end
cvx_end
%%CHECK sum constraint
sum(sum(K))
%CHECK isometry
count = 0;
for ii = 1:nn_con
            t = K(i(ii),i(ii)) -2*K(i(ii),j(ii)) + K(j(ii),j(ii)) <= pow_abs(vecnorm(X(i(ii),:) - X(j(ii),:)),2);       
            count = count + t;
end
nn_con
count
% % all verified
% % sum(sum(k) us 3e-7, close enough  (it is 'failed')

%% other option for code, DID NOT USE
% %% CVX LOOP, SDP version
% % cvx_begin sdp
% %     variable L(m,m) symmetric;
%     dual variable LL;
%     dual variables z{nn_con};
%     expression K(n,n)
%     % QLQ^T, not sure if I need this as expression - CHECK
%     K = Q*L*transpose(Q);% QLQ^T
%     % objective
%     maximize(trace(K)); 
%     subject to
%         % 1)
%         %L == semidefinite(m);
%         L >= 0 : LL;
%         % sum constraint 2)
%         sum(sum(K)) == 0;
%         % all the 3) NN constraints for isometry
%         for ii = 1:length(i) % is same length as j, and they are pairs, only one for loop
%             %K(i(ii),i(ii)) -2*K(i(ii),j(ii)) + K(j(ii),j(ii)) <= pow_abs(vecnorm(X(i(ii),:) - X(j(ii),:)),2); % I HOPE THIS WORK
%             % might have to re-write the norm^2 in terms of indiviual element multiplication
%             %ii
%             % re-writing as element-wise multiplication
%             (K(i(ii),i(ii)) -2*K(i(ii),j(ii)) + K(j(ii),j(ii))) <= ((X(i(ii),1) - X(j(ii),1))*(X(i(ii),1) - X(j(ii),1)) +(X(i(ii),2) - X(j(ii),2))*(X(i(ii),2) - X(j(ii),2)) + (X(i(ii),3) - X(j(ii),3))*(X(i(ii),3) - X(j(ii),3))) : z{ii};
%         end
% cvx_end
