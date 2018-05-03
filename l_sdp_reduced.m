function [L] = l_sdp_reduced(X,Q,n,m,k) %[L,time] = l_sdp_reduced(X,Q,n,m,k)
% attempting to look only at isometry constraints of NN first, then
% checking which ones are satisfied, adding, and re-running
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
clear mu
%% finding i,j where i,j <= m (first m are landmarks)
im = find(i <= m);
jm = find(j <= m);

% starting the reduced isomtery constraints
idx_red = union(im,jm); % to not repeat if landmarks are NN of each other
ired = i(idx_red);
jred = j(idx_red);
%[ired jred] % to check they all are landmark NN
nn_red = length(ired);
% maybe keep idx_red? don't think so
clear im jm idx_red
% so 477 landmarks NN
%% CVX SPECS - make sure to have run cvx_setup prior to this
cvx_solver sdpt3
% cvx_solver sedumi % check both, oe works sometimes for ertain data wile the other doesnt
cvx_solver
cvx_precision high
clear k
%% CVX loop, non-sdp (the one we weill probably use)
isometry = 1;

while isometry % keep going until all isometry constraints are satisfied
    cvx_begin
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
            % double-sum is less memory, sum(K(:)), more than elemnt-wise
%             sum(sum(K)) == 0;
            % element-wise is slower, but it'll work memory-wise
            %for large N
                    for l = 1:n
                        for h = 1:n
                            sumK = sumK + K(l,h);
                        end
                    end
                    sumK == 0;

            % REDUCED 3)
            for ii = 1:nn_red
                % A little more memory, not an issue at all
                K(ired(ii),ired(ii)) -2*K(ired(ii),jred(ii)) + K(jred(ii),jred(ii)) <= pow_abs(vecnorm(X(ired(ii),:) - X(jred(ii),:)),2);
                % re-writing as element-wise multiplication
                %(K(ired(ii),ired(ii)) -2*K(ired(ii),jred(ii)) + K(jred(ii),jred(ii))) <= ((X(ired(ii),1) - X(jred(ii),1))*(X(ired(ii),1) - X(jred(ii),1)) +(X(ired(ii),2) - X(jred(ii),2))*(X(ired(ii),2) - X(jred(ii),2)) + (X(ired(ii),3) - X(jred(ii),3))*(X(ired(ii),3) - X(jred(ii),3)));
            end
    cvx_end
 
    %%CHECK isometry, add those who fail
    count = 0;
    for ii = 1:nn_con
        t = K(i(ii),i(ii)) -2*K(i(ii),j(ii)) + K(j(ii),j(ii)) <= pow_abs(vecnorm(X(i(ii),:) - X(j(ii),:)),2);
        if t == 0 % not satisfied
            isometry = 1;
            if ~ismember([i(ii),j(ii)],[ired jred],'rows') % if it is not in ired, jred yet
                ired = [ired; i(ii)];
                jred = [jred; j(ii)];
            end
        else
            count = count + 1;
        end
        
    end
    % new length of reduced constraints, if added, if not, it stays same
    nn_red = length(ired);
    % switch if all isometry satisfied
    if count == nn_con
        isometry = 0;
        disp(['**ALL constraints satisfied. total used:', num2str(nn_red)])
    else
        disp(['new # of isometry constraints:', num2str(nn_red)])
    end
end
    %%CHECK sum constraint
     sum(sum(K))
% all verified
% sum(sum(k) us 4.5e-6, close enough  (it is 'failed')
