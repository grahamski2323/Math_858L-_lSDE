function [X,c] = swiss(n,m,noise)
rng(22);


% h - range/magnitude of Swiss roll
% omega - tightness/ number of of curls
if n < 1000 % make 500 for testing another curl0
    h = 11*rand(n,1);
    omega = 0.75;
    %omega = 1.0;
elseif n < 10000
    h = 13*rand(n,1);
    omega = 1.0;
else
    h = 15*rand(n,1);
    omega = 2.0;
end

% if we want some noise
if noise
    level = 0.25;% noise magnitude, change if desired
else
    level = 0.0;
end

t = 3*pi/2  * (1 + 2*rand(n,1));
% creating the data
X = [t.*cos(omega*t), h, t.*sin(omega*t)] + level*rand(n,3);
% rng('shuffle') % adding noise dimensions - DOES NOT WORK
% X = [t.*cos(omega*t), h, t.*sin(omega*t) level*rand(n,1) level*rand(n,1)]
% creating colors for plots
c = kmeans(X,5);

% partitioning X,c, for Q matrix, landmarks 1st
% REPEATABLE FOR SDP, FOR ONLY TESTING W,Q, uncomment next 2 lines
%rng('shuffle');
%rng(43);
% rng(55);
if n == 500
    rng(33); % for n = 500 good
    marks = randperm(n,40);
else
    count=0;
    
    for ii = 19:4:178
        c(ii) = 6;
        count = count+1;
        marks(count) = ii;
    end
    
    clear count
end

c(marks) = 6;
temp = 1:n;
temp(marks) = [];
temp = [marks temp];

X = X(temp,:);
c = c(temp);


% plotting
figure(1);
scatter3(X(:,1),X(:,2),X(:,3),[],c,'fill','MarkerEdgeColor','k');
view(-20,5);
tit = ['Original Swiss Roll: n=', num2str(n), ' m=' num2str(m)];
% figure(11);
% scatter3(X(:,1),X(:,2),X(:,4),[],c,'fill','MarkerEdgeColor','k');
% view(-20,5);
% tit = ['Original Swiss Roll w/noise dimension, n=', num2str(n)];

title(tit)

