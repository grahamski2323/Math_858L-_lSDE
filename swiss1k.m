function [X,c] = swiss1k(n)
h = 11*linspace(0,1,n).';
rng(20);
t =  3*pi/2 * (1 + 2*rand(n,1));

X = [t.*cos(0.75*t), h, t.*sin(0.75*t)];
%X = [t.*cos(1.0*t), h, t.*sin(1.0*t)]; % trying more folds

c = kmeans(X,5);
rng(55);
%rng(23);
marks = randperm(n,40);
c(marks) = 6;
% count = 0;
% for ii = 15:25:990
%     count = count + 1;
%     c(ii) = 6;
%     marks(count) = ii;
% end
temp = 1:n;
temp(marks) = [];
temp = [marks temp];

X = X(temp,:);
c = c(temp);


% plotting
figure(1);
scatter3(X(:,1),X(:,2),X(:,3),[],c,'fill','MarkerEdgeColor','k');
view(-20,5);
tit = ['Original Swiss Roll, n=', num2str(n)];
title(tit);