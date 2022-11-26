figure(1)
pcolor(X,Y,T)
hold on
pcolor(X,-Y,T)
if strcmp(shape, 'linear')
    axis([0 l -h1/2 h1/2 -1 1]) %limits of x & y axes
elseif strcmp(shape, 'quadratic')
    axis([0 l -h1/2*1.15 h1/2*1.15 -1 1]) %limits of x & y axes
else
    axis([0 l -h1/2*1.45 h1/2*1.45 -1 1]) %limits of x & y axes
end
set(gcf, 'Position',[70,150,550,550]);

figure(2)
contour(X,Y,T)
hold on
contour(X,-Y,T)
set(gcf, 'Position',[635,150,550,550]);

figure(3)
surf(X,Y,T); hold on; surf(X,-Y,T)
set(gcf, 'Position',[1200,150,550,550]);
view(45,24);            %view angle