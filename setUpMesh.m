function [X,Y] = setUpMesh(dimY, dimX, l, formfunction)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 1/(dimX-1);
dy = 1/(dimY-1);
[X,Y] = meshgrid(0:dx:1,1:-dy:0);
for  i=1:size(X,2)
    Y(:,i) = Y(:,i) * formfunction(X(1,i));
end
X = X.*l;
% pcolor(X,Y,zeros(dimY,dimX))
end