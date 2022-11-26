function T = solveFVM(dimY, dimX, X, Y, boundary, TD, lamda, alpha, Tinf, dt, tend, timeintegrationType, theta, simulationType, T)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File solveFVM.m
%
% This routine set up the linear system and solve it
%
% input
% T         Spatial Matrix T
% X         Matrix x coordinates 
% Y         Matrix y coordinates
% boundary  String vector. Boundary types.
% TD        Temperature for each boundary (if Dirichlet)
% alpha     Convective heat transfer coefficient
% Tinf      Temperature of the surrouding fluid 
% dt        Timestep
% tend      Time end
% theta     Theta scheme parameter
%
% output
% T         Temperature field (matrix) for steady case
%           Temperature field at every timestep for transient case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Index maps the node position to the correct linear equation

index = @(ii, jj) ii + (jj-1) * dimY;


% B is the right-hand side of the linear system
B = zeros(1,dimY*dimX);

%% Set boundary conditions
% in case of Neumann BC we have qdot=0 so there is nothing to be added to B

% North
if strcmp(boundary.north, 'Dirichlet')
    for i = 1:dimX
        if i == 1 && strcmp(boundary.west, 'Dirichlet')
            B(index(1,i)) = B(index(1,i)) + TD.north/2;
        elseif i == dimX && strcmp(boundary.east, 'Dirichlet')
            B(index(1,i)) = B(index(1,i)) + TD.north/2;
        else
            B(index(1,i)) = TD.north;
        end
    end
end

% South
if strcmp(boundary.south, 'Dirichlet')
    for i = 1:dimX
        if i== 1 && strcmp(boundary.west, 'Dirichlet')
            B(index(dimY,i)) = B(index(dimY,i)) + TD.south/2;
        elseif i == dimX && strcmp(boundary.east, 'Dirichlet')
            B(index(dimY,i)) = B(index(dimY,i)) + TD.south/2;
        else
            B(index(dimY,i)) = TD.south;
        end
    end
end

% East
if strcmp(boundary.east, 'Dirichlet')
    for i = 1:dimY
        if i == 1 && strcmp(boundary.north, 'Dirichlet')
            B(index(i,dimX)) = B(index(i,dimX)) + TD.east/2;
        elseif i == dimY && strcmp(boundary.south, 'Dirichlet')
            B(index(i,dimX)) = B(index(i,dimX)) + TD.east/2;
        else
            B(index(i,dimX)) = TD.east;
        end
    end
end

% West
if strcmp(boundary.west, 'Dirichlet')
    for i = 1:dimY
        if i == 1 && strcmp(boundary.north, 'Dirichlet')
            B(index(i,1)) = B(index(i,1)) + TD.west/2;
        elseif i == dimY && strcmp(boundary.south, 'Dirichlet')
            B(index(i,1)) = B(index(i,1)) + TD.west/2;
        else
            B(index(i,1)) = TD.west;
        end
    end
end

%% Steady Case
% Set up the system matrix A
A = zeros(dimY*dimX);

for i = 1:dimY
    for j = 1:dimX
        % Fill the system matrix and the right-hand side for node (i,j)
        [A(index(i,j), :)] =  stamp(i, j, X, Y, lamda, alpha, Tinf, boundary,TD);
    end
end


% Solve the linear system
T1(:) = A \ B';

% Convert solution vector into matrix
T = zeros(dimY,dimX);
for j = 1:dimX
    for i = 1:dimY
        T(i,j) = T(i,j) + T1(index(i,j));
    end
end
