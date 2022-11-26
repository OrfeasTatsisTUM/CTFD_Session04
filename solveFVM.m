function T = solveFVM(dimY, dimX, X, Y, boundary, TD, lamda, alpha, Tinf, dt, tend, TimeIntegrType, theta, simulationType)


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
T(:,:,1) = reshape(T1, [dimY,dimX]);

%% Transient case
if strcmp(simulationType, 'unsteady')

    T(:,:,1) = reshape(B, [dimY,dimX]); %initial value

    %Explicit & Implicit are subcases of Theta Scheme
    if strcmp(TimeIntegrType, 'Explicit')
        theta = 0;
    elseif strcmp(TimeIntegrType, 'Implicit')
        theta = 1;
    end

    % Theta Scheme
    if ~strcmp(TimeIntegrType, 'RungeKutta4')

        Astar = eye(dimX*dimY) - A*dt*theta;                            % A*
        for t = 1:(tend/dt-1)
            Tvec = reshape(T(:,:,t), [dimX*dimY, 1]);                   % From matrix to vector (initial)
            Bstar = (eye(dimX*dimY) + dt*(1-theta)*A) * Tvec - dt*B';   % B*
            Ttr_vec(:) = Astar \ Bstar;                                 % Linear computarion: A* T = B*
            T(:,:,t+1) = reshape(Ttr_vec, [dimY,dimX]);                 % From vector to matrix (output)
        end

    % Runge Kutta 4 Scheme
    else
        
        for t = 1:(tend/dt-1)
            Tvec = reshape(T(:,:,t), [dimX*dimY, 1]);                   % From matrix to vector (initial)
            Ti   = Tvec + dt*(A*Tvec - B')/2;                            % Predictor for t = (n + 1/2)Δt
            Tii  = Tvec + dt*(A*Ti   - B')/2;                            % Corrector for t = (n + 1/2)Δt
            Tiii = Tvec + dt*(A*Tii  - B');                              % Predictor for t = (n + 1)Δt
            Ttr_vec(:) = Tvec + dt*(A*Tvec + 2*A*Ti + 2*A*Tii ...
                + A*Tiii)/6 - dt*B';                                     % Corrector for t = (n + 1)Δt
            T(:,:,t+1) = reshape(Ttr_vec, [dimY,dimX]);                 % From vector to matrix (output)
        end

    end
end