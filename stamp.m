function [stencil] = stamp(i, j, X, Y, lamda, alpha, Tinf, boundary,TD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stencil calculates the linear equation for node (i, j)

%  Input:
%      i         node number in x direction
%      j         node number in y direction
%      X         x position of the nodes
%      Y         y position of the nodes
%      b         right-hand side value for node (i,j)
%      alpha     alpha
%      Tinf      Tinf for Robin BC
%      boundary  defines the boundary conditions
%      verbose   verbosity level
%
%  Output:
%      stencil   linear equation for node (i,j)
%      b         new right-hand side value for node (i,j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init

n = size(X, 1);     %dimY
m = size(X, 2);     %dimX
stencil = zeros(1, n*m);
index=@(ii, jj) ii + (jj-1)*n;

% Determine the node positon

if (i==1 || i==n) && j==m
    nodePosition = 'EastCorner';
elseif (i==1 || i==n) && j==1
    nodePosition = 'WestCorner';
elseif j == 1
    nodePosition = 'West';
elseif j == m
    nodePosition = 'East';
elseif i== n
    nodePosition = 'South';
elseif i == 1
    nodePosition = 'North';
else
    nodePosition = 'inner Node';
end
% WestCorner = (i==1 || i==n) && j==1;

% Calculate the equation for the correct node position
switch nodePosition
    %% Inner
    case 'inner Node'

        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        data_inner
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        %$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$
        build_inner
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        % P
        stencil(index(i, j))     = lamda(i,j)      * D0;

        % East
        stencil(index(i, j+1))   = lamda(i,j+1)    * D3;

        % West
        stencil(index(i, j-1))   = lamda(i,j-1)    * D_3;

        % South
        stencil(index(i+1, j))   = lamda(i+1,j)    * D1;

        % North
        stencil(index(i-1, j))   = lamda(i-1,j)    * D_1;

        % NW
        stencil(index(i-1, j-1)) = lamda(i-1,j-1)  * D_4;

        % NE
        stencil(index(i-1, j+1)) = lamda(i-1,j+1)  * D2;

        % SW
        stencil(index(i+1, j-1)) = lamda(i+1, j-1) * D_2;

        % SE
        stencil(index(i+1, j+1)) = lamda(i+1, j+1) * D4;

        %% South
    case 'South'

        if strcmp(boundary.south, 'Dirichlet')
            stencil(index(i, j))     = 1;
        else

            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            data_south
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            bc_ctrl = strcmp(boundary.south, 'Robin'); % factor that includes T_P in 3.16 (A.14)

            %$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$
            build_south
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            % P
            stencil(index(i, j))     = lamda(i,j)      * D0;

            % East
            stencil(index(i, j+1))   = lamda(i,j+1)    * D3;

            % West
            stencil(index(i, j-1))   = lamda(i,j-1)    * D_3;

            % North
            stencil(index(i-1, j))   = lamda(i-1,j)    * D_1;

            % NW
            stencil(index(i-1, j-1)) = lamda(i-1,j-1)  * D_4;

            % NE
            stencil(index(i-1, j+1)) = lamda(i-1,j+1)  * D2;
        end

        %% North
    case 'North'

        if strcmp(boundary.north, 'Dirichlet')
            stencil(index(i, j))     = 1;
        else

            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            data_north
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            bc_ctrl = strcmp(boundary.north, 'Robin'); % factor that includes T_P in 3.16 (A.14)

            %$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$
            build_north
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            % P
            stencil(index(i, j))     = lamda(i,j)      * D0;

            % East
            stencil(index(i, j+1))   = lamda(i,j+1)    * D3;

            % West
            stencil(index(i, j-1))   = lamda(i,j-1)    * D_3;

            % South
            stencil(index(i+1, j))   = lamda(i+1,j)    * D1;

            % SW
            stencil(index(i+1, j-1)) = lamda(i+1, j-1) * D_2;

            % SE
            stencil(index(i+1, j+1)) = lamda(i+1, j+1) * D4;
        end

        %% East
    case 'East'
        if strcmp(boundary.east, 'Dirichlet')
            stencil(index(i, j))     = 1;
        else
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            data_east
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            bc_ctrl = strcmp(boundary.east, 'Robin'); % factor that includes T_P in 3.16 (A.14)

            %$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$
            build_east
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            % P
            stencil(index(i, j))     = lamda(i,j)      * D0;

            % West
            stencil(index(i, j-1))   = lamda(i,j-1)    * D_3;

            % North
            stencil(index(i-1, j))   = lamda(i-1,j)    * D_1;

            % NW
            stencil(index(i-1, j-1)) = lamda(i-1, j-1) * D_4;

            % South
            stencil(index(i+1, j))   = lamda(i+1,j)    * D1;

            % SW
            stencil(index(i+1, j-1)) = lamda(i+1, j-1) * D_2;
        end

        %% West
    case 'West'
        if strcmp(boundary.west, 'Dirichlet')
            stencil(index(i, j))     = 1;
        end

        %% Corners
    case 'EastCorner'
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        data_eastcorner
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        bc_ctrl_east = strcmp(boundary.east, 'Robin'); % factor that includes T_P in 3.16
        bc_ctrl_n = (i==1)*strcmp(boundary.north, 'Robin');
        bc_ctrl_s = (i==n)*strcmp(boundary.south, 'Robin');

        %$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$
        build_eastcorner
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        % P
        stencil(index(i, j))     = lamda(i,j)      * D0;

        % West
        stencil(index(i, j-1))   = lamda(i,j-1)    * D_3;

        if i ~= 1
            % North
            stencil(index(i-1, j))   = lamda(i-1,j)    * D_1;

            % NW
            stencil(index(i-1, j-1)) = lamda(i-1, j-1) * D_4;
        end

        if i ~= n
            % South
            stencil(index(i+1, j))   = lamda(i+1,j)    * D1;

            % SW
            stencil(index(i+1, j-1)) = lamda(i+1, j-1) * D_2;
        end

    case 'WestCorner'
        if strcmp(boundary.west, 'Dirichlet')
            % P
            stencil(index(i, j)) =1;
        else
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            data_westcorner
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            bc_ctrl_west = strcmp(boundary.west, 'Robin'); % factor that includes T_P in [3.16]
            bc_ctrl_n = (i==1)*strcmp(boundary.north, 'Robin');
            bc_ctrl_s = (i==n)*strcmp(boundary.south, 'Robin');

            %$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$
            build_westcorner
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            % P
            stencil(index(i, j))     = lamda(i,j)      * D0;

            % East
            stencil(index(i, j+1))   = lamda(i,j+1)    * D3;

            if i ~= 1
                % North
                stencil(index(i-1, j))   = lamda(i-1,j)    * D_1;

                % NE
                stencil(index(i-1, j+1)) = lamda(i-1,j+1)  * D2;
            end

            if i ~= n
                % South
                stencil(index(i+1, j))   = lamda(i+1,j)    * D1;

                % SE
                stencil(index(i+1, j+1)) = lamda(i+1, j+1) * D4;
            end
        end
end

