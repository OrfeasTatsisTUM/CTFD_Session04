% Nomenclature:
            % NorthEast
            %
		    %   W(i, j-1) - - w - - - - ω  - - - P (i,j)
		    %                 |         |
		    %       |         |         |        |
		    %                 |         |
		    %       ηW - - -  ηw ------ ηω ----- η
		    %                 |         |
		    %       |         |         |        |
		    %                 |         |
		    %       sW - - -  sw ------ sω ----- s
		    %
		    %                 |         |
		    %
		    %   SW(i+1,j-1) - Sw    -   Sω   -   S(i+1,j)
		    %
		    % Indexing of stecil: 
    
		    %    D_3 - D_0
		    %     |     | 
		    %    D_2 -  D1

            % SouthEast
            %
            %   NW(i-1,j-1) - Nw    -   Nω   -   N(i-1,j)
    		%
    		%                 |         |
    		%
    		%       nW - - -  nw  ----- nω ----- n
    		%                 |         |
    		%       |         |         |        |
    		%                 |         |
    		%       ηW - - -  ηw  ----- ηω ----- η
    		%                 |         |
    		%       |         |         |        |
    		%                 |         |
    		%   W(i, j-1) - - w - - - - ω  - - - P (i,j)
    		%
    		% Indexing of stecil:

    		%    D_4 - D_1
    		%     |     |
    		%    D_3 - D_0

            % Principal node coordinates
            y_P     = Y(i,j);       x_P     = X(i,j);
            y_W     = Y(i,j-1);     x_W     = X(i,j-1);

            % Auxiliary node coordinates
            y_w     = (y_P+y_W)/2;  x_w     = (x_P+x_W)/2;

            y_omega = (y_w+y_P)/2;  x_omega = (x_w+x_P)/2;

            % Inter node distances
            dy_w_P   = y_P - y_w;   dx_w_P   = x_P - x_w;

            
            if i==1 %NE
                % Principal node coordinates
                y_S     = Y(i+1,j);     x_S     = X(i+1,j);
                y_SW    = Y(i+1,j-1);   x_SW    = X(i+1,j-1);

                % Auxiliary node coordinates
                y_Sw = (y_SW + y_S)/2;  x_Sw = (x_SW + x_S)/2;
                y_sW = (y_SW + y_W)/2;  x_sW = (x_SW + x_W)/2;

                y_s  = (y_P + y_S)/2;   x_s  = (x_P + x_S)/2;
                y_sw = (y_Sw + y_w)/2;  x_sw = (x_Sw + x_w)/2;

                y_somega = (y_sw + y_s)/2;  x_somega = (x_sw + x_s)/2;
                y_Somega = (y_Sw + y_S)/2;  x_Somega = (x_Sw + x_S)/2;

                % Inter node distances

                % Around sω
                dy_S_P   = y_P - y_S;   dx_S_P   = x_P - x_S;
                dy_Sw_S  = y_S - y_Sw;  dx_Sw_S  = x_S - x_Sw;
                dy_w_Sw  = y_Sw - y_w;  dx_w_Sw  = x_Sw - x_w;
                dy_P_w   = y_w - y_P;   dx_P_w   = x_w - x_P;

                % Around ηw2
                dy_W_sW  = y_sW - y_W;  dx_W_sW  = x_sW - x_W;
                dy_sW_s  = y_s - y_sW;  dx_sW_s  = x_s - x_sW;
                dy_s_P   = y_P - y_s;   dx_s_P   = x_P - x_s;
                dy_P_W   = y_W - y_P;   dx_P_W   = x_W - x_P;

                % Around P
                dy_w_sw  = y_sw - y_w;  dx_w_sw  = x_sw - x_w;
                dy_sw_s  = y_s - y_sw;  dx_sw_s  = x_s - x_sw;

                % Areas
                S_somega    = abs((x_S*y_Sw  - x_Sw*y_S)  + (x_Sw*y_w   - x_w*y_Sw)   + (x_w*y_P  - x_P*y_w)  + (x_P*y_S - x_S*y_P))/2;
                S_etaomega2 = abs((x_s*y_sw  - x_sw*y_s)  + (x_sw*y_w   - x_w*y_sw)   + (x_w*y_P  - x_P*y_w)  + (x_P*y_s - x_s*y_P))/2;
                S_etaw2     = abs((x_s*y_sW  - x_sW*y_s)  + (x_sW*y_W   - x_W*y_sW)   + (x_W*y_P  - x_P*y_W)  + (x_P*y_s - x_s*y_P))/2;
            
            else %SE

                % Principal node coordinates
                y_N     = Y(i-1,j);     x_N     = X(i-1,j);
                y_NW    = Y(i-1,j-1);   x_NW    = X(i-1,j-1);

                % Auxiliary node coordinates
                y_nW = (y_NW + y_W)/2;  x_nW = (x_NW + x_W)/2;
                y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2;

                y_n  = (y_P + y_N)/2;   x_n  = (x_P + x_N)/2;
                y_nw = (y_Nw + y_w)/2;  x_nw = (x_Nw + x_w)/2;

                y_Nomega = (y_Nw + y_N)/2;  x_Nomega = (x_Nw + x_N)/2;
                y_nomega = (y_nw + y_n)/2;  x_nomega = (x_nw + x_n)/2;

                % Inter node distances

                % Around nω
                dy_P_N   = y_N - y_P;   dx_P_N   = x_N - x_P;
                dy_w_P   = y_P - y_w;   dx_w_P   = x_P - x_w;
                dy_Nw_w  = y_w - y_Nw;  dx_Nw_w  = x_w - x_Nw;
                dy_N_Nw  = y_Nw - y_N;  dx_N_Nw  = x_Nw - x_N;

                % Around ηw1
                dy_nW_W  = y_W - y_nW;  dx_nW_W  = x_W - x_nW;
                dy_n_nW  = y_nW - y_n;  dx_n_nW  = x_nW - x_n;
                dy_P_n   = y_n - y_P;   dx_P_n   = x_n - x_P;
                dy_W_P   = y_P - y_W;   dx_W_P   = x_P - x_W;

                % Around P
                dy_nw_w  = y_w - y_nw;  dx_nw_w  = x_w - x_nw;
                dy_n_nw  = y_nw - y_n;  dx_n_nw  = x_nw - x_n; 

                % Areas
                S_nomega    = abs((x_P*y_w  - x_w*y_P)  + (x_w*y_Nw   - x_Nw*y_w)   + (x_Nw*y_N  - x_N*y_Nw)  + (x_N*y_P - x_P*y_N))/2;
                S_etaomega1 = abs((x_P*y_w  - x_w*y_P)  + (x_w*y_nw   - x_nw*y_w)   + (x_nw*y_n  - x_n*y_nw)  + (x_n*y_P - x_P*y_n))/2;
                S_etaw1     = abs((x_P*y_W  - x_W*y_P)  + (x_W*y_nW   - x_nW*y_W)   + (x_nW*y_n  - x_n*y_nW)  + (x_n*y_P - x_P*y_n))/2;
            end