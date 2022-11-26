% Nomenclature:
            % NorthWest
		    %
		    %   P (i,j) - - ω  - - e  - - E (i,j+1)
		    %      |                         |
		    %      |        |       |        |      |
		    %      |                         |
		    %      η ------ ηω ---- ηe - -  ηE
		    %      |                         |
		    %      |        |       |        |      |
		    %      |                         |
		    %      s ------ sω ---- se - -  sE
		    %
		    %      |        |       |        |
		    %
		    %   S(i+1,j)  - Sω  -  Se  - -  SE(i+1,j+1)
    
		    %  Indexing of stecil: 
		    %     D_0 - D3
		    %      |     | 
		    %     D1 - D4

            % SouthWest
		    %
		    %   N(i-1,j)- - ω - - -  Ne - - NE(i-1,j+1)
		    %
		    %      |        |        |        |
		    %
		    %      n ------ nω ----- ne - -  nE
		    %      |                          |
		    %      |        |        |        |
		    %      |                          |
		    %      η ------ ηω ----- ηe - -  ηE
		    %      |                          |
		    %      |        |        |        |
		    %      |                          |
		    %   P (i,j) - - ω - - -  e - - E (i,j+1)
    
		    %  Indexing of stecil: 
		    %     D_1 - D2
		    %      |     | 
		    %     D_0 - D3

            % Principal node coordinates
            y_P     = Y(i,j);       x_P     = X(i,j);
            y_E     = Y(i,j+1);     x_E     = X(i,j+1);

            % Auxiliary node coordinates
            y_e     = (y_P+y_E)/2;  x_e     = (x_P+x_E)/2;

            y_omega = (y_e+y_P)/2;  x_omega = (x_e+x_P)/2;

            % Inter node distances
            dy_e_P   = y_P - y_e;   dx_e_P   = x_P - x_e;
            dy_P_e   = y_e - y_P;   dx_P_e   = x_e - x_P;
            
            if i==1 %NW
                % Principal node coordinates
                y_S     = Y(i+1,j);     x_S     = X(i+1,j);
                y_SE    = Y(i+1,j+1);   x_SE    = X(i+1,j+1);

                % Auxiliary node coordinates
                y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
                y_sE = (y_SE + y_E)/2;  x_sE = (x_SE + x_E)/2;

                y_s  = (y_P + y_S)/2;   x_s  = (x_P + x_S)/2;
                y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;

                y_somega = (y_se + y_s)/2;  x_somega = (x_se + x_s)/2;
                y_Somega = (y_Se + y_S)/2;  x_Somega = (x_Se + x_S)/2;

                % Inter node distances

                % Around sω
                dy_P_S = y_S - y_P;     dx_P_S  = x_S - x_P;
                dy_S_Se = y_Se - y_S;   dx_S_Se = x_Se - x_S;
                dy_Se_e = y_e - y_Se;   dx_Se_e = x_e - x_Se;
                dy_e_P = y_P - y_e;     dx_e_P  = x_P - x_e;

                % Around ηe2
                dy_P_s  = y_s - y_P;    dx_P_s  = x_s - x_P;
                dy_s_sE = y_sE - y_s;   dx_s_sE = x_sE - x_s;
                dy_sE_E = y_E - y_sE;   dx_sE_E = x_E - x_sE;
                dy_E_P  = y_P - y_E;    dx_E_P  = x_P - x_E;

                % Around P
                dy_se_e = y_e - y_se;   dx_se_e = x_e - x_se;
                dy_s_se = y_se - y_s;   dx_s_se = x_se - x_s;

                % Areas
                S_somega    = abs((x_Se*y_S  - x_S*y_Se)  + (x_S*y_P   - x_P*y_S)   + (x_P*y_e  - x_e*y_P)  + (x_e*y_Se - x_Se*y_e))/2;
                S_etaomega2 = abs((x_se*y_s  - x_s*y_se)  + (x_s*y_P   - x_P*y_s)   + (x_P*y_e  - x_e*y_P)  + (x_e*y_se - x_se*y_e))/2;
                S_etae2     = abs((x_sE*y_s  - x_s*y_sE)  + (x_s*y_P   - x_P*y_s)   + (x_P*y_E  - x_E*y_P)  + (x_E*y_sE - x_sE*y_E))/2;
            
            else %SW

                % Principal node coordinates
                y_N     = Y(i-1,j);     x_N     = X(i-1,j);
                y_NE    = Y(i-1,j+1);   x_NE    = X(i-1,j+1);

                % Auxiliary node coordinates
                y_nE = (y_NE + y_E)/2;  x_nE = (x_NE + x_E)/2;
                y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;

                y_n  = (y_P + y_N)/2;   x_n  = (x_P + x_N)/2;
                y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;

                y_Nomega = (y_Ne + y_N)/2;  x_Nomega = (x_Ne + x_N)/2;
                y_nomega = (y_ne + y_n)/2;  x_nomega = (x_ne + x_n)/2;

                % Inter node distances

                % Around nω
                dy_N_P  = y_P - y_N;     dx_N_P  = x_P - x_N;
                dy_P_e  = y_e - y_P;     dx_P_e  = x_e - x_P;
                dy_e_Ne = y_Ne - y_e;    dx_e_Ne = x_Ne - x_e;
                dy_Ne_N = y_N - y_Ne;    dx_Ne_N = x_N - x_Ne;

                % Around ηe1
                dy_n_P  = y_P - y_n;    dx_n_P  = x_P - x_n;
                dy_P_E  = y_E - y_P;    dx_P_E  = x_E - x_P;
                dy_E_nE = y_nE - y_E;   dx_E_nE = x_nE - x_E;
                dy_nE_n = y_n - y_nE;   dx_nE_n = x_n - x_nE;

                % Around P
                dy_e_ne = y_ne - y_e;   dx_e_ne = x_ne - x_e;
                dy_ne_n = y_n - y_ne;   dx_ne_n = x_n - x_ne;

                % Areas
                S_nomega    = abs((x_e*y_P  - x_P*y_e)  + (x_P*y_N   - x_N*y_P)   + (x_N*y_Ne  - x_Ne*y_N)  + (x_Ne*y_e - x_e*y_Ne))/2;
                S_etaomega1 = abs((x_e*y_P  - x_P*y_e)  + (x_P*y_n   - x_n*y_P)   + (x_n*y_ne  - x_ne*y_n)  + (x_ne*y_e - x_e*y_ne))/2;
                S_etae1     = abs((x_E*y_P  - x_P*y_E)  + (x_P*y_n   - x_n*y_P)   + (x_n*y_nE  - x_nE*y_n)  + (x_nE*y_E - x_E*y_nE))/2;
            end