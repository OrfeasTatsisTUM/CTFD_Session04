
            % Nomenclature:
            %
            %   NW(i-1,j-1) - Nw -  N(i-1,j) -  Ne   -   NE(i-1,j+1)
            %
            %                 |                 |
            %
            %      nW  - - - nw ------ n ------ ne - - - nE
            %                 |                 |
            %       |         |        |        |        |
            %                 |                 |
            %      ηW  - - - ηw ------ η ------ ηe - - - ηΕ
            %                 |                 |
            %       |         |        |        |        |
            %                 |                 |
            %   W(i, j-1) - - w - - P (i,j) - - e - - -  E (i,j+1)
            %
            % Indexing of stecil:

            %    D_4 - D_1 - D2
            %     |     |     |
            %    D_3 - D_0 - D3


        % Principal node coordinates
            y_W     = Y(i,j-1);     x_W     = X(i,j-1);
            y_P     = Y(i,j);       x_P     = X(i,j);
            y_E     = Y(i,j+1);     x_E     = X(i,j+1);
            y_NW    = Y(i-1,j-1);   x_NW    = X(i-1,j-1);
            y_N     = Y(i-1,j);     x_N     = X(i-1,j);
            y_NE    = Y(i-1,j+1);   x_NE    = X(i-1,j+1);

        % Auxiliary node coordinates
            y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2;
            y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;
            y_nE = (y_NE + y_E)/2;  x_nE = (x_NE + x_E)/2;
            y_nW = (y_NW + y_W)/2;  x_nW = (x_NW + x_W)/2;

            y_e  = (y_P + y_E)/2;   x_e  = (x_P + x_E)/2;
            y_n  = (y_P + y_N)/2;   x_n  = (x_P + x_N)/2;
            y_w  = (y_P + y_W)/2;   x_w  = (x_P + x_W)/2;
            y_nw = (y_Nw + y_w)/2;  x_nw = (x_Nw + x_w)/2;
            y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;

            y_etaW = (y_nW + y_W)/2;  x_etaW = (x_nW + x_W)/2;
            y_etaw = (y_nw + y_w)/2;  x_etaw = (x_nw + x_w)/2;
            y_eta  = (y_n  + y_P)/2;  x_eta  = (x_n  + x_P)/2;
            y_etae = (y_ne + y_e)/2;  x_etae = (x_ne + x_e)/2;
            y_etaE = (y_nE + y_E)/2;  x_etaE = (x_nE + x_E)/2;

        % Inter node distances

            % Around ηe
            dy_E_nE  = y_nE - y_E;    dx_E_nE  = x_nE - x_E;
            dy_P_E   = y_E  - y_P;    dx_P_E   = x_E  - x_P;
            dy_n_P   = y_P  - y_n;    dx_n_P   = x_P  - x_n;
            dy_nE_n  = y_n  - y_nE;   dx_nE_n  = x_n  - x_nE;

            % Around n    
            dy_w_e   = y_e - y_w;     dx_w_e   = x_e - x_w;
            dy_e_Ne  = y_Ne - y_e;    dx_e_Ne  = x_Ne - x_e;
            dy_Ne_Nw = y_Nw - y_Ne;   dx_Ne_Nw = x_Nw - x_Ne;
            dy_Nw_w  = y_w - y_Nw;    dx_Nw_w  = x_w - x_Nw;

            % Around ηw
            dy_P_n   = y_n - y_P;     dx_P_n   = x_n - x_P;
            dy_W_P   = y_P  - y_W;    dx_W_P   = x_P  - x_W;
            dy_nW_W  = y_W  - y_nW;   dx_nW_W  = x_W  - x_nW;
            dy_n_nW  = y_nW  - y_n;   dx_n_nW  = x_nW  - x_n;

            % Around P (w - e is common edge with control volume around n)
            dy_e_ne  = y_ne - y_e;    dx_e_ne  = x_ne - x_e;
            dy_ne_nw = y_nw - y_ne;   dx_ne_nw = x_nw - x_ne;
            dy_nw_w  = y_w - y_nw;   dx_nw_w   = x_w - x_nw;

        % Areas
            S_eta  = abs((x_ne*y_e  - x_e*y_ne)  + (x_e*y_w   - x_w*y_e)   + (x_w*y_nw  - x_nw*y_w)  + (x_nw*y_ne - x_ne*y_nw))/2;
            S_etae = abs((x_nE*y_E  - x_E*y_nE)  + (x_E*y_P   - x_P*y_E)   + (x_P*y_n  - x_n*y_P)  + (x_n*y_nE - x_nE*y_n))/2;
            S_n    = abs((x_Ne*y_e  - x_e*y_Ne)  + (x_e*y_w   - x_w*y_e)   + (x_w*y_Nw  - x_Nw*y_w)  + (x_Nw*y_Ne - x_Ne*y_Nw))/2;
            S_etaw = abs((x_n*y_P  - x_P*y_n)  + (x_P*y_W   - x_W*y_P)   + (x_W*y_nW  - x_nW*y_W)  + (x_nW*y_n - x_n*y_nW))/2;