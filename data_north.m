% Nomenclature:
        %
        %   W(i, j-1) - - w - - P (i,j) - - e  - - - E (i,j+1)
        %                 |                 |
        %       |         |        |        |        |
        %                 |                 |
        %       ηW - - - ηw ------ η ------ ηe - - - ηΕ
        %                 |                 |
        %       |         |        |        |        |
        %                 |                 |
        %       sW - - -  sw ----- s ------ se - - - sE
        %
        %                 |                 |
        %
        %   SW(i+1,j-1) - Sw  -  S(i+1,j)  - Se  - SE(i+1,j+1)
        %
        % Indexing of stecil:

        %    D_3 - D_0 - D3
        %     |     |     |
        %    D_2 -  D1 - D4

        % Principal node coordinates
            y_W     = Y(i,j-1);     x_W     = X(i,j-1);
            y_P     = Y(i,j);       x_P     = X(i,j);
            y_E     = Y(i,j+1);     x_E     = X(i,j+1);
            y_SE    = Y(i+1,j+1);   x_SE    = X(i+1,j+1);
            y_S     = Y(i+1,j);     x_S     = X(i+1,j);
            y_SW    = Y(i+1,j-1);   x_SW    = X(i+1,j-1);

        % Auxiliary node coordinates
            y_Sw = (y_SW + y_S)/2;  x_Sw = (x_SW + x_S)/2;
            y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
            y_sW = (y_SW + y_W)/2;  x_sW = (x_SW + x_W)/2;
            y_sE = (y_SE + y_E)/2;  x_sE = (x_SE + x_E)/2;

            y_e  = (y_P + y_E)/2;   x_e  = (x_P + x_E)/2;
            y_s  = (y_P + y_S)/2;   x_s  = (x_P + x_S)/2;
            y_w  = (y_P + y_W)/2;   x_w  = (x_P + x_W)/2;
            y_sw = (y_Sw + y_w)/2;  x_sw = (x_Sw + x_w)/2;
            y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;

            y_etaW = (y_sW + y_W)/2;  x_etaW = (x_sW + x_W)/2;
            y_etaw = (y_sw + y_w)/2;  x_etaw = (x_sw + x_w)/2;
            y_eta  = (y_s  + y_P)/2;  x_eta  = (x_s  + x_P)/2;
            y_etae = (y_se + y_e)/2;  x_etae = (x_se + x_e)/2;
            y_etaE = (y_sE + y_E)/2;  x_etaE = (x_sE + x_E)/2;

        % Inter node distances

            % Around ηe
            dy_s_sE  = y_sE - y_s;  dx_s_sE  = x_sE - x_s;
            dy_P_s   = y_s - y_P;   dx_P_s   = x_s - x_P;
            dy_E_P   = y_P - y_E;   dx_E_P   = x_P - x_E;
            dy_sE_E  = y_E - y_sE;  dx_sE_E  = x_E - x_sE;

            % Around s
            dy_Se_e  = y_e - y_Se;  dx_Se_e  = x_e - x_Se;
            dy_Sw_Se = y_Se - y_Sw; dx_Sw_Se = x_Se - x_Sw;
            dy_w_Sw  = y_Sw - y_w;  dx_w_Sw  = x_Sw - x_w;
            dy_e_w   = y_w - y_e;   dx_e_w   = x_w - x_e;

            % Around ηw
            dy_s_P  = y_P - y_s;    dx_s_P  = x_P - x_s;
            dy_sW_s = y_s - y_sW;   dx_sW_s = x_s - x_sW;
            dy_W_sW = y_sW - y_W;   dx_W_sW = x_sW - x_W;
            dy_P_W  = y_W - y_P;    dx_P_W  = x_W - x_P;

            % Around P (e -w is common edge with control volume around s)
            dy_se_e  = y_e - y_se;  dx_se_e  = x_e - x_se;
            dy_sw_se = y_se - y_sw; dx_sw_se = x_se - x_sw;
            dy_w_sw  = y_sw - y_w;  dx_w_sw  = x_sw - x_w;

        % Areas
            S_eta  = abs((x_e*y_se  - x_se*y_e)  + (x_se*y_sw   - x_sw*y_se)   + (x_sw*y_w  - x_w*y_sw)  + (x_w*y_e - x_e*y_w))/2;
            S_etae = abs((x_E*y_sE  - x_sE*y_E)  + (x_sE*y_s   - x_s*y_sE)   + (x_s*y_P  - x_P*y_s)  + (x_P*y_E - x_E*y_P))/2;
            S_s    = abs((x_e*y_Se  - x_Se*y_e)  + (x_Se*y_Sw   - x_Sw*y_Se)   + (x_Sw*y_w  - x_w*y_Sw)  + (x_w*y_e - x_e*y_w))/2;
            S_etaw = abs((x_P*y_s  - x_s*y_P)  + (x_s*y_sW   - x_sW*y_s)   + (x_sW*y_W  - x_W*y_sW)  + (x_W*y_P - x_P*y_W))/2;