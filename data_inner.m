            
% Nomenclature:
%
%   NW(i-1,j-1) - Nw -  N(i-1,j) -  Ne  - NE(i-1,j+1)
%
%                 |                 |
%
%       nW - - -  nw ----- n ------ ne - - nE
%                 |                 |
%       |         |        |        |      |
%                 |                 |
%   W(i, j-1) - - w - - P (i,j) - - e  - - E (i,j+1)
%                 |                 |
%       |         |        |        |      |
%                 |                 |
%       sW - - -  sw ----- s ------ se - - sE
%
%                 |                 |
%
%   SW(i+1,j-1) - Sw  -  S(i+1,j)  - Se  - SE(i+1,j+1)
%
% Indexing of stecil: 

%    D_4 - D_1 - D2
%     |     |     | 
%    D_3 - D_0 - D3
%     |     |     | 
%    D_2 -  D1 - D4
        
% Principal node coordinates
        y_NW    = Y(i-1,j-1);   x_NW    = X(i-1,j-1);
        y_N     = Y(i-1,j);     x_N     = X(i-1,j);
        y_NE    = Y(i-1,j+1);   x_NE    = X(i-1,j+1);
        y_E     = Y(i,j+1);     x_E     = X(i,j+1);
        y_SE    = Y(i+1,j+1);   x_SE    = X(i+1,j+1);
        y_S     = Y(i+1,j);     x_S     = X(i+1,j);
        y_SW    = Y(i+1,j-1);   x_SW    = X(i+1,j-1);
        y_W     = Y(i,j-1);     x_W     = X(i,j-1);
        y_P     = Y(i,j);       x_P     = X(i,j);
     
% Auxiliary node coordinates    
        y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2;
        y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;
        y_nE = (y_NE + y_E)/2;  x_nE = (x_NE + x_E)/2;
        y_sE = (y_SE + y_E)/2;  x_sE = (x_SE + x_E)/2;
        y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;
        y_Sw = (y_SW + y_S)/2;  x_Sw = (x_SW + x_S)/2;
        y_sW = (y_SW + y_W)/2;  x_sW = (x_SW + x_W)/2;
        y_nW = (y_NW + y_W)/2;  x_nW = (x_NW + x_W)/2;

        y_s  = (y_P + y_S)/2;   x_s  = (x_P + x_S)/2;
        y_e  = (y_P + y_E)/2;   x_e  = (x_P + x_E)/2;
        y_n  = (y_P + y_N)/2;   x_n  = (x_P + x_N)/2;
        y_w  = (y_P + y_W)/2;   x_w  = (x_P + x_W)/2;        
        y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;
        y_sw = (y_Sw + y_w)/2;  x_sw = (x_Sw + x_w)/2;
        y_nw = (y_Nw + y_w)/2;  x_nw = (x_Nw + x_w)/2;
        y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;

        

% Inter node distances

        % Around s 
        dy_Sw_Se = y_Se - y_Sw;   dx_Sw_Se = x_Se - x_Sw;
        dy_w_Sw  = y_Sw - y_w;    dx_w_Sw  = x_Sw - x_w;
        dy_e_w   = y_w - y_e;     dx_e_w   = x_w - x_e;
        dy_Se_e  = y_e - y_Se;    dx_Se_e  = x_e - x_Se;
  
        % Around e
        dy_s_sE  = y_sE - y_s;    dx_s_sE  = x_sE - x_s;
        dy_sE_nE = y_nE - y_sE;   dx_sE_nE = x_nE - x_sE;
        dy_nE_n  = y_n - y_nE;    dx_nE_n  = x_n - x_nE; 
        dy_n_s   = y_s - y_n;     dx_n_s   = x_s - x_n;

        % Around n        
        dy_w_e   = y_e - y_w;     dx_w_e   = x_e - x_w;
        dy_e_Ne  = y_Ne - y_e;    dx_e_Ne  = x_Ne - x_e;
        dy_Ne_Nw = y_Nw - y_Ne;   dx_Ne_Nw = x_Nw - x_Ne;
        dy_Nw_w  = y_w - y_Nw;    dx_Nw_w  = x_w - x_Nw;

        % Around w
        dy_sW_s  = y_s - y_sW;    dx_sW_s  = x_s - x_sW;
        dy_s_n   = y_n - y_s;     dx_s_n   = x_n - x_s;
        dy_n_nW  = y_nW - y_n;    dx_n_nW  = x_nW - x_n;
        dy_nW_sW = y_sW - y_nW;   dx_nW_sW = x_sW - x_nW;

        % Around P
        dy_sw_se = y_se - y_sw;  dx_sw_se = x_se - x_sw;
        dy_se_ne = y_ne - y_se;  dx_se_ne = x_ne - x_se;
        dy_ne_nw = y_nw - y_ne;  dx_ne_nw = x_nw - x_ne;
        dy_nw_sw = y_sw - y_nw;  dx_nw_sw = x_sw - x_nw;

% Areas
        
        S_P = abs((x_ne*y_se - x_se*y_ne) + (x_se*y_sw - x_sw*y_se) + (x_sw*y_nw - x_nw*y_sw) + (x_nw*y_ne - x_ne*y_nw))/2;
        S_s = abs((x_e*y_Se  - x_Se*y_e)  + (x_Se*y_Sw - x_Sw*y_Se) + (x_Sw*y_w  - x_w*y_Sw)  + (x_w*y_e   - x_e*y_w))/2;
        S_e = abs((x_nE*y_sE - x_sE*y_nE) + (x_sE*y_s  - x_s*y_sE)  + (x_s*y_n   - x_n*y_s)   + (x_n*y_nE  - x_nE*y_n))/2;
        S_n = abs((x_Ne*y_e  - x_e*y_Ne)  + (x_e*y_w   - x_w*y_e)   + (x_w*y_Nw  - x_Nw*y_w)  + (x_Nw*y_Ne - x_Ne*y_Nw))/2;
        S_w = abs((x_n*y_s   - x_s*y_n)   + (x_s*y_sW  - x_sW*y_s)  + (x_sW*y_nW - x_nW*y_sW) + (x_nW*y_n  - x_n*y_nW))/2;