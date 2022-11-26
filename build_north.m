

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

% Stecil 

% East 
D3=((dx_se_e*(dx_E_P/2 + (3*dx_sE_E)/4 + dx_s_sE/4))/S_etae + (dy_se_e*(dy_E_P/2 + (3*dy_sE_E)/4 + dy_s_sE/4))/S_etae + (dx_Se_e*dx_sw_se)/(4*S_s) + (dy_Se_e*dy_sw_se)/(4*S_s))/S_eta; 

% West 
D_3=((dx_w_sw*(dx_P_W/2 + (3*dx_W_sW)/4 + dx_sW_s/4))/S_etaw + (dy_w_sw*(dy_P_W/2 + (3*dy_W_sW)/4 + dy_sW_s/4))/S_etaw + (dx_w_Sw*dx_sw_se)/(4*S_s) + (dy_w_Sw*dy_sw_se)/(4*S_s))/S_eta; 

% South 
D1=((dx_se_e*(dx_P_s/4 + dx_s_sE/4))/S_etae + (dx_w_sw*(dx_s_P/4 + dx_sW_s/4))/S_etaw + (dy_se_e*(dy_P_s/4 + dy_s_sE/4))/S_etae + (dy_w_sw*(dy_s_P/4 + dy_sW_s/4))/S_etaw + (dx_sw_se*(dx_Se_e/4 + dx_w_Sw/4 + dx_Sw_Se))/S_s + (dy_sw_se*(dy_Se_e/4 + dy_w_Sw/4 + dy_Sw_Se))/S_s)/S_eta; 

% SW 
D_2=((dx_w_sw*(dx_W_sW/4 + dx_sW_s/4))/S_etaw + (dy_w_sw*(dy_W_sW/4 + dy_sW_s/4))/S_etaw + (dx_w_Sw*dx_sw_se)/(4*S_s) + (dy_w_Sw*dy_sw_se)/(4*S_s))/S_eta; 

% SE 
D4=((dx_se_e*(dx_sE_E/4 + dx_s_sE/4))/S_etae + (dy_se_e*(dy_sE_E/4 + dy_s_sE/4))/S_etae + (dx_Se_e*dx_sw_se)/(4*S_s) + (dy_Se_e*dy_sw_se)/(4*S_s))/S_eta; 

% P 
D0=((dx_se_e*(dx_E_P/2 + (3*dx_P_s)/4 + dx_s_sE/4))/S_etae + (dx_w_sw*(dx_P_W/2 + (3*dx_s_P)/4 + dx_sW_s/4))/S_etaw + (dy_se_e*(dy_E_P/2 + (3*dy_P_s)/4 + dy_s_sE/4))/S_etae + (dy_w_sw*(dy_P_W/2 + (3*dy_s_P)/4 + dy_sW_s/4))/S_etaw + (dx_sw_se*(dx_e_w + dx_Se_e/4 + dx_w_Sw/4))/S_s + (dy_sw_se*(dy_e_w + dy_Se_e/4 + dy_w_Sw/4))/S_s + (alpha*bc_ctrl*(dx_e_w - dy_e_w))/lamda(i,j))/S_eta; 

