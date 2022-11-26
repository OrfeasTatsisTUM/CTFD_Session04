
if i == 1 %ΝE corner
		% Nomenclature:
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

		% Stecil 

		% South 
		D1=((dx_sw_s*(dx_S_P/2 + (3*dx_Sw_S)/4 + dx_w_Sw/4))/S_somega + (dy_sw_s*(dy_S_P/2 + (3*dy_Sw_S)/4 + dy_w_Sw/4))/S_somega + (dx_w_sw*(dx_s_P/4 + dx_sW_s/4))/S_etaw2 + (dy_w_sw*(dy_s_P/4 + dy_sW_s/4))/S_etaw2 + (alpha*bc_ctrl_east*(dx_s_P - dy_s_P))/(4*lamda(i,j)))/S_etaomega2; 

		% West 
		D_3=((dx_w_sw*(dx_P_W/2 + (3*dx_W_sW)/4 + dx_sW_s/4))/S_etaw2 + (dy_w_sw*(dy_P_W/2 + (3*dy_W_sW)/4 + dy_sW_s/4))/S_etaw2 + (dx_sw_s*(dx_P_w/4 + dx_w_Sw/4))/S_somega + (dy_sw_s*(dy_P_w/4 + dy_w_Sw/4))/S_somega + (alpha*bc_ctrl_n*(dx_P_w - dy_P_w))/(4*lamda(i,j)))/S_etaomega2; 

		% SW 
		D_2=((dx_w_sw*(dx_W_sW/4 + dx_sW_s/4))/S_etaw2 + (dx_sw_s*(dx_Sw_S/4 + dx_w_Sw/4))/S_somega + (dy_w_sw*(dy_W_sW/4 + dy_sW_s/4))/S_etaw2 + (dy_sw_s*(dy_Sw_S/4 + dy_w_Sw/4))/S_somega)/S_etaomega2; 

		% P 
		D0=((dx_w_sw*(dx_P_W/2 + (3*dx_s_P)/4 + dx_sW_s/4))/S_etaw2 + (dx_sw_s*(dx_S_P/2 + (3*dx_P_w)/4 + dx_w_Sw/4))/S_somega + (dy_w_sw*(dy_P_W/2 + (3*dy_s_P)/4 + dy_sW_s/4))/S_etaw2 + (dy_sw_s*(dy_S_P/2 + (3*dy_P_w)/4 + dy_w_Sw/4))/S_somega + (3*alpha*bc_ctrl_n*(dx_P_w - dy_P_w))/(4*lamda(i,j)) + (3*alpha*bc_ctrl_east*(dx_s_P - dy_s_P))/(4*lamda(i,j)))/S_etaomega2; 


else %SE corner
		% Nomenclature:
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
		% Stecil 

		% North 
		D_1=((dx_n_nw*(dx_P_N/2 + (3*dx_N_Nw)/4 + dx_Nw_w/4))/S_nomega + (dy_n_nw*(dy_P_N/2 + (3*dy_N_Nw)/4 + dy_Nw_w/4))/S_nomega + (dx_nw_w*(dx_P_n/4 + dx_n_nW/4))/S_etaw1 + (dy_nw_w*(dy_P_n/4 + dy_n_nW/4))/S_etaw1 + (alpha*bc_ctrl_east*(dx_P_n - dy_P_n))/(4*lamda(i,j)))/S_etaomega1; 

		% West 
		D_3=((dx_nw_w*(dx_W_P/2 + (3*dx_nW_W)/4 + dx_n_nW/4))/S_etaw1 + (dy_nw_w*(dy_W_P/2 + (3*dy_nW_W)/4 + dy_n_nW/4))/S_etaw1 + (dx_n_nw*(dx_w_P/4 + dx_Nw_w/4))/S_nomega + (dy_n_nw*(dy_w_P/4 + dy_Nw_w/4))/S_nomega + (alpha*bc_ctrl_s*(dx_w_P - dy_w_P))/(4*lamda(i,j)))/S_etaomega1; 

		% NW 
		D_4=((dx_nw_w*(dx_nW_W/4 + dx_n_nW/4))/S_etaw1 + (dx_n_nw*(dx_N_Nw/4 + dx_Nw_w/4))/S_nomega + (dy_nw_w*(dy_nW_W/4 + dy_n_nW/4))/S_etaw1 + (dy_n_nw*(dy_N_Nw/4 + dy_Nw_w/4))/S_nomega)/S_etaomega1; 

		% P 
		D0=((dx_nw_w*(dx_W_P/2 + (3*dx_P_n)/4 + dx_n_nW/4))/S_etaw1 + (dx_n_nw*(dx_P_N/2 + (3*dx_w_P)/4 + dx_Nw_w/4))/S_nomega + (dy_nw_w*(dy_W_P/2 + (3*dy_P_n)/4 + dy_n_nW/4))/S_etaw1 + (dy_n_nw*(dy_P_N/2 + (3*dy_w_P)/4 + dy_Nw_w/4))/S_nomega + (3*alpha*bc_ctrl_s*(dx_w_P - dy_w_P))/(4*lamda(i,j)) + (3*alpha*bc_ctrl_east*(dx_P_n - dy_P_n))/(4*lamda(i,j)))/S_etaomega1; 

end