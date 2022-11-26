
if i == 1 %ΝW corner
		% Nomenclature:
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

		% Stecil 

		% South 
		D1=((dx_s_se*(dx_P_S/2 + (3*dx_S_Se)/4 + dx_Se_e/4))/S_somega + (dy_s_se*(dy_P_S/2 + (3*dy_S_Se)/4 + dy_Se_e/4))/S_somega + (dx_se_e*(dx_P_s/4 + dx_s_sE/4))/S_etae2 + (dy_se_e*(dy_P_s/4 + dy_s_sE/4))/S_etae2 + (alpha*bc_ctrl_west*(dx_P_s - dy_P_s))/(4*lamda(i,j)))/S_etaomega2; 

		% East 
		D3=((dx_se_e*(dx_E_P/2 + (3*dx_sE_E)/4 + dx_s_sE/4))/S_etae2 + (dy_se_e*(dy_E_P/2 + (3*dy_sE_E)/4 + dy_s_sE/4))/S_etae2 + (dx_s_se*(dx_e_P/4 + dx_Se_e/4))/S_somega + (dy_s_se*(dy_e_P/4 + dy_Se_e/4))/S_somega + (alpha*bc_ctrl_n*(dx_e_P - dy_e_P))/(4*lamda(i,j)))/S_etaomega2; 

		% SE 
		D4=((dx_se_e*(dx_sE_E/4 + dx_s_sE/4))/S_etae2 + (dx_s_se*(dx_S_Se/4 + dx_Se_e/4))/S_somega + (dy_se_e*(dy_sE_E/4 + dy_s_sE/4))/S_etae2 + (dy_s_se*(dy_S_Se/4 + dy_Se_e/4))/S_somega)/S_etaomega2; 

		% P 
		D0=((dx_se_e*(dx_E_P/2 + (3*dx_P_s)/4 + dx_s_sE/4))/S_etae2 + (dx_s_se*(dx_P_S/2 + (3*dx_e_P)/4 + dx_Se_e/4))/S_somega + (dy_se_e*(dy_E_P/2 + (3*dy_P_s)/4 + dy_s_sE/4))/S_etae2 + (dy_s_se*(dy_P_S/2 + (3*dy_e_P)/4 + dy_Se_e/4))/S_somega + (3*alpha*bc_ctrl_n*(dx_e_P - dy_e_P))/(4*lamda(i,j)) + (3*alpha*bc_ctrl_west*(dx_P_s - dy_P_s))/(4*lamda(i,j)))/S_etaomega2; 


else %SW corner
		% Nomenclature:
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

		% Stecil 

		% North 
		D_1=((dx_ne_n*(dx_N_P/2 + (3*dx_Ne_N)/4 + dx_e_Ne/4))/S_nomega + (dy_ne_n*(dy_N_P/2 + (3*dy_Ne_N)/4 + dy_e_Ne/4))/S_nomega + (dx_e_ne*(dx_n_P/4 + dx_nE_n/4))/S_etae1 + (dy_e_ne*(dy_n_P/4 + dy_nE_n/4))/S_etae1 + (alpha*bc_ctrl_west*(dx_n_P - dy_n_P))/(4*lamda(i,j)))/S_etaomega1; 

		% East 
		D3=((dx_e_ne*(dx_P_E/2 + (3*dx_E_nE)/4 + dx_nE_n/4))/S_etae1 + (dy_e_ne*(dy_P_E/2 + (3*dy_E_nE)/4 + dy_nE_n/4))/S_etae1 + (dx_ne_n*(dx_P_e/4 + dx_e_Ne/4))/S_nomega + (dy_ne_n*(dy_P_e/4 + dy_e_Ne/4))/S_nomega + (alpha*bc_ctrl_s*(dx_P_e - dy_P_e))/(4*lamda(i,j)))/S_etaomega1; 

		% NE 
		D2=((dx_e_ne*(dx_E_nE/4 + dx_nE_n/4))/S_etae1 + (dx_ne_n*(dx_Ne_N/4 + dx_e_Ne/4))/S_nomega + (dy_e_ne*(dy_E_nE/4 + dy_nE_n/4))/S_etae1 + (dy_ne_n*(dy_Ne_N/4 + dy_e_Ne/4))/S_nomega)/S_etaomega1; 

		% P 
		D0=((dx_e_ne*(dx_P_E/2 + (3*dx_n_P)/4 + dx_nE_n/4))/S_etae1 + (dx_ne_n*(dx_N_P/2 + (3*dx_P_e)/4 + dx_e_Ne/4))/S_nomega + (dy_e_ne*(dy_P_E/2 + (3*dy_n_P)/4 + dy_nE_n/4))/S_etae1 + (dy_ne_n*(dy_N_P/2 + (3*dy_P_e)/4 + dy_e_Ne/4))/S_nomega + (3*alpha*bc_ctrl_s*(dx_P_e - dy_P_e))/(4*lamda(i,j)) + (3*alpha*bc_ctrl_west*(dx_n_P - dy_n_P))/(4*lamda(i,j)))/S_etaomega1; 

end