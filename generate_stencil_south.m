% This code generates a 2-D FVM stencil for south nodes.
% The faces of the cells do not have to be aligned to a 
% cartesian grid on any side, i.e. a cell can be any convex quadrilateral

% The code is divided in two parts. The first part generates the stencil
% and prints it to a file given in variable 'target file' after the 
% expression 'Start_stecil' (case sensivtive). The second part
% replaces the scalar expressions with matrix expression and allows to
% calculate the stencil for all inner nodes at once. This accelerates the
% calculation of the system Matrix a lot. 

% Either way all distances and surfaces (areas) have to be provided by you.



clear; clc;

% Print results to file 
target_file = 'build_south.m';

fclose(fopen(target_file, 'w'));

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around etae
syms dy_P_E dy_E_nE dy_nE_n dy_n_P real
syms dx_P_E dx_E_nE dx_nE_n dx_n_P real

% Around n
syms dy_w_e dy_e_Ne dy_Ne_Nw dy_Nw_w  real
syms dx_w_e dx_e_Ne dx_Ne_Nw dx_Nw_w real

% Around etaw
syms dy_W_P dy_P_n dy_n_nW dy_nW_W real
syms dx_W_P dx_P_n dx_n_nW dx_nW_W real 

% Around P
syms dx_e_ne dx_ne_nw dx_nw_w real
syms dy_e_ne dy_ne_nw dy_nw_w real

% Areas
syms S_etae S_etaw S_eta S_n real

% Temperatures
syms  T_E T_W T_N T_NW T_NE T_P real

% inner Temperatures
syms T_e T_w T_ne T_nw T_eta T_etaE T_etaW T_nE T_n T_nW real

% T_P in A.14
syms bc_ctrl alpha lamda Tinf

% Define inner Temperatures as interpolation of outer Temperatures (3.33 - 3.36)
T_e   =(T_P  + T_E)/2;
T_w   =(T_P  + T_W)/2;
T_ne  =(T_NE + T_N  + T_P  + T_E)/4;
T_nw  =(T_NW + T_N  + T_P  + T_W)/4;

T_n   =(T_P  + T_N)/2;
T_nE  =(T_E  + T_NE)/2;
T_nW  =(T_W  + T_NW)/2;

T_eta =(3*T_P+T_N)/4;
T_etaE=(T_NE + 3*T_E)/4;
T_etaW=(T_NW + 3*T_W)/4;

% Gradients (Greens theorem) (A.15 -A.20)
dTdx_etae =   (dy_P_E*T_e + dy_E_nE*T_etaE + dy_nE_n*T_ne + dy_n_P*T_eta) /S_etae;
dTdy_etae =  -(dx_P_E*T_e + dx_E_nE*T_etaE + dx_nE_n*T_ne + dx_n_P*T_eta) /S_etae;

dTdx_n    =   (dy_w_e*T_P + dy_e_Ne*T_ne + dy_Ne_Nw*T_N + dy_Nw_w*T_nw) /S_n;
dTdy_n    =  -(dx_w_e*T_P + dx_e_Ne*T_ne + dx_Ne_Nw*T_N + dx_Nw_w*T_nw) /S_n;

dTdx_etaw =   (dy_W_P*T_w + dy_P_n*T_eta + dy_n_nW*T_nw + dy_nW_W*T_etaW) /S_etaw;
dTdy_etaw =  -(dx_W_P*T_w + dx_P_n*T_eta + dx_n_nW*T_nw + dx_nW_W*T_etaW) /S_etaw;


% Build whole stecil acounting for quadratic lamda like in Helmholtz (A.14)

 DDT= (-bc_ctrl*(dy_w_e-dx_w_e)*alpha/lamda*(T_P-Tinf)...
     + dy_e_ne*dTdx_etae - dx_e_ne*dTdy_etae...
     + dy_ne_nw*dTdx_n   - dx_ne_nw*dTdy_n...
     + dy_nw_w*dTdx_etaw - dx_nw_w*dTdy_etaw ) /S_eta;...
%      / (bc_ctrl * (dy_w_e-dx_w_e)/S_eta + 1 * ~bc_ctrl);
 
% Make Temperature vector 
T=[T_E; T_W; T_N; T_NW; T_NE; T_P];


% NOTE: The Jacobian function is a "smart" way of factorizing the
% expresssion DDT with respect the temperature vector

stecil=jacobian(DDT,T);  % 3.27

% Find position in file
fileID2 = fopen(target_file, 'r+');

        fprintf(fileID2,'\n\n');
        fprintf(fileID2,'%% Nomenclature:\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%   NW(i-1,j-1) - Nw -  N(i-1,j) -  Ne   -   NE(i-1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%      nW  - - - nw ------ n ------ ne - - - nE\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |        |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%      ηW  - - - ηw ------ η ------ ηe - - - ηΕ\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |        |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%   W(i, j-1) - - w - - P (i,j) - - e - - -  E (i,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%% Indexing of stecil: \n\n');
        fprintf(fileID2,'%%    D_4 - D_1 - D2\n');
        fprintf(fileID2,'%%     |     |     | \n');
        fprintf(fileID2,'%%    D_3 - D_0 - D3\n\n');
        
        %fprintf(fileID2,'lambda=boundary.lambda; \n\n');
      
        fprintf(fileID2,'%% Stecil \n\n');
        fprintf(fileID2,'%% East \n');
        fprintf(fileID2,'D3=%s; \n\n',replace(char(stecil(1)), 'lamda', 'lamda(i,j)'));
   
        fprintf(fileID2,'%% West \n');
        fprintf(fileID2,'D_3=%s; \n\n',replace(char(stecil(2)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% North \n');
        fprintf(fileID2,'D_1=%s; \n\n',replace(char(stecil(3)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% NW \n');
        fprintf(fileID2,'D_4=%s; \n\n',replace(char(stecil(4)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% NE \n');
        fprintf(fileID2,'D2=%s; \n\n',replace(char(stecil(5)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% P \n');
        fprintf(fileID2,'D0=%s; \n\n',replace(char(stecil(6)), 'lamda', 'lamda(i,j)'));

fclose(fileID2);
