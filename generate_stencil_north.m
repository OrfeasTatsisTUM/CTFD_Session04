% This code generates a 2-D FVM stencil for north nodes.
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
target_file = 'build_north.m';

fclose(fopen(target_file, 'w'));

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around etae
syms dy_s_sE dy_sE_E dy_E_P dy_P_s real
syms dx_s_sE dx_sE_E dx_E_P dx_P_s real

% Around s
syms dy_Sw_Se dy_Se_e dy_e_w dy_w_Sw  real
syms dx_Sw_Se dx_Se_e dx_e_w dx_w_Sw real

% Around etaw
syms dy_sW_s dy_s_P dy_P_W dy_W_sW real
syms dx_sW_s dx_s_P dx_P_W dx_W_sW real 

% Around P
syms dy_se_e dy_sw_se dy_w_sw real
syms dx_se_e dx_sw_se dx_w_sw real

% Areas
syms S_etae S_etaw S_eta S_s real

% Temperatures
syms  T_E T_W T_S T_SW T_SE T_P real

% inner Temperatures
syms T_e T_w T_se T_sw T_eta T_etaE T_etaW T_sE T_s T_sW real

% T_P in A.14
syms bc_ctrl alpha lamda Tinf

% Define inner Temperatures as interpolation of outer Temperatures (3.33 - 3.36)
T_e   =(T_P  + T_E)/2;
T_w   =(T_P  + T_W)/2;
T_se  =(T_SE + T_S  + T_P  + T_E)/4;
T_sw  =(T_SW + T_S  + T_P  + T_W)/4;

T_s   =(T_P  + T_S)/2;
T_sE  =(T_E  + T_SE)/2;
T_sW  =(T_W  + T_SW)/2;

T_eta =(3*T_P+T_S)/4;
T_etaE=(T_SE + 3*T_E)/4;
T_etaW=(T_SW + 3*T_W)/4;

% Gradients (Greens theorem) (A.15 -A.20)
dTdx_etae =   (dy_s_sE*T_se + dy_sE_E*T_etaE + dy_E_P*T_e + dy_P_s*T_eta) /S_etae;
dTdy_etae =  -(dx_s_sE*T_se + dx_sE_E*T_etaE + dx_E_P*T_e + dx_P_s*T_eta) /S_etae;

dTdx_s    =   (dy_Sw_Se*T_S + dy_Se_e*T_se + dy_e_w*T_P + dy_w_Sw*T_sw) /S_s;
dTdy_s    =  -(dx_Sw_Se*T_S + dx_Se_e*T_se + dx_e_w*T_P + dx_w_Sw*T_sw) /S_s;

dTdx_etaw =   (dy_sW_s*T_sw + dy_s_P*T_eta + dy_P_W*T_w + dy_W_sW*T_etaW) /S_etaw;
dTdy_etaw =  -(dx_sW_s*T_sw + dx_s_P*T_eta + dx_P_W*T_w + dx_W_sW*T_etaW) /S_etaw;


% Build whole stecil acounting for quadratic lamda like in Helmholtz (A.14)

 DDT= (-bc_ctrl*(dy_e_w-dx_e_w)*alpha/lamda*(T_P-Tinf)...
     + dy_se_e*dTdx_etae - dx_se_e*dTdy_etae...
     + dy_sw_se*dTdx_s   - dx_sw_se*dTdy_s...
     + dy_w_sw*dTdx_etaw - dx_w_sw*dTdy_etaw ) /S_eta;...
 
% Make Temperature vector 
T=[T_E; T_W; T_S; T_SW; T_SE; T_P];


% NOTE: The Jacobian function is a "smart" way of factorizing the
% expresssion DDT with respect the temperature vector

stecil=jacobian(DDT,T);  % 3.27

% Find position in file
fileID2 = fopen(target_file, 'r+');

        fprintf(fileID2,'\n\n');
        fprintf(fileID2,'%% Nomenclature:\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%   W(i, j-1) - - w - - P (i,j) - - e  - - - E (i,j+1)\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |        |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       ηW - - - ηw ------ η ------ ηe - - - ηΕ\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |        |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       sW - - -  sw ----- s ------ se - - - sE\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%   SW(i+1,j-1) - Sw  -  S(i+1,j)  - Se  - SE(i+1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%% Indexing of stecil: \n\n');
        fprintf(fileID2,'%%    D_3 - D_0 - D3\n');
        fprintf(fileID2,'%%     |     |     | \n');
        fprintf(fileID2,'%%    D_2 -  D1 - D4\n\n');
        
        %fprintf(fileID2,'lambda=boundary.lambda; \n\n');
      
        fprintf(fileID2,'%% Stecil \n\n');
        fprintf(fileID2,'%% East \n');
        fprintf(fileID2,'D3=%s; \n\n', replace(char(stecil(1)), 'lamda', 'lamda(i,j)'));
   
        fprintf(fileID2,'%% West \n');
        fprintf(fileID2,'D_3=%s; \n\n',replace(char(stecil(2)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% South \n');
        fprintf(fileID2,'D1=%s; \n\n',replace(char(stecil(3)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% SW \n');
        fprintf(fileID2,'D_2=%s; \n\n',replace(char(stecil(4)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% SE \n');
        fprintf(fileID2,'D4=%s; \n\n',replace(char(stecil(5)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% P \n');
        fprintf(fileID2,'D0=%s; \n\n',replace(char(stecil(6)), 'lamda', 'lamda(i,j)'));

fclose(fileID2);
