% This code generates a 2-D FVM stencil for east nodes.
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
target_file = 'build_east.m';

fclose(fopen(target_file, 'w'));

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around somega
syms dy_w_Sw dy_Sw_S dy_S_P dy_P_w real
syms dx_w_Sw dx_Sw_S dx_S_P dx_P_w real

% Around w
syms dy_nW_sW dy_sW_s dy_s_n dy_n_nW  real
syms dx_nW_sW dx_sW_s dx_s_n dx_n_nW real

% Around nomega
syms dy_Nw_w dy_w_P dy_P_N dy_N_Nw real
syms dx_Nw_w dx_w_P dx_P_N dx_N_Nw real 

% Around P
syms dx_s_n dx_sw_s dx_nw_sw dx_n_nw real
syms dy_s_n dy_sw_s dy_nw_sw dy_n_nw dl_s_n real

% Areas
syms S_somega S_nomega S_omega S_w real

% Temperatures
syms  T_S T_W T_N T_NW T_SW T_P real

% inner Temperatures
syms T_n T_s T_sw T_nw T_omega T_Nomega T_Somega T_Nw T_w T_Sw real

% T_P in A.14
syms bc_ctrl alpha lamda Tinf

% Define inner Temperatures as interpolation of outer Temperatures (3.33 - 3.36)
T_n   =(T_P  + T_N)/2;
T_s   =(T_P  + T_S)/2;
T_sw  =(T_SW + T_S  + T_P  + T_W)/4;
T_nw  =(T_NW + T_N  + T_P  + T_W)/4;

T_w   =(T_P  + T_W)/2;
T_Nw  =(T_N  + T_NW)/2;
T_Sw  =(T_S  + T_SW)/2;

T_omega = (3*T_P+T_W)/4;
T_Nomega= (T_NW + 3*T_N)/4;
T_Somega= (T_SW + 3*T_S)/4;

% Gradients (Greens theorem) (A.15 -A.20)
dTdx_somega =   (dy_w_Sw*T_sw + dy_Sw_S*T_Somega + dy_S_P*T_s + dy_P_w*T_omega) /S_somega;
dTdy_somega =  -(dx_w_Sw*T_sw + dx_Sw_S*T_Somega + dx_S_P*T_s + dx_P_w*T_omega) /S_somega;

dTdx_w    =   (dy_nW_sW*T_W + dy_sW_s*T_sw + dy_s_n*T_P + dy_n_nW*T_nw) /S_w;
dTdy_w    =  -(dx_nW_sW*T_W + dx_sW_s*T_sw + dx_s_n*T_P + dx_n_nW*T_nw) /S_w;

dTdx_nomega =   (dy_Nw_w*T_nw + dy_w_P*T_omega + dy_P_N*T_n + dy_N_Nw*T_Nomega) /S_nomega;
dTdy_nomega =  -(dx_Nw_w*T_nw + dx_w_P*T_omega + dx_P_N*T_n + dx_N_Nw*T_Nomega) /S_nomega;


% Build whole stecil acounting for quadratic lamda like in Helmholtz (A.14)

 DDT= (-bc_ctrl*dl_s_n*alpha/lamda*(T_P-Tinf)...
     + dy_sw_s*dTdx_somega - dx_sw_s*dTdy_somega...
     + dy_nw_sw*dTdx_w   - dx_nw_sw*dTdy_w...
     + dy_n_nw*dTdx_nomega - dx_n_nw*dTdy_nomega ) /S_omega;
 
% Make Temperature vector 
T=[T_S; T_W; T_N; T_NW; T_SW; T_P];


% NOTE: The Jacobian function is a "smart" way of factorizing the
% expresssion DDT with respect the temperature vector

stecil=jacobian(DDT,T);  % 3.27

% Find position in file
fileID2 = fopen(target_file, 'r+');

        fprintf(fileID2,'\n\n');
        fprintf(fileID2,'%% Nomenclature:\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%   NW(i-1,j-1) - Nw    -   Nω   -   N(i-1,j)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                 |         |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%       nW - - -  nw  ----- nω ----- n\n');
        fprintf(fileID2,'%%                 |         |\n');
        fprintf(fileID2,'%%       |         |         |        |\n');
        fprintf(fileID2,'%%                 |         |\n');
        fprintf(fileID2,'%%   W(i, j-1) - - w - - - - ω  - - - P (i,j)\n');
        fprintf(fileID2,'%%                 |         |\n');
        fprintf(fileID2,'%%       |         |         |        |\n');
        fprintf(fileID2,'%%                 |         |\n');
        fprintf(fileID2,'%%       sW - - -  sw ------ sω ----- s\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                 |         |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%   SW(i+1,j-1) - Sw    -   Sω   -   S(i+1,j)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%% Indexing of stecil: \n\n');
        fprintf(fileID2,'%%    D_4 - D_1\n');
        fprintf(fileID2,'%%     |     | \n');
        fprintf(fileID2,'%%    D_3 - D_0\n');
        fprintf(fileID2,'%%     |     | \n');
        fprintf(fileID2,'%%    D_2 -  D1\n\n');
              
        fprintf(fileID2,'%% Stecil \n\n');
        fprintf(fileID2,'%% South \n');
        fprintf(fileID2,'D1=%s; \n\n', replace(char(stecil(1)), 'lamda', 'lamda(i,j)'));
   
        fprintf(fileID2,'%% West \n');
        fprintf(fileID2,'D_3=%s; \n\n',replace(char(stecil(2)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% North \n');
        fprintf(fileID2,'D_1=%s; \n\n', replace(char(stecil(3)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% NW \n');
        fprintf(fileID2,'D_4=%s; \n\n',replace(char(stecil(4)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% SW \n');
        fprintf(fileID2,'D_2=%s; \n\n', replace(char(stecil(5)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'%% P \n');
        fprintf(fileID2,'D0=%s; \n\n',replace(char(stecil(6)), 'lamda', 'lamda(i,j)'));

fclose(fileID2);
