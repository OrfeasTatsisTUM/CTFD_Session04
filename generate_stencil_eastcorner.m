% This code generates a 2-D FVM stencil for the two east corner nodes.
% The faces of the cells do not have to be aligned to a 
% cartesian grid on any side, i.e. a cell can be any convex quadrilateral

% The code is divided in two parts. The first part generates the stencil
% and prints it to a file given in variable 'target file' after the 
% expression 'Start_stecil' (case sensivtive). The second part
% replaces the scalar expressions with matrix expression and allows to
% calculate the stencil for all inner nodes at once. This accelerates the
% calculation of the system Matrix a lot. 

% Either way all distances and surfaces (areas) have to be provided by you.


% 1 = SE Corner; 2 = NE Corner (of the whole matrix)

clear; clc;

% Print results to file 
target_file = 'build_eastcorner.m';

fclose(fopen(target_file, 'w'));

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around sω
syms dy_w_Sw dy_Sw_S dy_S_P dy_P_w real
syms dx_w_Sw dx_Sw_S dx_S_P dx_P_w real

% Around ηw1
syms dy_nW_W dy_W_P dy_P_n dy_n_nW  real
syms dx_nW_W dx_W_P dx_P_n dx_n_nW real

% Around ηw2
syms dy_W_sW dy_sW_s dy_s_P dy_P_W  real
syms dx_W_sW dx_sW_s dx_s_P dx_P_W real

% Around nω
syms dy_Nw_w dy_w_P dy_P_N dy_N_Nw real
syms dx_Nw_w dx_w_P dx_P_N dx_N_Nw real 

% Around P
syms dx_nw_w dx_w_sw dx_n_nw dx_sw_s real
syms dy_nw_w dy_w_sw dy_n_nw dy_sw_s real

% Areas
syms S_somega S_nomega S_etaomega1 S_etaomega2 S_etaw1 S_etaw2 real

% Temperatures
syms  T_S T_W T_N T_NW T_SW T_P real

% inner Temperatures
syms T_n T_s T_sw T_nw T_omega T_Nomega T_Somega T_Nw T_nW T_sW T_w T_Sw real

% T_P in A.14
syms bc_ctrl_east bc_ctrl_n bc_ctrl_s alpha lamda Tinf

% Define inner Temperatures as interpolation of outer Temperatures (3.33 - 3.36)
T_n   =(T_P  + T_N)/2;
T_s   =(T_P  + T_S)/2;
T_sw  =(T_SW + T_S  + T_P  + T_W)/4;
T_nw  =(T_NW + T_N  + T_P  + T_W)/4;
T_nW  =(T_W  + T_NW)/2;    T_sW  =(T_W + T_SW)/2;

T_w   =(T_P  + T_W)/2;
T_Nw  =(T_N  + T_NW)/2;
T_Sw  =(T_S  + T_SW)/2;

%1 = SE Corner; 2 = NE Corner
T_eta1  =(T_N + 3*T_P)/4;     T_eta2  =(T_S  + 3*T_P)/4;
T_etaW1 =(T_NW + 3*T_W)/4;    T_etaW2 =(T_SW + 3*T_W)/4;
T_omega =(3*T_P  + T_W)/4;
T_Nomega=(T_NW + 3*T_N)/4;
T_Somega=(T_SW + 3*T_S)/4;

% Gradients (Greens theorem) (A.15 -A.20)
dTdx_somega =   (dy_w_Sw*T_sw + dy_Sw_S*T_Somega + dy_S_P*T_s + dy_P_w*T_omega) /S_somega;
dTdy_somega =  -(dx_w_Sw*T_sw + dx_Sw_S*T_Somega + dx_S_P*T_s + dx_P_w*T_omega) /S_somega;

dTdx_etaw1  =   (dy_nW_W*T_etaW1 + dy_W_P*T_w + dy_P_n*T_eta1 + dy_n_nW*T_nw) /S_etaw1;
dTdy_etaw1  =  -(dx_nW_W*T_etaW1 + dx_W_P*T_w + dx_P_n*T_eta1 + dx_n_nW*T_nw) /S_etaw1;

dTdx_etaw2  =   (dy_W_sW*T_etaW2 + dy_sW_s*T_sw + dy_s_P*T_eta2 + dy_P_W*T_w) /S_etaw2;
dTdy_etaw2  =  -(dx_W_sW*T_etaW2 + dx_sW_s*T_sw + dx_s_P*T_eta2 + dx_P_W*T_w) /S_etaw2;

dTdx_nomega =   (dy_Nw_w*T_nw + dy_w_P*T_omega + dy_P_N*T_n + dy_N_Nw*T_Nomega) /S_nomega;
dTdy_nomega =  -(dx_Nw_w*T_nw + dx_w_P*T_omega + dx_P_N*T_n + dx_N_Nw*T_Nomega) /S_nomega;


% Build whole stecil acounting for quadratic lamda like in Helmholtz (A.14)
% bc_ctrl_east, n & s are =1 when we have Robin on the respective edge...
% otherwise they are =0
 DDT1= ((dTdx_etaw1*dy_nw_w - dTdy_etaw1*dx_nw_w...
     - bc_ctrl_s*(dy_w_P - dx_w_P)*alpha/lamda*(T_omega-Tinf) ...
     - bc_ctrl_east*(dy_P_n - dx_P_n)*alpha/lamda*(T_eta1-Tinf)...
     + dTdx_nomega*dy_n_nw -dTdy_nomega*dx_n_nw) /S_etaomega1);

  DDT2= ((dTdx_etaw2*dy_w_sw - dTdy_etaw2*dx_w_sw...
     - bc_ctrl_n*(dy_P_w - dx_P_w)*alpha/lamda*(T_omega-Tinf) ...
     - bc_ctrl_east*(dy_s_P - dx_s_P)*alpha/lamda*(T_eta2-Tinf)...
     + dTdx_somega*dy_sw_s -dTdy_somega*dx_sw_s) /S_etaomega2);
 
% Make Temperature vector 
T1=[T_N; T_W; T_NW; T_P];
T2=[T_S; T_W; T_SW; T_P];

% NOTE: The Jacobian function is a "smart" way of factorizing the
% expresssion DDT with respect the temperature vector

stecil1=jacobian(DDT1,T1);  % 3.27
stecil2=jacobian(DDT2,T2);  % 3.27

% Find position in file
fileID2 = fopen(target_file, 'r+');

        fprintf(fileID2,'\n');
        fprintf(fileID2,"if i == 1 %%ΝE corner\n");
        fprintf(fileID2,'\t\t%% Nomenclature:\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%   W(i, j-1) - - w - - - - ω  - - - P (i,j)\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%       |         |         |        |\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%       ηW - - -  ηw ------ ηω ----- η\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%       |         |         |        |\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%       sW - - -  sw ------ sω ----- s\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%   SW(i+1,j-1) - Sw    -   Sω   -   S(i+1,j)\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%% Indexing of stecil: \n\n');
        fprintf(fileID2,'\t\t%%    D_3 - D_0\n');
        fprintf(fileID2,'\t\t%%     |     | \n');
        fprintf(fileID2,'\t\t%%    D_2 -  D1\n\n');
              
        fprintf(fileID2,'\t\t%% Stecil \n\n');
        fprintf(fileID2,'\t\t%% South \n');
        fprintf(fileID2,'\t\tD1=%s; \n\n', replace(char(stecil2(1)), 'lamda', 'lamda(i,j)'));
   
        fprintf(fileID2,'\t\t%% West \n');
        fprintf(fileID2,'\t\tD_3=%s; \n\n',replace(char(stecil2(2)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'\t\t%% SW \n');
        fprintf(fileID2,'\t\tD_2=%s; \n\n', replace(char(stecil2(3)), 'lamda', 'lamda(i,j)'));
        
        fprintf(fileID2,'\t\t%% P \n');
        fprintf(fileID2,'\t\tD0=%s; \n\n',replace(char(stecil2(4)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'\n');
        fprintf(fileID2,"else %%SE corner\n");
        fprintf(fileID2,'\t\t%% Nomenclature:\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%   NW(i-1,j-1) - Nw    -   Nω   -   N(i-1,j)\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%       nW - - -  nw  ----- nω ----- n\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%       |         |         |        |\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%       ηW - - -  ηw  ----- ηω ----- η\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%       |         |         |        |\n');
        fprintf(fileID2,'\t\t%%                 |         |\n');
        fprintf(fileID2,'\t\t%%   W(i, j-1) - - w - - - - ω  - - - P (i,j)\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%% Indexing of stecil: \n\n');
        fprintf(fileID2,'\t\t%%    D_4 - D_1\n');
        fprintf(fileID2,'\t\t%%     |     | \n');
        fprintf(fileID2,'\t\t%%    D_3 - D_0\n');

        fprintf(fileID2,'\t\t%% Stecil \n\n');
        fprintf(fileID2,'\t\t%% North \n');
        fprintf(fileID2,'\t\tD_1=%s; \n\n', replace(char(stecil1(1)), 'lamda', 'lamda(i,j)'));
   
        fprintf(fileID2,'\t\t%% West \n');
        fprintf(fileID2,'\t\tD_3=%s; \n\n',replace(char(stecil1(2)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'\t\t%% NW \n');
        fprintf(fileID2,'\t\tD_4=%s; \n\n', replace(char(stecil1(3)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'\t\t%% P \n');
        fprintf(fileID2,'\t\tD0=%s; \n\n',replace(char(stecil1(4)), 'lamda', 'lamda(i,j)'));
        fprintf(fileID2,"end");

fclose(fileID2);
