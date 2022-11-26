% This code generates a 2-D FVM stencil for the two west corner nodes.
% The faces of the cells do not have to be aligned to a 
% cartesian grid on any side, i.e. a cell can be any convex quadrilateral

% The code is divided in two parts. The first part generates the stencil
% and prints it to a file given in variable 'target file' after the 
% expression 'Start_stecil' (case sensivtive). The second part
% replaces the scalar expressions with matrix expression and allows to
% calculate the stencil for all inner nodes at once. This accelerates the
% calculation of the system Matrix a lot. 

% Either way all distances and surfaces (areas) have to be provided by you.


% 1 = SW Corner; 2 = NW Corner (of the whole matrix)

clear; clc;

% Print results to file 
target_file = 'build_westcorner.m';

fclose(fopen(target_file, 'w'));

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around sω
syms dy_P_S dy_S_Se dy_Se_e dy_e_P real
syms dx_P_S dx_S_Se dx_Se_e dx_e_P real

% Around ηe1
syms dy_n_P dy_P_E dy_E_nE dy_nE_n  real
syms dx_n_P dx_P_E dx_E_nE dx_nE_n real

% Around ηe2
syms dy_P_s dy_s_sE dy_sE_E dy_E_P  real
syms dx_P_s dx_s_sE dx_sE_E dx_E_P real

% Around nω
syms dy_N_P dy_P_e dy_e_Ne dy_Ne_N real
syms dx_N_P dx_P_e dx_e_Ne dx_Ne_N real 

% Around P
syms dy_e_ne dy_ne_n dy_se_e dy_s_se real
syms dx_e_ne dx_ne_n dx_se_e dx_s_se real

% Areas
syms S_somega S_nomega S_etaomega1 S_etaomega2 S_etae1 S_etae2 real

% Temperatures
syms  T_S T_E T_N T_NE T_SE T_P real

% inner Temperatures
syms T_n T_s T_se T_ne T_omega T_Nomega T_Somega T_Ne T_nE T_sE T_e T_Se real

% T_P in A.14
syms bc_ctrl_west bc_ctrl_n bc_ctrl_s alpha lamda Tinf

% Define inner Temperatures as interpolation of outer Temperatures (3.33 - 3.36)
T_n   =(T_P  + T_N)/2;
T_s   =(T_P  + T_S)/2;
T_se  =(T_SE + T_S  + T_P  + T_E)/4;
T_ne  =(T_NE + T_N  + T_P  + T_E)/4;
T_nE  =(T_E  + T_NE)/2;    T_sE  =(T_E + T_SE)/2;

T_e   =(T_P  + T_E)/2;
T_Ne  =(T_N  + T_NE)/2;
T_Se  =(T_S  + T_SE)/2;

%1 = SW Corner; 2 = NW Corner
T_eta1  =(T_N + 3*T_P)/4;     T_eta2  =(T_S  + 3*T_P)/4;
T_etaE1 =(T_NE + 3*T_E)/4;    T_etaE2 =(T_SE + 3*T_E)/4;
T_omega =(3*T_P  + T_E)/4;
T_Nomega=(T_NE + 3*T_N)/4;
T_Somega=(T_SE + 3*T_S)/4;

% Gradients (Greens theorem) (A.15 -A.20)
dTdx_somega =   (dy_P_S*T_s + dy_S_Se*T_Somega + dy_Se_e*T_se + dy_e_P*T_omega) /S_somega;
dTdy_somega =  -(dx_P_S*T_s + dx_S_Se*T_Somega + dx_Se_e*T_se + dx_e_P*T_omega) /S_somega;

dTdx_etae1  =   (dy_n_P*T_eta1 + dy_P_E*T_e + dy_E_nE*T_etaE1 + dy_nE_n*T_ne) /S_etae1;
dTdy_etae1  =  -(dx_n_P*T_eta1 + dx_P_E*T_e + dx_E_nE*T_etaE1 + dx_nE_n*T_ne) /S_etae1;

dTdx_etae2  =   (dy_P_s*T_eta2 + dy_s_sE*T_se + dy_sE_E*T_etaE2 + dy_E_P*T_e) /S_etae2;
dTdy_etae2  =  -(dx_P_s*T_eta2 + dx_s_sE*T_se + dx_sE_E*T_etaE2 + dx_E_P*T_e) /S_etae2;

dTdx_nomega =   (dy_N_P*T_n + dy_P_e*T_omega + dy_e_Ne*T_ne + dy_Ne_N*T_Nomega) /S_nomega;
dTdy_nomega =  -(dx_N_P*T_n + dx_P_e*T_omega + dx_e_Ne*T_ne + dx_Ne_N*T_Nomega) /S_nomega;


% Build whole stecil acounting for quadratic lamda like in Helmholtz (A.14)
% bc_ctrl_west, n & s are =1 when we have Robin on the respective edge...
% otherwise they are =0
 DDT1= (dTdx_etae1*dy_e_ne - dTdy_etae1*dx_e_ne...
     - bc_ctrl_s*(dy_P_e - dx_P_e)*alpha/lamda*(T_omega-Tinf) ...
     - bc_ctrl_west*(dy_n_P - dx_n_P)*alpha/lamda*(T_eta1-Tinf)...
     + dTdx_nomega*dy_ne_n - dTdy_nomega*dx_ne_n) /S_etaomega1;
    
  DDT2= (dTdx_etae2*dy_se_e - dTdy_etae2*dx_se_e...
     - bc_ctrl_n*(dy_e_P - dx_e_P)*alpha/lamda*(T_omega-Tinf) ...
     - bc_ctrl_west*(dy_P_s - dx_P_s)*alpha/lamda*(T_eta2-Tinf)...
     + dTdx_somega*dy_s_se -dTdy_somega*dx_s_se) /S_etaomega2;
 
% Make Temperature vector 
T1=[T_N; T_E; T_NE; T_P];
T2=[T_S; T_E; T_SE; T_P];

% NOTE: The Jacobian function is a "smart" way of factorizing the
% expresssion DDT with respect the temperature vector

stecil1=jacobian(DDT1,T1);  % 3.27
stecil2=jacobian(DDT2,T2);  % 3.27

% Find position in file
fileID2 = fopen(target_file, 'r+');

        fprintf(fileID2,'\n');
        fprintf(fileID2,"if i == 1 %%ΝW corner\n");
        fprintf(fileID2,'\t\t%% Nomenclature:\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%   P (i,j) - - ω  - - e  - - E (i,j+1)\n');
        fprintf(fileID2,'\t\t%%      |                         |\n');
        fprintf(fileID2,'\t\t%%      |        |       |        |      |\n');
        fprintf(fileID2,'\t\t%%      |                         |\n');
        fprintf(fileID2,'\t\t%%      η ------ ηω ---- ηe - -  ηE\n');
        fprintf(fileID2,'\t\t%%      |                         |\n');
        fprintf(fileID2,'\t\t%%      |        |       |        |      |\n');
        fprintf(fileID2,'\t\t%%      |                         |\n');
        fprintf(fileID2,'\t\t%%      s ------ sω ---- se - -  sE\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%      |        |       |        |\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%   S(i+1,j)  - Sω  -  Se  - -  SE(i+1,j+1)\n');
        fprintf(fileID2,'\n');
        fprintf(fileID2,'\t\t%%  Indexing of stecil: \n');
        fprintf(fileID2,'\t\t%%     D_0 - D3\n');
        fprintf(fileID2,'\t\t%%      |     | \n');
        fprintf(fileID2,'\t\t%%     D1 - D4\n\n');
              
        fprintf(fileID2,'\t\t%% Stecil \n\n');
        fprintf(fileID2,'\t\t%% South \n');
        fprintf(fileID2,'\t\tD1=%s; \n\n', replace(char(stecil2(1)), 'lamda', 'lamda(i,j)'));
   
        fprintf(fileID2,'\t\t%% East \n');
        fprintf(fileID2,'\t\tD3=%s; \n\n',replace(char(stecil2(2)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'\t\t%% SE \n');
        fprintf(fileID2,'\t\tD4=%s; \n\n', replace(char(stecil2(3)), 'lamda', 'lamda(i,j)'));
        
        fprintf(fileID2,'\t\t%% P \n');
        fprintf(fileID2,'\t\tD0=%s; \n\n',replace(char(stecil2(4)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'\n');
        fprintf(fileID2,"else %%SW corner\n");
        fprintf(fileID2,'\t\t%% Nomenclature:\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%   N(i-1,j)- - ω - - -  Ne - - NE(i-1,j+1)\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%      |        |        |        |\n');
        fprintf(fileID2,'\t\t%%\n');
        fprintf(fileID2,'\t\t%%      n ------ nω ----- ne - -  nE\n');
        fprintf(fileID2,'\t\t%%      |                          |\n');
        fprintf(fileID2,'\t\t%%      |        |        |        |\n');
        fprintf(fileID2,'\t\t%%      |                          |\n');
        fprintf(fileID2,'\t\t%%      η ------ ηω ----- ηe - -  ηE\n');
        fprintf(fileID2,'\t\t%%      |                          |\n');
        fprintf(fileID2,'\t\t%%      |        |        |        |\n');
        fprintf(fileID2,'\t\t%%      |                          |\n');
        fprintf(fileID2,'\t\t%%   P (i,j) - - ω - - -  e - - E (i,j+1)\n');
        fprintf(fileID2,'\n');
        fprintf(fileID2,'\t\t%%  Indexing of stecil: \n');
        fprintf(fileID2,'\t\t%%     D_1 - D2\n');
        fprintf(fileID2,'\t\t%%      |     | \n');
        fprintf(fileID2,'\t\t%%     D_0 - D3\n\n');

        fprintf(fileID2,'\t\t%% Stecil \n\n');
        fprintf(fileID2,'\t\t%% North \n');
        fprintf(fileID2,'\t\tD_1=%s; \n\n', replace(char(stecil1(1)), 'lamda', 'lamda(i,j)'));
   
        fprintf(fileID2,'\t\t%% East \n');
        fprintf(fileID2,'\t\tD3=%s; \n\n',replace(char(stecil1(2)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'\t\t%% NE \n');
        fprintf(fileID2,'\t\tD2=%s; \n\n', replace(char(stecil1(3)), 'lamda', 'lamda(i,j)'));

        fprintf(fileID2,'\t\t%% P \n');
        fprintf(fileID2,'\t\tD0=%s; \n\n',replace(char(stecil1(4)), 'lamda', 'lamda(i,j)'));
        fprintf(fileID2,"end");

fclose(fileID2);
