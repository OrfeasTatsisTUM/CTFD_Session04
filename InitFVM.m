%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Geometry:
%
%    |---
%    |   ---
%    |      ----
%    |          |
% h1 |----------|  h2  <- symmetry axis
%    |          |
%    |      ----
%    |   ---
%    |---
%
%    |<--  l -->|

% Shape of the Cooling Fin
% 1) 'linear'
% 2) 'quadratic'
% 3) 'crazy'
shape = 'quadratic';

% Define dimension of the trapezoidal domain
% h2 <= h1 !
h1 = 10;
hm = 4;  % only necessary for quatratic option 
h2 = 3;
l  = 10;

% Number of degrees of freedom (number of nodes per length)
dimX = 35;
dimY = 30;

switch shape
    
    case 'linear'
        formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2;

    case 'quadratic'
        c1 = h2+2*h1/2-2*hm;
        c2 = 2*hm - 3*h1/2 - h2/2;
        c3 = h1/2;
        formfunction = @(xnorm) c1*xnorm.^2 +c2*xnorm + c3;
        
    case 'crazy'
        d1 = 3;
        d2 = 4;
        formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2+ (sin(2*pi*d1*xnorm)).*(1-(1-1/d2)*xnorm);

    otherwise
        error(false, 'false shape specified: %s', shape);
    
end

%% Parameter for Conjugated Heat Transfer (For Session 04)
alpha = 5;
Tinf = 0;

%% Boundary conditions (Only Dirichlet applied in Session 03) 
% Type: 1) 'Dirichlet'    2) 'Neumann'    3) 'Robin'
boundary.south = 'Neumann'; % q=0 for mirroring
boundary.north = 'Robin';
boundary.east  = 'Robin';
boundary.west  = 'Dirichlet'; % only Dirichlet can be applied

% Values for Dirichlet BC
TD.north = 10;
TD.south = 50;
TD.west  = 100;
TD.east  = 10;

%% Thermal Conductivity Parameters:
% Thermal conductivity Coefficient 0.05(ice) - 400(pure copper) [W/(m*K)]
% 1) homgenous
% 2) non_homogenous (region with different K)
% 3) random
% 4) linear (changing through x)
heat_conduc = 'homogenous';

% Define the Heat conductivity coefficient values
% for non_homogenous & linear cases it has been assumed that lamda changes on x axis
minlamda = 1;                       % minimum lamda value
deltalamda = 100;                   % lamda difference from side to side (x axis)

maxlamda = minlamda + deltalamda;   % maximum lamda value
switch heat_conduc

    case 'homogenous'
        lamda(1:dimY,1:dimX) = minlamda;

    case 'non_homogenous'
        lamda(1:dimY,1:round(dimX/2)) = minlamda;
        lamda(1:dimY,(round(dimX/2)+1):dimX) = maxlamda;

    case 'random'
        lamda = (deltalamda).*rand(dimY,dimX) + minlamda;

    case 'linear'
        lamda = zeros (dimY,dimX);
        dlamda = deltalamda/l;
        for i=1:dimX
            lamda(1:dimY,i) = minlamda + dlamda * (i-1)*l/(dimX-1);
        end
        
    otherwise
        error(false, 'false lamda type specified: %s', heat_conduc);
        
end
clear minlamda maxlamda deltalamda dlamda

%% Time Dicretization (Session 04)

% Steady/unsteady simulation
% 1) 'steady'
% 2) 'unsteady'
simulationType = 'unsteady';

% Type of time integration
% 1) 'Explicit'
% 2) 'Implicit'
% 3) 'Theta'
% 4) 'RungeKutta4'
TimeIntegrType = 'Theta';

% Parameter for theta scheme  (0 <= θ <= 1)
%  θ = 0: Explicit       θ = 0.5: Crank-Nicolson        θ = 1: Implicit
theta = 0;

% Timestep size and Endtime for unsteady case
dt = 0.003;  % timestep  (when explicit or Runge Kutta 4 it should be around 0.001)
tend = 12; % end-time

% Critical Δt
if strcmp(simulationType, 'unsteady')
    if (theta >= 0 && theta < 0.5)
        D = 60; % Thermal conductivity (λ/(ρ*c_p)) [mm^2/s]
        % https://en.wikipedia.org/wiki/Thermal_diffusivity
        % Material	                                Thermal diffusivity (mm2/s)
        % Pyrolytic graphite, parallel to layers	1220
        % Carbon/carbon composite at 25 °C	        216.5
        % Helium (300 K, 1 atm)	                    190
        % Silver, pure (99.9%)	                    165.63
        % Hydrogen (300 K, 1 atm)	                160
        % Gold	                                    127
        % Copper at 25 °C	                        111
        % Aluminium	                                97
        % Silicon	                                88
        % Al-10Si-Mn-Mg (Silafont 36) at 20 °C	    74.2
        % Aluminium 6061-T6 Alloy	                64
        % Molybdenum (99.95%) at 25 °C	            54.3
        % Al-5Mg-2Si-Mn (Magsimal-59) at 20 °C	    44
        % Tin	                                    40
        % Water vapor (1 atm, 400 K)	            23.38
        % Iron	                                    23
        % Argon (300 K, 1 atm)	                    22
        % Nitrogen (300 K, 1 atm)	                22
        % Air (300 K)	                            19
        % Steel, AISI 1010 (0.1% carbon)	        18.8
        % Aluminium oxide (polycrystalline)	        12
        % Steel, 1% carbon	                        11.72
        % Si3N4 with CNTs 26 °C	                    9.142
        % Si3N4 without CNTs 26 °C	                8.605
        % Steel, stainless 304A at 27 °C	        4.2
        % Pyrolytic graphite, normal to layers	    3.6
        % Steel, stainless 310 at 25 °C	            3.352
        % Inconel 600 at 25 °C	                    3.428
        % Quartz	                                1.4
        % Sandstone	                                1.15
        % Ice at 0 °C	                            1.02
        % Silicon dioxide (polycrystalline)	        0.83
        % Brick, common	                            0.52
        % Glass, window	                            0.34
        % Brick, adobe	                            0.27
        % PC (polycarbonate) at 25 °C	            0.144
        % Water at 25 °C	                        0.143
        % PTFE (Polytetrafluorethylene) at 25 °C	0.124
        % PP (polypropylene) at 25 °C	            0.096
        % Nylon	                                    0.09
        % Rubber	                                0.089 - 0.13
        % Wood (yellow pine)	                    0.082
        % Paraffin at 25 °C	                        0.081
        % PVC (polyvinyl chloride)	                0.08
        % Oil, engine (saturated liquid, 100 °C)	0.0738
        % Alcohol	                                0.07
        c = 0.9; % a constant smaller than 1 (calibration)
        AR = 1;  % Aspect Ratio (dimX/dimY) AR=1->square

        % find average h
        h_avg = 0;
        [x_cur,y_cur] = meshgrid(0:(1/29):1,1:(-1/29):0);
        for i = 1:30
            h_cur = h1 * formfunction(x_cur(1,i));
            h_avg = h_avg + h_cur;
        end
        h_avg = h_avg/30;

        i = 1;
        for n=300:4:10000
            dt_crit(1,i) = n;
            %dt_crit(2,i) = c * h_avg*l/(n - sqrt(n*(1+AR^2)/AR) + 1) / (4*D);
            dx = l / (sqrt(n/AR)-1);
            dy = h_avg / (sqrt(AR*n)-1);
            dt_crit(2,i) = dx^2 * dy^2 / (dx^2 + dy^2) / (2*D*(1-2*theta));
            i = i + 1;
        end
        plot(dt_crit(1,:),dt_crit(2,:));
    end
end