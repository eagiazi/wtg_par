clear
clc

d=11.56*1e-3; % Diameter of conductor in meters
a=d/2;
mu = 4*pi*1e-7; % Permeability of earth equals to permeability of vacuum
epsilon_r = 10; % Primitivity of earth
epsilon_o = 8.854*1e-12;
epsilon =epsilon_r*epsilon_o;
r_ring = 6.5/2;
theta = 2*pi/20;

% ----------------------USER INPUT----------------------

rho = 1000;      % Resistivity of soil
l = 1;       % Length in meters
d = 0.55;        % Depth in meters
choice = 2;     % 1 for horizontal and 2 for vertical

% ----------------------USER INPUT----------------------

if choice == 1
    % 1. Horizontal DC Resistance
    R = (rho / (pi * l)) * (log((2 * l) / sqrt(2 * a * d)) - 1)
elseif choice == 2;
    % 2. Vertical DC Resistance
    R = (rho / (2 * pi * l)) * (log((4 * l) / a) - 1)
end

% Inductance (mH)
L_mH = ((mu * l) / (2 * pi)) * (log((2 * l) / a) - 1) * 1e3

% Capacitance (uF)
C_uF = ((rho * epsilon) / R) * 1e6

% M=4.1870 e-06
% L=3.1805e-08
% Leq= -1.7530e-11/(6.3610e-08-8.3740e-06) % Leq=(L^2-M^2)/(2L-2M)

fprintf('--------------NEXT MODEL--------------\n');


% %Model 2
r=a;
rho_cu = 1.77*1e-8;
% r_galvanized = 2.5*1e-7;
% rho_cu=r_galvanized;
%l,rho is the same as above


R = rho_cu*l/(pi*r^2)

log_segment = log( (2*l)/sqrt(2*r*d) ) - 1;

R_t = 2*(rho/(pi*l))*(log_segment)/2

L_mH_2 = ((mu*l/(2*pi))*(log_segment)) * 1e3

C_uF_2 = ((pi*epsilon*l)/log_segment) * 1e6



% R_ring = (rho/(pi*theta*r_ring))*log(8*r_ring/sqrt(2*r*d))

% R_ring = (2*rho/(pi*l))*log(1.27*l/sqrt(2*r*d))
% R_ring = (2*rho/(pi*(r_ring*theta)))*log(1.27*(r_ring*theta)/sqrt(2*r*d))
% 
% % r_ring= 13.6/2;
% % theta = 2*pi/43;
% % R_ring3 = (2*rho/(pi*(r_ring*2*pi/43)))*log(1.27*(r_ring*2*pi/43)/sqrt(2*r*d))
% 
% L_ring = (mu*(r_ring*theta)/(2*pi))*log(1.27*(r_ring*theta)/sqrt(2*r*d));
% L_ring = L_ring * 1e3
% C_ring = ((pi*epsilon*(r_ring*theta))/log(1.27*(r_ring*theta)/sqrt(2*r*d))) * 1e6

%1.0324e-07*1e3

%3.5860e-07*1e3
