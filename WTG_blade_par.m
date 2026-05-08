clear; clc;
% Το αιολικό πάρκο Άγιος Γεώργιος περιλαμβάνει 23 Α/Γ συνολικής ισχύος 73,2
% ΜW. Ο τύπος Α/Γ είναι VESTAS, V90 3,0 ΜW (9), και  VESTAS, V112 3,3 ΜW
% (14). Υποθέτουμε στην εργασία μας ότι όλες οι Α/Γ έχουν πύργο ύψους 80μ
% και το μήκος των πτερύγιων 52μ. Σύμφωνα με το βιβλίο lighting protection
% ο αγωγός καθόδου πρέπει να έχει ελάχιστη διάμετρο 10mm και κατά προτίμηση
% 12mm στην εργασία διαλέγετε διάμετρος 16mm με ακτίνα 8mm. Υποθέτουμε ότι
% την στιγμή του κεραυνικού πλήγματος το φτερό βρίσκευται στην στην μέγιστη
% απόσταση από το έδαφος δηλαδή στα 171m ή 143m

% Μπορει για την 3.3MW (14 Α/Γ) να έχει φτερό με μήκος 52μ και πύργο 119μ,
% στην περίπτωση των 3.0MW (9 Α/Γ) έχει φτερό μήκους 38μ και πύργο 105μ. 


%% Constants and user parameters
% Constants for blade
c = 3e8;                % Velocity of light (m/s)
mu_0 = 4*pi*10^-7;      % Permeability of free space (H/m)
mu_r = 1.002;               % Relative permeability (typically 1 for copper/aluminum)
mu = mu_r * mu_0;       % Total material permeability
eps0 = 8.854e-12;

% Material Properties for tower
rho_steel = 1.8e-7; % Resistivity of steel (Ohm-m)
thickness = 0.04;   % Tower shell thickness (m)
mu_r_steel = 300;   % Relative permeability of steel (approx)
mu_tower = mu_r_steel * mu_0;
r1 = 1.5;   % Radius at top (m)
r2 = 2.2;   % Radius at middle (m)
r3 = 3.0;   % Radius at bottom (m)
h1 = 59.5;  % Height of upper to middle (m)
h2 = 59.5;  % Height of middle to lower (m)

%   User-defined parameters blade #############################################
f_m = 5e6;              % Max frequency likely to affect transient (Hz)
rho = 69e-8;          % Material resistivity (Ohm-m) (e.g.,  Stainless steel) blade
r0 = 0.008;             % Radius of the down conductor (m) blade
% H_tower = 119;
H_blade = 52;

% For shielding  RG58A/U coaxial cable
R_shield_km = 21;       % DC resistance given in Ohm/km
D12 = 0.050;            % Distance between conductors or outer diameter (m)
r_core = 0.005;         % Radius of the inner conductor (m)
mu_shield = mu_0;       % Permeability (H/m)

% For control line  RG58A/U coaxial cable

%%   Calculate Segment Length
gamma = 1/rho;          % Conductivity (S/m)
ls_max = (1/10) * (c / f_m);
ls = ls_max/2;      % Choosing a value slightly below the limit
if ls >= ls_max
    warning('Segment length ls does not satisfy the propagation condition!');
end

fprintf('Maximum Segment Length (ls_max): %.4f m\n', ls_max);
fprintf('Chosen Segment Length (ls): %.4f m\n', ls);



%% Blade parameters
%   Calculate Segment Resistance (Rj)
sigma = 1 / sqrt(pi * f_m * mu * gamma);
term1 = exp(-r0 / sigma);
numerator = rho * ls;
denominator = pi * sigma * (1 - term1) * (2 * r0 - sigma * (1 - term1));

Rj = numerator / denominator;  % For 1m 0.0744Ω
fprintf('\n');
fprintf('Resistance of segment j (Rj): %.6f Ohms\n', Rj);

%   Calculate Cj and Lj for the entire Blade Height
z_bounds = H_tower : ls : (H_tower + H_blade);
if z_bounds(end) < (H_tower + H_blade)
    z_bounds = [z_bounds, (H_tower + H_blade)];
end

num_segments = length(z_bounds) - 1;
C_total = zeros(num_segments, 1);
L_total = zeros(num_segments, 1);

for j = 1:num_segments
    z1 = z_bounds(j);
    z2 = z_bounds(j+1);

    
    % Re-calculate Lambda for each specific height
    L_j = r0 * (2 + WTG_f_function((z2 - z1)/r0) + WTG_f_function(-(z2 - z1)/r0));
    
    % Image integral (changes as you go higher)
    L_j_p = r0 * (WTG_f_function(-2*z1/r0) + WTG_f_function(-2*z2/r0) - 2*WTG_f_function(-(z1+z2)/r0));
    
    % Calculate Cj
    num_C = 4 * pi * eps0 * ((z2 - z1)^2);
    den_C = L_j - L_j_p;
    
    C_total(j) = (num_C / den_C) * 1e6;     %   Is in uF
    L_total(j) = (mu_0/(4*pi))*(L_j + L_j_p) * 1e3;       %   Is in mH
    
end

% %% 8. Slip Ring Impedance (Carbon Brush + Holder)
% %Από το https://www.bgbinnovation.com/catalogue/wind-specific/bgb-branded-brushes/product/generator-brush-as-65-v90-spc100432-01
% %Εχουμε Brush Size (Width mm) 40, Brush Radius (mm) 93, Brush Size (Thick
% %mm) 12, Material Type Silver Graphite
% %  
% 
% % Brush Dimensions 
% w = 0.032;      % Width (m)
% d_cb = 0.020;   % Depth (m)
% sigma_c = 3*1e5;  % Conductivity of Carbon (S/m)
% 
% % Calculate Brush Impedance
% delta_c = 1 / sqrt(pi * f_m * mu * sigma_c); % Skin depth in carbon
% R_cb = 2 / (w * delta_c * sigma_c);
% L_cb = (mu * d_cb) / w;
% 
% % Brush Holder Dimensions (Example values)
% a = 0.025; % Inner radius
% b = 0.040; % Outer radius
% t_bh = 0.050; % Thickness/Length
% 
% % Calculate Holder Impedance
% R_bh = (2 * log(b/a)) / (sigma * pi * t_bh);
% L_bh = (2 * mu * log(b/a)) / (pi * t_bh);
% 
% % Total Slip Ring Impedance
% R_total_SR = R_cb + R_bh;
% L_total_SR = L_cb + L_bh;
% 
% fprintf('--- Slip Ring Parameters ---\n');
% fprintf('Carbon Brush: R = %.4f Ohm, L = %.4e mH\n', R_cb, L_cb*1e3);
% fprintf('Brush Holder: R = %.4f Ohm, L = %.4e mH\n', R_bh, L_bh*1e3);

%% Tower Parameters
% Tower Dimensions (Example for a 119m Tower)

% 1. Calculate Average Tower Radius (rtower) - Eq (3)
rtower = (r1*h2 + r2*H_tower + r3*h1) / (2 * H_tower);

% 2. Calculate Tower Inductance (Ltower) - Eq (2)
% c_ratio: inner radius / outer radius
r_inner = rtower - thickness;
c_ratio = r_inner / rtower;

term_L = log((2 * H_tower) / rtower) - 1 - (mu_tower / mu_0) * log(c_ratio);
Ltower = (mu_0 * H_tower / (2 * pi)) * term_L;

% 3. Calculate Tower Resistance (Rtower) - Eq (4)
% Conducting Area A (Annulus)
Area = pi * (rtower^2 - r_inner^2);
Rtower = rho_steel * (H_tower / Area);

% 4. Calculate Tower Capacitance (C0) - image_da3541
C0 = (2 * pi * eps0 * H_tower) / log((2 * H_tower) / rtower);

N = H_tower/ls; % Number of parts
fprintf('\n');
fprintf('--- Tower Parameters per part ---\n');
fprintf('Average Radius: %.4f m\n', rtower);
fprintf('Segmented parts for tower: %.4f parts of %.4f m length\n', N, ls);
fprintf('Total Resistance: %.10f Ohms\n', Rtower/N);
fprintf('Total Inductance: %.10f mH\n', (Ltower/N) * 1e3);
fprintf('Total Capacitance: %.10f uF\n', (C0/N) * 1e6);

%% Shield Parameters
% Convert Ohm/km to Ohm for the specific segment length
R_shield = (R_shield_km / 1000) * ls;

% 2. Calculate Shield Inductance (L_shield) - Eq (5)
% L = (mu / 2*pi) * ln( D12 / (0.7788 * r) )
% Note: 0.7788*r is the Geometric Mean Radius (GMR) for a solid conductor
L_shield_per_meter = (mu_shield / (2 * pi)) * log(D12 / (0.7788 * r_core));
L_shield = L_shield_per_meter * ls;

% Outputs for ATPDraw
fprintf('\n');
fprintf('--- Shielding Layer Parameters ---\n');
fprintf('Resistance R_shield: %.4f Ohms\n', R_shield);
fprintf('Inductance L_shield: %.4e H (%.4f mH)\n', L_shield, L_shield * 1e3);

%% Control Line Parameters
% Geometric Inputs
r_tower_local = rtower;  % Using average radius calculated previously (m)
r_cable = 0.012;         % Outer radius of the cable insulation (m)
D23 = 0.5;               % Distance from tower center to cable center (m)
eps_r = 2.3;             % Relative permittivity of insulation (e.g., Polyethylene)
eps = eps0 * eps_r;      % Total permittivity

% Calculate C12 per meter - image_db133a
numerator_C12 = 2 * pi * eps;
argument_acosh = (r_tower_local^2 + r_cable^2 - D23^2) / (2 * r_tower_local * r_cable);

% Use acosh for inverse hyperbolic cosine
C12_per_meter = numerator_C12 / acosh(argument_acosh);

% Total capacitance for a specific segment length
C12_total = C12_per_meter * ls;
fprintf('\n');
fprintf('--- Cable-to-Tower Capacitance ---\n');
fprintf('C12 per meter: %.4e F/m\n', C12_per_meter);
fprintf('Total C12 (for %.1fm): %.4f pF\n', ls, C12_total * 1e12);


% Calculation: (uH -> mH is /1000) and (1/ft -> 1/m is /0.3048)
L_C = 0.065; % μH/ft
ft_to_m = 0.3048;  % 1 foot = 0.3048 meters
L_C_mH_m = (L_C / 1000) / ft_to_m;
R_C = 8.8; % Ohm/ft
ft_to_m = 0.3048;
R_C_m = R_C / ft_to_m; 


% Geometric Inputs
d_cond = 0.005;          % Conductor radius (m)
D_insul = 0.012;         % Outer radius of insulation (m)
eps_r_jacket = 2.3;      % Relative permittivity of insulation

% Calculate C23 in uF/km - Eq (6)
% Note: MATLAB's log10() is used for base-10 logarithm
C23_per_km = (0.02413 * eps_r_jacket) / log10(D_insul / d_cond);

% Convert to Total Farads for SI consistency
% (uF/km) * (1e-6 F/uF) * (1 km / 1000 m) * length
C23_total_F = (C23_per_km * 1e-6 / 1000) * ls;
fprintf('\n');
fprintf('--- Core-to-Shield Capacitance (C23) ---\n');
fprintf('C23 per km: %.4f uF/km\n', C23_per_km);
fprintf('Total C23 for ATPDraw: %.8f uF\n', C23_per_km * (ls / 1000));


%% Three-Phase Cable Capacitance (C1 and C2)
% Inputs
eps_r_three = 6.5;       % Relative permittivity (using 6~7 range)
h_tm = ls;     % Segment length (m)
D23_shield = 0.50;      % Distance between shield and core center (m)
r1_core = 0.008;         % Radius of the cable core (m)
S_dist = 0.040;          % Center-to-center distance between two cores (m)
r2_outer = D23_shield;

% 1. Core-to-Shield Capacitance (C1) - Eq (14)
% Note: 'lg' usually refers to log10 in these papers
C1_total = (0.02413 * eps0 * eps_r_three * h_tm) / log10(D23_shield / r1_core);

% 2. Mutual Capacitance between cores (C2) - Eq (15)
num_C2 = pi * eps0 * eps_r_three * h_tm;
den_C2 = log((S_dist + sqrt(S_dist^2 - 4 * r1_core^2)) / (2 * r1_core));
C2_total = num_C2 / den_C2;

% 2. Core Self-Inductance (L2) - Eq (16)
term_L2 = log(r2_outer / (0.7788 * r1_core));
L2 = (mu_0 * h_tm / (2 * pi)) * term_L2;

fprintf('\n');
fprintf('--- Three-Phase Cable Capacitance ---\n');
fprintf('C1 (Core-Shield): %.4e F (%.10f μF)\n', C1_total, C1_total * 1e6);
fprintf('C2 (Core-Core): %.4e F (%.10f μF)\n', C2_total, C2_total * 1e6);
fprintf('L2 (Core Inductance): %.4e H (%.6f mH)\n', L2, L2 * 1e3);