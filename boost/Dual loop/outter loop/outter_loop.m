%% Boost Converter Double-Loop PI Controller Design (Exact Method)
%
% Method: At target crossover frequencies, satisfy open-loop conditions:
%   Magnitude = 1, Phase = -180 degrees + PM
%
% Plant 1 (Inner Loop): Duty cycle to inductor current Gid(s)
% Plant 2 (Outer Loop): Inductor current to output voltage Gvi(s)
%
% Feature: Uses exact closed-loop transfer function of the inner loop 
%          via the 'feedback' command instead of a 1st-order approximation.

clear; clc; close all;

%% =========================================================
%  1. Parameter Setup
% =========================================================
% --- Circuit Parameters ---
Vin = 12;           % Input voltage (V)
Vo  = 24;           % Steady-state output voltage (V)
L   = 24e-6;        % Inductance (H)
C   = 100e-6;       % Output capacitance (F)
R   = 6;            % Full-load resistance (ohm) -> P = Vo^2/R = 96 W
fsw = 100e3;        % Switching frequency (Hz)

% --- Sensor Gains ---
Hi  = 1;            % Current sensor gain (V/A)
Hv  = 1;            % Voltage sensor gain (V/V), e.g., 1 for direct feedback

% --- Inner Loop (Current) Specifications ---
fci   = fsw / 10;   % Target crossover frequency (Hz), set to 10 kHz
phi_mi = 60;        % Target phase margin (degrees)

% --- RHPZ Calculation & Outer Loop (Voltage) Specifications ---
D = 1 - Vin / Vo;   % Steady-state duty cycle
% Calculate Right-Half-Plane Zero (RHPZ) frequency
f_RHPZ = (R * (1-D)^2) / (2 * pi * L); 

% The voltage loop bandwidth MUST be lower than f_RHPZ / 3 and fci / 5
fcv   = min(fci / 5, f_RHPZ / 4); % Automatically set safe voltage loop bandwidth
phi_mv = 60;        % Target phase margin (degrees)


%% =========================================================
%  2. Inner Loop (Current) Controller Design
% =========================================================
wci = 2 * pi * fci;

% 1. Calculate Gid(j*wci) exactly
num_real_i = 1;
num_imag_i = wci * R * C / 2;
den_real_i = 1 - wci^2 * L * C / (1 - D)^2;
den_imag_i = wci * L / (R * (1 - D)^2);

Gid_DC = 2 * Vo / (R * (1 - D)^2);
Gid_jwci = Gid_DC * complex(num_real_i, num_imag_i) / complex(den_real_i, den_imag_i);

Gid_mag = abs(Gid_jwci);
theta_plant_i = angle(Gid_jwci) * (180 / pi);

% 2. Inner PI parameters
phi_req_i = -180 + phi_mi - theta_plant_i;
phi_req_rad_i = phi_req_i * pi / 180;

Kp_i = cos(phi_req_rad_i) / (Gid_mag * Hi);
Ki_i = Kp_i * wci * tan(-phi_req_rad_i);

%% =========================================================
%  3. Outer Loop (Voltage) Controller Design (Exact Method)
% =========================================================
s = tf('s');

% 1. Construct Inner Loop Transfer Functions
Gid = Gid_DC * (1 + s * R * C / 2) / (1 + s * L / (R * (1-D)^2) + s^2 * L * C / (1-D)^2);
Gpi_i = Kp_i + Ki_i / s;

% 2. EXPLICITLY calculate exact closed-loop inner transfer function
% Gid_close = (iL / vc) = (Gpi_i * Gid) / (1 + Gpi_i * Gid * Hi)
Gid_close = feedback(Gpi_i * Gid, Hi); 

% 3. Construct Gvi(s): Inductor current to Output voltage
Gvi_DC = R * (1 - D) / 2;
Gvi = Gvi_DC * (1 - s * L / (R * (1 - D)^2)) / (1 + s * R * C / 2);

% 4. Total Voltage Plant = Closed Inner Loop * Gvi
Gv_plant = Gid_close * Gvi;

% 5. Evaluate Total Voltage Plant at Outer Loop Crossover Frequency
wcv = 2 * pi * fcv;
[mag_v, phase_v] = bode(Gv_plant, wcv);
mag_v = squeeze(mag_v);         % Reduce to scalar
phase_v = squeeze(phase_v);     % Reduce to scalar

% 6. Outer PI parameters
phi_req_v = -180 + phi_mv - phase_v;
phi_req_rad_v = phi_req_v * pi / 180;

Kp_v = cos(phi_req_rad_v) / (mag_v * Hv);
Ki_v = Kp_v * wcv * tan(-phi_req_rad_v);

Gpi_v = Kp_v + Ki_v / s;


%% =========================================================
%  4. Results Output
% =========================================================
fprintf('============ System Basics ============\n');
fprintf('Duty cycle               D         = %.4f\n',  D);
fprintf('Right-Half-Plane Zero    f_RHPZ    = %.1f Hz\n', f_RHPZ);

fprintf('\n============ INNER LOOP (Current) ============\n');
fprintf('Target bandwidth         fci       = %.1f Hz\n', fci);
fprintf('Kp_i = %.6g\n', Kp_i);
fprintf('Ki_i = %.6g\n', Ki_i);

fprintf('\n============ OUTER LOOP (Voltage) ============\n');
fprintf('Target bandwidth         fcv       = %.1f Hz\n', fcv);
fprintf('Kp_v = %.6g\n', Kp_v);
fprintf('Ki_v = %.6g\n', Ki_v);


%% =========================================================
%  5. Bode Plot Verification (Both Loops)
% =========================================================
% --- Inner Loop ---
T_open_i = Gpi_i * Gid * Hi;
figure('Name', 'Inner Current Loop', 'NumberTitle', 'off', 'Color', 'w');
margin(T_open_i);
grid on;
title(sprintf('Current Loop Open-Loop Bode (fc = %.0f Hz, PM = %.0f deg)', fci, phi_mi));

% --- Outer Loop ---
T_open_v = Gpi_v * Gv_plant * Hv;
figure('Name', 'Outer Voltage Loop', 'NumberTitle', 'off', 'Color', 'w');
margin(T_open_v);
grid on;
title(sprintf('Voltage Loop Open-Loop Bode (fc = %.0f Hz, PM = %.0f deg)', fcv, phi_mv));

% Extract actual outer loop margins
[Gm_v, Pm_v, ~, Wpm_v] = margin(T_open_v);
fprintf('\n============ Outer Loop Verification ============\n');
fprintf('Actual fcv = %.1f Hz (Target: %.1f Hz)\n', Wpm_v / (2*pi), fcv);
fprintf('Actual PM  = %.1f deg (Target: %.1f deg)\n', Pm_v, phi_mv);


%% =========================================================
%  6. Full Closed-Loop Step Response (Reference to Output Voltage)
% =========================================================
% Total closed-loop transfer function: (Vo / Vref)
T_closed_v = feedback(Gpi_v * Gv_plant, Hv);

figure('Name', 'Full System Output Voltage Step Response', 'NumberTitle', 'off', 'Color', 'w');
step(T_closed_v);
grid on;
title(sprintf('Output Voltage Step Response (fc = %.0f Hz, PM = %.0f deg)', fcv, phi_mv));
xlabel('Time (s)');
ylabel('Output Voltage (V)');

% Extract step response characteristics
info_v = stepinfo(T_closed_v);
fprintf('\n============ Closed-Loop Step Response ============\n');
fprintf('Settling time (2%%)          ts        = %.2e s\n', info_v.SettlingTime);
fprintf('Rise time                   tr        = %.2e s\n', info_v.RiseTime);
fprintf('Overshoot                   OS        = %.2f %%\n', info_v.Overshoot);
fprintf('Peak value                  peak      = %.4f V\n', info_v.Peak);