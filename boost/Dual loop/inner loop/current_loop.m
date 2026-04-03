%% Boost Converter Current Loop PI Controller Design
%
% Method: At target crossover frequency fc, satisfy open-loop conditions:
%   Magnitude: |Gpi(jwc)| * |Gid(jwc)| * Hi = 1
%   Phase: angle Gpi(jwc) + angle Gid(jwc) = -180 degrees + PM
%
% Plant (duty cycle to inductor current):
%   Gid(s) = [2Vo / (R(1-D)^2)] * (1 + s*RC/2)
%            / [1 + s*L/(R(1-D)^2) + s^2*LC/(1-D)^2]
%
% PI controller: Gpi(s) = Kp + Ki/s

clear; clc; close all;

%% =========================================================
%  1. Parameter Setup
% =========================================================

% --- Circuit Parameters ---
Vin = 12;           % Input voltage (V)
Vo  = 24;           % Steady-state output voltage (V)
L   = 24e-6;        % Inductance (H)
C   = 100e-6;       % Output capacitance (F)
R   = 6;            % Full-load resistance (ohm)  ->  P = Vo^2/R = 96 W
fsw = 100e3;        % Switching frequency (Hz)
Hi  = 1;          % Current sensor gain (V/A), e.g., 0.1 ohm sense resistor

% --- Controller Design Specifications ---
fc    = fsw / 10;   % Target crossover frequency (Hz), set to fsw/10
phi_m = 60;         % Target phase margin (degrees)

%% =========================================================
%  2. Core Calculations
% =========================================================

% 1. Steady-state duty cycle
D = 1 - Vin / Vo;

% 2. Target crossover angular frequency
wc = 2 * pi * fc;

% 3. Calculate Gid(j*wc) at s = j*wc term by term
%    Numerator: 1 + j*wc * RC/2
num_real =  1;
num_imag =  wc * R * C / 2;

%    Denominator: (1 - wc^2 * LC/(1-D)^2) + j * wc * L/(R*(1-D)^2)
den_real =  1 - wc^2 * L * C / (1 - D)^2;
den_imag =  wc * L / (R * (1 - D)^2);

Gid_DC  = 2 * Vo / (R * (1 - D)^2);            % DC gain coefficient
Gid_jwc = Gid_DC * complex(num_real, num_imag) ...
                 / complex(den_real, den_imag);

Gid_mag     = abs(Gid_jwc);                     % |Gid(jwc)|
theta_plant = angle(Gid_jwc) * (180 / pi);      % angle Gid(jwc), unit: degrees

% 4. Phase compensation required by PI controller
phi_req = -180 + phi_m - theta_plant;            % unit: degrees

% 5. PI controller parameters
phi_req_rad = phi_req * pi / 180;

Kp = cos(phi_req_rad) / (Gid_mag * Hi);
Ki = Kp * wc * tan(-phi_req_rad);

%% =========================================================
%  3. Results Output
% =========================================================

fprintf('============ Intermediate Variables ============\n');
fprintf('Steady-state duty cycle    D        = %.4f\n',  D);
fprintf('Crossover angular frequency wc       = %.2f rad/s\n', wc);
fprintf('Plant magnitude            |Gid|    = %.4f\n',  Gid_mag);
fprintf('Plant phase                theta_plant = %.2f degrees\n', theta_plant);
fprintf('PI phase compensation      phi_req    = %.2f degrees\n', phi_req);

fprintf('\n============ PI Controller Parameters ============\n');
fprintf('Kp = %.6g\n', Kp);
fprintf('Ki = %.6g\n', Ki);
fprintf('Controller zero frequency  fz = Ki/(2*pi*Kp) = %.1f Hz\n', Ki / (2 * pi * Kp));

%% =========================================================
%  4. Bode Plot Verification
% =========================================================

s = tf('s');

% Plant transfer function
Gid = Gid_DC * (1 + s * R * C / 2) ...
    / (1 + s * L / (R * (1-D)^2) + s^2 * L * C / (1-D)^2);

% PI controller
Gpi = Kp + Ki / s;

% Open-loop transfer function
T_open = Gpi * Gid * Hi;

figure('Name', 'Current Loop Open-Loop Bode Plot', 'NumberTitle', 'off');
margin(T_open);
grid on;
title(sprintf('Current Loop Open-Loop Bode Plot  (Target fc = %.0f Hz, PM = %.0f degrees)', fc, phi_m));

% Extract actual crossover frequency and phase margin
[Gm_lin, Pm_act, ~, Wpm_act] = margin(T_open);

fprintf('\n============ Open-Loop Verification ============\n');
fprintf('Actual crossover frequency  fc_actual = %.1f Hz\n',  Wpm_act / (2*pi));
fprintf('Actual phase margin         PM        = %.1f degrees\n',   Pm_act);
fprintf('Gain margin                 GM        = %.2f dB\n',  20*log10(Gm_lin));

%% =========================================================
%  5. Closed-Loop Step Response
% =========================================================

% Closed-loop transfer function (from reference to inductor current)
T_closed = T_open / (1 + T_open);

% Step response
figure('Name', 'Current Loop Closed-Loop Step Response', 'NumberTitle', 'off');
step(T_closed);
grid on;
title(sprintf('Current Loop Closed-Loop Step Response  (fc = %.0f Hz, PM = %.0f degrees)', fc, phi_m));
xlabel('Time (s)');
ylabel('Inductor Current (A)');

% Extract step response characteristics
info = stepinfo(T_closed);
fprintf('\n============ Closed-Loop Step Response ============\n');
fprintf('Settling time (2%%)          ts        = %.2e s\n', info.SettlingTime);
fprintf('Rise time                   tr        = %.2e s\n', info.RiseTime);
fprintf('Overshoot                   OS        = %.2f %%\n', info.Overshoot);
fprintf('Peak value                  peak      = %.4f A\n', info.Peak);
