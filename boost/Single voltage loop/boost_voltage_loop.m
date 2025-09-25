%% Strategy

% The book first obtains a discrete model for the converter. Then
% approximates it back to a continuous time domain using the bilinear
% trnasformation (Tustin), and the controller is designed in the continuous
% time. The most accurate way would be designing the controller in discrete
% time directly. However, tools for continuous time are usually more
% familiar.

% Since here we start with a continuous time model for the plant, an 
% equivalent continuous time model of the PID is used to design the
% controller entirely in the continuous time domain.


%% Parameters

Vin = 5;
Vo = 12;
L = 23e-6;
C = 3*150e-6 + 20e-6;
Rl = 0.05;
Rc = 0.05;
R = [4 24];

fs = 100e3;
Ts = 1/fs;

T_aa = 2*10e3/(2+10)*470e-12;   % anti aliasing filter

% derived parameters
D = 1-Vin/Vo;

Gdo = Vin/(1-D)^2;
wz1 = 1/(Rc*C);
wrhpz = (1-D)^2*(R-Rl)/L;
w0 = 1/sqrt(L*C)*sqrt((Rl+(1-D)^2*R)./R);
Q = w0./(Rl./L + 1./(C*(R+Rc)));

s=tf('s');

figure

for k=1:length(R)
    Gvd(k) = Gdo*(1+s/wz1)*(1-s/wrhpz(k))/(1 + s/(w0(k)*Q(k)) + s^2/w0(k)^2);
    margin(Gvd(k));
    hold on
    [~,pm] = margin(Gvd(k));
    Pms(k) = pm;
end

legend({'4 Ohm','24 Ohm'})

figure 
pzmap(Gvd(1))
title('Poles and zeros - Open loop')

% Add delays due to PWM (0.5) and calculation (1)
% The controller is designed for the highest load case (4 Ohm)
% The antialiasing filter dynamics is also added

Gvd_delay = Gvd(1)*exp(-3/2*Ts*s)/(1 + s*T_aa);

figure
margin(Gvd_delay)


% The uncompensated loop gain is simply the plant, since the gain of the
% acquisition system is 1 (voltage divider gain is comepnsated in software).
Tu = Gvd_delay;

%% Controller design

% Desired dynamics
w_c = 2*pi*1.2e3;   % crossing frequency (slightly bigger than asked, because we can)
PM = (pi/180)*50;   % Phase margin

% The discrete PID controller has the following form

% G_PID_d(z^-1) = Kp + Ki*(1-z^-1) + Kd*(1-z^-1)

% A continuous time approximation can be obtained with the bilinear
% transform

% G_PID = Kp + Ki/Ts*(1 + s/w_p)/s + Kd*Ts*s/(1 + s/w_p)

% where 
w_p = 2/Ts;

% Its multiplicative form is

% G_PID = G_PI_inf*(1 + w_PI/s)*G_PD_0*(1 + s/w_PD)/(1 + s/w_p)

% The gains will be found in this multiplicative form and converted back
% into PID gains in the end

% Corrected crossing frequency due to errors in the bilinear mapping
% (prewarping)
% (this makes perfect sense in the strategy of the book, here is
% questinable) -  effect is minimum, anyway

w_cp = 2/Ts*tan(w_c*Ts/2);

% The integral part (PI) is neglected now (it mostly influences the low
% frequency response), and the PD part of the controller is obtained by
% forcing the open loop response to yield the desired w_cp and PM. Two
% equations, one for the crossing frequency and other for the PM, allows 
% the determination of the two variables - P and D.

% Original gain at the crossing frequency
Gvd_delay_at_wc = evalfr(Gvd_delay,1i*w_cp);
Gvd_delay_abs = abs(Gvd_delay_at_wc);
Gvd_delay_ang = angle(Gvd_delay_at_wc);

% Necessary PD controller to bring it to the desired values
w_PD = w_cp/tan(PM - pi - Gvd_delay_ang + atan(w_cp/w_p));
G_PD_0 = 1/Gvd_delay_abs*sqrt(1+(w_cp/w_p)^2)/sqrt(1+(w_cp/w_PD)^2);

% PD controller
G_PD = G_PD_0*(1 + s/w_PD)/(1 + s/w_p);

% Open loop response with PD only
G_OL_PD = Gvd_delay*G_PD;

figure 
margin(G_OL_PD);
hold on

% Adding the PI part should not alter the gain/phase of G_OL_PD at the
% crossing frequency. For that, w_PI must be much smaller then w_cp. The
% gain at high frequncies must also be unitary.

w_PI = w_cp/5;
G_PI_inf = 1;

% The final transfer function is
G_PID = G_PI_inf*(1 + w_PI/s)*G_PD_0*(1 + s/w_PD)/(1 + s/w_p);

% Final open loop transfer function
G_OL_PID = Gvd_delay*G_PID;

margin(G_OL_PID);
legend({'OL with PD','OL with PID'})

% Closed loop response
G_CL_PID = feedback(G_OL_PID,1)
figure
step(G_CL_PID)


% The original gains are calculated with

K_P = G_PI_inf*G_PD_0*(1 + w_PI/w_PD - 2*w_PI/w_p)
K_I = 2*G_PI_inf*G_PD_0*w_PI/w_p
K_D = 0.5*G_PI_inf*G_PD_0*(1 - w_PI/w_p)*(w_p/w_PD - 1)

% These gains should be used directly, as in the PID transf. function shown
% above. No need to multiply divide by Ts.




