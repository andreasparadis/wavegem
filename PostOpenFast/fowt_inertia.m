clc; clear all;

rho = 1025
g = 9.81

% Masses
%m_rnb = 1E+5 + 1.25E+5 + 0.36E+5; % Mass of rotor, nacelle and blades
%m_tow = 2.1E+5; % Tower mass
%m_float = 1.032001E+6; % Floating platform mass (exc. ballast)
%m_bal = 8.872140E+6; % Total ballast mass

Disp = 1.006431E+04; % Displacement
m_rnb = 230540.344;  % Mass of rotor, nacelle and blades
m_tow = 256154.031;  % Tower mass
m_float = 1.032001E+6; % Floating platform mass (exc. ballast)
%m_bal = 8.872140E+6; % Total ballast mass
%m_moor = 6.7294e+04;
m_moor = 0;

mTOT = rho*Disp;
m_turb = m_rnb + m_tow;
m_ptfm = mTOT - m_turb - m_moor
m_bal = m_ptfm - m_float;

%m_ptfm = m_float + m_bal + m_moor;
%mTOT = m_turb + m_ptfm;

%Disp = mTOT/rho

% CoG (relative to SWL)
CG_rnb = 87;
CG_tow = 44.5;
CG_float = -14.982;
CG_bal = -20.677;
CG_moor = -78.265;

CG_turb = (m_rnb*CG_rnb + m_tow*CG_tow)/m_turb;
CG_ptfm = (m_float*CG_float + m_bal*CG_bal + m_moor*CG_moor)/m_ptfm;
CGTOT = (m_turb*CG_turb + m_ptfm*CG_ptfm)/mTOT;

% Inertia
Ixx_turb = m_turb*CG_turb^2 + m_turb*(CG_turb-CGTOT)^2;
Ixx_ptfm = m_ptfm*CG_ptfm^2 + m_ptfm*(CG_ptfm-CGTOT)^2;
%Ixx = mTOT*CGTOT^2;
Ixx = Ixx_turb + Ixx_ptfm;

Ix1 = 1/12 * 9.8646e4 * (6*2.5^2+40^2);
Ix2 = 1/12 * 3.42178e5 * (6*10^2+30^2);
Ix31 = 1/12 * 6.0445e4 * (6*17.5^2+3.5^2);
Ix32 = 1/4 * 1.51034e5 * 17.5^2;
Ix33 = 1/4 * 1.01736e5 * 17.5^2;
Ix4 = 1/12 * (72.28363+1.27170)*1e5 * (3*(2.5^2+10^2)+9.5^2);
Ix5 = 1/12 * (16.43777+1.50805)*1e5 * (3*(10^2+17.5^2)+3.3^2);

Ixx1 = Ix1 + 9.8646e4 * (-5-CG_float)^2;
Ixx2 = Ix2 + 3.42178e5 * (-10-CG_float)^2;
Ixx31 = Ix31 + 6.0445e4 * (-23.25-CG_float)^2;
Ixx32 = Ix32 + 1.51034e5 * (-25-CG_float)^2;
Ixx33 = Ix33 + 1.01736e5 * (-21.5-CG_float)^2;
Ixx4 = Ix4 + (72.28363+1.27170)*1e5 * (-20.25-CG_float)^2;
Ixx5 = Ix5 + (16.43777+1.50805)*1e5 * (-23.5-CG_float)^2;

Ixx_ptfm_alt = Ixx1 + Ixx2 + Ixx31 + Ixx32 + Ixx33 + Ixx4 + Ixx5
