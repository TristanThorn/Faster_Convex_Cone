
addpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1');
%% Load data and spectrum
load('spectrum_HbO2_Cope.mat'); % Molar extinction coefficient of oxyhemoglobin
load('spectrum_Hb_Cope.mat'); % Molar extinction coefficient of deoxyhemoglobin

mu_a = (spectrum_HbO2 .* sBO2 + spectrum_Hb .* (1 - sBO2));


%% Local Fluence Correction
% First block: Calculate mu_a & mu_s'
