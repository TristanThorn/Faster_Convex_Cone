%% Single-Coefficient Angular Distance Compensation (in-vivo example)
% Uses the generic compensation function for human data.
clear; close all;

select_wv = 3:2:21;

%% Measurement points (pixels)
points.Ni    = [760, 454];    % NiSO4
points.Cu    = [407, 484];    % CuSO4
points.Others = [628,659; 387,676; 443,803; 620,884];

%% Load reference spectra
load('Substance_spectra/spectrumNiSO4_extin.mat');
load('Substance_spectra/spectrumCuSO4_extin.mat');
refSpectra.Ni = spectrum_extin_NiSO4(select_wv) * 12.7;
refSpectra.Cu = spectrum_extin_CuSO4(select_wv);

% Spectra for linear unmixing
load('Substance_spectra/spectrum_Hb_Cope.mat');
load('Substance_spectra/spectrum_HbO2_Cope.mat');

% Path pattern to the reconstructed images
dataPattern = 'human/PA_ml70_recon%d.mat';

[opt_comp, corrected, so2] = bi_compensate_core(points, refSpectra, dataPattern, ...
    'SelectWv', select_wv, ...
    'UnmixSpectra', {spectrum_HbO2, spectrum_Hb});

% Extract individual corrected spectra
corrected_Vessel1 = corrected.Others(1,:);
corrected_Vessel2 = corrected.Others(2,:);
corrected_Vessel3 = corrected.Others(3,:);
corrected_Vessel4 = corrected.Others(4,:);

save('human_deep','corrected_Vessel1','corrected_Vessel2','corrected_Vessel3','corrected_Vessel4','so2');
