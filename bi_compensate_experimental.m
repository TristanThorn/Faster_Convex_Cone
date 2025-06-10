%% Bi-spectrum compensation for human experiment using helper function
clear; close all;
select_wv = 3:2:21;

% Measurement points for the human dataset
points.Ni = [760,454];
points.Cu = [407,484];
points.extras = [628,659; 387,676; 443,803; 620,884];

% Custom offsets for local averaging
offsets.Ni = [5 10];
offsets.Cu = [5 10];
offsets.extras = [8 15; 8 15; 3 7; 3 7];

load('Substance_spectra/spectrumNiSO4_extin.mat');
load('Substance_spectra/spectrumCuSO4_extin.mat');
load('Substance_spectra/spectrum_Hb_Cope.mat');
load('Substance_spectra/spectrum_HbO2_Cope.mat');
ref.Ni = spectrum_extin_NiSO4(select_wv)*12.7;
ref.Cu = spectrum_extin_CuSO4(select_wv);

[opt_compensation, corrected_signals, so2] = bi_compensate_func(points, ref, 'human/PA_ml70_recon%d.mat', 'SelectWV', select_wv, 'Offsets', offsets, 'LinSpectra', {spectrum_HbO2, spectrum_Hb});

save('human_deep','corrected_signals','so2');

