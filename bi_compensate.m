%% Single-Coefficient Angular Distance Compensation (phantom example)
% This script demonstrates usage of the generic compensation function.
clear; close all;

select_wv = 3:2:21;

%% Measurement points (pixels)
points.Ni    = [341, 388];    % CuSO4
points.Cu    = [687, 397];    % NiSO4
points.Others = [320,760;     % 60%
                  578,736;     % 25%
                  820,751];    % 80%

%% Load reference spectra
load('Substance_spectra/spectrumNiSO4_extin.mat');
load('Substance_spectra/spectrumCuSO4_extin.mat');
refSpectra.Ni = spectrum_extin_NiSO4(select_wv) * 12.7;
refSpectra.Cu = spectrum_extin_CuSO4(select_wv);

% Path pattern to the reconstructed images
dataPattern = '0520/PA_deep_recon_%d.mat';

[opt_comp, corrected, so2] = bi_compensate_core(points, refSpectra, dataPattern, ...
    'SelectWv', select_wv, ...
    'UnmixSpectra', {refSpectra.Ni, refSpectra.Cu});

% Example access to corrected spectra
corrected_60 = corrected.Others(1,:);
corrected_25 = corrected.Others(2,:);
corrected_80 = corrected.Others(3,:);

save('phantom_deep_0520','corrected_25','corrected_60','corrected_80','so2');
