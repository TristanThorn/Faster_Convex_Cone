%% Bi-spectrum compensation using helper function
clear; close all;
select_wv = 3:2:21;

% Measurement points for 0520 phantom data
points.Ni = [341,388];   % CuSO4 in image
points.Cu = [687,397];   % NiSO4 in image
points.extras = [320,760;   % 60%
                  578,736;   % 25%
                  820,751];  % 80%

load('Substance_spectra/spectrumNiSO4_extin.mat');
load('Substance_spectra/spectrumCuSO4_extin.mat');
ref.Ni = spectrum_extin_NiSO4(select_wv)*12.7;
ref.Cu = spectrum_extin_CuSO4(select_wv);

[opt_compensation, corrected_signals, ~] = bi_compensate_func(points, ref, '0520/PA_deep_recon_%d.mat', 'SelectWV', select_wv);

save('phantom_deep_0520','corrected_signals');

