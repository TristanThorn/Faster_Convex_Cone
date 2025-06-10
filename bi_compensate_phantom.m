% Example phantom compensation using bi_compensate_func
clear; close all;
select_wv = 3:2:21;

points.Ni = [341,388];
points.Cu = [687,397];
points.extras = [320,760; 578,736; 820,751];

load('Substance_spectra/spectrumNiSO4_extin.mat');
load('Substance_spectra/spectrumCuSO4_extin.mat');
ref.Ni = spectrum_extin_NiSO4(select_wv)*12.7;
ref.Cu = spectrum_extin_CuSO4(select_wv);

[opt_compensation, corrected_signals] = bi_compensate_func(points, ref, '0520/PA_deep_recon_%d.mat', 'SelectWV', select_wv);

save('phantom_deep_0520','corrected_signals');
