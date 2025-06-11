addpath('Toolbox');
load('Substance_Spectra/spectrum_Hb_Cope.mat');
load('Substance_Spectra/spectrum_HbO2_Cope.mat');
PA = 0.6 * spectrum_Hb + 0.4 * spectrum_HbO2;  % simulated measurement
F_base = repmat(linspace(0.8, 1.2, numel(spectrum_Hb))', 1, numel(spectrum_Hb));
[SO2, ~, ~, ~] = convexConeSO2(PA, spectrum_Hb, spectrum_HbO2, F_base);
fprintf('Estimated sO2 (convex cone) = %.2f\n', SO2);
