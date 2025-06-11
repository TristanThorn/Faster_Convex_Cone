addpath('Toolbox');
load('Substance_Spectra/spectrum_Hb_Cope.mat');
load('Substance_Spectra/spectrum_HbO2_Cope.mat');
PA = 0.6 * spectrum_HbO2 + 0.4 * spectrum_Hb;  % simulated measurement
SO2 = linearUnmixing(PA, spectrum_HbO2, spectrum_Hb);
fprintf('Estimated sO2 (linear) = %.2f\n', SO2);
