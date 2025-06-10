%% Single-Coefficient Angular Distance Compensation
clear all; close all;
select_wv = [3:2:21];
%% 1. Define Measurement Points
% 0503 deeper
point1 = [407, 484];   % CuSO4
point2 = [760, 454];   % NiSO4

point3 = [628, 659];
point4 = [387, 676];
point5 = [443, 803];
point6 = [620, 884];




%% 2. Load Reference Spectra
% spectrum_NiSO4 = [1.646392,1.644816,1.634696,1.628864,1.613904,1.601186,1.587817,1.555868,1.498385,...
%     1.40218, 1.277022, 1.126735,0.966332,0.821277,0.698008,0.593562,0.517698,0.463545,0.430178,0.418456,0.425192];
load( ['Substance_spectra\spectrumNiSO4_extin.mat' ] );
true_spectrum_Ni = spectrum_extin_NiSO4(select_wv)*12.7;

load( ['Substance_spectra\spectrumCuSO4_extin.mat' ] );
true_spectrum_Cu = spectrum_extin_CuSO4(select_wv);

load( ['Substance_spectra/spectrum_Hb_Cope.mat' ] );
load( ['Substance_spectra/spectrum_HbO2_Cope.mat' ] );

%% 3. Initialize

wv_all = 700:10:900;
wv_values = wv_all(select_wv);
num_wavelengths = length(wv_values);
snr_Ni = zeros(1, num_wavelengths);
snr_Cu = zeros(1, num_wavelengths);
snr_Target = zeros(1, num_wavelengths);

%% 4. Calculate SNR
for idx = 1:num_wavelengths
    img = load(sprintf('human/PA_ml70_recon%d.mat', wv_values(idx))).img;
    imshow(img, []);

    % Define neighborhood size
    row_offset = 5; % up/down ±2
    col_offset = 10; % left/right ±3

    row_offset1 = 8; % up/down ±2
    col_offset1 = 15; % left/right ±3

    row_offset2 = 3; % up/down ±2
    col_offset2 = 7; % left/right ±3

    % Noise region
    noise_data = img(720:820, 480:580);
    noise_mean = mean2(noise_data);
    noise_std = std2(noise_data(:));


    % Get signals using local average
    region_Ni = img(point2(2)-row_offset:point2(2)+row_offset, point2(1)-col_offset:point2(1)+col_offset);
    signal_Ni = mean(region_Ni(:));

    region_Cu = img(point1(2)-row_offset:point1(2)+row_offset, point1(1)-col_offset:point1(1)+col_offset);
    signal_Cu = mean(region_Cu(:));

    % Compute SNR
    snr_Ni(idx) = signal_Ni;
    snr_Cu(idx) = signal_Cu;

    region_vessel1 = img(point3(2)-row_offset1:point3(2)+row_offset1, point3(1)-col_offset1:point3(1)+col_offset1);
    signal_vessel1(idx) = mean(region_vessel1(:));

    region_vessel2 = img(point4(2)-row_offset1:point4(2)+row_offset1, point4(1)-col_offset1:point4(1)+col_offset1);
    signal_vessel2(idx) = mean(region_vessel2(:));

    region_vessel3 = img(point5(2)-row_offset2:point5(2)+row_offset2, point5(1)-col_offset2:point5(1)+col_offset2);
    signal_vessel3(idx) = mean(region_vessel3(:));

    region_vessel4 = img(point6(2)-row_offset2:point6(2)+row_offset2, point6(1)-col_offset2:point6(1)+col_offset2);
    signal_vessel4(idx) = mean(region_vessel4(:));

    pause(0.3);
end


%% 5. Define Optimization Problem with Regularization

%% Optimize compensation vector to match normalized spectra
normalize = @(x) x / norm(x);  % L2 norm
% normalize = @(x) (x - min(x)) / norm(x - min(x));


% Objective: fit spectrum shape with smoothness regularization
objective = @(comp) ...
    norm(normalize(snr_Ni ./ comp) - normalize(true_spectrum_Ni))^2 + ...
    norm(normalize(snr_Cu ./ comp) - normalize(true_spectrum_Cu))^2 + ...
    0.001 * norm(diff(comp))^2;

% Initial compensation vector
init_comp = ones(size(snr_Ni));

% Bounds to avoid zero-division or extreme values
lb = 0.1 * ones(size(snr_Ni));
ub = 10  * ones(size(snr_Ni));

% Optimization settings
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 5000, ...
    'MaxIterations', 1000);

% Run optimization
opt_compensation = fmincon(objective, init_comp, [], [], [], [], lb, ub, [], options);

%% 6. Apply Compensation
corrected_Ni = snr_Ni ./ opt_compensation;
corrected_Cu = snr_Cu ./ opt_compensation;

corrected_Vessel1 = signal_vessel1./opt_compensation;
corrected_Vessel2 = signal_vessel2./opt_compensation;
corrected_Vessel3 = signal_vessel3./opt_compensation;
corrected_Vessel4 = signal_vessel4./opt_compensation;




%% 7. Normalize All Spectra for Fair Shape Comparison
normalize = @(x) x / norm(x);  % L2-norm normalization

snr_Ni_norm       = normalize(snr_Ni);
snr_Cu_norm       = normalize(snr_Cu);
corrected_Ni_norm = normalize(corrected_Ni);
corrected_Cu_norm = normalize(corrected_Cu);
true_Ni_norm      = normalize(true_spectrum_Ni);
true_Cu_norm      = normalize(true_spectrum_Cu);

%% 8. Plot Normed Comparison
figure;

subplot(2,1,1);
plot(wv_values, snr_Ni_norm, '--c', 'LineWidth', 1); hold on;
plot(wv_values, corrected_Ni_norm, '-bo', 'LineWidth', 1.5);
plot(wv_values, true_Ni_norm, '-k', 'LineWidth', 2);
legend('Original SNR_{Ni} (norm)', 'Corrected Ni (norm)', 'True Ni (norm)', 'Location', 'Best');
xlabel('Wavelength (nm)'); ylabel('Normalized Amplitude');
title('Normalized NiSO4 Spectrum Comparison');
grid on;

subplot(2,1,2);
plot(wv_values, snr_Cu_norm, '--m', 'LineWidth', 1); hold on;
plot(wv_values, corrected_Cu_norm, '-ro', 'LineWidth', 1.5);
plot(wv_values, true_Cu_norm, '-k', 'LineWidth', 2);
legend('Original SNR_{Cu} (norm)', 'Corrected Cu (norm)', 'True Cu (norm)', 'Location', 'Best');
xlabel('Wavelength (nm)'); ylabel('Normalized Amplitude');
title('Normalized CuSO4 Spectrum Comparison');
grid on;

%% 9. Plot Compensation Curve
figure;
plot(wv_values, opt_compensation, 'g-o', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Compensation Factor');
title('Optimized Compensation Vector');
grid on;


%% 8. Plot Compensation Curve
figure;
plot(wv_values, opt_compensation, 'g-o', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Compensation Factor');
title('Optimized Compensation Vector');
grid on;



%% Linear unmixing
%load absroption spectrum

addpath('toolbox');
% SO2_tube_linear = linearUnmixing(spectrum2(:), spectrum_HbO2(5:2:end), spectrum_Hb(5:2:end));
SO2_tube_linear = linearUnmixing(corrected_Vessel1, spectrum_HbO2(select_wv),spectrum_Hb(select_wv))
SO2_tube_linear = linearUnmixing(corrected_Vessel2, spectrum_HbO2(select_wv),spectrum_Hb(select_wv))
SO2_tube_linear = linearUnmixing(corrected_Vessel3, spectrum_HbO2(select_wv),spectrum_Hb(select_wv))
SO2_tube_linear = linearUnmixing(corrected_Vessel4, spectrum_HbO2(select_wv),spectrum_Hb(select_wv))

save('human_deep','corrected_Vessel1','corrected_Vessel2','corrected_Vessel3','corrected_Vessel4')
