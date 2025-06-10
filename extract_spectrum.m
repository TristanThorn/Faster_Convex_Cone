clear all;
close all;
% Define parameters
wv_values = 720:20:900; % 21 wavelength values (700, 710, ..., 900)
snr1 = zeros(size(wv_values));
snr2 = zeros(size(wv_values));
snr3 = zeros(size(wv_values));

% Coordinates of interest (row, column)
% point1 = [405, 862];   % First measurement point
% point2 = [572, 879];   % Second measurement point
% point3 = [724, 875];   % Third measurement point

% point1 = [552, 869];   % First measurement point
% point2 = [360, 862];   % Second measurement point
% point3 = [745, 876];   % Third measurement point

% point1 = [497, 571];   % Second measurement point
% point2 = [314, 568];   % First measurement point  0422_1
% point3 = [706, 575];   % Third measurement point

% point1 = [504, 676];   % Second measurement point
% point2 = [318, 674];   % First measurement point  0422_2
% point3 = [710, 687];   % Third measurement point

% point1 = [368, 704];   % First measurement point  0423_new
% point2 = [517, 712];   % Second measurement point
% point3 = [667, 717];   % Third measurement point

point1 = [280, 561];   % First measurement point  0424
point2 = [545 710];   % Third measurement point
point3 = [677, 519];   % Second measurement point 
% reference
load( ['Substance_spectra\spectrumNiSO4_extin.mat' ] );
true_spectrum = spectrum_extin_NiSO4*12.7;
true_spectrum = true_spectrum(3:2:end);

% load( 'Substance_spectra/spectrumCuSO4_extin.mat');
% true_spectrum = spectrum_extin_CuSO4(1:2:end);
% spectrum_CuSO4 = [0.743217,0.83517,0.921447,1.004543,1.075969,1.13385,1.180269,1.214376,1.24127,...
%     1.259913, 1.272579,1.271596,1.271596,1.255983,1.239902,1.239902,1.197211,1.1716,1.147071,1.122401,1.09586];
% spectrum_CuSO4 = [1.164756,1.18318,0.964564,1.051407,1.12858,1.195505,1.244105,1.272197,1.293193,1.34025,1.350587,...
%     1.331975,1.349069,1.325291,1.308238,1.282679,1.271985,1.245274,1.214709,1.179103,1.164756];
% true_spectrum = spectrum_CuSO4(3:2:end);

% calculate compensation
for idx = 1:length(wv_values)
    % Read image
    img_name = sprintf('0426/PA_deep_recon%d.mat', wv_values(idx));
    load(img_name)
    imshow(img,[]);
    % Extract intensity values at three points
    refer = img(point3(2), point3(1));
    
    % Calculate noise standard deviation
    noise_data = img(730:830,500:600); % 0419
    noise_std = std2(noise_data(:));
    noise_mean = mean2(noise_data);
    
    % Compute SNR
    refer_spectrum(idx) = (refer-noise_mean)./noise_std;
end
compensation = refer_spectrum./true_spectrum;


% Process each wavelength image
for idx = 1:length(wv_values)
    % Read image
    % img_name = sprintf('0419_raw/PA_recon%d.mat', wv_values(idx));
    img_name = sprintf('0426/PA_deep_recon%d.mat', wv_values(idx));
    load(img_name)
    imshow(img,[]);
    % Extract intensity values at three points
    val1 = img(point1(2), point1(1));
    val2 = img(point2(2), point2(1));
    val3 = img(point3(2), point3(1));
    
    % Calculate noise standard deviation
    noise_data = img(730:830,500:600); % 0419
    % noise_data = img(720:820,480:580);
    noise_mean = mean2(noise_data);
    noise_std = std2(noise_data(:));
    
    % Compute SNR
    snr1(idx) = (val1-noise_mean)./noise_std;
    snr2(idx) = (val2-noise_mean)./noise_std;
    snr3(idx) = (val3-noise_mean)./noise_std;
end
spectrum1 = snr1./compensation;
spectrum2 = snr2./compensation;
spectrum3 = snr3./compensation;
% Plot results
figure;
plot(wv_values, spectrum1, 'r-o', 'LineWidth', 1.5, 'DisplayName', '0% Ni');
hold on;
plot(wv_values, spectrum2, 'b-s', 'LineWidth', 1.5, 'DisplayName', '50% Ni');
plot(wv_values, spectrum3, 'g-^', 'LineWidth', 1.5, 'DisplayName', '100% Ni');

% figure;
% plot(wv_values, snr1, 'r-o', 'LineWidth', 1.5, 'DisplayName', '0% Ni');
% hold on;
% plot(wv_values, snr2, 'b-s', 'LineWidth', 1.5, 'DisplayName', '50% Ni');
% plot(wv_values, snr3, 'g-^', 'LineWidth', 1.5, 'DisplayName', '100% Ni');


set(gcf, 'Position', [200, 400, 600, 400]);
% Format plot
xlabel('Wavelength (nm)');
ylabel('SNR');
title('SNR Spectrum with compensation using NiSO4');
% title('Raw signal');
legend('Location', 'best');
grid on;
hold off;


%% Linear unmixing
%load absroption spectrum
load( ['Substance_spectra\spectrumCuSO4_extin.mat' ] );
load( ['Substance_spectra\spectrumNiSO4_extin.mat' ] );
load( ['Substance_spectra\spectrum_H2O.mat' ] );
% Plot results
figure;
plot(720:20:900, spectrum_extin_CuSO4(3:2:end), 'r-o', 'LineWidth', 1.5, 'DisplayName', '0% Ni');
hold on;
plot(720:20:900, spectrum_extin_CuSO4(3:2:end)*0.5+spectrum_extin_NiSO4(3:2:end)*0.5*12.7, 'b-s', 'LineWidth', 1.5, 'DisplayName', '50% Ni');
plot(720:20:900, spectrum_extin_NiSO4(3:2:end)*12.7, 'g-^', 'LineWidth', 1.5, 'DisplayName', '100% Ni');

set(gcf, 'Position', [200, 400, 600, 400]);
% Format plot
xlabel('Wavelength (nm)');
ylabel('SNR');
title('Reference spectrum');
legend('Location', 'best');
grid on;
hold off;

spectrum_Hb   = spectrum_extin_CuSO4;
spectrum_HbO2 = spectrum_extin_NiSO4*12.7;

addpath('toolbox');
SO2_tube_linear = linearUnmixing(spectrum2(:)', spectrum_HbO2(3:2:end), spectrum_Hb(3:2:end));
%SO2_tube_linear = linearUnmixing(spectrum2(1:end), spectrum1(1:end), spectrum3(1:end));