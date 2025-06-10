clear all;
close all;
% Define parameters
wv_values = 720:20:900; % 21 wavelength values (700, 710, ..., 900)
snr1 = zeros(size(wv_values));
snr2 = zeros(size(wv_values));
snr3 = zeros(size(wv_values));

% Coordinates of interest (row, column)
% point1 = [405, 862];   % First measurement point  0419
% point2 = [572, 879];   % Second measurement point
% point3 = [724, 875];   % Third measurement point

% point1 = [552, 869];   % First measurement point  0421_1
% point2 = [360, 862];   % Second measurement point
% point3 = [745, 876];   % Third measurement point

% point1 = [522, 872];   % First measurement point  0421_2
% point2 = [360, 864];   % Second measurement point
% point3 = [748, 878];   % Third measurement point

% point1 = [314, 568];   % First measurement point  0422_1
% point2 = [497, 571];   % Second measurement point
% point3 = [706, 575];   % Third measurement point

% point1 = [318, 674];   % First measurement point  0422_2
% point2 = [504, 676];   % Second measurement point
% point3 = [710, 687];   % Third measurement point
 
% point1 = [373, 675];   % First measurement point  0423_1
% point2 = [555, 678];   % Second measurement point
% point3 = [762, 680];   % Third measurement point

% point1 = [367, 704];   % First measurement point  0423_new
% point2 = [518, 711];   % Second measurement point
% point3 = [668, 715];   % Third measurement point

% point1 = [310, 702];   % First measurement point  0423_new2
% point2 = [454, 711];   % Second measurement point
% point3 = [612, 703];   % Third measurement point

% point1 = [306, 406];   % First measurement point  0424
% point2 = [771, 420];   % Second measurement point
% point3 = [604, 566];   % Third measurement point

point1 = [280, 561];   % First measurement point  0424
point2 = [677, 519];   % Second measurement point 
point3 = [392, 630];   % Third measurement point %75
% [641 636] 25%
% [449 709] 50%
% [545 710] 0%

% reference
% load( ['C:\Users\Shunyao Zhang\Desktop\Spectral_unmixing\Convex-Cone-Approach-main\Substance_spectra\spectrumNiSO4_extin.mat' ] );
% true_spectrum = spectrum_extin_NiSO4*12.7;

% load( 'Substance_spectra/spectrumCuSO4_extin.mat');
% true_spectrum = spectrum_extin_CuSO4;
spectrum_CuSO4 = [0.743217,0.83517,0.921447,1.004543,1.075969,1.13385,1.180269,1.214376,1.24127,...
    1.259913, 1.272579,1.271596,1.271596,1.255983,1.239902,1.239902,1.197211,1.1716,1.147071,1.122401,1.09586];

spectrum_NiSO4 = [1.646392,1.644816,1.634696,1.628864,1.613904,1.601186,1.587817,1.555868,1.498385,...
    1.40218, 1.277022, 1.126735,0.966332,0.821277,0.698008,0.593562,0.517698,0.463545,0.430178,0.418456,0.425192];
true_spectrum = spectrum_CuSO4;

% Process each wavelength image
for idx = 1:length(wv_values)
% for idx = 5
    % Read image
    % img_name = sprintf('0419_raw/PA_recon%d.mat', wv_values(idx));
    img_name = sprintf('0426/PA_deeper_recon%d.mat', wv_values(idx))
    load(img_name)
    imshow(img,[]);
    % Extract intensity values at three points
    val1 = img(point1(2), point1(1));
    val2 = img(point2(2), point2(1));
    val3 = img(point3(2), point3(1));
    
    % Calculate noise standard deviation
    % noise_data = img(730:830,500:600); % 0419
    noise_data = img(720:820,480:580);
    noise_mean = mean2(noise_data);
    noise_std = std2(noise_data(:));
    
    % Compute SNR
    snr1(idx) = (val1-noise_mean)./noise_std;
    snr2(idx) = (val2-noise_mean)./noise_std;
    snr3(idx) = (val3-noise_mean)./noise_std;
    % snr1(idx) = val1;
    % snr2(idx) = val2;
    % snr3(idx) = val3;
end
spectrum1 = snr1;
spectrum2 = snr2;
spectrum3 = snr3;

f_experimental = spectrum1./spectrum3;
f_experimental = f_experimental./max(f_experimental(:));
f_experimental = (f_experimental-min(f_experimental(:)))/(max(f_experimental(:))-min(f_experimental(:)));


load( ['Substance_spectra\spectrumCuSO4_extin.mat' ] );
load( ['Substance_spectra\spectrumNiSO4_extin.mat' ] );
load( ['Substance_spectra\spectrum_H2O.mat' ] );


% f_spectrophotometer = spectrum_extin_CuSO4./(spectrum_extin_NiSO4*12.7);
f_spectrophotometer = spectrum_CuSO4./spectrum_NiSO4;
f_spectrophotometer = f_spectrophotometer(3:2:end);
f_spectrophotometer = f_spectrophotometer./(max(f_spectrophotometer(:)));
f_spectrophotometer = (f_spectrophotometer-min(f_spectrophotometer(:)))/(max(f_spectrophotometer(:))-min(f_spectrophotometer(:)));

figure;
plot(wv_values, f_experimental(), 'r-o', 'LineWidth', 1.5, 'DisplayName', 'PA');
hold on;
plot(wv_values, f_spectrophotometer, 'b-s', 'LineWidth', 1.5, 'DisplayName', 'spectrophotometer');

set(gcf, 'Position', [200, 400, 600, 400]);
% Format plot
xlabel('Wavelength (nm)');
ylabel('SNR');
title('SNR Spectrum with compensation using CuSO4');
legend('Location', 'best');
grid on;
hold off;

