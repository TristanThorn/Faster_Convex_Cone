%% Single-Coefficient Angular Distance Compensation
clear all; close all;
select_wv = [3:2:21];
%% 1. Define Measurement Points
% 0426 deep
% point1 = [280, 561];   % CuSO4
% point2 = [677, 519];   % NiSO4
% 
% point3 = [392, 630]; % 75%
% point4 = [641 636]; % 25%
% point5 = [449 709]; % 50%
% point6 = [545 710]; % 0%

% 0426 deeper
% point1 = [349, 445];   % CuSO4
% point2 = [629, 433];   % NiSO4
% 
% point3 = [390, 557]; % 75%
% point4 = [632 564]; % 25%
% point5 = [448 640]; % 50%
% point6 = [544 639]; % 0%

% 0427 deep
% point1 = [386, 503];   % CuSO4
% point2 = [816, 511];   % NiSO4
% 
% point3 = [422, 626]; % 75%
% point4 = [656 620]; % 25%
% point5 = [602 721]; % 50%
% point6 = [713 723]; % 0%

% % 0427 deeper
% point1 = [323, 370];   % CuSO4
% point2 = [773, 373];   % NiSO4
% 
% point3 = [331, 594]; % 75%
% point4 = [553 589]; % 25%
% point5 = [491 688]; % 50%
% point6 = [614 690]; % 0%

% % 0501 deeper
% point1 = [326, 597];   % CuSO4
% point2 = [597, 608];   % NiSO4
% 
% point3 = [340, 568]; % 75%
% point4 = [584 760]; % 25%
% point5 = [727 792]; % 50%
% point6 = [616 820]; % 0%

% % 0507 deeper
% point1 = [364, 448];   % CuSO4
% point2 = [769, 464];   % NiSO4
% 
% point3 = [391, 533]; % 50%
% point4 = [747 545]; % 80%
% point5 = [543 619]; % 20%

% % 0507 deeper
% point1 = [473, 382];   % CuSO4
% point2 = [673, 390];   % NiSO4
% 
% point3 = [437, 525]; % 20%
% point4 = [701 512]; % 50%
% point5 = [399 589]; % 25%
% point6 = [793 601]; % 80%

% 0520
point2 = [341, 388];   % CuSO4
point1 = [687, 397];   % NiSO4

point3 = [320, 760]; % 60%
point4 = [578 736]; % 25%
point5 = [820 751]; % 80%


%% 2. Load Reference Spectra
% spectrum_NiSO4 = [1.646392,1.644816,1.634696,1.628864,1.613904,1.601186,1.587817,1.555868,1.498385,...
%     1.40218, 1.277022, 1.126735,0.966332,0.821277,0.698008,0.593562,0.517698,0.463545,0.430178,0.418456,0.425192];
load( ['Substance_spectra\spectrumNiSO4_extin.mat' ] );
true_spectrum_Ni = spectrum_extin_NiSO4(select_wv)*12.7;

% CuSO4 spectrum (example values)
% spectrum_CuSO4 = [0.743,0.835,0.921,1.005,1.076,1.134,1.180,1.214,1.241,...
%                  1.260,1.273,1.272,1.272,1.256,1.240,1.240,1.197,1.172,...
%                  1.147,1.122,1.096];
% spectrum_CuSO4 = [1.164756,1.18318,0.964564,1.051407,1.12858,1.195505,1.244105,1.272197,1.293193,1.34025,1.350587,...
%     1.331975,1.349069,1.325291,1.308238,1.282679,1.271985,1.245274,1.214709,1.179103,1.164756];

load( ['Substance_spectra\spectrumCuSO4_extin.mat' ] );
true_spectrum_Cu = spectrum_extin_CuSO4(select_wv);

%% 3. Initialize
wv_all = 700:10:900;
wv_values = wv_all(select_wv);
num_wavelengths = length(wv_values);
snr_Ni = zeros(1, num_wavelengths);
snr_Cu = zeros(1, num_wavelengths);
snr_Target = zeros(1, num_wavelengths);

%% 4. Calculate SNR
for idx = 1:num_wavelengths
    img = load(sprintf('0520/PA_deep_recon_%d.mat', wv_values(idx))).img;
    imshow(img, []);

    % Define neighborhood size
    row_offset = 4; % up/down ±2
    col_offset = 8; % left/right ±3

    % Get signals using local average
    region_Ni = img(point2(2)-row_offset:point2(2)+row_offset, point2(1)-col_offset:point2(1)+col_offset);
    signal_Ni = mean(region_Ni(:));

    region_Cu = img(point1(2)-row_offset:point1(2)+row_offset, point1(1)-col_offset:point1(1)+col_offset);
    signal_Cu = mean(region_Cu(:));

    % Compute SNR
    snr_Ni(idx) = signal_Ni;
    snr_Cu(idx) = signal_Cu;

    % Other points using local average
    region_60 = img(point3(2)-row_offset:point3(2)+row_offset, point3(1)-col_offset:point3(1)+col_offset);
    snr_60(idx) = mean(region_60(:));

    region_25 = img(point4(2)-row_offset:point4(2)+row_offset, point4(1)-col_offset:point4(1)+col_offset);
    snr_25(idx) = mean(region_25(:));

    region_80 = img(point5(2)-row_offset:point5(2)+row_offset, point5(1)-col_offset:point5(1)+col_offset);
    snr_80(idx) = mean(region_80(:));

    % region_80 = img(point6(2)-row_offset:point6(2)+row_offset, point6(1)-col_offset:point6(1)+col_offset);
    % snr_80(idx) = mean(region_80(:));
end


%% 5. Define Optimization Problem with Regularization

%% 优化补偿向量以拟合目标归一化谱形
normalize = @(x) x / norm(x);  % L2 norm
% normalize = @(x) (x - min(x)) / norm(x - min(x));


% 目标函数：拟合谱形 + 平滑惩罚项
objective = @(comp) ...
    norm(normalize(snr_Ni ./ comp) - normalize(true_spectrum_Ni))^2 + ...
    norm(normalize(snr_Cu ./ comp) - normalize(true_spectrum_Cu))^2 + ...
    0.001 * norm(diff(comp))^2;

% 初始补偿向量
init_comp = ones(size(snr_Ni));

% 上下界（防止除以0或过大）
lb = 0.1 * ones(size(snr_Ni));
ub = 10  * ones(size(snr_Ni));

% 优化设置
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 5000, ...
    'MaxIterations', 1000);

% 执行优化
opt_compensation = fmincon(objective, init_comp, [], [], [], [], lb, ub, [], options);

%% 6. Apply Compensation
corrected_Ni = snr_Ni ./ opt_compensation;
corrected_Cu = snr_Cu ./ opt_compensation;
corrected_Target = snr_Target./opt_compensation;

% corrected_0 = snr_0 ./ opt_compensation;
corrected_60 = snr_60 ./ opt_compensation;
corrected_25 = snr_25 ./ opt_compensation;
corrected_80 = snr_80 ./ opt_compensation;
% corrected_80 = snr_80 ./ opt_compensation;


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
measured_compensation = [6.69
8.22
9.193333333
9.823333333
9.706666667
10.08
9.973333333
9.9
9.843333333
9.573333333
];
addpath('toolbox');
% SO2_tube_linear = linearUnmixing(spectrum2(:), spectrum_HbO2(5:2:end), spectrum_Hb(5:2:end));
% SO2_tube_linear = linearUnmixing(corrected_0, true_spectrum_Ni,true_spectrum_Cu);
SO2_tube_linear = linearUnmixing(corrected_60, true_spectrum_Ni,true_spectrum_Cu)
SO2_tube_linear = linearUnmixing(corrected_25, true_spectrum_Ni,true_spectrum_Cu)
SO2_tube_linear = linearUnmixing(corrected_80, true_spectrum_Ni,true_spectrum_Cu)
% SO2_tube_linear = linearUnmixing(corrected_80, true_spectrum_Ni,true_spectrum_Cu)

save('phantom_deep_0520','corrected_25','corrected_60',"corrected_80");