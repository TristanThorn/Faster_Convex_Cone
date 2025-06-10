% 1. 加载参考谱线（假设你已经load了这两个变量）
load( ['Substance_spectra\spectrumNiSO4_extin.mat' ] );
load( 'Substance_spectra/spectrumCuSO4_extin.mat');

% 2. 自定义测量谱线
spectrum_CuSO4 = [1.164756,1.18318,0.964564,1.051407,1.12858,1.195505,1.244105,1.272197,1.293193,1.34025,1.350587,...
    1.331975,1.349069,1.325291,1.308238,1.282679,1.271985,1.245274,1.214709,1.179103,1.164756];

spectrum_NiSO4 = [1.763306,1.747571,1.746373,1.724271,1.69942,1.706706,1.6967,1.649732,1.652714,...
    1.534028,1.416467,1.278659,1.111068,0.943837,0.80629,0.696053,0.612729,0.554713,0.516803,0.502615,0.508012];

% 3. 波长（如果参考谱线是等间隔取样，默认起始720nm）
wavelengths = 700:10:900;

% 4. 定义L2范数归一化函数
normalize = @(x) x / norm(x);

% 5. 归一化所有谱线
spectrum_CuSO4_norm       = normalize(spectrum_CuSO4(3:end));
spectrum_extin_CuSO4_norm = normalize(spectrum_extin_CuSO4(3:end));

spectrum_NiSO4_norm       = normalize(spectrum_NiSO4(3:end));
spectrum_extin_NiSO4_norm = normalize(spectrum_extin_NiSO4(3:end));

% 6. 绘图
figure;

% NiSO4
subplot(2,1,1);
plot(wavelengths(3:end), spectrum_NiSO4_norm, '-bo', 'LineWidth', 2); hold on;
plot(wavelengths(3:end), spectrum_extin_NiSO4_norm, '--k*', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Normalized Amplitude');
title('NiSO₄ Spectrum Comparison (L2 Norm)');
legend('Measured NiSO₄ (norm)', 'Reference NiSO₄ (norm)', 'Location', 'Best');
grid on;

% CuSO4
subplot(2,1,2);
plot(wavelengths(3:end), spectrum_CuSO4_norm, '-ro', 'LineWidth', 2); hold on;
plot(wavelengths(3:end), spectrum_extin_CuSO4_norm, '--k*', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Normalized Amplitude');
title('CuSO₄ Spectrum Comparison (L2 Norm)');
legend('Measured CuSO₄ (norm)', 'Reference CuSO₄ (norm)', 'Location', 'Best');
grid on;
