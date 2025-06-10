load('Substance_spectra/spectrum_HbO2_Cope.mat');
load('Substance_spectra/spectrum_Hb_Cope.mat');

sele_wavlen = 1:21;
wavelengths = 700 : 10 : 900;

figure;
plot(wavelengths, spectrum_Hb, 'r', 'LineWidth', 2.5); hold on;   % Hb in red
plot(wavelengths, spectrum_HbO2, 'b', 'LineWidth', 2.5);          % HbO₂ in blue

% Format axes and labels
xlabel('Wavelength (nm)', 'FontSize', 12);                    % X-axis label
ylabel('Absorption coefficient', 'FontSize', 12);             % Y-axis label
set(gca, 'YTick', []);                                        % Remove Y-axis tick values
set(gcf, 'Color', 'w');                                       % White background
legend({'Hb', 'HbO₂'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
axis tight;

% Create a colormap from blue to red
cmap = [linspace(0, 1, 5); zeros(1, 5); linspace(1, 0, 5)]'; % Blue to Red
index = 1;
for sBO2 = 0.75
    % Calculate blood absorption coefficient (mu_a)
    blood_mu_a = (spectrum_HbO2 .* sBO2 + spectrum_Hb .* (1 - sBO2)) ...
               .* (150 / 64500) * log(10);

    fluence_spec = 1./((spectrum_HbO2 .* 0.2 + spectrum_Hb .* (1 - 0.2)) ...
               .* (150 / 64500) * log(10) .* (spectrum_HbO2 .* 0.9 + spectrum_Hb .* (1 - 0.9)) ...
               .* (150 / 64500) * log(10)).^0.2;

    figure('Color', 'w', 'Position', [100, 100, 600, 300]);
    plot(wavelengths, fluence_spec, 'LineWidth', 5, 'Color', cmap(index, :)); % Set color based on index
    axis tight;
    axis off; % Hide axes

    index = index + 1;
end

%% Phantom
% Load spectra
% Load spectra
load('Substance_spectra\spectrumNiSO4_extin.mat');
load('Substance_spectra\spectrumCuSO4_extin.mat');

% Define wavelength range (700–900 nm, step 10 nm)
wavelengths = 700:10:900;

% Select wavelengths (adjust as needed)
select_wv = 1:length(wavelengths);

% Extract and scale spectra
true_spectrum_Ni = spectrum_extin_NiSO4(select_wv) * 12.7;
true_spectrum_Cu = spectrum_extin_CuSO4(select_wv);

% Plot
figure;
plot(wavelengths, true_spectrum_Ni, 'r', 'LineWidth', 3); hold on;
plot(wavelengths, true_spectrum_Cu, 'b', 'LineWidth', 3);

% Format appearance
xlabel('Wavelength (nm)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Molar extinction coefficient (M^{-1}·cm^{-1})', ...
       'FontSize', 20, 'FontWeight', 'bold');

set(gca, 'FontSize', 20, 'FontWeight', 'bold');  % Axis ticks bold
legend({'NiSO₄ × 12.7', 'CuSO₄'}, ...
       'Location', 'northeast', 'Box', 'off', ...
       'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'tex');

set(gcf, 'Color', 'w');  % White background
axis tight;
