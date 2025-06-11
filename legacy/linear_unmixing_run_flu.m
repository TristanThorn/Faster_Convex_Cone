clear all;
close all;
cd '/home/zz111/Spectral_unmixing/'
i = 0;

load( ['spectrum_Hb_Cope.mat' ] );
load( ['spectrum_HbO2_Cope.mat' ] );
sele_wavlen = [1:21];
sTO2 = 0.8;

wavelengths = linspace(700, 900, 21);
wave_combinations = nchoosek(1:21, 2);
num_combinations = size(wave_combinations, 1); % 计算组合总数
output_folder = 'lu_result';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

for num = 1
%     a = wave_combinations(num, 1);
%     b = wave_combinations(num, 2);
    SO2_linear_unmixing_diff = zeros(11,10);
    i = 0;
    for sBO2 = 0:0.1:1
        i = i+1;
    
        %filename = sprintf('fluence_sO2_%.1f.mat', sO2);
        %filename = sprintf('fluence_1e9_top_sTO2_%.1f_sO2_%.1f.mat', sTO2, sO2); 
        filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_fluence_1e9_sTO2%.1f_sBO2%.1f.mat', sTO2, 0.8); 
        blood_mu_a = (spectrum_HbO2 .* sBO2 + spectrum_Hb .* (1 - sBO2)) ...
            .* (150/64500) * log(10);
        load(filename)
        SO2_linear_unmixing_diff(i ,1 ) = linearUnmixing( fluence_spectrum1(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
        SO2_linear_unmixing_diff(i ,2 ) = linearUnmixing( fluence_spectrum2(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
        SO2_linear_unmixing_diff(i ,3 ) = linearUnmixing( fluence_spectrum3(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
        SO2_linear_unmixing_diff(i ,4 ) = linearUnmixing( fluence_spectrum4(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
%         SO2_linear_unmixing_diff(i ,5 ) = linearUnmixing( fluence_spectrum5(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
%         SO2_linear_unmixing_diff(i ,6 ) = linearUnmixing( fluence_spectrum6(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
%         SO2_linear_unmixing_diff(i ,7 ) = linearUnmixing( fluence_spectrum7(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
%         SO2_linear_unmixing_diff(i ,8 ) = linearUnmixing( fluence_spectrum8(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
%         SO2_linear_unmixing_diff(i ,9 ) = linearUnmixing( fluence_spectrum9(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
%         SO2_linear_unmixing_diff(i ,10 ) = linearUnmixing( fluence_spectrum10(sele_wavlen).*blood_mu_a(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    % figure;
    % plot(wavelengths, fluence_spectrum4, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);hold on;
    % xlabel('Wavelength (nm)', 'FontSize', 12);
    % ylabel('Residual fluence', 'FontSize', 12);
    end
    
    
    sO2s = 0:0.1:1; % Generate 10 points
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'}; % 10 different shapes
    %markers = {'o', 's', 'd', '^', 'v', '>', '<'}; % 10 different shapes
    % Create a scatter plot
    figure; hold on;
    for i = 1:4
        i
        scatter(sO2s, SO2_linear_unmixing_diff(:, i) ,60, markers{i});
    end
    
    % Plot the reference line y = x
    plot(sO2s, sO2s, 'k--', 'LineWidth', 2);
    
    % Add labels and title
    xlabel('True sO2');
    ylabel('Predicted sO2 (Linear Unmixing)');
    %title(sprintf('Linear Unmixing sO2, sTO2 = %.1f, wavelength = %d,%d',sTO2,wavelengths(a),wavelengths(b)));
    title(sprintf('Linear Unmixing sO2, sTO2 = %.1f, 21 wavelength ',sTO2));
    
    % Add legend
    legend({'Vessel 1 - 2mm', 'Vessel 2 - 3mm', 'Vessel 3 - 4mm', 'Vessel 4 - 5.4mm', ...
            'Ground Truth'}, 'Location', 'Best');

%     legend({'Vessel 1 - 1mm', 'Vessel 2 - 5mm', 'Vessel 3 - 10mm', 'Vessel 4 - 15mm', 'Vessel 5 - 20mm','Vessel 6 - 25mm' 'Ground Truth'}, 'Location', 'Best');
    
    set(gcf, 'Color', 'w');
    hold off;
%     filename = sprintf('LU_sTO2_%.1f_wavelength_%d_%d.png', sTO2, wavelengths(a), wavelengths(b));
%     filepath = fullfile(output_folder, filename);
% 
%     saveas(gcf, filepath)
end

%% Compute mean absolute error for each sO2 level
error_matrix = abs(SO2_linear_unmixing_diff(:,:) - sO2s');
mean_error = mean(error_matrix, 2);
std_error = std(error_matrix, 0, 2);

% Create a bar plot with error bars
figure; 
b = bar(sO2s, mean_error, 'FaceColor', [1, 0, 0]); % Set bar color to red
hold on;

% Add error bars
errorbar(sO2s, mean_error, std_error, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set x-axis ticks and labels
xticks(sO2s);
xticklabels(arrayfun(@num2str, sO2s, 'UniformOutput', false));

% Add labels and title
xlabel('True sO2');
ylabel('Mean Prediction Error (Across 10 Vessels)');
title(sprintf('Mean Error of Linear Unmixing Prediction with Standard Deviation, sTO2 = %.1f',sTO2));

grid on;
set(gcf, 'Color', 'w');
hold off;

%% Compute mean absolute error for each depth level
depths = 1:10;
error_matrix = abs(SO2_linear_unmixing_diff(:,:) - sO2s');
mean_error = mean(error_matrix, 1);
std_error = std(error_matrix, 0, 1);

% Create a bar plot with error bars
figure; 
b = bar(depths, mean_error, 'FaceColor', [0, 1, 0]); % Set bar color to red
hold on;

% Add error bars
% errorbar(depths, mean_error, std_error, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set x-axis ticks and labels
xticks(depths);
xticklabels(arrayfun(@num2str, depths, 'UniformOutput', false));

% Add labels and title
xlabel('Depth(mm)');
ylabel('Mean Prediction Error (Across 10 Vessels)');
title(sprintf('Mean Error of Linear Unmixing Prediction with Standard Deviation, sTO2 = %.1f',sTO2));

grid on;
set(gcf, 'Color', 'w');
hold off;

%% show
%vol = cfg.vol;
%volumeViewer(vol);
% load('fluence_top_sO2_0.8.mat')
wavelengths = linspace(700, 900, 21); 
sBO2 = 0.8;
figure;
plot(wavelengths, fluence_spectrum4, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);hold on;
xlabel('Wavelength (nm)', 'FontSize', 12);
ylabel('Residual fluence', 'FontSize', 12);
title(sprintf('Residual fluence ratio vs. Wavelength, sTO_2 = %.1f(MC-v4)', sTO2), 'FontSize', 14);
grid on;
xticks(min(xlim):10:max(xlim));
xline(730, '--k', 'LineWidth', 1 ,'Color','r');
xline(760, '--k', 'LineWidth', 1 ,'Color','r');
xline(800, '--k', 'LineWidth', 1 ,'Color','r');
set(gca, 'FontSize', 12,'Color', 'w');
set(gcf, 'Color', 'w');
% saveas(gcf, filename);