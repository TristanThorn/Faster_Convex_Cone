%% 
clear all
close all;
%add the toolbox path

addpath('/home/zz111/Spectral_unmixing/generate_cc/cc_set/')
addpath('/home/zz111/Spectral_unmixing/Convex-Cone-Approach-main/Toolbox')
%Load Colorbase
load( 'Human_experiment/MC_light_model/ColorBase_s/ColorBase_s.mat' );
% load( 'Simulation/MC_light_model/ColorBase_s/ColorBase_s.mat' );
% load('/home/zz111/Spectral_unmixing/generate_cc/My_ColorBase_s_large.mat')
% load('/home/zz111/Spectral_unmixing/cc_diffusion_theory/subset_cita_distance0.005.mat')
% load('/home/zz111/Spectral_unmixing/cc_diffusion_theory/My_ColorBase_large_shrink_comp_200_200_3dt.mat')
% load('/home/zz111/Spectral_unmixing/cc_diffusion_theory/My_ColorBase_large_shrink_spe_200_200_1.mat');
load( ['Substance_spectra/spectrum_Hb_Cope.mat' ] );
load( ['Substance_spectra/spectrum_HbO2_Cope.mat' ] );
load( ['Substance_spectra/spectrum_H2O.mat' ] );

sele_wavlen = [1:21];
wavelengths   = 700 : 10 : 900;
figure
plot(wavelengths,spectrum_Hb/norm(spectrum_Hb));
hold on
plot(wavelengths,spectrum_HbO2/norm(spectrum_HbO2));
hold on
plot(wavelengths,spectrum_H2O/norm(spectrum_H2O));
legend('Hb', 'Hb02', 'H2O');

sTO2 = 0.8;
sVO2 = 0.8;

My_ColorBase_s = ColorBase_s(:,:);

% My_ColorBase_s = double(subsetData(:,:));
% fluence_spectrum1 = My_ColorBase_s(1,:);
% % 
% fluence_spectrum2 = My_ColorBase_s(2,:);
% % figure;plot(wavelengths,fluence_spectrum2);
% fluence_spectrum3 = My_ColorBase_s(3,:);
% fluence_spectrum4 = My_ColorBase_s(4,:);
% fluence_spectrum5 = My_ColorBase_s(5,:);
My_ColorBase_s = double(full(My_ColorBase_s(1:end,:)));
%Load PA spectrum
i = 0;
for sBO2 = 0:0.1:1
    i = i+1;
%     filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_deep_fluence_sTO2%.1f_sBO2%.1f_dt.mat', sTO2, sVO2);
    filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_deep_fluence_sTO2%.1f_sBO2%.1f_dt.mat', sTO2, sVO2);
%     %filename = sprintf('fluence_1e9_top_sTO2_%.1f_sO2_%.1f.mat', sTO2, sO2); 
% %     filename = sprintf('Simulation_size5/fluence_1e9_top_sTO2_%.1f_sBO2_0.6.mat',sTO2); 
   load(filename)
    blood_mu_a = (spectrum_HbO2 .* sBO2 + spectrum_Hb .* (1 - sBO2)) ...
            .* (150/64500) * log(10);
    SO2_cc_diff(i ,1 ) = convexConeSO2(fluence_spectrum1(sele_wavlen).*blood_mu_a(sele_wavlen),spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,2 ) = convexConeSO2(fluence_spectrum2(sele_wavlen).*blood_mu_a(sele_wavlen),spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
     SO2_cc_diff(i ,3 ) = convexConeSO2(fluence_spectrum3(sele_wavlen).*blood_mu_a(sele_wavlen),spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
     SO2_cc_diff(i ,4 ) = convexConeSO2(fluence_spectrum4(sele_wavlen).*blood_mu_a(sele_wavlen),spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,5 ) = convexConeSO2(fluence_spectrum5(sele_wavlen).*blood_mu_a(sele_wavlen),spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,6 ) = convexConeSO2(fluence_spectrum6(sele_wavlen).*blood_mu_a(sele_wavlen),spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,7 ) = convexConeSO2(fluence_spectrum7(sele_wavlen).*blood_mu_a(sele_wavlen),spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    
%     SO2_cc_diff(i ,1 ) = convexConeSO2(fluence_spectrum1(sele_wavlen),ones(1,21) ,ones(1,21), My_ColorBase_s);
%     SO2_cc_diff(i ,2 ) = convexConeSO2(fluence_spectrum2(sele_wavlen),ones(1,21) ,ones(1,21), My_ColorBase_s);
%      SO2_cc_diff(i ,3 ) = convexConeSO2(fluence_spectrum3(sele_wavlen),ones(1,21) ,ones(1,21), My_ColorBase_s);
%      SO2_cc_diff(i ,4 ) = convexConeSO2(fluence_spectrum4(sele_wavlen),ones(1,21) ,ones(1,21), My_ColorBase_s);
%     SO2_cc_diff(i ,5 ) = convexConeSO2(fluence_spectrum5(sele_wavlen),ones(1,21) ,ones(1,21), My_ColorBase_s);


end

%% Show
sO2s = 0:0.1:1; % Generate 10 points
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'}; % 10 different shapes
% Create a scatter plot
figure; hold on;
for i = 1:7
    scatter(sO2s, SO2_cc_diff(:, i), 60, markers{i});
end

% Plot the reference line y = x
plot(sO2s, sO2s, 'k--', 'LineWidth', 2);

% Add labels and title
xlabel('True sO2');
ylabel('Predicted sO2 (Convex cone)');
title(sprintf('CC method flu sO2, sTO2 = %.1f, 21 wavelength',sTO2));

% Add legend
legend({'Vessel 1 - 1.5mm', 'Vessel 2 - 2.5mm', 'Vessel 3 - 3.5mm', 'Vessel 4 - 4.5mm', 'Vessel 5 - 6.5'}, 'Location', 'Best');

set(gcf, 'Color', 'w');
hold off;
save(sprintf('/home/zz111/Spectral_unmixing/Convex-Cone-Approach-main/sO2_predictionDT_deep_theirHcc_direct_%.1f_%.1f.mat',sTO2,sVO2),"SO2_cc_diff");
%% Compute mean absolute error for each sO2 level
error_matrix = abs(SO2_cc_diff(:,:) - sO2s');
mean_error = mean(error_matrix, 2);
std_error = std(error_matrix, 0, 2);

% Create a bar plot with error bars
figure('Position', [100, 100, 800, 400]);
b = bar(sO2s, mean_error, 'FaceColor', [1, 0, 0]); % Set bar color to red
hold on;

% Add error bars
errorbar(sO2s, mean_error, std_error, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set x-axis ticks and labels
xticks(sO2s);
xticklabels(arrayfun(@num2str, sO2s, 'UniformOutput', false));

% Add labels and title
xlabel('True sO2');
ylabel('Mean Prediction Error (Across 5 Vessels)');
title(sprintf('Mean Error of Linear Unmixing Prediction with Standard Deviation, sTO2 = %.1f',sTO2));

grid on;
set(gcf, 'Color', 'w');
hold off;