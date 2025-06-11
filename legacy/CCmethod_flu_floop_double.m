clear;
close all;
%add the toolbox path

addpath('/Users/cassandrayang/Documents/GitHub/Faster_Convex_Cone/Legacy/CC_Set')
addpath('/Users/cassandrayang/Documents/GitHub/Faster_Convex_Cone/Toolbox')
load( ['/Users/cassandrayang/Documents/GitHub/Faster_Convex_Cone/Substance_Spectra/spectrum_Hb_Cope.mat' ] );
load( ['/Users/cassandrayang/Documents/GitHub/Faster_Convex_Cone/Substance_Spectra/spectrum_HbO2_Cope.mat' ] );
load( ['/Users/cassandrayang/Documents/GitHub/Faster_Convex_Cone/Substance_Spectra/spectrum_H2O.mat' ] );

sele_wavlen = [1:21];
wavelengths   = 700 : 20 : 900;

spectrum_Hb = spectrum_Hb(sele_wavlen);
spectrum_HbO2 = spectrum_HbO2(sele_wavlen);
spectrum_H2O = spectrum_H2O(sele_wavlen);

save_num = 0;
% % for save_num = 2
for phantom_num = 1:8
%Load Colorbase
load(sprintf('/home/zz111/Spectral_unmixing/cc_diffusion_theory/My_ColorBase_large_deep%d_VesselDiff_200_200', phantom_num), "My_ColorBase_s");
for para_num = 1:8
    save_num = save_num+1
% figure
% plot(wavelengths,spectrum_Hb/norm(spectrum_Hb));
% hold on
% plot(wavelengths,spectrum_HbO2/norm(spectrum_HbO2));
% hold on
% plot(wavelengths,spectrum_H2O/norm(spectrum_H2O));
% legend('Hb', 'Hb02', 'H2O');


My_ColorBase_s = double(full(My_ColorBase_s(1:end,sele_wavlen)));
%Load PA spectrum
i = 0;
for sBO2 = 0:0.1:1
    i = i+1;
%     filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_deep_nonuni_fluence_num%d_sBO20.6_dt.mat', save_num);
   filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_deepmore_nonuni_VesselDiff_fluence_num%d_dt.mat', save_num);
%     
%     %filename = sprintf('fluence_1e9_top_sTO2_%.1f_sO2_%.1f.mat', sTO2, sO2); 
% %     filename = sprintf('Simulation_size5/fluence_1e9_top_sTO2_%.1f_sBO2_0.6.mat',sTO2); 
   load(filename)
    blood_mu_a = (spectrum_HbO2 .* sBO2 + spectrum_Hb .* (1 - sBO2)) ...
            .* (150/64500) * log(10);
    SO2_cc_diff(i ,1 ) = convexConeSO2(fluence_spectrum1(sele_wavlen).*blood_mu_a,spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,2 ) = convexConeSO2(fluence_spectrum2(sele_wavlen).*blood_mu_a,spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,3 ) = convexConeSO2(fluence_spectrum3(sele_wavlen).*blood_mu_a,spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,4 ) = convexConeSO2(fluence_spectrum4(sele_wavlen).*blood_mu_a,spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,5 ) = convexConeSO2(fluence_spectrum5(sele_wavlen).*blood_mu_a,spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,6 ) = convexConeSO2(fluence_spectrum6(sele_wavlen).*blood_mu_a,spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,7 ) = convexConeSO2(fluence_spectrum7(sele_wavlen).*blood_mu_a,spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);
    SO2_cc_diff(i ,8 ) = convexConeSO2(fluence_spectrum8(sele_wavlen).*blood_mu_a,spectrum_Hb ,spectrum_HbO2, My_ColorBase_s);



end
sO2s = 0:0.1:1;
error_matrix = abs(SO2_cc_diff(:,:) - sO2s');
mean(error_matrix(:))
%% Compute mean absolute error for each depth level
error_matrix = abs(SO2_cc_diff(:,:) - sO2s');
mean_error = mean(error_matrix, 1);
std_error = std(error_matrix, 0, 1);
depths = [2 4 6 8 10 12 15 20];
% Create a bar plot with error bars
figure('Position', [100, 100, 800, 400]);
b = bar(depths, mean_error, 'FaceColor', [1, 0, 0]); % Set bar color to red
hold on;

% Add error bars
errorbar(depths, mean_error, std_error, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set x-axis ticks and labels
xticks(depths);
xticklabels(arrayfun(@num2str, depths, 'UniformOutput', false));

% Add labels and title
xlabel('Depth');
ylabel('Mean Prediction Error (Across 5 Vessels)');
title('Mean Error of Linear Unmixing Prediction with Standard Deviation non-uni');

grid on;
set(gcf, 'Color', 'w');
hold off;
save(sprintf('/home/zz111/Spectral_unmixing/Convex-Cone-Approach-main/sO2_predictionDT_deepmore_VesselDiff_MylDEEPcc_direct_savenum%d.mat',save_num),"SO2_cc_diff");

end
end