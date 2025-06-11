% %% 
% clear all
% %add the toolbox path
% addpath('/home/zz111/Spectral_unmixing/Convex-Cone-Approach-main')
% addpath('/home/zz111/Spectral_unmixing/Convex-Cone-Approach-main/Toolbox')
% %Load Colorbase
% load( 'Simulation/MC_light_model/ColorBase_s/ColorBase_s.mat' );
% load( ['Substance_spectra/spectrum_Hb_Cope.mat' ] );
% load( ['Substance_spectra/spectrum_HbO2_Cope.mat' ] );
% load( ['Substance_spectra/spectrum_H2O.mat' ] );
% 
% wavelengths   = 700 : 10 : 900;
% figure
% plot(wavelengths,spectrum_Hb/norm(spectrum_Hb));
% hold on
% plot(wavelengths,spectrum_HbO2/norm(spectrum_HbO2));
% hold on
% plot(wavelengths,spectrum_H2O/norm(spectrum_H2O));
% legend('Hb', 'Hb02', 'H2O');
% 
% sTO2 = 0.4;
% %Load PA spectrum
% i = 0;
% for sBO2 = 0:0.1:1
%     i = i+1
%     %filename = sprintf('fluence_sO2_%.1f.mat', sO2);
%     %filename = sprintf('fluence_1e9_top_sTO2_%.1f_sO2_%.1f.mat', sTO2, sO2); 
%     filename = sprintf('/home/zz111/Spectral_unmixing/pa_spectrum/PA_1e9_top_sTO2_%.1f_sBO2_%.1f.mat', sTO2, sBO2); 
%     load(filename)
%     SO2_cc_diff(i ,1 ) = convexConeSO2(gather(pa_spectrum1),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
%     SO2_cc_diff(i ,2 ) = convexConeSO2(gather(pa_spectrum2),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
%     SO2_cc_diff(i ,3 ) = convexConeSO2(gather(pa_spectrum3),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
%     SO2_cc_diff(i ,4 ) = convexConeSO2(gather(pa_spectrum4),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
%     SO2_cc_diff(i ,5 ) = convexConeSO2(gather(pa_spectrum5),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
%     SO2_cc_diff(i ,6 ) = convexConeSO2(gather(pa_spectrum6),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
%     SO2_cc_diff(i ,7 ) = convexConeSO2(gather(pa_spectrum7),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
%     SO2_cc_diff(i ,8 ) = convexConeSO2(gather(pa_spectrum8),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
%     SO2_cc_diff(i ,9 ) = convexConeSO2(gather(pa_spectrum9),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
%     SO2_cc_diff(i ,10 ) = convexConeSO2(gather(pa_spectrum10),spectrum_Hb ,spectrum_HbO2, ColorBase_s);
% 
% end


sO2s = 0:0.1:1; % Generate 10 points
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'}; % 10 different shapes
% Create a scatter plot
figure; hold on;
for i = 5:10
    i
    scatter(sO2s, SO2_cc_diff(:, i), 60, markers{i});
end

% Plot the reference line y = x
plot(sO2s, sO2s, 'k--', 'LineWidth', 2);

% Add labels and title
xlabel('True sO2');
ylabel('Predicted sO2 (Convex cone)');
title(sprintf('CC method PA sO2, sTO2 = %.1f, 21 wavelength',sTO2));

% Add legend
% legend({'Vessel 1 - 1mm', 'Vessel 2 - 2mm', 'Vessel 3 - 3mm', 'Vessel 4 - 4mm', 'Ground Truth',...
     legend({'Vessel 5 - 5mm', 'Vessel 6 - 6mm', 'Vessel 7 - 7mm', 'Vessel 8 - 8mm', ...
        'Vessel 9 - 9mm', 'Vessel 10 - 10mm', 'Ground Truth'}, 'Location', 'Best');
% legend({'Vessel 1 - 1mm', 'Vessel 2 - 4mm', 'Vessel 3 - 7mm', 'Vessel 4 - 10mm', 'Ground Truth'}, 'Location', 'Best');

set(gcf, 'Color', 'w');
hold off;

