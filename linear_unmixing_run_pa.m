clear all;
% close all;
i = 0;
load( ['spectrum_Hb_Cope.mat' ] );
load( ['spectrum_HbO2_Cope.mat' ] );

sele_wavlen = [1:21];
sTO2 = 0.6;
% Fluence_spectrum_array = ColorBase_s;
% fluence_spectrum1 = Fluence_spectrum_array(100,:) ;
% fluence_spectrum2 = Fluence_spectrum_array(2000,:) ;
% fluence_spectrum3 = Fluence_spectrum_array(300,:) ;
% fluence_spectrum4 = Fluence_spectrum_array(40,:) ;
% fluence_spectrum5 = Fluence_spectrum_array(500,:) ;
for sBO2 = 0:0.1:1
    i = i+1
    %filename = sprintf('fluence_sO2_%.1f.mat', sO2);
    %filename = sprintf('fluence_1e9_top_sTO2_%.1f_sO2_%.1f.mat', sTO2, sO2); 
    filename = sprintf('pa_spectrum/PA_1e9_top_sTO2_%.1f_sBO2_%.1f.mat', sTO2, sBO2); 
    load(filename)
    SO2_linear_unmixing_diff(i ,1 ) = linearUnmixing( pa_spectrum1(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    SO2_linear_unmixing_diff(i ,2 ) = linearUnmixing( pa_spectrum2(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    SO2_linear_unmixing_diff(i ,3 ) = linearUnmixing( pa_spectrum3(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    SO2_linear_unmixing_diff(i ,4 ) = linearUnmixing( pa_spectrum4(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    SO2_linear_unmixing_diff(i ,5 ) = linearUnmixing( pa_spectrum5(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    SO2_linear_unmixing_diff(i ,6 ) = linearUnmixing( pa_spectrum6(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    SO2_linear_unmixing_diff(i ,7 ) = linearUnmixing( pa_spectrum7(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    SO2_linear_unmixing_diff(i ,8 ) = linearUnmixing( pa_spectrum8(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    SO2_linear_unmixing_diff(i ,9 ) = linearUnmixing( pa_spectrum9(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));
    SO2_linear_unmixing_diff(i ,10 ) = linearUnmixing( pa_spectrum10(sele_wavlen), spectrum_HbO2(sele_wavlen), spectrum_Hb(sele_wavlen));

end


sO2s = 0:0.1:1; % Generate 10 points
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'}; % 10 different shapes
% Create a scatter plot
figure; hold on;
for i = 1:1
    i
    scatter(sO2s, SO2_linear_unmixing_diff(:, i), 60, markers{i});
end

% Plot the reference line y = x
plot(sO2s, sO2s, 'k--', 'LineWidth', 2);

% Add labels and title
xlabel('True sO2');
ylabel('Predicted sO2 (Linear Unmixing)');
title(sprintf('Linear Unmixing PA sO2, sTO2 = %.1f, 21 wavelength',sTO2));

% Add legend
legend({'Vessel 1 - 1mm', 'Vessel 2 - 2mm', 'Vessel 3 - 3mm', 'Vessel 4 - 4mm', ...
        'Vessel 5 - 5mm', 'Vessel 6 - 6mm', 'Vessel 7 - 7mm', 'Vessel 8 - 8mm', ...
        'Vessel 9 - 9mm', 'Vessel 10 - 10mm', 'Ground Truth'}, 'Location', 'Best');
% legend({'Vessel 1 - 1mm', 'Vessel 2 - 4mm', 'Vessel 3 - 7mm', 'Vessel 4 - 10mm', 'Ground Truth'}, 'Location', 'Best');

set(gcf, 'Color', 'w');
hold off;

