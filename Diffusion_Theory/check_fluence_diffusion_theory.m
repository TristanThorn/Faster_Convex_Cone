clear all
close all

wavelengths = linspace(700, 900, 21);
% cd '/home/zz111/Spectral_unmixing/Simulation_size5'
sTO2 = 0.6;
% for sBO2 = 0.5
sBO2 = 0.5;
num = 1;
for i = 1:length(wavelengths)
    i
    wavelength = wavelengths(i);
    file_name = sprintf('nirfast_result/nirfast_wavelength%d_num%d',wavelength,num);
    load(file_name);
    fluence_spectrum4(i) = phi_amplitude(10);


end
% filename = sprintf('TestPhantom_fluence_1e9_sTO2%.1f_sBO2%.1f.mat', sTO2, sBO2); 
% % save(filename, 'fluence_spectrum1','fluence_spectrum2','fluence_spectrum3','fluence_spectrum4','fluence_spectrum5', ...
% %     'fluence_spectrum6','fluence_spectrum7','fluence_spectrum8','fluence_spectrum9','fluence_spectrum10');
% 
% save(filename, 'fluence_spectrum1','fluence_spectrum2','fluence_spectrum3','fluence_spectrum4');

% end


%% show fluence
wavelengths = linspace(700, 900, 21); 
% load('fluence_top_sO2_0.8.mat')
figure;
set(gcf, 'Position', [100, 100, 600, 400]);
plot(wavelengths, fluence_spectrum4, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);hold on;
xlabel('Wavelength (nm)', 'FontSize', 12);
ylabel('Residual fluence', 'FontSize', 12);
% title(sprintf('Residual fluence ratio vs. Wavelength, sTO_2 = %.1f(MC-10 mm)', sTO2), 'FontSize', 8);
title(sprintf('Residual fluence ratio(DT)'), 'FontSize', 8);
grid on;
xticks(min(xlim):10:max(xlim));
xline(730, '--k', 'LineWidth', 1 ,'Color','r');
xline(760, '--k', 'LineWidth', 1 ,'Color','r');
xline(800, '--k', 'LineWidth', 1 ,'Color','r');
set(gca, 'FontSize', 12,'Color', 'w');
set(gcf, 'Color', 'w');
% saveas(gcf, filename);