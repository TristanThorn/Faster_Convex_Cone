% clear all
% close all

wavelengths = linspace(700, 900, 21);
% cd '/home/zz111/Spectral_unmixing/Simulation_size5'
sTO2 = 0.9;
% for sBO2 = 0.5
sBO2 = 0.9;
num = 100;
for i = 1:length(wavelengths)
    i
    wavelength = wavelengths(i);
%     filename = sprintf('/home/zz111/Spectral_unmixing/generate_cc/cc_set/CC_MC1e8_0.05_wavelength%d_num_%d.mat', wavelength, num);
    seed = 0;
    filename = sprintf('/home/zz111/Spectral_unmixing/mcx_result/TestPhantom_spe_sTO2%.1f_sBO2%.1f_wavelength%d_seed%d.mat',sTO2, sBO2, wavelength, seed);
    load(filename);

%     seed = 1;
%     filename = sprintf('/home/zz111/Spectral_unmixing/Simulation_size5/check_noise1e9_wavelength%d_seed%d.mat', wavelength, seed);
%     load(filename);
%     seed1 = squeeze(mcxdata(45, :,:));
%     dim = size(mcxdata);
%     error = abs((seed1-seed0)./seed0);
%     yi = ((1:dim(2)-25)) * 0.2;
%     zi = (1:dim(3)-25) * 0.2;
%     error(error>0.1) = 1;
%     imagesc(yi, zi, (abs(squeeze(error(:, 25:end))))');

%      figure;imagesc((log(abs(squeeze(mcxdata(2, :, :)))))');
%     imagesc(log(abs(squeeze(mcxdata(400, :, :))))');
 
% % origional
%      fluence_spectrum1(i) = mcxdata(2,36,35);
%     fluence_spectrum2(i) = mcxdata(2,65,40);
%     fluence_spectrum3(i) = mcxdata(2,95,45);
%     fluence_spectrum4(i) = mcxdata(2,114,52);

 % spe
     fluence_spectrum1(i) = mcxdata(2,35,33);
    fluence_spectrum2(i) = mcxdata(2,95,38);
    fluence_spectrum3(i) = mcxdata(2,75,43);
    fluence_spectrum4(i) = mcxdata(2,115,48);
    fluence_spectrum5(i) = mcxdata(2,50,58); 



% 4cm
%     fluence_spectrum1(i) = mcxdata(400,200,120);
%     fluence_spectrum2(i) = mcxdata(400,540,200);
%     fluence_spectrum3(i) = mcxdata(400,360,300);
%     fluence_spectrum4(i) = mcxdata(400,600,400);
%     fluence_spectrum5(i) = mcxdata(400,520,500);
%     fluence_spectrum6(i) = mcxdata(400,240,600);
%     fluence_spectrum7(i) = mcxdata(400,400,700);

% mimic CC
%     fluence_spectrum1(i) = mcxdata(200,240,140);
%     fluence_spectrum2(i) = mcxdata(200,480,180);
%     fluence_spectrum3(i) = mcxdata(200,360,160);
%     fluence_spectrum4(i) = mcxdata(200,560,200);

end
filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_spe_fluence_1e9_sTO2%.1f_sBO2%.1f.mat', sTO2, sBO2)
% save(filename, 'fluence_spectrum1','fluence_spectrum2','fluence_spectrum3','fluence_spectrum4','fluence_spectrum5', ...
%     'fluence_spectrum6','fluence_spectrum7','fluence_spectrum8','fluence_spectrum9','fluence_spectrum10');

save(filename, 'fluence_spectrum1','fluence_spectrum2','fluence_spectrum3','fluence_spectrum4','fluence_spectrum5');

% end


%% show fluence
wavelengths = linspace(700, 900, 21); 
figure;
set(gcf, 'Position', [100, 100, 600, 400]);
fluence_spectrum4 = fluence_spectrum4/(max(fluence_spectrum4(:)));
plot(wavelengths, fluence_spectrum4, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);hold on;
xlabel('Wavelength (nm)', 'FontSize', 12);
ylabel('Residual fluence', 'FontSize', 12);
% title(sprintf('Residual fluence ratio vs. Wavelength, sTO_2 = %.1f(MC-10 mm)', sTO2), 'FontSize', 8);
% title(sprintf('Residual fluence ratio(MC-5.4 mm)'), 'FontSize', 8);
title(sprintf('Test'), 'FontSize', 8);
grid on;
xticks(min(xlim):10:max(xlim));
xline(730, '--k', 'LineWidth', 1 ,'Color','r');
xline(760, '--k', 'LineWidth', 1 ,'Color','r');
xline(800, '--k', 'LineWidth', 1 ,'Color','r');
set(gca, 'FontSize', 12,'Color', 'w');
set(gcf, 'Color', 'w');
% saveas(gcf, filename);