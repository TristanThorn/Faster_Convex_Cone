clear all;
close all;
addpath('pa_result')
wavelengths = linspace(700, 900, 21);
cd '/home/zz111/Spectral_unmixing'
sTO2 = 0.4;
for sBO2 = 0:0.1:1
for i = 1:length(wavelengths)
    wavelength = wavelengths(i);
    filename = sprintf('pa_result/PA1e9_4cm_sTO2_%.1f_wavelength%d_sBO2_%.1f.mat', sTO2, wavelength, sBO2);
    
    % filename = sprintf('MC_wavelength%d_sO2_%.1f.mat', wavelength, sO2);
    load(filename);
%      figure; imshow(p0_recon,[]);

    pa_spectrum1(i) = p0_recon(110,250);
    pa_spectrum2(i) = p0_recon(130,590);
    pa_spectrum3(i) = p0_recon(150,410);
    pa_spectrum4(i) = p0_recon(170,650);
    pa_spectrum5(i) = p0_recon(190,570);
    pa_spectrum6(i) = p0_recon(210,290);
    pa_spectrum7(i) = p0_recon(230,450);
    pa_spectrum8(i) = p0_recon(250,250);
    pa_spectrum9(i) = p0_recon(270,570);
    pa_spectrum10(i) = p0_recon(290,490);
%     fluence_spectrum4(i) = mcxdata(100,260,125);
%     fluence_spectrum5(i) = mcxdata(100,180,130);

        %     fluence_spectrum1(i) = mcxdata(100,170,75);
        % fluence_spectrum2(i) = mcxdata(100,50,100);
        % fluence_spectrum3(i) = mcxdata(100,260,115);
        % fluence_spectrum4(i) = mcxdata(100,210,140);
        % fluence_spectrum5(i) = mcxdata(100,130,150);
end
filename = sprintf('pa_spectrum/PA_1e9_top_sTO2_%.1f_sBO2_%.1f.mat', sTO2, sBO2); 
save(filename, 'pa_spectrum1', 'pa_spectrum2', 'pa_spectrum3','pa_spectrum4', 'pa_spectrum5', ...
    'pa_spectrum6','pa_spectrum7', 'pa_spectrum8', 'pa_spectrum9','pa_spectrum10');
end