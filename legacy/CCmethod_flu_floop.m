clear all
close all;
%add the toolbox path

addpath('/home/zz111/Spectral_unmixing/generate_cc/cc_set/')
addpath('/home/zz111/Spectral_unmixing/Convex-Cone-Approach-main/Toolbox')
load( ['Substance_spectra/spectrum_Hb_Cope.mat' ] );
load( ['Substance_spectra/spectrum_HbO2_Cope.mat' ] );
load( ['Substance_spectra/spectrum_H2O.mat' ] );

for seed = 1:32
%Load Colorbase
load(sprintf('/home/zz111/Spectral_unmixing/cc_diffusion_theory/My_ColorBase_large_shrink_comp_200_200_%d.mat',seed));

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


My_ColorBase_s = double(full(My_ColorBase_s(1:end,:)));
%Load PA spectrum
i = 0;
for sBO2 = 0:0.1:1
    i = i+1;
    filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_deep_fluence_sTO2%.1f_sBO2%.1f_dt.mat', sTO2, sVO2);
%     
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
    


end
sO2s = 0:0.1:1;
error_matrix = abs(SO2_cc_diff(:,:) - sO2s');
mean(error_matrix(:))
save(sprintf('/home/zz111/Spectral_unmixing/Convex-Cone-Approach-main/sO2_predictionDT_deep_MylCOMcc_direct%d_%.1f_%.1f.mat',seed,sTO2,sVO2),"SO2_cc_diff");

end