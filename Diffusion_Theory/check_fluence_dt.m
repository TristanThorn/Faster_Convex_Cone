clear all;
% close all;
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
cd '/home/zz111/Spectral_unmixing/cc_diffusion_theory'
rng = 42;%42
My_ColorBase_dt = zeros(10000,21);

filename = sprintf('nirfast_result/nirfast_0.8_spe_wavelength%d_num%d',700,0);
load(filename);

mesh = load_mesh('digital_vessel_shrink_spe_mimicCC_0.2_flip_nirfast_mesh');
% find(mesh.region==5);
num_points = 1;
selected_indices1 = find(mesh.nodes(:,1) >=4.5 & mesh.nodes(:,1) <= 5.5 &  ...
    mesh.nodes(:,2) >=22.8  & mesh.nodes(:,2) <= 23.2 & ...
    mesh.nodes(:,3) >=10.3 & mesh.nodes(:,3) <=10.7);
num_samples = min(num_points, length(selected_indices1));
point_indices1 = selected_indices1(randperm(length(selected_indices1), num_samples));
% mesh.nodes(point_indices1,:)
% mesh.region(selected_indices1)

selected_indices2 = find(mesh.nodes(:,1) >=4.5 & mesh.nodes(:,1) <= 5.5 &  ...
    mesh.nodes(:,2) >=10.8  & mesh.nodes(:,2) <= 11.2 & ...
    mesh.nodes(:,3) >=9.3 & mesh.nodes(:,3) <=9.7);
num_samples = min(num_points, length(selected_indices2));
point_indices2 = selected_indices2(randperm(length(selected_indices1), num_samples));
% mesh.nodes(point_indices2,:)
% mesh.region(selected_indices2)

selected_indices3 = find(mesh.nodes(:,1) >=4.5 & mesh.nodes(:,1) <= 5.5 &  ...
    mesh.nodes(:,2) >=14.8  & mesh.nodes(:,2) <= 15.2 & ...
    mesh.nodes(:,3) >=8.3 & mesh.nodes(:,3) <=8.7);
num_samples = min(num_points, length(selected_indices3));
point_indices3 = selected_indices3(randperm(length(selected_indices3), num_samples));
% mesh.nodes(point_indices3,:)
% mesh.region(selected_indices3)

selected_indices4 = find(mesh.nodes(:,1) >=4.5 & mesh.nodes(:,1) <= 5.5 &  ...
    mesh.nodes(:,2) >=6.8  & mesh.nodes(:,2) <= 7.2 & ...
    mesh.nodes(:,3) >=7.3 & mesh.nodes(:,3) <=7.7);
num_samples = min(num_points, length(selected_indices4));
point_indices4 = selected_indices4(randperm(length(selected_indices4), num_samples));
% mesh.nodes(point_indices4,:)
% mesh.region(selected_indices4)

selected_indices5 = find(mesh.nodes(:,1) >=4.5 & mesh.nodes(:,1) <= 5.5 &  ...
    mesh.nodes(:,2) >=19.8  & mesh.nodes(:,2) <= 20.2 & ...
    mesh.nodes(:,3) >=5.3 & mesh.nodes(:,3) <=5.7);
num_samples = min(num_points, length(selected_indices5));
point_indices5 = selected_indices5(randperm(length(selected_indices5), num_samples));
% mesh.nodes(point_indices4,:)
% mesh.region(selected_indices5)

wavelengths = linspace(700, 900, 21); 

num =0;

for i = 1:21
    wavelength = wavelengths(i);
    filename = sprintf('nirfast_result/nirfast_0.8_spe_wavelength%d_num%d',wavelength,num); 
    load(filename);
    interested_area = phi_amplitude;
%     sampled_values = interested_area(point_indices1);
    fluence_spectrum1(i) = interested_area(point_indices1);
    fluence_spectrum2(i) = interested_area(point_indices2);
    fluence_spectrum3(i) = interested_area(point_indices3);
    fluence_spectrum4(i) = interested_area(point_indices4);
    fluence_spectrum5(i) = interested_area(point_indices5);
    
end


% % Show fluence curve
% for ii = 1
%     wavelengths = linspace(700, 900, 21); 
%     figure;
%     set(gcf, 'Position', [100, 100, 600, 400]);
%     plot(wavelengths,My_ColorBase_s(ii,:)/max(My_ColorBase_s(ii,:),[],2));
%     xlabel('Wavelength (nm)', 'FontSize', 12);
%     ylabel('Residual fluence', 'FontSize', 12);
%     % title(sprintf('Residual fluence ratio vs. Wavelength, sTO_2 = %.1f(MC-10 mm)', sTO2), 'FontSize', 8);
%     % title(sprintf('Residual fluence ratio(MC-5.4 mm)'), 'FontSize', 8);
%     title(sprintf('Test'), 'FontSize', 8);
%     grid on;
%     xticks(min(xlim):10:max(xlim));
%     xline(730, '--k', 'LineWidth', 1 ,'Color','r');
%     xline(760, '--k', 'LineWidth', 1 ,'Color','r');
%     xline(800, '--k', 'LineWidth', 1 ,'Color','r');
%     set(gca, 'FontSize', 12,'Color', 'w');
%     set(gcf, 'Color', 'w');
%     
% end

filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_spe_fluence_1e9_sTO2%.1f_sBO2%.1f_dt.mat', 0.8, 0.8);
save(filename, 'fluence_spectrum1','fluence_spectrum2','fluence_spectrum3','fluence_spectrum4','fluence_spectrum5');