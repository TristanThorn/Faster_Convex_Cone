clear all;
% close all;
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
cd '/home/zz111/Spectral_unmixing/cc_diffusion_theory'
rng = 42;%42
My_ColorBase_dt = zeros(10000,21);
sTO2 = 0.6;
sBO2 = 0.6;
filename = sprintf('nirfast_result/nirfas_deep_nonuni_VesselDiff_wavelength%d_num%d.mat',700, 1);
load(filename);

mesh = load_mesh('digital_VesselDiff_deep_mimicCC_0.2_flip_nirfast_mesh');
% find(mesh.region==5);
num_points = 1;
selected_indices1 = find(mesh.nodes(:,1) >=4 & mesh.nodes(:,1) <= 5 &  ...
    mesh.nodes(:,2) >=22.5  & mesh.nodes(:,2) <= 23.5 & ...
    mesh.nodes(:,3) >=22.5 & mesh.nodes(:,3) <=23.5);
num_samples = min(num_points, length(selected_indices1));
point_indices1 = selected_indices1(randperm(length(selected_indices1), num_samples));
% mesh.nodes(point_indices1,:)
mesh.region(selected_indices1)

selected_indices2 = find(mesh.nodes(:,1) >=4 & mesh.nodes(:,1) <= 5 &  ...
    mesh.nodes(:,2) >=10.5  & mesh.nodes(:,2) <= 11.5 & ...
    mesh.nodes(:,3) >=20.5 & mesh.nodes(:,3) <=21.5);
num_samples = min(num_points, length(selected_indices2));
point_indices2 = selected_indices2(randperm(length(selected_indices1), num_samples));
% mesh.nodes(point_indices2,:)
mesh.region(selected_indices2)

selected_indices3 = find(mesh.nodes(:,1) >=4 & mesh.nodes(:,1) <= 5 &  ...
    mesh.nodes(:,2) >=6.5  & mesh.nodes(:,2) <= 7.5 & ...
    mesh.nodes(:,3) >=18.5 & mesh.nodes(:,3) <=19.5);
num_samples = min(num_points, length(selected_indices3));
point_indices3 = selected_indices3(randperm(length(selected_indices3), num_samples));
% mesh.nodes(point_indices3,:)
mesh.region(selected_indices3)

selected_indices4 = find(mesh.nodes(:,1) >=4 & mesh.nodes(:,1) <= 5 &  ...
    mesh.nodes(:,2) >=21.5  & mesh.nodes(:,2) <= 22.5 & ...
    mesh.nodes(:,3) >=16.5 & mesh.nodes(:,3) <=17.5);
num_samples = min(num_points, length(selected_indices4));
point_indices4 = selected_indices4(randperm(length(selected_indices4), num_samples));
% mesh.nodes(point_indices4,:)
mesh.region(selected_indices4)

selected_indices5 = find(mesh.nodes(:,1) >=4 & mesh.nodes(:,1) <= 5 &  ...
    mesh.nodes(:,2) >=14.5  & mesh.nodes(:,2) <= 15.5 & ...
    mesh.nodes(:,3) >=12.5 & mesh.nodes(:,3) <=13.5);
num_samples = min(num_points, length(selected_indices5));
point_indices5 = selected_indices5(randperm(length(selected_indices5), num_samples));
% mesh.nodes(point_indices4,:)
mesh.region(selected_indices5)

selected_indices6 = find(mesh.nodes(:,1) >=4 & mesh.nodes(:,1) <= 5 &  ...
    mesh.nodes(:,2) >=24.5  & mesh.nodes(:,2) <= 25.5 & ...
    mesh.nodes(:,3) >=9.5 & mesh.nodes(:,3) <=10.5);
num_samples = min(num_points, length(selected_indices5));
point_indices6 = selected_indices6(randperm(length(selected_indices5), num_samples));
% mesh.nodes(point_indices4,:)
mesh.region(selected_indices6)

selected_indices7 = find(mesh.nodes(:,1) >=4 & mesh.nodes(:,1) <= 5 &  ...
    mesh.nodes(:,2) >=5.5  & mesh.nodes(:,2) <= 6.5 & ...
    mesh.nodes(:,3) >=4.5 & mesh.nodes(:,3) <=5.5);
num_samples = min(num_points, length(selected_indices5));
point_indices7 = selected_indices7(randperm(length(selected_indices5), num_samples));
% mesh.nodes(point_indices4,:)
mesh.region(selected_indices7)

wavelengths = linspace(700, 900, 21); 

num =0;
for save_num = 1:32
    for i = 1:21
        wavelength = wavelengths(i);
    %     filename = sprintf('nirfast_result/nirfast_0.8_spe_wavelength%d_num%d',wavelength,num); 
    %     filename = sprintf('nirfast_result/nirfast_large_deep_%.1f_%.1f_wavelength%d_num%d',sTO2, sBO2, wavelength, num);
        filename = sprintf('nirfast_result/nirfas_deep_nonuni_VesselDiff_wavelength%d_num%d.mat',wavelength,save_num);
        load(filename);
        interested_area = phi_amplitude;
    %     sampled_values = interested_area(point_indices1);
        fluence_spectrum1(i) = interested_area(point_indices1);
        fluence_spectrum2(i) = interested_area(point_indices2);
        fluence_spectrum3(i) = interested_area(point_indices3);
        fluence_spectrum4(i) = interested_area(point_indices4);
        fluence_spectrum5(i) = interested_area(point_indices5);
        fluence_spectrum6(i) = interested_area(point_indices6);
        fluence_spectrum7(i) = interested_area(point_indices7);
        
    end
    
    
    % Show fluence curve
    for ii = 1
        wavelengths = linspace(700, 900, 21); 
%         figure;
        set(gcf, 'Position', [100, 100, 600, 400]);
        plot(wavelengths,fluence_spectrum4);
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
        
    end

filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_deep_nonuni_VesselDiff_fluence_num%d_dt.mat', save_num);
save(filename, 'fluence_spectrum1','fluence_spectrum2','fluence_spectrum3', ...
    'fluence_spectrum4','fluence_spectrum5',"fluence_spectrum6","fluence_spectrum7");

end