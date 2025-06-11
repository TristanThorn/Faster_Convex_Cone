clear all;
% close all;
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
cd '/home/zz111/Spectral_unmixing/cc_diffusion_theory'
rng = 42;

filename = sprintf('nirfast_result/nirfast_mlarge_spe_wavelength%d_num%d',700,1);
load(filename);

mesh = load_mesh('digital_vessel_shrink_spe_mimicCC_0.2_flip_nirfast_mesh');
depths = [10.4 9.4 8.4 7.4 5.4];
for v = 1:5
    depth = depths(v);
    % find(mesh.region==5);
    % selected_indices = find(mesh.nodes(:,1) >=8.5 & mesh.nodes(:,1) <= 9.5 &  ...
    %     mesh.nodes(:,2) >=15  & mesh.nodes(:,2) <= 16 & ...
    %     mesh.nodes(:,3) >= 11.9 & mesh.nodes(:,3) <=12.1);
    selected_indices = find(mesh.nodes(:,3) >= depth-1.5 & mesh.nodes(:,3) <=depth+1.5);
    num_points = 50;
    num_samples = min(num_points, length(selected_indices));
    random_indices = selected_indices(randperm(length(selected_indices), num_samples));
    mesh.nodes(random_indices,:)
    
    wavelengths = linspace(700, 900, 21); 
    
    for num =1:100
        num
        for i = 1:21
            wavelength = wavelengths(i);
            filename = sprintf('nirfast_result/nirfast_mlarge_spe_wavelength%d_num%d',wavelength,num); 
            load(filename);
            interested_area = phi_amplitude;
            sampled_values = interested_area(random_indices);
            num_start = 1+(num-1)*num_samples;
            num_end = num*num_samples;
            My_ColorBase_sub(num_start:num_end,i) = sampled_values;
    %         My_ColorBase_s(:,i) = sampled_values;
        end
    end
    My_ColorBase_s{v} = My_ColorBase_sub;
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

save('My_ColorBase_mlarge_shrink_spe_dt_loc', "My_ColorBase_s");