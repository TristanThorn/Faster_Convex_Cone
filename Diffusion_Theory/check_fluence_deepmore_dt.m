clear all;
% close all;
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
cd '/home/zz111/Spectral_unmixing/cc_diffusion_theory'
rng = 42;%42
i=0;
for phantom_num = 1:8
    for para_num = 1:8
        i = i+1;
My_ColorBase_dt = zeros(40000,21);
filename = sprintf('nirfast_result/nirfas_deepmore_nonuni_VesselDiff_wavelength%d_num%d.mat',700,i);
load(filename);

mesh = load_mesh(sprintf('digital_VesselDiff_deep_mimicCC_0.2_num%d_nirfast_mesh',phantom_num));
% find(mesh.region==5);
num_points = 1;

z_offsets = [2, 4, 6, 8, 10, 12, 15, 20];
selected_nodes = cell(1, 8);
for depth = 1:length(z_offsets)
    z_min = 25 - z_offsets(depth)-0.5;
    z_max = z_min+1;
    
    selected_indices = find(...
        mesh.nodes(:,1) >=4 & mesh.nodes(:,1) <=5 & ...
        mesh.nodes(:,3) >=z_min & mesh.nodes(:,3) <=z_max);
    
    valid_indices = selected_indices(mesh.region(selected_indices) == 4+depth);

    point_indices{depth} = valid_indices(randperm(length(valid_indices), 1));
    selected_nodes = mesh.nodes(point_indices{depth}, :); 
    disp(mesh.region(point_indices{depth}));
end

wavelengths = linspace(700, 900, 21); 

num =0;
    for wv = 1:21
        wavelength = wavelengths(wv);
    %     filename = sprintf('nirfast_result/nirfast_0.8_spe_wavelength%d_num%d',wavelength,num); 
    %     filename = sprintf('nirfast_result/nirfast_large_deep_%.1f_%.1f_wavelength%d_num%d',sTO2, sBO2, wavelength, num);
        filename = sprintf('nirfast_result/nirfas_deepmore_nonuni_VesselDiff_wavelength%d_num%d.mat',wavelength,i);
        load(filename);
        interested_area = phi_amplitude;
    %     sampled_values = interested_area(point_indices1);
        fluence_spectrum1(wv) = interested_area(point_indices{1});
        fluence_spectrum2(wv) = interested_area(point_indices{2});
        fluence_spectrum3(wv) = interested_area(point_indices{3});
        fluence_spectrum4(wv) = interested_area(point_indices{4});
        fluence_spectrum5(wv) = interested_area(point_indices{5});
        fluence_spectrum6(wv) = interested_area(point_indices{6});
        fluence_spectrum7(wv) = interested_area(point_indices{7});
        fluence_spectrum8(wv) = interested_area(point_indices{8});
        
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

filename = sprintf('/home/zz111/Spectral_unmixing/TestPhantom_deepmore_nonuni_VesselDiff_fluence_num%d_dt.mat', i);
save(filename, 'fluence_spectrum1','fluence_spectrum2','fluence_spectrum3', ...
    'fluence_spectrum4','fluence_spectrum5',"fluence_spectrum6","fluence_spectrum7","fluence_spectrum8");

end
end