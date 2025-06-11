clear all;
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
cd '/home/zz111/Spectral_unmixing/cc_diffusion_theory'

current_seed = 0;
% Define 8 random seeds
for phantom_num = 1:8
    mesh = load_mesh(sprintf('digital_VesselDiff_deep_mimicCC_0.2_num%d_nirfast_mesh',phantom_num));
% Load mesh data (single load for efficiency)
% mesh = load_mesh('digital_VesselDiff_deep_mimicCC_0.2_flip_nirfast_mesh');

% Find region of interest coordinates
selected_indices = find(mesh.nodes(:,1) >=3 & mesh.nodes(:,1) <=6 & ...
                    mesh.nodes(:,3) >=1 & mesh.nodes(:,3) <=24);

wavelengths = linspace(700, 900, 21);

% Main processing loop for each seed
    rng(current_seed);  % Initialize random generator with current seed
    
    % Random spatial points selection
    num_spatial_samples = 200;
    random_indices = selected_indices(randperm(length(selected_indices), num_spatial_samples));
    
    % Generate unique numeric samples using full permutation
    % 1. Shuffle all 800 numbers
    % 2. Select first 200 elements
    num_samples_num = randperm(800, 200);  % Direct method (MATLAB 2011b+)
    
    % Initialize data container (200 samples × 200 points × 21 wavelengths)
    My_ColorBase_s = zeros(40000, 21);  % 200*200 = 40000 rows
    
    % Process each numeric sample
    for num_idx = 1:200
        current_num = num_samples_num(num_idx);
        
        % Process all wavelengths for current num
        for i = 1:21
            wavelength = wavelengths(i);
            % Load precomputed data file
%             filename = sprintf('nirfast_result/nirfast_large_deep_VesselDiff_wavelength%d_num%d.mat',wavelength,current_num); 
            filename = sprintf('nirfast_result/nirfast_large_deep%d_VesselDiff_wavelength%d_num%d.mat',phantom_num,wavelength,num_idx);
            load(filename);
            
            % Extract amplitude values at selected spatial points
            sampled_values = phi_amplitude(random_indices);
            
            % Calculate storage indices
            num_start = 1 + (num_idx-1)*num_spatial_samples;
            num_end = num_idx*num_spatial_samples;
            
            % Store wavelength-specific data
            My_ColorBase_s(num_start:num_end, i) = sampled_values;
        end
    end
    
    % Save seed-specific results
    save(sprintf('My_ColorBase_large_deep%d_VesselDiff_200_200', phantom_num), "My_ColorBase_s");
    fprintf('Phantom %d processing completed\n', phantom_num);
end

