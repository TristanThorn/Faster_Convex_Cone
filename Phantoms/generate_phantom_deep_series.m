% Generate a series of deep phantoms with varying vessel configurations.
clear all;
close all;

% Basic parameters setup
pixel_size = 0.2;
fac = 0.1 / pixel_size; % Conversion factor: pixel size = 0.1/fac mm
nx = 120 * fac; % Size in x-direction
ny = 300 * fac; % Size in y-direction
nz = 250 * fac; % Size in z-direction (axial direction)

% Base cylinder parameters (reference configuration)
base_cylinders = {
    [round(70*fac), round(31*fac), round(11*fac), 5];  % Cylinder 1
    [round(190*fac), round(53*fac), round(13*fac), 6];  % Cylinder 2  
    [round(240*fac), round(78*fac), round(18*fac), 7];  % Cylinder 3
    [round(80*fac), round(98*fac), round(18*fac), 8];   % Cylinder 4
    [round(200*fac), round(116*fac), round(16*fac), 9]; % Cylinder 5
    [round(150*fac), round(140*fac), round(20*fac), 10]; % Cylinder 6
    [round(50*fac), round(165*fac), round(15*fac), 11]; % Cylinder 7
    [round(240*fac), round(214*fac), round(14*fac), 12]}; % Cylinder 8

% Generate 8 different phantom configurations
for iter = 1:8
    % ------------------------
    % 1. Create layered tissue structure
    % ------------------------
    rng(iter); % Set different random seed for each iteration
    phantom = ones(nx, ny, nz); % Initialize 3D matrix
    
    % Generate random layer thicknesses
    layer1_thickness = round((5+5*rand) * fac); 
    layer2_thickness = round((20+20*rand) * fac); 
    layer3_thickness = round((35+15*rand) * fac); 
    layer4_thickness = nz - layer1_thickness - layer2_thickness - layer3_thickness;
    
    % Validate thickness values
    if layer4_thickness <= 0
        error('Invalid thickness values: total thickness exceeds nz.');
    end
    
    % Assign values to each layer
    phantom(:, :, 1:layer1_thickness) = 1;          % Layer 1 (value 1)
    phantom(:, :, layer1_thickness+1:layer1_thickness+layer2_thickness) = 2; % Layer 2 (value 2)
    phantom(:, :, layer1_thickness+layer2_thickness+1:layer1_thickness+layer2_thickness+layer3_thickness) = 3; % Layer 3 (value 3)
    phantom(:, :, layer1_thickness+layer2_thickness+layer3_thickness+1:end) = 4; % Layer 4 (value 4)

    % ------------------------
    % 2. Adjust cylinder parameters
    % ------------------------
    if iter == 1
        cylinders = base_cylinders; % Use original parameters for first iteration
    else
        % Initialize cell array for modified cylinders
        ref_cylinders = cylinders;
        
        % Shift Center_y up by one position
        for i = 1:length(base_cylinders)
            if i == 1
                cylinders{i}(1) = ref_cylinders{end}(1); % First element takes last
            else
                cylinders{i}(1) = ref_cylinders{i-1}(1); % Others take previous
            end
        end
        
        % Shift Radius down by one position
        for i = 1:length(ref_cylinders)
            if i == length(ref_cylinders)
                cylinders{i}(3) = ref_cylinders{1}(3); % Last element takes first
            else
                cylinders{i}(3) = ref_cylinders{i+1}(3); % Others take next
            end
        end
        
        % Adjust Center_z to maintain (Center_z - Radius) difference
        for i = 1:length(ref_cylinders)
            original_diff = ref_cylinders{i}(2) - ref_cylinders{i}(3);
            cylinders{i}(2) = cylinders{i}(3) + original_diff;
        end
        
        % Keep Value column unchanged
        for i = 1:length(ref_cylinders)
            cylinders{i}(4) = ref_cylinders{i}(4);
        end
    end
    
    % ------------------------
    % 3. Embed cylinders into phantom matrix
    % ------------------------
    for i = 1:length(cylinders)
        cy = cylinders{i}(1);       % y-center coordinate
        cz = cylinders{i}(2);       % z-center coordinate
        radius = cylinders{i}(3);   % Radius
        value = cylinders{i}(4);    % Intensity value
        
        % Rasterize cylinder into 3D space
        for x = 1:nx
            for z = 1:nz
                for y = 1:ny
                    % Calculate distance from cylinder axis
                    dist = sqrt((z - cz)^2 + (y - cy)^2);
                    if dist <= radius
                        phantom(x, y, z) = value; % Assign value if inside cylinder
                    end
                end
            end
        end
    end
    
    % ------------------------
    % 4. Save results
    % ------------------------
    % Reorient volume for proper visualization
    vol = flip(flip(flip(phantom, 1), 2), 3);
    vol = uint8(vol); % Convert to 8-bit unsigned integer
    
    % Save as NIfTI file
    % filename = sprintf('digital_VesselDiff_deep_mimicCC_0.2_num%d.nii', iter);
    % niftiwrite(vol, filename);
    
    % Visualize first iteration (optional)

    figure;
    imshow(squeeze(max(vol,[],1)), []);
    % title(sprintf('Iteration %d - Maximum intensity projection', iter));
    % colorbar;

end
