clear all;
close all;

pixel_size = 0.2;
fac=0.1/pixel_size;% pixel size = 0.1/fac mm
% Initialize the size of the 3D matrix
nx = 150 * fac; % Size in the x-direction %180
ny = 200 * fac; % Size in the y-direction
nz = 130 * fac; % Size in the z-direction (axial direction) % 200
phantom = ones(nx, ny, nz); % Initialize the matrix with zeros

% ------------------------
% 1. Layered structure (Tissue layers)
% ------------------------

% Specify the thickness of the first two layers
layer1_thickness = 40 * fac; 
layer2_thickness = 20 * fac; 
layer3_thickness = 40 * fac; 

% Automatically calculate the thickness of the third layer
layer4_thickness = nz - layer1_thickness - layer2_thickness - layer3_thickness;

% Check if the thickness values are valid
if layer4_thickness <= 0
    error('Invalid thickness values: the total thickness exceeds nz.');
end

% Assign values to each layer
phantom(:, :, 1:layer1_thickness) = 1; % First layer with value 1
phantom(:, :, layer1_thickness+1:layer1_thickness+layer2_thickness) = 2; % Second layer with value 2
phantom(:, :, layer1_thickness + layer2_thickness+1 : layer1_thickness + layer2_thickness + layer3_thickness) = 3; % Second layer with value 2
phantom(:, :, layer1_thickness+layer2_thickness+layer3_thickness+1:end) = 4; % Third layer with value 3

% ------------------------
% 2. Cylinders (Blood vessels)
% ------------------------

cylinders = {
    % Center_y, Center_z, Radius, Value
    [round(65*fac), round(37*fac), round(3*fac), 5];
    [round(135*fac), round(37*fac), round(3*fac), 6]; 
    [round(95*fac), round(70*fac), round(3*fac), 7];
    [round(125*fac), round(70*fac), round(3*fac), 8];};


% Iterate through each cylinder and embed it into the phantom matrix
for i = 1:length(cylinders)
    cy = cylinders{i}(1);       % Cylinder center in y-direction
    cz = cylinders{i}(2);       % Cylinder center in z-direction
    radius = cylinders{i}(3);   % Radius of the cylinder
    value = cylinders{i}(4);    % Value assigned to the cylinder

    % Fill the cylinder into the phantom matrix
    for x = 1:nx
        for z = 1:nz
            for y = 1:ny
                % Calculate the distance from the current point to the cylinder center
                dist = sqrt((z - cz)^2 + (y - cy)^2);
                if dist <= radius
                    phantom(x, y, z) = value; % Assign the value if within the cylinder radius
                end
            end
        end
    end
end

% ------------------------
% 3. Visualization
% ------------------------

% Open volume viewer for the phantom matrix
vol = phantom;
vol = flip(flip(flip(vol, 1), 2), 3);
vol = uint8(vol);
phantom = uint8(phantom);
figure;imshow(squeeze(max(vol,[],1)),[]);
figure;imshow(squeeze(max(phantom,[],1)),[]);
volumeViewer(vol);

% ------------------------
% 4. Save to a .mat file
% ------------------------

% Save the resized phantom matrix to a .mat file
% niftiwrite(vol, 'digital_phantom_large_tube.nii')
