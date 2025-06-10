clear;
close all;

pixel_size = 0.05;
fac=0.1/pixel_size;% pixel size = 0.1/fac mm
% Initialize the size of the 3D matrix
nx = 180 * fac; % Size in the x-direction %180
ny = 300 * fac; % Size in the y-direction
nz = 200 * fac; % Size in the z-direction (axial direction) % 200
phantom = ones(nx, ny, nz); % Initialize the matrix with zeros

% ------------------------
% 1. Layered structure (Tissue layers)
% ------------------------

% Specify the thickness of the first two layers
layer1_thickness = 10 * fac; 
layer2_thickness = 20 * fac; 
layer3_thickness = 30 * fac; 

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

% Define cylinder parameters: center, radius, and value
% cylinders = {
%     % Center_y, Center_z, Radius, Value
%     [100*fac, 13*fac, 3*fac, 5];   % Cylinder 1: radius = 30, value = 4
%     [270*fac, 25*fac, 5*fac, 5];   % Cylinder 1: radius = 30, value = 4
%     [180*fac, 36*fac, 6*fac, 5]; % Cylinder 2: radius = 20, value = 
%     [300*fac, 46*fac, 6*fac, 5]; % Cylinder 2: radius = 20, value = 5
%     [260*fac, 54*fac, 4*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [120*fac, 67*fac, 7*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [200*fac, 75*fac, 5*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [100*fac, 84*fac, 4*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [260*fac, 95*fac, 5*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [220*fac, 107*fac, 7*fac, 5];  % Cylinder 4: radius = 13, value = 6
% };

% cylinders = {
%     % Center_y, Center_z, Radius, Value
%     [100*fac, 15*fac, 5*fac, 5];   % Cylinder 1: radius = 30, value = 4
%     [270*fac, 25*fac, 5*fac, 5];   % Cylinder 1: radius = 30, value = 4
%     [180*fac, 35*fac, 5*fac, 5]; % Cylinder 2: radius = 20, value = 
%     [300*fac, 45*fac, 5*fac, 5]; % Cylinder 2: radius = 20, value = 5
%     [260*fac, 55*fac, 5*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [120*fac, 65*fac, 5*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [200*fac, 75*fac, 5*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [100*fac, 85*fac, 5*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [260*fac, 95*fac, 5*fac, 5];  % Cylinder 3: radius = 15, value = 6
%     [220*fac, 105*fac, 5*fac, 5];  % Cylinder 4: radius = 13, value = 6
% };

% cylinders = {
%     % Center_y, Center_z, Radius, Value
%     [100*fac, 15*fac, 5*fac, 5];
%     [270*fac, 55*fac, 5*fac, 5]; 
%     [180*fac, 105*fac, 5*fac, 5];  
%     [300*fac, 155*fac, 5*fac, 5]; 
%     [260*fac, 205*fac, 5*fac, 5];  
%     [120*fac, 255*fac, 5*fac, 5];  
%     [200*fac, 305*fac, 5*fac, 5]; };

cylinders = {
    % Center_y, Center_z, Radius, Value
    [round(70*fac), round(31*fac), round(11*fac), 5];
    [round(190*fac), round(53*fac), round(13*fac), 5]; 
    [round(130*fac), round(47*fac), round(17*fac), 5];  
    [round(230*fac), round(68*fac), round(18*fac), 5]; };


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
% niftiwrite(vol, 'digital_vessel_com_mimicCC_0.2_flip.nii')
% save('digital_vessel_new_mimicCC_0.2.mat', 'phantom');
% 
% 
% close all;
% clear all;
% load('digital_vessel_mimicCC_0.1.mat');
% vol = uint8(vol);
% volumeViewer(vol);
% niftiwrite(vol, 'digital_vessel_mimicCC_0.1.nii')