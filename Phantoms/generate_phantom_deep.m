% Generate a deeper layered phantom with multiple vessels.
clear all;
close all;

pixel_size = 0.05;
fac=0.1/pixel_size;% pixel size = 0.1/fac mm
% Initialize the size of the 3D matrix
nx = 120 * fac; % Size in the x-direction %180
ny = 300 * fac; % Size in the y-direction
nz = 250 * fac; % Size in the z-direction (axial direction) % 200
phantom = ones(nx, ny, nz); % Initialize the matrix with zeros


% ------------------------
% 1. Layered structure (Tissue layers)
% ------------------------
rng(1);
% Specify the thickness of the first two layers
layer1_thickness = round((5+5*rand) * fac); 
layer2_thickness = round((20+20*rand) * fac); 
layer3_thickness = round((35+15*rand) * fac); 

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
% 1
cylinders = {
    % Center_y, Center_z, Radius, Value
    [round(70*fac), round(31*fac), round(11*fac), 5];
    [round(190*fac), round(53*fac), round(13*fac), 6];   
    [round(240*fac), round(78*fac), round(18*fac), 7];
    [round(80*fac), round(98*fac), round(18*fac), 8];
    [round(200*fac), round(116*fac), round(16*fac), 8];
    [round(150*fac), round(140*fac), round(20*fac), 9];
    [round(50*fac), round(165*fac), round(15*fac), 10];
    [round(240*fac), round(214*fac), round(14*fac), 11];};



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
% A = rand(4, 4, 4);
% sO2_map = imresize3(A,[nx,ny,nz],'cubic');
% sO2_map = (sO2_map-min(sO2_map(:)))/(max(sO2_map(:))-min(sO2_map(:)))*11;
% sO2_map(phantom==5) = phantom(phantom==5);
% sO2_map(phantom==6) = phantom(phantom==6);
% sO2_map(phantom==7) = phantom(phantom==7);
% sO2_map(phantom==8) = phantom(phantom==8);
% sO2_map(phantom==9) = phantom(phantom==9);
% sO2_map(phantom==10) = phantom(phantom==10);
% sO2_map(phantom==11) = phantom(phantom==11);
% imshow(squeeze(max(sO2_map,[],1)),[]);
% volumeViewer(vol);

% ------------------------
% 4. Save to a .mat file
% ------------------------

% Save the resized phantom matrix to a .mat file
niftiwrite(vol, fullfile('data','phantoms','digital_VesselDiff_deep_mimicCC_0.2_flip2.nii'))
% save('digital_vessel_new_mimicCC_0.2.mat', 'phantom');

% close all;
% clear all;
% load('digital_vessel_mimicCC_0.1.mat');
% vol = uint8(vol);
% volumeViewer(vol);
% niftiwrite(vol, 'digital_vessel_mimicCC_0.1.nii')
