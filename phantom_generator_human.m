clear all;
close all;

pixel_size = 0.2;
fac=0.1/pixel_size;% pixel size = 0.1/fac mm
% Initialize the size of the 3D matrix
nx = 120 * fac; % Size in the x-direction %180
ny = 200 * fac; % Size in the y-direction
nz = 200 * fac; % Size in the z-direction (axial direction) % 200
phantom = ones(nx, ny, nz); % Initialize the matrix with zeros


% ------------------------
% 1. Layered structure (Tissue layers)
% ------------------------
rng(1);
% Specify the thickness of the first two layers
layer1_thickness = round(5 * fac); % skin
layer2_thickness = round(15 * fac); % skin
layer3_thickness = round(38 * fac); % subcu
layer4_thickness = round(36 * fac); % muscle

% Automatically calculate the thickness of the third layer
layer5_thickness = nz - layer1_thickness - layer2_thickness - layer3_thickness - layer4_thickness; % subcu

% Check if the thickness values are valid
if layer4_thickness <= 0
    error('Invalid thickness values: the total thickness exceeds nz.');
end

% Assign values to each layer
phantom(:, :, 1:layer1_thickness) = 1; % First layer with value 1 (epidermis)
phantom(:, :, layer1_thickness+1:layer1_thickness+layer2_thickness) = 2; % Second layer with value 2 (dermis)
phantom(:, :, layer1_thickness + layer2_thickness+1 : layer1_thickness + layer2_thickness + layer3_thickness) = 3; % value 3 (subcutenous)
phantom(:, :, layer1_thickness+layer2_thickness+layer3_thickness+1:layer1_thickness+layer2_thickness+layer3_thickness+layer4_thickness) = 4; %value 4 (muscle)
phantom(:, :, layer1_thickness+layer2_thickness+layer3_thickness+layer4_thickness+1:end) = 3; % value 3 (subcutenous)

% ------------------------
% 2. Cylinders (Blood vessels)
% ------------------------
% 1
% Define elliptical cylinders as [Center_y, Center_z, Radius_y, Radius_z, Value]
cylinders = {
    [round(81*fac), round(46*fac), round(9*fac), round(4*fac), 5];
    [round(30*fac), round(51*fac), round(9*fac), round(3*fac), 6];   
    [round(150*fac), round(65*fac), round(9*fac), round(9*fac), 7];
    [round(74*fac), round(74*fac), round(10*fac), round(11*fac), 8];
    [round(135*fac), round(91*fac), round(20*fac), round(17*fac), 8];
    [round(105*fac), round(98*fac), round(10*fac), round(8*fac), 9];
    [round(115*fac), round(147*fac), round(15*fac), round(10*fac), 10];};

% Iterate through each elliptical cylinder and embed it into the phantom matrix
for i = 1:length(cylinders)
    cy = cylinders{i}(1);       % Center y
    cz = cylinders{i}(2);       % Center z
    ry = cylinders{i}(3);       % Radius in y-direction
    rz = cylinders{i}(4);       % Radius in z-direction
    value = cylinders{i}(5);    % Assigned value

    for x = 1:nx
        for z = 1:nz
            for y = 1:ny
                dy = y - cy;
                dz = z - cz;
                if (dy^2)/(ry^2) + (dz^2)/(rz^2) <= 1
                    phantom(x, y, z) = value;
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
niftiwrite(vol, 'digital_phantom_human.nii')
% save('digital_vessel_new_mimicCC_0.2.mat', 'phantom');

% close all;
% clear all;
% load('digital_vessel_mimicCC_0.1.mat');
% vol = uint8(vol);
% volumeViewer(vol);
% niftiwrite(vol, 'digital_vessel_mimicCC_0.1.nii')
