clear al
clearvars
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
% In this example, we will use a 2D case to demonstrate how we can generate
% simulated forward data, and use tomographic recontruction from these data
% to recover location of anamalies
%% Load the 2D multiwavelength wavelength mesh
% This mesh is 43 mm in diameter and has 16 source/detectors on the
% periphery
% The wavelengths being modelled are 661   735   761   785   808   826
% 849 nm
mesh = load_mesh('circle2000_86_spec'); 
% Set anomaly size and properties
% Anomaly 1
blob.x=-14; % x location
blob.y=14; % y location
blob.r=10; % radius (mm)
blob.region=2;
blob.sa=1;
blob.sp=1;
blob.HbO=0.02;
blob.deoxyHb=0.01;
blob.Water=0.4;
mesh_anom = add_blob(mesh,blob);
% Anomaly 2
blob.x=14; % x location
blob.y=14; % y location
blob.r=10; % radius (mm)
blob.region=3;
blob.sa=1;
blob.sp=1;
blob.HbO=0.01;
blob.deoxyHb=0.02;
blob.Water=0.4;
mesh_anom = add_blob(mesh_anom,blob);
% Anomaly 3
blob.x=-14; % x location
blob.y=-14; % y location
blob.r=10; % radius (mm)
blob.region=4;
blob.sa=1;
blob.sp=1;
blob.HbO=0.01;
blob.deoxyHb=0.01;
blob.Water=0.8;
mesh_anom = add_blob(mesh_anom,blob);
% Anomaly 4
blob.x=14; % x location
blob.y=-14; % y location
blob.r=10; % radius (mm)
blob.region=5;
blob.sa=0.5;
blob.sp=1;
blob.HbO=0.01;
blob.deoxyHb=0.01;
blob.Water=0.4;
mesh_anom = add_blob(mesh_anom,blob);
% Anomaly 5
blob.x=14; % x location
blob.y=-14; % y location
blob.r=10; % radius (mm)
blob.region=6;
blob.sa=1;
blob.sp=0.5;
blob.HbO=0.01;
blob.deoxyHb=0.01;
blob.Water=0.4;
mesh_anom = add_blob(mesh_anom,blob);
% And view the mesh
plotmesh(mesh_anom,1);

%% Generate forward data using Frequency Domain model at 100 MHz
data_anom = femdata(mesh_anom,100);
% and add random noise (1% amplitude, and 2 degrees in Phase)
data_anom_noise = add_noise(data_anom,1,2);
% and view data
plot_data(data_anom_noise);

%% Reconstruct
lambda.type='Automatic';
lambda.value=1;
[mesh_recon,pj] = reconstruct(mesh,[30 30],100,data_anom_noise,40,lambda,'spectral_example',0);
plotmesh(mesh_recon);

%% Reconstruct assumin some information about location
% in this case, we assume we know the region of the blobs and will use this
% as priori knowledge
% Copy region label into mesh for reconstructing
mesh.region = mesh_anom.region;
lambda=1;
[mesh_recon_spatial,pj] = reconstruct_spectral_spatial(mesh,[30 30],100,data_anom_noise,40,lambda,'spectral_spatial_example',0);
plotmesh(mesh_recon_spatial);