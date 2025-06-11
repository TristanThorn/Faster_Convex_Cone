clear all
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
%%mesh loading
% we can create a mesh from a segmented mha file:
%image2mesh_gui
% and then follow the guide in powerpoint 'Head Model Tutorial.pptx'
% OR use one made for you already...NOTE we will use a spectral mesh
mesh = load_mesh('digital_vessel_mimicCC_0.2_nirfast_mesh');


mesh.source.coord

mesh.source.distributed = 1; % All source are activated togather

% detectors, located at
mesh.meas.coord 
rng(1);
global rand_num;
for rand_num = 1:100
    %% parameter assignment
    % assign tissue properties: Note we are assigning HbO, Hb, Water and
    % scattering...see powerpoint 'Head Model Tutorial.pptx' for values
    % epidermis (region 1)
    values.sa = 10;
    values.sp = 1.4;
    values.HbO =  0;
    values.deoxyHb = 0;
    values.Water = 0.5;
    mesh = set_mesh(mesh,1,values);
    % dermis (region 2)
    sO2 = 0.5 + (0.7 - 0.5) * rand();
    Cblood = 0.03;
    values.sa = 20;
    values.sp = 1.4;
    values.HbO =  sO2 * Cblood;
    values.deoxyHb = (1-sO2) * Cblood;
    values.Water = 0.7;
    mesh = set_mesh(mesh,2,values);
    % subcutaneous (region 3)
    sO2 = 0.5 + (0.7 - 0.5) * rand();
    Cblood = 0.03;
    values.sa = 15;
    values.sp = 0.35;
    values.HbO =  sO2 * Cblood;
    values.deoxyHb = (1-sO2) * Cblood;
    values.Water = 0.7;
    mesh = set_mesh(mesh,3,values);
    % muscle (regions 4)
    sO2 = 0.5 + (0.7 - 0.5) * rand();
    Cblood = 0.03;
    values.sa = 5;
    values.sp = 1;
    values.HbO =  sO2 * Cblood;
    values.deoxyHb = (1-sO2) * Cblood;
    values.Water = 0.5;
    mesh = set_mesh(mesh,4,values);
    % blood (regions 5)
    sO2 = 0.5 + (0.7 - 0.5) * rand();
    Cblood = 1;
    values.sa = 16.13;
    values.sp = 0.66;
    values.HbO =  sO2 * Cblood;
    values.deoxyHb = (1-sO2) * Cblood;
    values.Water = 0;
    mesh = set_mesh(mesh,5,values);
    clear values;
    
    % let us look at the cross section of HbO
    % To do this, we will interpolate values onto a uniform grid, at a cross
    % section
    [x,y,z]=meshgrid(min(mesh.nodes(:,1)):1:max(mesh.nodes(:,1)),100,min(mesh.nodes(:,3)):1:max(mesh.nodes(:,3)));
    F = TriScatteredInterp(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),mesh.conc(:,1));
    val = F(x,y,z);
    figure; surf(squeeze(x),squeeze(y),squeeze(z),squeeze(val));
    view(0,0); shading interp; colorbar; axis tight
    
    tic
    data = femdata_spectral(mesh,0); % frequency = 0
    toc
end