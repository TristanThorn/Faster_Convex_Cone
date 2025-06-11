clear all
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
addpath('/home/zz111/Spectral_unmixing/generate_cc')
load('spectrum_HbO2_Cope.mat'); % Molar extinction coefficient of oxyhemoglobin
load('spectrum_Hb_Cope.mat'); % Molar extinction coefficient of deoxyhemoglobin
load('spectrum_H2O.mat'); % Absorption coefficient of water
load('spectrum_fat.mat');
load('spectrum_collagen.mat');
load('spectrum_melanin.mat');

%%mesh loading
mesh = load_mesh('digital_VesselDiff_deep_mimicCC_0.2_num8_nirfast_mesh');

% mesh = load_mesh('digital_vessel_deep_mimicCC_0.2_flip_nirfast_mesh');
mesh.source.coord

mesh.source.distributed = 1; % All source are activated togather

% detectors, located at
mesh.meas.coord 

rng(1);
for save_num = 1:200
    % Initialization
    mua_spec = {};
    musr_spec = {};
    kappa_spec = {};

    mesh.mua = zeros(size(mesh.sa));
    mesh.mus = zeros(size(mesh.sa));
    mesh.kappa = zeros(size(mesh.sa));

    %% parameter assignment
    % assign tissue properties: Note we are assigning HbO, Hb, Water and
    % scattering...see powerpoint 'Head Model Tutorial.pptx' for values
    % epidermis (region 1)
    region = 1;
ME = 0.2*rand();
[mua_spec{region}, musr_spec{region}] = epidermis(ME,spectrum_H2O,spectrum_melanin);
musr_spec{region} = musr_spec{region}*0.1;
mua_spec{region} = musr_spec{region}*0.1;
kappa_spec{region} = 1./(3*(mua_spec{region}+musr_spec{region}));


% dermis (region 2)
region = 2;
sO2 = 0.6+0.4*rand;
[mua_spec{region}, musr_spec{region}] = dermis(sO2,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O);
musr_spec{region} = musr_spec{region}*0.1;
mua_spec{region} = mua_spec{region}*0.1;
kappa_spec{region} = 1./(3*(mua_spec{region}+musr_spec{region}));

% subcutaneous (region 3)
region = 3;
sO2 = 0.6+0.4*rand;
[mua_spec{region}, musr_spec{region}]=subcutaneous(sO2,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O,spectrum_fat);
musr_spec{region} = musr_spec{region}*0.1;
mua_spec{region} = mua_spec{region}*0.1;
kappa_spec{region} = 1./(3*(mua_spec{region}+musr_spec{region}));

% muscle (regions 4)
region = 4;
sO2 = 0.6+0.4*rand;
[mua_spec{region}, musr_spec{region}] = muscle(sO2,spectrum_HbO2,spectrum_Hb,spectrum_H2O);
musr_spec{region} = musr_spec{region}*0.1;
mua_spec{region} = mua_spec{region}*0.1;
kappa_spec{region} = 1./(3*(mua_spec{region}+musr_spec{region}));

% blood (regions 5-12)
for region = 5:12
    mu_s0_redu = 16.13;
    b = 0.66;
    sO2 = 0.5+rand*0.5;
    wavelengths = linspace(700, 900, 21);
    mua_spec{region} = (spectrum_HbO2 .* sO2 + spectrum_Hb .* (1 - sO2)) ...
        .* (150/64500) * log(10);
    musr_spec{region} = mu_s0_redu * (wavelengths / 800) .^ (-b);
    musr_spec{region} = musr_spec{region}*0.1;
    mua_spec{region} = mua_spec{region}*0.1;
    kappa_spec{region} = 1./(3*(mua_spec{region}+musr_spec{region}));
end
    

clear mua mus_r kappa ME sO2

mesh.mua_spec = mua_spec;
mesh.musr_spec = musr_spec;
mesh.kappa_spec = kappa_spec;
    
    % let us look at the cross section of HbO
    % To do this, we will interpolate values onto a uniform grid, at a cross
    % section
    [x,y,z]=meshgrid(min(mesh.nodes(:,1)):1:max(mesh.nodes(:,1)),100,min(mesh.nodes(:,3)):1:max(mesh.nodes(:,3)));
    F = TriScatteredInterp(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),mesh.conc(:,1));
    val = F(x,y,z);
%     figure; surf(squeeze(x),squeeze(y),squeeze(z),squeeze(val));
%     view(0,0); shading interp; colorbar; axis tight
    
    tic
    data = My_femdata_spectral(mesh,0,save_num,true); % frequency = 0
    toc
end

%% functions
function [Epid_mu_a, Epid_mu_s_redu] = epidermis(ME,spectrum_H2O,spectrum_melanin)
    % Mua ------------------------------
    SW = 0.15 + 0.1*rand; % Water volume fraction  
    wavelengths = linspace(700, 900, 21);  
    Epid_mu_a = spectrum_H2O .* SW + spectrum_melanin .* ME;
    
    % Mus' -------------------------------
    % Define parameters
    lambda0 = 800; % Reference wavelength
    mu_s0_redu = 7.44 + 40*rand;   % Scattering coefficient at lambda0
    b = 1.29 + 0.2*rand;       
    Epid_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);
end

function [der_mu_a,der_mu_s_redu]=dermis(sO2,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O)

    % Mua ------------------------------
    SW = 0.525 + 0.35*rand; % Water volume fraction
    Cblood = 0.05*rand; % Blood concentration
    %Cblood = 0.03;
    COLL = 0.5*rand; % Collagen volume fraction    
    % Define wavelengths (example: 700 nm to 900 nm with 21 points)
    wavelengths = linspace(700, 900, 21); 
    
% Mua ------------------------------
    
    % Compute absorption coefficient mu_a(λ)
    der_mu_a = (spectrum_HbO2 .* sO2 + spectrum_Hb .* (1 - sO2)) ...
           .* Cblood .* (150/64500) * log(10) + spectrum_collagen .* COLL...
           + spectrum_H2O .* SW;
    
    % Mus' -------------------------------
    % Define parameters
    lambda0 = 800; % Reference wavelength
    mu_s0_redu = 7.44 + 40*rand;   % Scattering coefficient at lambda0
    b = 1.29 + 0.2*rand;       % Scattering power
    
    % Compute reduced scattering coefficient mu_s'
    der_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);
end

function [subcu_mu_a, subcu_mu_s_redu]=subcutaneous(sO2,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O,spectrum_fat)
    % Mua ------------------------------
    % Define parameters
    SW = 0.525 + 0.35*rand; % Water volume fraction
    Cblood = 0.05*rand; % Blood concentration
    %Cblood = 0.03;
    FA = rand; % Fat volume fraction
    COLL = 0.5 * rand; % Collagen volume fraction

    wavelengths = linspace(700, 900, 21); 
    
    % Compute absorption coefficient mu_a(λ)
    subcu_mu_a = (spectrum_HbO2 .* sO2 + spectrum_Hb .* (1 - sO2)) ...
           .* Cblood .* (150/64500) * log(10) + spectrum_collagen .* COLL...
           + spectrum_H2O .* SW + spectrum_fat .* FA;
    
    
    % Mus' -------------------------------
    % Define parameters
    lambda0 = 800; % Reference wavelength
    mu_s0_redu = 4.4 + 24*rand;   % Scattering coefficient at lambda0
    b = 0.285 + 0.2*rand;       % Scattering power
    
    % Compute reduced scattering coefficient mu_s'
    subcu_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);

end


function [mus_mu_a,mus_mu_s_redu] = muscle(sO2,spectrum_HbO2,spectrum_Hb,spectrum_H2O)
    % Mua ------------------------------
    % Define parameters
    SW = 0.375 + rand*0.25; % Water content (50%)
    Cblood = 0.05*rand; % Blood concentration (3%)
    %Cblood = 0.03;
    wavelengths = linspace(700, 900, 21); 
    mus_mu_a = (spectrum_HbO2 .* sO2 + spectrum_Hb .* (1 - sO2)) ...
           .* Cblood .* (150/64500) * log(10) ...
           + spectrum_H2O .* SW;
    
    % Mus' -------------------------------
    % Define parameters
    lambda0 = 800; % Reference wavelength
    mu_s0_redu = 2.58 + 14*rand;   % Scattering coefficient at lambda0
    b = 0.95 + 0.2*rand;       % Scattering power    

    mus_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);
end