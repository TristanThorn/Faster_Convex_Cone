clear all
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
addpath('/home/zz111/Spectral_unmixing/generate_cc')

save_num = 1;
for deep_num = 1:8

load('spectrum_HbO2_Cope.mat'); % Molar extinction coefficient of oxyhemoglobin
load('spectrum_Hb_Cope.mat'); % Molar extinction coefficient of deoxyhemoglobin
load('spectrum_H2O.mat'); % Absorption coefficient of water
load('spectrum_fat.mat');
load('spectrum_collagen.mat');
load('spectrum_melanin.mat');

%%mesh loading
% mesh = load_mesh('digital_vessel_mimicCC_0.2_flip_nirfast_mesh');
% mesh = load_mesh('digital_VesselDiff_deep_mimicCC_0.2_flip_nirfast_mesh');
mesh = load_mesh(sprintf('digital_VesselDiff_deep_mimicCC_0.2_num%d_nirfast_mesh',deep_num));

mesh.source.coord

mesh.source.distributed = 1; % All source are activated togather

% detectors, located at
mesh.meas.coord 

rng(1);

for distri_num = 1:8

    % Initialization
    mua_spec = [];
    musr_spec = [];
    kappa_spec = [];
    
    mesh.mua = zeros(size(mesh.sa));
    mesh.mus = zeros(size(mesh.sa));
    mesh.kappa = zeros(size(mesh.sa));
    
    
    %% parameter assignment
    
    % Creating non-uni map
    % for sO2
    X_mm = max(squeeze(mesh.nodes(:,1))); % mm
    Y_mm = max(squeeze(mesh.nodes(:,2)));
    Z_mm = max(squeeze(mesh.nodes(:,3)));
    X_pixel = ceil(X_mm/0.1)+1; % pixel size 0.1 mm
    Y_pixel = ceil(Y_mm/0.1)+1;
    Z_pixel = ceil(Z_mm/0.1)+1;
    A = rand(4, 4, 4);
    sO2_map = imresize3(A,[X_pixel,Y_pixel,Z_pixel],'cubic');
    sO2_map = (sO2_map-min(sO2_map(:)))/(max(sO2_map(:))-min(sO2_map(:)));
    % imshow(squeeze(max(sO2_map,[],1)),[]);
    
    % for Cblood
    B = rand(4, 4, 4);
    Cblood_map = imresize3(B,[X_pixel,Y_pixel,Z_pixel],'cubic');
    Cblood_map = (Cblood_map-min(Cblood_map(:)))/(max(Cblood_map(:))-min(Cblood_map(:)));
    % imshow(squeeze(max(Cblood_map,[],1)),[]);
    
    
    % epidermis (region 1)
    region = 1;
    idxs = find(mesh.region==region);
    ME = 0.01;
    for i = 1:length(idxs)
        idx = idxs(i);
        [mua_spec(idx,:), musr_spec(idx,:)] = epidermis(ME,spectrum_H2O,spectrum_melanin);
        musr_spec(idx,:) = musr_spec(idx,:)*0.1; % change unit to mm
        mua_spec(idx,:) = mua_spec(idx,:)*0.1;
        kappa_spec(idx,:) = 1./(3*(mua_spec(idx,:)+musr_spec(idx,:)));
    end
    
    
    % dermis (region 2)
    region = 2;
    idxs = find(mesh.region==region);
    for i = 1:length(idxs)
        idx = idxs(i);
        co_pixel = round(mesh.nodes(idx,:)/0.1)+1;
    
        sO2_fac = sO2_map(co_pixel(1),co_pixel(2),co_pixel(3));
        sO2 = 0.6 + 0.4 * sO2_fac;
    
        Cblood_fac = Cblood_map(co_pixel(1),co_pixel(2),co_pixel(3));
        Cblood = 0.05 * Cblood_fac;
    
        [mua_spec(idx,:), musr_spec(idx,:)] = dermis(sO2,Cblood,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O);
        musr_spec(idx,:) = musr_spec(idx,:)*0.1; % change unit to mm
        mua_spec(idx,:) = mua_spec(idx,:)*0.1;
        kappa_spec(idx,:) = 1./(3*(mua_spec(idx,:)+musr_spec(idx,:)));
    end
    
    % subcutaneous (region 3)
    region = 3;
    idxs = find(mesh.region==region);
    for i = 1:length(idxs)
        idx = idxs(i);
        co_pixel = round(mesh.nodes(idx,:)/0.1)+1;
    
        sO2_fac = sO2_map(co_pixel(1),co_pixel(2),co_pixel(3));
        sO2 = 0.6 + 0.4 * sO2_fac;
    
        Cblood_fac = Cblood_map(co_pixel(1),co_pixel(2),co_pixel(3));
        Cblood = 0.05 * Cblood_fac;
    
        [mua_spec(idx,:), musr_spec(idx,:)] = subcutaneous(sO2,Cblood,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O,spectrum_fat);
        musr_spec(idx,:) = musr_spec(idx,:)*0.1; % change unit to mm
        mua_spec(idx,:) = mua_spec(idx,:)*0.1;
        kappa_spec(idx,:) = 1./(3*(mua_spec(idx,:)+musr_spec(idx,:)));
    end
    
    % muscle (regions 4)
    region = 4;
    idxs = find(mesh.region==region);
    for i = 1:length(idxs)
        idx = idxs(i);
        co_pixel = round(mesh.nodes(idx,:)/0.1)+1;
    
        sO2_fac = sO2_map(co_pixel(1),co_pixel(2),co_pixel(3));
        sO2 = 0.6 + 0.4 * sO2_fac;
    
        Cblood_fac = Cblood_map(co_pixel(1),co_pixel(2),co_pixel(3));
        Cblood = 0.05 * Cblood_fac;
    
        [mua_spec(idx,:), musr_spec(idx,:)] = muscle(sO2,Cblood,spectrum_HbO2,spectrum_Hb,spectrum_H2O);
        musr_spec(idx,:) = musr_spec(idx,:)*0.1; % change unit to mm
        mua_spec(idx,:) = mua_spec(idx,:)*0.1;
        kappa_spec(idx,:) = 1./(3*(mua_spec(idx,:)+musr_spec(idx,:)));
    end
    
    
    % blood (regions 5-11)
    for region = 5:12
        idxs = find(mesh.region==region);
        mu_s0_redu = 16.13;
        b = 0.66;
        sO2 = 0.5 + 0.5 * rand;
        wavelengths = linspace(700, 900, 21);
        for i = 1:length(idxs)
            idx = idxs(i);
            mua_spec(idx,:) = (spectrum_HbO2 .* sO2 + spectrum_Hb .* (1 - sO2)) ...
                .* (150/64500) * log(10);
            musr_spec(idx,:) = mu_s0_redu * (wavelengths / 800) .^ (-b);
            musr_spec(idx,:) = musr_spec(idx,:)*0.1;
            mua_spec(idx,:) = mua_spec(idx,:)*0.1;
            kappa_spec(idx,:) = 1./(3*(mua_spec(idx,:)+musr_spec(idx,:)));
        end
    end
    
    % for i = 1:5
    %     figure;plot(mua_spec{i});
    % %     figure;plot(musr_spec{i});
    % end
    
    clear mua mus_r kappa ME sO2
    
    mesh.mua_spec = mua_spec;
    mesh.musr_spec = musr_spec;
    mesh.kappa_spec = kappa_spec;
    
    % % let us look at the cross section of HbO
    % % To do this, we will interpolate values onto a uniform grid, at a cross
    % % section
    % [x,y,z]=meshgrid(min(mesh.nodes(:,1)):1:max(mesh.nodes(:,1)),100,min(mesh.nodes(:,3)):1:max(mesh.nodes(:,3)));
    % F = TriScatteredInterp(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),mesh.conc(:,1));
    % val = F(x,y,z);
    % figure; surf(squeeze(x),squeeze(y),squeeze(z),squeeze(val));
    % view(0,0); shading interp; colorbar; axis tight
    
    tic
    data = My_femdata_spectral_nonuni(mesh,0,save_num,true); % frequency = 0
    toc
save_num = save_num + 1;
end
end


%% functions
function [Epid_mu_a, Epid_mu_s_redu] = epidermis(ME,spectrum_H2O,spectrum_melanin)
    % Mua ------------------------------
    SW = 0.2; % Water volume fraction  
    wavelengths = linspace(700, 900, 21);  
    Epid_mu_a = spectrum_H2O .* SW + spectrum_melanin .* ME;
    
    % Mus' -------------------------------
    % Define parameters
    lambda0 = 800; % Reference wavelength
    mu_s0_redu = 10.0;   % Scattering coefficient at lambda0
    b = 1.4;       
    Epid_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);
end

function [der_mu_a,der_mu_s_redu]=dermis(sO2,Cblood,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O)

    % Mua ------------------------------
    SW = 0.7; % Water volume fraction
    COLL = 0.25; % Collagen volume fraction    
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
    mu_s0_redu = 20.0;   % Scattering coefficient at lambda0
    b = 1.4;       % Scattering power
    
    % Compute reduced scattering coefficient mu_s'
    der_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);
end

function [subcu_mu_a, subcu_mu_s_redu]=subcutaneous(sO2,Cblood,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O,spectrum_fat)
    % Mua ------------------------------
    % Define parameters
    SW = 0.7; % Water volume fraction
    FA = 0.3; % Fat volume fraction
    COLL = 0.25; % Collagen volume fraction

    wavelengths = linspace(700, 900, 21); 
    
    % Compute absorption coefficient mu_a(λ)
    subcu_mu_a = (spectrum_HbO2 .* sO2 + spectrum_Hb .* (1 - sO2)) ...
           .* Cblood .* (150/64500) * log(10) + spectrum_collagen .* COLL...
           + spectrum_H2O .* SW + spectrum_fat .* FA;
    
    
    % Mus' -------------------------------
    % Define parameters
    lambda0 = 800; % Reference wavelength
    mu_s0_redu = 15.0;   % Scattering coefficient at lambda0
    b = 0.35;       % Scattering power
    
    % Compute reduced scattering coefficient mu_s'
    subcu_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);

end


function [mus_mu_a,mus_mu_s_redu] = muscle(sO2,Cblood,spectrum_HbO2,spectrum_Hb,spectrum_H2O)
    % Mua ------------------------------
    % Define parameters
    SW = 0.5; % Water content (50%)

    wavelengths = linspace(700, 900, 21); 
    mus_mu_a = (spectrum_HbO2 .* sO2 + spectrum_Hb .* (1 - sO2)) ...
           .* Cblood .* (150/64500) * log(10) ...
           + spectrum_H2O .* SW;
    
    % Mus' -------------------------------
    % Define parameters
    lambda0 = 800; % Reference wavelength
    mu_s0_redu = 5.0;   % Scattering coefficient at lambda0
    b = 1;       % Scattering power    

    mus_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);
end