clear all
close all


load('spectrum_HbO2_Cope.mat'); % Molar extinction coefficient of oxyhemoglobin
load('spectrum_Hb_Cope.mat'); % Molar extinction coefficient of deoxyhemoglobin
load('spectrum_H2O.mat'); % Absorption coefficient of water
load('spectrum_fat.mat');
load('spectrum_collagen.mat');
load('spectrum_melanin.mat');

rng(1);
gpuDevice(1);

addpath('/home/zz111/Spectral_unmixing/mcxlab');

ME = 0.01;
sTO2 = 0.9;
sDO2 = sTO2;
sSO2 = sTO2;
sMO2 = sTO2;
sBO2 = 0.9;
[Epid_mu_a, Epid_mu_s_redu] = epidermis(ME,spectrum_H2O,spectrum_melanin);
[der_mu_a,der_mu_s_redu]=dermis(sDO2,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O);
[subcu_mu_a, subcu_mu_s_redu]=subcutaneous(sSO2,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O,spectrum_fat);
[mus_mu_a,mus_mu_s_redu] = muscle(sMO2,spectrum_HbO2,spectrum_Hb,spectrum_H2O);
wavelengths = linspace(700, 900, 21); 

% figure;plot(der_mu_a);

lambda0 = 800; % Reference wavelength
mu_s0_redu = 16.13;   % Scattering coefficient at lambda0
b = 0.66;       % Scattering power

    blood_mu_a = (spectrum_HbO2 .* sBO2 + spectrum_Hb .* (1 - sBO2)) ...
        .* (150/64500) * log(10);
    % Compute reduced scattering coefficient mu_s'
    blood_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);
    for i = 1:21
        i
        tic
        cfg.seed = 0;
        wavelength = wavelengths(i);
        load("digital_vessel_spe_mimicCC_0.2.mat")
        vol = phantom;
        [vol_x, vol_y, vol_z] = size(vol);
%         imshow(squeeze(max(vol,[],1)),[]);
        vol_background = ones(vol_x,vol_y,vol_z+25); 
        vol_background(:, :, 26:end) = vol(:,:,1:end) +1;
        cfg.vol = vol_background;
        cfg.unitinmm = 0.2; % define the pixel size in terms of mm

        cfg.prop = [0         0    1    1   % medium 0: the environment
                    0         0    1    1

%                     der_mu_a(i)*0.1      der_mu_s_redu(i)    0.9    1.3700
%                    der_mu_a(i)*0.1      der_mu_s_redu(i)    0.9    1.3700
%                    subcu_mu_a(i)*0.1    subcu_mu_s_redu(i)    0.9    1.3700
%                    mus_mu_a(i)*0.1      mus_mu_s_redu(i)    0.9    1.3700
%                    der_mu_a(i)*0.1      der_mu_s_redu(i)    0.9    1.3700];

                   Epid_mu_a(i)*0.1     Epid_mu_s_redu(i)  0.9   1.3700
                   der_mu_a(i)*0.1      der_mu_s_redu(i)    0.9    1.3700
                   subcu_mu_a(i)*0.1    subcu_mu_s_redu(i)    0.9    1.3700
                   mus_mu_a(i)*0.1      mus_mu_s_redu(i)    0.9    1.3700
                   blood_mu_a(i)*0.1    blood_mu_s_redu(i)           0.9    1.3700];

% Simulation  
        [Nx,Ny,Nz] = size(cfg.vol);
        cfg.nphoton = 1e9;
        cfg.issrcfrom0 = 1;
        cfg.srcpos = [Nx/2 5 20]; % strat point
        cfg.tstart = 0;
        cfg.tend = 5e-8;%5e-8
        cfg.tstep = 5e-8;%5e-8
        cfg.srcdir = [0 0 1];
        cfg.gpuid=3;
        cfg.srctype = 'line';
        %cfg.srcparam1 = [0.3 / cfg.unitinmm 0 0 0];
        cfg.srcparam1 = [0,Ny-10,0]; %size of light source
        cfg.isreflect = 0;
        cfg.autopilot = 1;
        cfg.debuglevel = 'P';
        
        % cfg.outputtype='energy';
        cfg.outputtype = 'flux';
        flux = mcxlab(cfg);
        
        if (strcmp(cfg.outputtype, 'energy'))
            mcxdata = flux.data / ((cfg.unitinmm / 10)^3);
        else
            mcxdata = flux.data * 100;
        end
        
        if (strcmp(cfg.outputtype, 'flux'))
            mcxdata = mcxdata * cfg.tstep;
        end
        
%             imagesc(log10(abs(squeeze(mcxdata(Nx/2, :, :))))');
        if (strcmp(cfg.outputtype, 'energy'))
            set(gca, 'clim', [-2.4429 4.7581]);
        else
            set(gca, 'clim', [0.5 2.8]);
        
        toc
        end
        toc
        mcxdata = mcxdata(Nx/2-1:Nx/2+1,:,:);
%         interested_area = mcxdata(:,201:600,101:250);
%         imagesc(log10(abs(squeeze(interested_area(3, :, :))))');
    filename = sprintf('/home/zz111/Spectral_unmixing/mcx_result/TestPhantom_spe_sTO2%.1f_sBO2%.1f_wavelength%d_seed%d.mat', sTO2, sBO2, wavelength, cfg.seed);    
    save(filename, "mcxdata");


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

function [der_mu_a,der_mu_s_redu]=dermis(sO2,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O)

    % Mua ------------------------------
    SW = 0.7; % Water volume fraction
    Cblood = 0.03; % Blood concentration
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

function [subcu_mu_a, subcu_mu_s_redu]=subcutaneous(sO2,spectrum_HbO2,spectrum_Hb,spectrum_collagen,spectrum_H2O,spectrum_fat)
    % Mua ------------------------------
    % Define parameters
    SW = 0.7; % Water volume fraction
    Cblood = 0.03; % Blood concentration
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


function [mus_mu_a,mus_mu_s_redu] = muscle(sO2,spectrum_HbO2,spectrum_Hb,spectrum_H2O)
    % Mua ------------------------------
    % Define parameters
    SW = 0.5; % Water content (50%)
    Cblood = 0.03; % Blood concentration (3%)

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