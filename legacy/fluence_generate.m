%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  mcxyz skinvessel benchmark
%
%  must change mcxyz maketissue.m boundaryflag variable from 2 to 1 to get
%  comparable absorption fraction (40%), otherwise, mcxyz obtains slightly
%  higher absorption (~42%) with boundaryflag=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gpuDevice(2);
clear all;
close all;
% only clear cfg to avoid accidentally clearing other useful data
clear cfg flux;
% cd '/home/zz111/Spectral_unmixing'
% addpath(genpath('/home/zz111/Spectral_unmixing/mcx'));

addpath('/home/zz111/Spectral_unmixing/spectrum');
addpath('/home/zz111/Spectral_unmixing/mcxlab');
% load mcxyz_skinvessel.mat
for sTO2 = 0.8
filename = sprintf('subcutaneous_coeff_%.1f.mat', sTO2);
load(filename);
filename = sprintf('muscle_coeff_%.1f.mat', sTO2);
load(filename);
filename = sprintf('dermis_coeff_%.1f.mat', sTO2);
load(filename);
filename = 'Epidermis_coeff.mat';
load(filename);
load('spectrum_H2O.mat')
wavelengths = linspace(700, 900, 21); 

% coefficient for blood
load('spectrum_HbO2_Cope.mat'); % Molar extinction coefficient of oxyhemoglobin
load('spectrum_Hb_Cope.mat'); % Molar extinction coefficient of deoxyhemoglobin
lambda0 = 800; % Reference wavelength
mu_s0_redu = 16.13;   % Scattering coefficient at lambda0
b = 0.66;       % Scattering power


for sBO2 = 0.8
    blood_mu_a = (spectrum_HbO2 .* sBO2 + spectrum_Hb .* (1 - sBO2)) ...
        .* (150/64500) * log(10);
    figure;plot(blood_mu_a);
    % Compute reduced scattering coefficient mu_s'
    blood_mu_s_redu = mu_s0_redu * (wavelengths / lambda0) .^ (-b);
    figure;plot(blood_mu_s_redu);
    for i = 1:21
        i
        tic
        cfg.seed = 1;
        wavelength = wavelengths(i);
        load("digital_vessel_mimicCC_0.2.mat")
        [vol_x, vol_y, vol_z] = size(vol);
%          vol = round(imresize3(vol,[800,1600,600],'Method','linear'));
        vol_background = ones(vol_x,vol_y,vol_z+25); 
        vol_background(:, :, 26:end) = vol(:,:,1:end) +1;
        cfg.vol = vol_background;
        cfg.unitinmm = 0.2; % define the pixel size in terms of mm
        cfg.prop = [0         0    1    1   % medium 0: the environment
                    0         0    1    1
                   Epid_mu_a(i)*0.1     Epid_mu_s_redu(i)  0.9   1.3700
                   subcu_mu_a(i)*0.1    subcu_mu_s_redu(i)    0.9    1.3700
                   mus_mu_a(i)*0.1      mus_mu_s_redu(i)    0.9    1.3700
                   blood_mu_a(i)*0.1    blood_mu_s_redu(i)           0.9    1.3700];
        
                   %% visulization simulation setup
                    % vol_a = vol_background;
                    % vol_a(vol_background == 1) = spectrum_H2O(i);
                    % vol_a(vol_background == 2) = Epid_mu_a(i);
                    % vol_a(vol_background == 3) = der_mu_a(i);
                    % vol_a(vol_background == 4) = subcu_mu_a(i);
                    % vol_a(vol_background == 5) = mus_mu_a(i);
                    % vol_a(vol_background == 6) = blood_mu_a(i);
                    % % imagesc(squeeze(max(vol_a,[],1)));
                    % figure;imagesc(squeeze(vol_background(45, :, :))');
                    % figure;imagesc(squeeze(vol_a(45, :, :))');
                    % 
                    % figure;imagesc(squeeze(vol_background(100, :, :))');

%% Simulation  
        [Nx,Ny,Nz] = size(cfg.vol);
        cfg.nphoton = 1e9;
        cfg.issrcfrom0 = 1;
        cfg.srcpos = [Nx/2 5 20]; % strat point
        cfg.tstart = 0;
        cfg.tend = 5e-8;%5e-8
        cfg.tstep = 5e-8;%5e-8
        cfg.srcdir = [0 0 1];
        cfg.gpuid=4;
        cfg.srctype = 'line';
        cfg.srcparam1 = [0,Ny-10,0]; %size of light source
        cfg.isreflect = 0;
        cfg.autopilot = 1;
        cfg.debuglevel = 'P';
        
        cfg.outputtype = 'flux';
        flux = mcxlab(cfg);
        
        % convert mcx solution to mcxyz's output
        % 'energy': mcx outputs normalized energy deposition, must convert
        % it to normalized energy density (1/cm^3) as in mcxyz
        % 'flux': cfg.tstep is used in mcx's fluence normalization, must
        % undo 100 converts 1/mm^2 from mcx output to 1/cm^2 as in mcxyz
        if (strcmp(cfg.outputtype, 'energy'))
            mcxdata = flux.data / ((cfg.unitinmm / 10)^3);
        else
            mcxdata = flux.data * 100;
        end
        
        if (strcmp(cfg.outputtype, 'flux'))
            mcxdata = mcxdata * cfg.tstep;
        end
        
%         figure;
%         dim = size(cfg.vol);
%         yi = ((1:dim(2)) - floor(dim(2) / 2)) * cfg.unitinmm;
%         zi = (1:dim(3)) * cfg.unitinmm;
% 
%             imagesc(abs(squeeze(mcxdata(Nx/2, :, 25:end)))');
%             mcxdata(Nx/2,75,125)
%             figure;plot(mcxdata(Nx/2,:,125));
%             mcxdata(Nx/2,75,175)
%             figure;plot(mcxdata(Nx/2,:,175));
        % axis equal;
        % colormap(jet);
        % colorbar;
        if (strcmp(cfg.outputtype, 'energy'))
            set(gca, 'clim', [-2.4429 4.7581]);
        else
            set(gca, 'clim', [0.5 2.8]);
        
        end
        
        mcxdata = mcxdata(Nx/2-1:Nx/2+1,:,:);
        filename = sprintf('/home/zz111/Spectral_unmixing/mcx_result/TestPhantom_sTO2%.1f_sBO2%.1f_wavelength%d_seed%d.mat', sTO2, sBO2, wavelength, cfg.seed);        
        save(filename, "mcxdata");
        
        toc


    end

end
end
