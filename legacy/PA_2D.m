clear all;
close all;

gpuDevice(1)
% =========================================================================
% SIMULATION
% =========================================================================
addpath('/scratch/k-Wave-linux/k-wave-toolbox-version-1.4/');

% create the computational grid
PML_size = 10;              % size of the PML in grid points
Nx = 500 - 2 * PML_size;     % number of grid points in the x direction   4cm:1000
Ny = 900;                    % number of grid points in the y direction
dx = 0.05e-3;                % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1508;	% [m/s]

%% Define initial pressure parameters
sTO2 = 0.4;
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


load("digital_vessel_5vessel_size5_0.05.mat");

[vol_x, vol_y, vol_z] = size(vol);
vol_background = ones(vol_x,vol_y,vol_z+100); 
vol_background(:, :, 101:end) = vol(:,:,1:end) +1;


%% %%%%%%%%%%%%%%%%%%%%%%% build transducer arcs %%%%%%%%%%%%%%%%%%%%%%%



% set the input arguements
% load('mimic_company_data_reduceFOV_f2_cut_5.mat')
% sensor_data = reshape(raw_data5, 6151, 96*292);



for sBO2 = 0:0.1:1
    
    blood_mu_a = (spectrum_HbO2 .* sBO2 + spectrum_Hb .* (1 - sBO2)) ...
                .* (150/64500) * log(10);
    for i = 1:21
        wavelength = wavelengths(i);

        p0 = zeros(Nx, Ny);
        filename = sprintf('Simulation_size5/MC1e9_0.05_sTO2_%.1f_wavelength%d_sBO2_%.1f_seed0.mat', sTO2, wavelength, sBO2);
        load(filename);

        vol_a = vol_background;
        vol_a(vol_background == 1) = spectrum_H2O(i);
        vol_a(vol_background == 2) = Epid_mu_a(i);
        vol_a(vol_background == 3) = der_mu_a(i);
        vol_a(vol_background == 4) = subcu_mu_a(i);
        vol_a(vol_background == 5) = mus_mu_a(i);
        vol_a(vol_background == 6) = blood_mu_a(i);
% 
%          figure;imagesc(log10(squeeze(mcxdata(200, :, :)))');
%          figure;imagesc(log10(squeeze(vol_a(200, :, :)))');

        p_3D = mcxdata.*vol_a;
        p_2D = squeeze(p_3D(200,:,101:end))';
        object_0 = imresize(p_2D, 1); % size 0.1mm

        [m, n] = size(object_0); % 获取矩阵的维度
        
        % 计算中心区域的起始和结束索引
        start_x = floor((Nx-m) / 2) + 1; % 中心的起始索引（第一维）
        end_x = start_x + m - 1; % 中心的结束索引（第一维）
        
        start_y = floor((Ny-n) / 2) + 1; % 中心的起始索引（第二维）
        end_y = start_y + n - 1; % 中心的结束索引（第二维）
        
        p0(start_x:end_x, start_y:end_y)=object_0;
        % p0 = smooth(kgrid, p0, true);
%         figure;imshow(p0,[]);
        source.p0 = p0;

        % assign to the sensor structure
center_freq = 10e6;      % [Hz]
bandwidth = 80;         % [%]
sensor.frequency_response = [center_freq, bandwidth];

sensor_mask = zeros(Nx,Ny);
sensor_mask(2,:) = 1;

sensor.mask = sensor_mask;
        
        
        % run the simulation
        input_args = {'PMLSize', PML_size, 'PMLInside', false, 'PlotPML', false, ...
    'Smooth', true,'DataCast', 'gpuArray-single',  'PlotSim', false, 'CartInterp', 'nearest'};

        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        % reset the initial pressure
        source.p0 = 0;
        
        % assign the time reversal data

        sensor.time_reversal_boundary_data = sensor_data;
        
        % run the time-reversal reconstruction
        p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        p0_recon = abs(hilbert2D(p0_recon));
        save_filename = sprintf('pa_result/PA1e9_4cm_sTO2_%.1f_wavelength%d_sBO2_%.1f.mat', sTO2, wavelength, sBO2);
        save(save_filename, "p0_recon");
        clear sensor
%          figure; imshow(p0_recon,[]);
    end
end




function H = hilbert2D(data)

    F = fft2(data);
    
    % generate hilbert filter
    [nx, ny] = size(data);
    u = ifftshift(-fix(nx/2):ceil(nx/2)-1);
    v = ifftshift(-fix(ny/2):ceil(ny/2)-1);
    [U, V] = meshgrid(u, v);
    H_filter = double(U > 0);
    
    H = ifft2(F .* H_filter');
end