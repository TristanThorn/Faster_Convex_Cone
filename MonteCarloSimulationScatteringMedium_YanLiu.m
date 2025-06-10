% clear all 
% format compact
% % Set random seed for reproducibility
% s = RandStream('mt19937ar','Seed','shuffle');
% RandStream.setGlobalStream(s);
% % Simulation parameters
% N = 1e7;             % Max number of photon packets
% g = 0.9;             % Anisotropy factor
% us = 100;            % Scattering coefficient [cm^-1]
% ua = 0.1;            % Absorption coefficient [cm^-1]
% dz = 1e-4;           % Grid resolution along z [cm]
% dr = 1e-5;           % Grid resolution along r [cm]
% Nr = 5001;           % Number of radial bins
% n1 = 1.33;           % Refractive index
% NA1 = 0.1;           % Numerical aperture relative to air
% NA = NA1 * n1;       % Absolute NA
% lambda = 0.57e-4;    % Wavelength [cm]
% ut = us + ua;                       % Total attenuation coefficient
% us_prime = us * (1 - g);            % Reduced scattering
% ut_p = us_prime + ua;              % Reduced total coefficient
% % Beam and depth setup
% r_spot = 0.3183 * lambda / NA;     % Gaussian beam 1/e^2 radius
% zf_list = 0.03:0.02:0.17;          % Focal depths [cm]
% colors = lines(length(zf_list));   % Color map for plotting
% fluence_curves = cell(length(zf_list), 1);  % To store results
% figure;
% hold on;
% for idx = 1:length(zf_list)
%     zf = zf_list(idx);
%     % Adjust photon number for shallow focus
%     if zf < 0.07
%         N_current = 1e6;
%     else
%         N_current = 5e6;
%     end
%     fprintf('Simulating zf = %.2f cm with N = %.0e photons...\n', zf, N_current);
%     tic;
%     Nz = ceil((zf + 0.02) / dz);   % Number of z grid points
%     % Preallocate matrices
%     A_zr = zeros(Nz, Nr);  % Absorbed energy
%     n_zr = zeros(Nz, Nr);  % Number of steps
%     S_zr = zeros(Nz, Nr);  % Scattering events
%     % Beam focusing parameters
%     z0 = pi * r_spot^2 / lambda;
%     k = zf / z0;
%     k1 = k^2 / (k^2 + 1);
%     k2 = 1 + k^2;
%     r_spot_sqrt_k2 = r_spot * sqrt(k2);
%     photon_progress_step = round(N_current / 10);  % Show progress every 10%
%     for n = 1:N_current
%         if mod(n, photon_progress_step) == 0
%             fprintf('  Photon %d / %d (%.0f%%)\n', n, N_current, 100 * n / N_current);
%         end
%         W = 1;  % Photon weight
%         r = sqrt(-log(rand)/2) * r_spot_sqrt_k2;
%         ag = 2*pi*rand;
%         x = r * cos(ag);
%         y = r * sin(ag);
%         z = 0;
%         % Direction initialization
%         if rand < 0.5
%             X = -k*(k*x + y)/k2;
%             Y = k*(x - k*y)/k2;
%         else
%             X = k*(y - k*x)/k2;
%             Y = -k*(k*y + x)/k2;
%         end
%         ulength = sqrt(k1 * r^2 + zf^2);
%         ux = X / ulength;
%         uy = Y / ulength;
%         uz = zf / ulength;
%         S1 = 0;
%         s_ = 0;
%         dead = false;
%         while ~dead
%             if s_ == 0
%                 s_ = -log(rand);
%             end
%             s = s_ / ut;
%             % Update position
%             if z + s * uz > 0
%                 x = x + ux * s;
%                 y = y + uy * s;
%                 z = z + uz * s;
%                 r = sqrt(x^2 + y^2);
%                 % Energy absorbed
%                 d_W = ua / ut * W;
%                 W = W - d_W;
%                 z_idx = ceil(z / dz);
%                 z_idx = min(max(z_idx, 1), Nz);
%                 r_idx = ceil(r / dr);
%                 r_idx = min(max(r_idx, 1), Nr);
%                 A_zr(z_idx, r_idx) = A_zr(z_idx, r_idx) + d_W;
%                 S_zr(z_idx, r_idx) = S_zr(z_idx, r_idx) + S1;
%                 n_zr(z_idx, r_idx) = n_zr(z_idx, r_idx) + 1;
%                 S1 = S1 + 1;
%                 s_ = 0;
%                 % Scattering update (Henyey-Greenstein)
%                 temp1 = (1 - g^2) / (1 - g + 2 * g * rand);
%                 cos_s = (1 + g^2 - temp1^2) / (2 * g);
%                 sin_s = sqrt(1 - cos_s^2);
%                 phi = 2 * pi * rand;
%                 cos_phi = cos(phi);
%                 sin_phi = sin(phi);
%                 if abs(uz) < 0.99999
%                     temp = sqrt(1 - uz^2);
%                     ux_p = sin_s * (ux * uz * cos_phi - uy * sin_phi) / temp + ux * cos_s;
%                     uy_p = sin_s * (uy * uz * cos_phi + ux * sin_phi) / temp + uy * cos_s;
%                     uz_p = -temp * sin_s * cos_phi + uz * cos_s;
%                 else
%                     ux_p = sin_s * cos_phi;
%                     uy_p = sin_s * sin_phi;
%                     uz_p = sign(uz) * cos_s;
%                 end
%                 ux = ux_p; uy = uy_p; uz = uz_p;
%             else
%                 dead = true;
%             end
%             % Russian roulette
%             if (W < 1e-4) && ~dead
%                 if rand < 0.1
%                     W = W * 10;
%                 else
%                     dead = true;
%                 end
%             end
%         end
%     end
%     % Normalize absorption to get fluence
%     A_zr1 = A_zr ./ (ones(Nz,1) * (0.5:1:(Nr-0.5))) / N_current / dz / dr^2 / (2 * pi);
%     z_idx = round(zf / dz);
%     fluence = A_zr1(z_idx,:);
%     fluence_sym = 0.5 * (fluence + fliplr(fluence));     % Symmetrize
%     fluence_norm = fluence_sym / max(fluence_sym);       % Normalize
%     % Generate symmetric radial profile [micron]
%     r_um_half = ((1:Nr)-0.5) * dr * 1e4;
%     r_um_full = [-fliplr(r_um_half(1:1000)), r_um_half(1:1000)];
%     fluence_full = [fliplr(fluence_norm(1:1000)), fluence_norm(1:1000)];
%     % Store curve
%     fluence_curves{idx} = struct( ...
%         'zf', zf, ...
%         'r_um', r_um_full, ...
%         'fluence', fluence_full ...
%     );
%     % Plot curve
%     plot(r_um_full, fluence_full, '-', 'LineWidth', 1.5, 'Color', colors(idx,:)); hold on;
%     elapsed = toc;
%     fprintf('Finished zf = %.2f cm in %.1f sec.\n\n', zf, elapsed);
%     filename = sprintf('fluence_zf_%03dum.mat', round(zf*1e4));  % e.g., zf = 0.01 → 0010
%     save(filename, 'A_zr', 'zf');
% end
% xlabel('Lateral position [\mum]');
% ylabel('Normalized fluence');
% title('Symmetrized lateral fluence at focal plane for different z_f');
% legend(arrayfun(@(zf) sprintf('z_f = %.2f cm', zf), zf_list, 'UniformOutput', false));
% grid on;

clear; close all;

%% Common parameters
dz = 1e-4;        % z resolution [cm]
dr = 1e-5;        % r resolution [cm]
Nr = 5001;        % number of radial bins
g = 0.9;          % anisotropy
us = 100;         % scattering coeff [cm^-1]
ua = 0.1;         % absorption coeff [cm^-1]
ut = us + ua;     % total coeff
dz = 1e-4; dr = 1e-5;

% beam
n1 = 1.33; NA1 = 0.1;
NA = NA1 * n1;
lambda = 0.57e-4;          
r_spot = 0.3183 * lambda / NA;

zf_list = 0.01:0.08:0.17;     % focal depths      
colors = lines(length(zf_list));
figure; hold on;

for idx = 1:length(zf_list)
    zf = zf_list(idx);
    % choose N depending on depth
    if zf < 0.07, Np = 1e6; else Np = 5e6; end

    Nz = ceil((zf + 0.02) / dz);
    A_zr = zeros(Nz, Nr);
    n_zr = zeros(Nz, Nr);
    S_zr = zeros(Nz, Nr);

    % focusing parameters
    z0 = pi * r_spot^2 / lambda;
    k = zf / z0; k1 = k^2/(k^2+1); k2 = 1 + k^2;
    r_mult = r_spot * sqrt(k2);

    for p = 1:Np
        % launch photon
        W = 1; S1 = 0; s_ = 0; dead = false;
        r0 = sqrt(-log(rand)/2) * r_mult;
        theta = 2*pi*rand; x = r0*cos(theta); y = r0*sin(theta); z = 0;
        if rand<0.5
            X = -k*(k*x+y)/k2; Y = k*(x-k*y)/k2;
        else
            X = k*(y-k*x)/k2; Y = -k*(k*y+x)/k2;
        end
        L = sqrt(k1*r0^2 + zf^2);
        ux = X/L; uy = Y/L; uz = zf/L;

        while ~dead
            if s_==0, s_ = -log(rand); end
            s = s_/ut;            
            if z + s*uz <= 0
                dead = true; break;
            end
            % move
            x = x + ux*s; y = y + uy*s; z = z + uz*s;
            r = sqrt(x^2+y^2);
            % absorb
            dW = ua/ut * W; W = W - dW;
            iz = min(max( ceil(z/dz),1), Nz);
            ir = min(max( ceil(r/dr),1), Nr);
            A_zr(iz,ir) = A_zr(iz,ir) + dW;
            S_zr(iz,ir) = S_zr(iz,ir) + S1;
            n_zr(iz,ir) = n_zr(iz,ir) + 1;
            S1 = S1 + 1; s_ = 0;
            % scatter
            tmp = (1-g^2)/(1-g+2*g*rand);
            cos_s = (1+g^2-tmp^2)/(2*g);
            sin_s = sqrt(1-cos_s^2);
            phi = 2*pi*rand; cphi=cos(phi); sphi=sin(phi);
            if abs(uz)<0.99999
                t = sqrt(1-uz^2);
                ux_p = sin_s*(ux*uz*cphi-uy*sphi)/t + ux*cos_s;
                uy_p = sin_s*(uy*uz*cphi+ux*sphi)/t + uy*cos_s;
                uz_p = -t*sin_s*cphi + uz*cos_s;
            else
                ux_p = sin_s*cphi; uy_p = sin_s*sphi; uz_p = sign(uz)*cos_s;
            end
            ux=ux_p; uy=uy_p; uz=uz_p;
            % roulette
            if W<1e-4
                if rand<0.1, W=W*10; else dead=true; end
            end
        end
    end

    % compute mean scatter events at focal plane
    iz = round(zf/dz);
    N_scatter = S_zr(iz,:) ./ max(n_zr(iz,:),1);  % avoid /0

    % symmetrize around first positive sample (r=dr)
    center = 2;
    left = N_scatter(1:center-1);   right = N_scatter(center:end);
    L = min(numel(left),numel(right));
    half = 0.5*(left(end:-1:end-L+1) + right(1:L));
    scatter_full = [ half(end:-1:1), right(1:L) ];
    ru = (-(L-1):0)*(dr*1e4);
    rpos = (0:(L-1))*(dr*1e4);
    r_full = [ru, rpos];

    % plot
    plot(r_full, scatter_full, 'LineWidth',1.5,'Color',colors(idx,:));
end

xlabel('Radial position [\mum]');
ylabel('Number of scattering events');
title('Scattering events vs radial position for various z_f');
legend(arrayfun(@(z) sprintf('z_f=%.2f cm',z),zf_list,'UniformOutput',false));
grid on;

