function [opt_compensation, corrected, so2] = bi_compensate_core(points, refSpectra, dataPattern, varargin)
%BI_COMPENSATE_CORE  Apply angular distance compensation to PA spectra.
%   [OPT_COMP, CORRECTED, SO2] = BI_COMPENSATE_CORE(POINTS, REFSPECTRA, PATTERN)
%   performs the compensation using measurement POINTS and reference spectra.
%
%   POINTS is a struct with fields:
%       Ni    - [x y] coordinates of the NiSO4 location
%       Cu    - [x y] coordinates of the CuSO4 location
%       Others - N-by-2 matrix of additional measurement points (optional)
%
%   REFSPECTRA is a struct with fields 'Ni' and 'Cu' that contain the
%   reference spectra sampled at the selected wavelengths.
%
%   DATAPATTERN is a sprintf pattern to load reconstructed images, e.g.
%   '0520/PA_deep_recon_%d.mat'. The placeholder will be replaced with the
%   wavelength value.
%
%   Name-Value pairs:
%       'SelectWv'     - wavelength indices (default 3:2:21)
%       'Window'       - [row col] averaging window size (default [4 8])
%       'UnmixSpectra' - {spec1, spec2} perform linear unmixing using these
%                        reference spectra (optional)
%       'ShowPlots'    - true/false to display plots (default true)
%
%   The function returns OPT_COMPENSATION, the corrected spectra structure
%   CORRECTED with fields Ni, Cu and Others, and optionally the SO2 values
%   from linear unmixing if 'UnmixSpectra' is provided.

% Parse inputs
p = inputParser;
addParameter(p, 'SelectWv', 3:2:21);
addParameter(p, 'Window', [4 8]);
addParameter(p, 'UnmixSpectra', []);
addParameter(p, 'ShowPlots', true);
parse(p, varargin{:});
opt = p.Results;

select_wv = opt.SelectWv;
wv_all = 700:10:900;
wv_values = wv_all(select_wv);
num_wavelengths = numel(wv_values);

% Preallocate
snr_Ni  = zeros(1, num_wavelengths);
snr_Cu  = zeros(1, num_wavelengths);
num_other = 0;
if isfield(points, 'Others') && ~isempty(points.Others)
    num_other = size(points.Others,1);
    snr_other = zeros(num_other, num_wavelengths);
else
    snr_other = [];
end

row_off = opt.Window(1);
col_off = opt.Window(2);

% Read images and accumulate mean signals
for idx = 1:num_wavelengths
    data = load(sprintf(dataPattern, wv_values(idx)));
    img = data.img;

    region_Ni = img(points.Ni(2)-row_off:points.Ni(2)+row_off, ...
                    points.Ni(1)-col_off:points.Ni(1)+col_off);
    snr_Ni(idx) = mean(region_Ni(:));

    region_Cu = img(points.Cu(2)-row_off:points.Cu(2)+row_off, ...
                    points.Cu(1)-col_off:points.Cu(1)+col_off);
    snr_Cu(idx) = mean(region_Cu(:));

    for pIdx = 1:num_other
        pt = points.Others(pIdx,:);
        region = img(pt(2)-row_off:pt(2)+row_off, pt(1)-col_off:pt(1)+col_off);
        snr_other(pIdx, idx) = mean(region(:));
    end
end

normalize = @(x) x / norm(x);
objective = @(comp) ...
    norm(normalize(snr_Ni ./ comp) - normalize(refSpectra.Ni))^2 + ...
    norm(normalize(snr_Cu ./ comp) - normalize(refSpectra.Cu))^2 + ...
    0.001 * norm(diff(comp))^2;

init_comp = ones(1, num_wavelengths);
lb = 0.1 * ones(size(init_comp));
ub = 10  * ones(size(init_comp));
options = optimoptions('fmincon', 'Display', 'none', ...
    'MaxFunctionEvaluations', 5000, 'MaxIterations', 1000);
opt_compensation = fmincon(objective, init_comp, [], [], [], [], lb, ub, [], options);

% Apply compensation
corrected.Ni = snr_Ni ./ opt_compensation;
corrected.Cu = snr_Cu ./ opt_compensation;
corrected.Others = snr_other ./ opt_compensation;

% Plot results if requested
if opt.ShowPlots
    plot_results(wv_values, snr_Ni, snr_Cu, corrected, refSpectra, opt_compensation);
end

% Optional linear unmixing
so2 = [];
if ~isempty(opt.UnmixSpectra)
    if exist('linearUnmixing', 'file') ~= 2
        addpath('Toolbox');
    end
    spectra1 = opt.UnmixSpectra{1};
    spectra2 = opt.UnmixSpectra{2};
    so2 = zeros(num_other,1);
    for pIdx = 1:num_other
        so2(pIdx) = linearUnmixing(corrected.Others(pIdx,:), ...
                                   spectra1(select_wv), spectra2(select_wv));
    end
end
end

function plot_results(wv, snr_Ni, snr_Cu, corrected, refSpectra, comp)
    normalize = @(x) x / norm(x);

    figure;
    subplot(2,1,1);
    plot(wv, normalize(snr_Ni), '--c', 'LineWidth', 1); hold on;
    plot(wv, normalize(corrected.Ni), '-bo', 'LineWidth', 1.5);
    plot(wv, normalize(refSpectra.Ni), '-k', 'LineWidth', 2);
    legend({'Original SNR_{Ni}', 'Corrected Ni', 'True Ni'}, 'Location', 'Best');
    xlabel('Wavelength (nm)'); ylabel('Normalized Amplitude');
    title('Normalized NiSO_4 Spectrum Comparison');
    grid on;

    subplot(2,1,2);
    plot(wv, normalize(snr_Cu), '--m', 'LineWidth', 1); hold on;
    plot(wv, normalize(corrected.Cu), '-ro', 'LineWidth', 1.5);
    plot(wv, normalize(refSpectra.Cu), '-k', 'LineWidth', 2);
    legend({'Original SNR_{Cu}', 'Corrected Cu', 'True Cu'}, 'Location', 'Best');
    xlabel('Wavelength (nm)'); ylabel('Normalized Amplitude');
    title('Normalized CuSO_4 Spectrum Comparison');
    grid on;

    figure;
    plot(wv, comp, 'g-o', 'LineWidth', 2);
    xlabel('Wavelength (nm)');
    ylabel('Compensation Factor');
    title('Optimized Compensation Vector');
    grid on;
end
