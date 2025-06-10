function [opt_compensation, corrected_signals, lin_unmix] = bi_compensate_func(points, ref_spectra, data_template, varargin)
%BI_COMPENSATE_FUNC  Estimate and apply compensation factor for two reference spectra.
%   [COMP, CORRECTED, LIN_UNMIX] = BI_COMPENSATE_FUNC(POINTS, REF_SPECTRA, TEMPLATE)
%   performs angular distance compensation using measurement POINTS, reference
%   spectra REF_SPECTRA and image files specified by TEMPLATE. POINTS is a
%   structure with fields Ni and Cu giving [x y] pixel positions.  Additional
%   points can be supplied in POINTS.extras as an N-by-2 matrix. REF_SPECTRA must
%   contain fields Ni and Cu with the reference spectral values. TEMPLATE is a
%   sprintf string that should contain one integer placeholder for the
%   wavelength. Optional name/value pairs are:
%       'SelectWV'    - indices of wavelengths (default 3:2:21)
%       'Offsets'     - struct with fields matching POINTS specifying [row col]
%                       neighbourhood half sizes (default [4 8] for all)
%       'LinSpectra'  - {HbO2, Hb} to perform linear unmixing on extras
%
%   The function returns the optimized compensation vector, the corrected
%   signals (structure with same fields as POINTS) and, if linear unmixing is
%   requested, the unmixing results for each extra point.

% Parse optional inputs
p = inputParser;
p.addParameter('SelectWV', 3:2:21);
p.addParameter('Offsets', struct());
p.addParameter('LinSpectra', {});
p.parse(varargin{:});
select_wv = p.Results.SelectWV;
offsets = p.Results.Offsets;
lin_spec = p.Results.LinSpectra;

wv_all = 700:10:900;
wv_values = wv_all(select_wv);
num_wv = numel(wv_values);

% Prepare offsets
if ~isfield(offsets, 'Ni'), offsets.Ni = [4 8]; end
if ~isfield(offsets, 'Cu'), offsets.Cu = [4 8]; end
if isfield(points, 'extras')
    if ~isfield(offsets, 'extras')
        offsets.extras = repmat(offsets.Ni, size(points.extras,1),1);
    end
end

snr_Ni = zeros(1, num_wv);
snr_Cu = zeros(1, num_wv);
extra_signals = [];
if isfield(points, 'extras')
    extra_signals = zeros(size(points.extras,1), num_wv);
end

for idx = 1:num_wv
    data = load(sprintf(data_template, wv_values(idx)));
    img = data.img;

    snr_Ni(idx) = localAverage(img, points.Ni, offsets.Ni);
    snr_Cu(idx) = localAverage(img, points.Cu, offsets.Cu);
    if isfield(points, 'extras')
        for e = 1:size(points.extras,1)
            extra_signals(e,idx) = localAverage(img, points.extras(e,:), offsets.extras(e,:));
        end
    end
end

% optimisation
normalize = @(x) x / norm(x);
objective = @(comp) norm(normalize(snr_Ni./comp) - normalize(ref_spectra.Ni))^2 + ...
                     norm(normalize(snr_Cu./comp) - normalize(ref_spectra.Cu))^2 + ...
                     0.001 * norm(diff(comp))^2;
init_comp = ones(1,num_wv);
opts = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',5000,'MaxIterations',1000);
opt_compensation = fmincon(objective, init_comp, [], [], [], [], 0.1*ones(1,num_wv), 10*ones(1,num_wv), [], opts);

% Apply compensation
corrected_signals.Ni = snr_Ni ./ opt_compensation;
corrected_signals.Cu = snr_Cu ./ opt_compensation;
if isfield(points,'extras')
    for e = 1:size(points.extras,1)
        corrected_signals.extras(e,:) = extra_signals(e,:) ./ opt_compensation;
    end
end

% Plot
plotNormalizedSpectra(wv_values, snr_Ni, snr_Cu, corrected_signals, ref_spectra);
plotCompensation(wv_values, opt_compensation);

% Optional linear unmixing
lin_unmix = [];
if ~isempty(lin_spec) && isfield(points,'extras')
    addpath('Toolbox');
    lin_unmix = zeros(size(points.extras,1),1);
    for e = 1:size(points.extras,1)
        lin_unmix(e) = linearUnmixing(corrected_signals.extras(e,:), lin_spec{1}(select_wv), lin_spec{2}(select_wv));
    end
end
end

function val = localAverage(img, pt, off)
%LOCALAVERAGE Mean of a rectangular region with boundary clipping.
rows = max(pt(2)-off(1),1):min(pt(2)+off(1), size(img,1));
cols = max(pt(1)-off(2),1):min(pt(1)+off(2), size(img,2));
region = img(rows, cols);
val = mean(region(:));
end

function plotNormalizedSpectra(wv, snr_Ni, snr_Cu, corr, ref)
normalize = @(x) x / norm(x);
figure;
subplot(2,1,1);
plot(wv, normalize(snr_Ni), '--c','LineWidth',1); hold on;
plot(wv, normalize(corr.Ni), '-bo','LineWidth',1.5);
plot(wv, normalize(ref.Ni), '-k','LineWidth',2);
legend('Original SNR_{Ni}','Corrected Ni','True Ni','Location','Best');
xlabel('Wavelength (nm)'); ylabel('Normalized Amplitude');
title('Normalized NiSO4 Spectrum Comparison'); grid on;

subplot(2,1,2);
plot(wv, normalize(snr_Cu), '--m','LineWidth',1); hold on;
plot(wv, normalize(corr.Cu), '-ro','LineWidth',1.5);
plot(wv, normalize(ref.Cu), '-k','LineWidth',2);
legend('Original SNR_{Cu}','Corrected Cu','True Cu','Location','Best');
xlabel('Wavelength (nm)'); ylabel('Normalized Amplitude');
title('Normalized CuSO4 Spectrum Comparison'); grid on;
end

function plotCompensation(wv, comp)
figure;
plot(wv, comp,'g-o','LineWidth',2);
xlabel('Wavelength (nm)'); ylabel('Compensation Factor');
title('Optimized Compensation Vector'); grid on;
end
