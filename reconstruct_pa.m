function img = reconstruct_pa(rf_data, params)
%RECONSTRUCT_PA Reconstruct a PA image from RF data using rekon\_OA\_freqdom.
%   IMG = RECONSTRUCT_PA(RF_DATA, PARAMS) reconstructs a photoacoustic image
%   using the specified parameters structure.  Missing parameters are
%   substituted with commonly used defaults.
%
%   Required fields of PARAMS (with defaults):
%       F         - Sampling frequency in MHz (default 40)
%       pitch     - Transducer pitch in cm (default 0.315)
%       c         - Sound speed (default 1.48)
%       delay     - Acquisition delay (default 0)
%       zeroX     - Zero pad in lateral direction (default 1)
%       zeroT     - Zero pad in axial direction (default 1)
%       coeffT    - Number of coefficients for interpolation (default 5)
%       samplingX - Image lateral oversampling factor (default 8)

    if nargin < 2, params = struct(); end

    if ~isfield(params, 'F'),         params.F = 40;       end
    if ~isfield(params, 'pitch'),     params.pitch = 0.315; end
    if ~isfield(params, 'c'),         params.c = 1.48;     end
    if ~isfield(params, 'delay'),     params.delay = 0;    end
    if ~isfield(params, 'zeroX'),     params.zeroX = 1;    end
    if ~isfield(params, 'zeroT'),     params.zeroT = 1;    end
    if ~isfield(params, 'coeffT'),    params.coeffT = 5;   end
    if ~isfield(params, 'samplingX'), params.samplingX = 8; end

    [img, ~] = rekon_OA_freqdom(rf_data, params.F, params.pitch, params.c, ...
        params.delay, params.zeroX, params.zeroT, params.coeffT, params.samplingX);
end
