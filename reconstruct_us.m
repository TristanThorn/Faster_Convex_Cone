function img = reconstruct_us(rf_data, params)
%RECONSTRUCT_US Reconstruct a US image from RF data using rekon\_US\_freqdom1.
%   IMG = RECONSTRUCT_US(RF_DATA, PARAMS) reconstructs an ultrasound image
%   using the specified parameters structure.  Missing parameters are
%   substituted with commonly used defaults.
%
%   Required fields of PARAMS (with defaults):
%       F         - Sampling frequency in MHz (default 20)
%       pitch     - Transducer pitch in cm (default 0.315)
%       c         - Sound speed (default 1.48)
%       delay     - Acquisition delay (default 0)
%       zeroX     - Zero pad in lateral direction (default 1)
%       zeroT     - Zero pad in axial direction (default 1)
%       coeffT    - Number of coefficients for interpolation (default 5)
%       samplingX - Image lateral oversampling factor (default 8)
%       alpha     - Steering angle (default 0)

    if nargin < 2, params = struct(); end

    if ~isfield(params, 'F'),         params.F = 20;       end
    if ~isfield(params, 'pitch'),     params.pitch = 0.315; end
    if ~isfield(params, 'c'),         params.c = 1.48;     end
    if ~isfield(params, 'delay'),     params.delay = 0;    end
    if ~isfield(params, 'zeroX'),     params.zeroX = 1;    end
    if ~isfield(params, 'zeroT'),     params.zeroT = 1;    end
    if ~isfield(params, 'coeffT'),    params.coeffT = 5;   end
    if ~isfield(params, 'samplingX'), params.samplingX = 8; end
    if ~isfield(params, 'alpha'),     params.alpha = 0;    end

    img = rekon_US_freqdom1(rf_data, params.F, params.pitch, params.c, ...
        params.delay, params.zeroX, params.zeroT, params.coeffT, ...
        params.samplingX, params.alpha);
    img = real(img);
    img(isnan(img) | isinf(img)) = 0;
end
