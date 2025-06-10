%Fast Fourier Reconstruction using inverse k-space interpolation

% ACUTAL VERSION FOR ANGLED TRANSMISSION ECHO ULTRASOUND RECONSTRUCTION
% 22/11/2012
% following new method, where signal data is shifted in time, and kt2 is
% differently calculated than in common echo ultrasound reconstruction.
% signal data shift in time according to arrival time of TX transient at
% respective elements makes a two-step reconstruction process as used in
% previous method unnecessary.
%
% ACTUAL VERSION: 15.8.2011
% Frequency domain reconstruction of 2D image from US (echo ultrasound)
% signals acquired with a linear array transducer.
% This algorithm is based on the paper:
% "Fourier reconstruction in optoacoustic imaging using complex k-space
% interpolation" by Jaeger et al., Inverse Problems 23 (2007)
% However, due to the different nature of signal generation, there are
% significant differences to the OA algorithm.
%
% bugs corrected 15.8.11:
% - back rotation of image at end
% - when samplingX>1, Xextent needs to be adapted, too
% - calculation of zero crossing point of TX transient with array was
% erroneous when samplingX>1. crossing point is now calculated at start
% of algorithm. Because it is a physical dimension, it does not afford
% adaptation with samplingX>1.

% ------------------------------------------------------------------------
% Parameter list (for further description see documentation):
% sig: US signal matrix, transducer element index corresponds to columns
% F: sampling frequency (e.g. 40 MHz) at which signal was acquired
% pitch: transducer element pitch (center of element distance)
% delay: time delay from laser irradiation to signal acquisition start
% zeroX: 1 = zero pad in lateral (x) direction; 0 = don't
% zeroT: 1 = zero pad in axial (t, time) direction; 0 = don't
% coeffT: defines from how many signal fourier coefficients a single
% image fourier coefficient is interpolated (typically 5)
% samplingX: defines how many image lines are reconstructed per transducer
% element. For value>1, the additional image lines are
% equidistantly placed between the transducer elements
% alpha: transmit angle
% ------------------------------------------------------------------------
% Output list:
% rekon: resulting image, has equal number of samples in axial (z)
% direction as the signal. The number of image lines is:
% (number_of_elements - 1)*samplingX + 1
% rekonuncut: if the zeroX flag was set, then this variable contains the
% full image corresponding to the full (virtual) size of the padded
% aperture, while rekon contains only the part of the image
% corresponding to the real aperture
% ------------------------------------------------------------------------
% typical usage (for Kontron transducer):
% rekon = rekon_us_freqdom1(rf_signal,50,0.245,1.54,1.33,0,5,1,0);
function rekon = rekon_us_freqdom1(sig,F,pitch,c,delay,zeroX,zeroT,coeffT,samplingX,alpha)

sig = sig'; % signal data is transposed in order to fit to the algorithm
            % which was initially implemented for element index
            % corresponding to matrix lines
            
% c = 1.54; % [mm/us]; sound speed in tissue, might be adapted for real tissue
zquell = 0; % [mm]; depth position where

% ------------------- DATA DIMENSIONS -------------------------------------
% dimension of the signal and image in number of samples
X = size(sig,1);
Z = size(sig,2);
T = Z;
% corresponding physical dimension of signal and image in [mm]
Xextent = X*pitch;
Zextent = Z*c/F/2;
Textent = T*c/F;
crossingpoint = Xextent/2;

% ------------------ DATA PRE-CONDITIONING --------------------------------
% lateral oversampling of the image, compared to real element pitch.
% this measure leads to finer images.
% for this purpose, zero signals are inserted between real signals
% and the element pitch is correspondingly reduced.
if samplingX>1
    sig2 = zeros(samplingX*(X-1)+1,Z); % a zero array is generated having
    % the appropriate size
    sig2([0:X-1]*samplingX+1,:) = sig(:,:); % signal lines corresponding to
    % real elements are filled
    pitch = pitch/samplingX; % the pitch is reduced
    X = samplingX*(X-1)+1; % the number of samples is increased
    Xextent = X*pitch; % the aperture is adapted, bug corrected 15.8.11
    sig = sig2;
end

% zero padding of the signal data in T-direction, by a factor of two.
% this measure reduces aliasing artifacts when a strong OA source is
% located either close to the start or to the end of the image frame.
% generally the influence is quite low, if coeffT is set to a large enough
% value (5).

if zeroT
    sig = [sig zeros(X,T)]; % the signal matrix is padded to double T-size
    Z = 2*Z; Zextent = 2*Zextent; % dimensions of image frame are doubled
    T = 2*T; Textent = 2*Textent; % dimensions of signal frame are doubled
end
% zero padding of the signal data in X-direction, by a factor of two.
% this measure reduces aliasing artifacts when a strong OA source is
% located close to either side of the image frame. if a source is located
% outside the image frame, it may be correctly reconstructed, leading to
% less artifacts.
if zeroX
    deltaX = round(X/2);
    deltaXextent = Xextent/2;
    sig = [zeros(deltaX,T); sig; zeros(deltaX,T)];
    X = X + 2*deltaX; Xextent = Xextent + 2*deltaXextent;
    crossingpoint = crossingpoint + deltaXextent;
end

% rearrange signal data according to transmission angle
xses = [0:X-1]*pitch - crossingpoint;
tses = xses/c * sin(alpha); % delay profile

for lineind = 1:X
    delayind = round(tses(lineind)*F);
    sig(lineind,:) = circshift(sig(lineind,:),[0 -delayind]);
end

% ---------------------- RECONSTRUCTION ALGORITHM -------------------------
% 2D arrays containing kx and kt values of k-vectors of the signal fourier space
[kx,kt] = ndgrid((-ceil((X-1)/2):floor((X-1)/2))/Xextent,...
                    (-ceil((Z-1)/2):floor((Z-1)/2))/Textent);

% 2D array containing kz values of k-vectors of the image fourier space
kz = ones(X,1) * (-ceil((Z-1)/2):floor((Z-1)/2))/Zextent;

% 1D array containing values of kt
kate = (-ceil((Z-1)/2):floor((Z-1)/2))/Textent;
% maximum value of kt
katemax = ceil((Z-1)/2)/Textent;

% projection of (kx,kz) onto (kx,kt) describing the signal formation in the
% fourier space
% ####### the projection law is the main difference to OA algorithm #######
warning off MATLAB:divideByZero

kk = sqrt(kx.^2+kz.^2);
gamma = asin(kx./kk);

kt2 = sqrt(kx.^2+kz.^2) .* cos(gamma - alpha)./(1 + cos(2*gamma-2*alpha));

kt2(isnan(kt2)) = 0;
kt2(isinf(kt2)) = 0;
kt2(cos(2*gamma - alpha)<=0) = 0; % these components can not be detected

% calculate the jacobiante of the projection of kz onto kt
kz(kt==0&kx==0) = 1;
jakobiante = 2*(0.5 - 0.5 * kx.^2./kz.^2);
jakobiante(kt==0&kx==0) = 1;
jakobiante(cos(2*gamma - alpha)<=0) = 0; % these components can not be detected
kz(kt==0&kx==0) = 0;

% % if kt2 out of bound of kt_prime, take the modulo value. Improvement 19/11/12
% samplfreq = T/Textent;
% kt2 = mod(kt2+samplfreq/2,samplfreq) - samplfreq/2;

% if kt2 out of bound of kt, take the modulo value
kt2(kt2*Textent+ceil((T-1)/2)+1<0.5) =...
    kt2(kt2*Textent+ceil((T-1)/2)+1<0.5) + T/Textent;
kt2(kt2*Textent+ceil((T-1)/2)+1>=T+0.5) =...
    kt2(kt2*Textent+ceil((T-1)/2)+1>=T+0.5) - T/Textent;

% #########################################################################

% 2D fourier spectrum of signal data
sigtrans = fft2(sig);
sigtrans = fftshift(sigtrans);
sigtrans(kt<0) = 0;

% initialize array for image fourier coefficients
ptrans = zeros(X,Z);

% numbers of signal fourier coefficients used for interpolation of image
% fourier coefficient, in up and down direction
nTup = ceil((coeffT-1)/2);
nTdo = floor((coeffT-1)/2);

% helper index vector covering the range of fourier amplitudes going
% into interpolation
%ktrange = [0:coeffT-1];
ktrange = -nTdo:nTup;

% loop for interpolation of image fourier coefficients from signal fourier
% coefficients
for xind = 1:X % for each image line do
    % calculate kt-index into kt-dimension of signal fourier space, from
    % the kt value resulting from the projection of kz to kt
    ktind = round(kt2(xind,:)*Textent)+ ceil((T-1)/2) + 1;
    ktind = ktind'*ones(1,coeffT) + ones(T,1)*ktrange;
    
    ktind(ktind>T) = ktind(ktind>T)-T;
    ktind(ktind<1) = ktind(ktind<1)+T;
    
    % V is a helper matrix, which contains for a given kx all the
    % kt-related signal fourier amplitudes which then will be used for
    % interpolation of all kz-related image fourier amplitudes
    V = sigtrans(xind + (ktind-1)*X);
    % Kt is a helper matrix, which contains for a given kx the kt values
    % corresponding to all the kt-related signal fourier amplitudes which
    % will be used for interpolation
    Kt = kt(xind + (ktind-1)*X);
    
    % calculate the complex interpolation weights for interpolation of
    % kz-related image fourier amplitudes from the corresponding multiples
    % of kt-related signal fourier amplitudes (based on paper by Jaeger et al.)
    deltakt = kt2(xind,:)'*ones(1,coeffT) - Kt; % the distance between kt2 and kt
    coeff = ones(T,coeffT); % prepare the coefficient matrix
    help = deltakt~=0;
    coeff(help) = (1 - exp(-2*pi*1i*deltakt(help)*Textent))./...
    (2*pi*1i*deltakt(help)*Textent);

    % interpolation as a weighted sum over multiple of kt-related signal
    % fourier amplitudes, resulting in the kz-related image fourier
    % amplitudes
    ptrans(xind,:) = sum(V.*coeff,2).'.*jakobiante(xind,:);
end

ptrans(kz<0) = 0;

% if the signal acquisition was started at a time different to the
% excitation time, fourier amplitudes have to be compensated for the
% corresponding phase
ptrans = ptrans.*exp(-2*pi*1i*kt2*(delay*c+zquell));
ptrans = ptrans.*exp(2*pi*1i*kz*0.5*(delay*c+zquell));

% inverse fourier transformation of image fourier amplitudes to get the
% image
p = real(ifft2(ifftshift(ptrans)));

% if zeropadding in t-direction was used, the image is cropped to the frame
% corresponding to the initial data
if zeroT
    Z = Z/2; Zextent = Zextent/2;
    p = p(:,1:Z);
    T = T/2; Textent = Textent/2;
end

% if zeropadding in x-direction was used, the image is cropped to the frame
% corresponding to the initial data
puncut = p;
if zeroX
    X = X - 2*deltaX; Xextent = Xextent - 2*deltaXextent;
    p = p(deltaX + (1:X),:);
end

% if lateral subsampling was used, the resulting image must be laterally
% filtered in order to reduce undersampling artifacts
if samplingX > 1
p = conditioner(p,samplingX/2,1,0);
puncut = conditioner(puncut,1,1,0);
end

rekon = p';
% figure,
rekon = hilbert(abs(rekon));
% rekon = rekon/max(max(rekon));
% imagesc((1:128)*0.3,(1:length(rekon))*(c/F/2),abs(rekon)),colormap(gray),colorbar,xlabel('x(mm)'),ylabel('z(mm)'),

% title('US Image');

colorbar;
% axis equal;
% axis tight;
% plot(rekon(:,33)./20+33,'--r')
% figure,imagesc(sig'),colormap(hot),colorbar
% figure(5),hold on
% plot(rekon(:,25),'--r')
% sig=sig';
% plot(rf_signal(:,33)','-r')
rekonuncut = puncut';

end
