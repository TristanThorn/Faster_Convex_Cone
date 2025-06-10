% 19/11/12: Improvements 

% ACTUAL VERSION: 15.9.2011 
% Frequency domain reconstruction of 2D image from OA signals acquired 
% with a linear array transducer. 
% This algorithm is based on the paper: 
% "Fourier reconstruction in optoacoustic imaging using complex k-space
% interpolation" by Jaeger et al., Inverse Problems 23 (2007) 
% 
% bugs corrected 15.8.11: 
% - when samplingX>1, Xextent needs to be adapted, too 

% ------------------------------------------------------------------------
% Parameter list (for further description see documentation): 
% sig: OA signal matrix, transducer element index corresponds to columns 
% F: sampling frequency (e.g. 40 MHz) at which signal was acquired  
% pitch: transducer element pitch (center of element distance) 
% delay: time delay from laser irradiation to signal acquisition start 
% zeroX: 1 = zero pad in lateral (x) direction; 0 = don't  
% zeroT: 1 = zero pad in axial (t, time) direction; 0 = don't 
% coeffT: defines from how many signal fourier coefficients a single 
%         image fourier coefficient is interpolated (typically 5) 
% samplingX: defines how many image lines are reconstructed per transducer
%         element. For value>1, the additional image lines are
%         equidistantly placed between the transducer elements 
% ------------------------------------------------------------------------
% Output list: 
% rekon: resulting image, has equal number of samples in axial (z)
%        direction as the signal. The number of image lines is: 
%               (number_of_elements - 1)*samplingX + 1 
% rekonuncut: if the zeroX flag was set, then this variable contains the
%        full image corresponding to the full (virtual) size of the padded
%        aperture, while rekon contains only the part of the image
%        corresponding to the real aperture 
% ------------------------------------------------------------------------
% typical usage (for Kontron transducer): 
% rekon = rekon_OA_freqdom(signal,40,0.3,0,0,0,5,3);  

function [rekon,rekonuncut] = rekon_OA_freqdom(sig,F,pitch,c,delay,zeroX,zeroT,coeffT,samplingX) 

 

sig = sig'; % signal data is transposed in order to fit to the algorithm 
            % which was initially implemented for element index
            % corresponding to matrix lines 

%c = 1.5;    % [mm/us]; sound speed in tissue, might be adapted for real tissue 

% ------------------- DATA DIMENSIONS -------------------------------------
% dimension of the signal and image in number of samples 
X = size(sig,1); 
Z = size(sig,2); 
T = Z; 
% corresponding physical dimension of signal and image in [mm]
Xextent = X*pitch; 
Zextent = Z*c/F; 
Textent = T*c/F; 

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
end 

% ---------------------- RECONSTRUCTION ALGORITHM -------------------------
% 2D arrays containing kx and kz values of k-vectors of the image fourier space 
[kx,kz] = ndgrid([-ceil((X-1)/2):floor((X-1)/2)]/Xextent,...
                 [-ceil((Z-1)/2):floor((Z-1)/2)]/Zextent); 

% 2D array containing kt values of k-vectors of the signal fourier space  
kt = kz; 
             
% 1D array containing values of kt 
kate = [-ceil((Z-1)/2):floor((Z-1)/2)]/Zextent; 
% maximum value of kt 
katemax = ceil((Z-1)/2)/Zextent; 

% projection of (kx,kz) onto (kx,kt) describing the signal formation in the
% fourier space 
kt2 = - sqrt(kz.^2+kx.^2); % projection of kz onto kt 

% calculate the jacobiante of the projection of kz onto kt % DEFINITION JACOBIANTE MATRIX A
kt2(kz==0&kx==0) = 1; 
jakobiante = kz./kt2; 
jakobiante(kz==0&kx==0) = 1; 
kt2(kz==0&kx==0) = 0; 

% if kt2 out of bound of kt, take the modulo value. Improvement 19/11/12  
samplfreq = T/Textent; 
kt2 = mod(kt2+samplfreq/2,samplfreq) - samplfreq/2; 

% #########################################################################

% 2D fourier spectrum of signal data 
sigtrans = fft2(sig); 
sigtrans = fftshift(sigtrans); 
sigtrans(kz>0) = 0; % only negatively propagating waves are detected, 19/11/12  

% initialize array for image fourier coefficients 
ptrans = zeros(X,Z); 

% numbers of signal fourier coefficients used for interpolation of image
% fourier coefficient, in up and down direction 
nTup = ceil((coeffT-1)/2); 
nTdo = floor((coeffT-1)/2); 

% helper index vector covering the range of fourier amplitudes going 
% into interpolation  
ktrange = [-nTdo:nTup]; 

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

ptrans(kt>0) = 0; % only negative kt are valid, 19/11/12 

% if the signal acquisition was started at a time different to the
% excitation time, fourier amplitudes have to be compensated for the
% corresponding phase 
% ptrans
ptrans = ptrans.*exp(-2*pi*1i*kt2*delay*c); 
ptrans = ptrans.*exp(2*pi*1i*kz*delay*c); 

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
    p = p(deltaX + [1:X],:); 
end 

% if lateral subsampling was used, the resulting image must be laterally
% filtered in order to reduce undersampling artifacts 
if samplingX > 1 
    p = conditioner(p,samplingX/2,1,0); 
    puncut = conditioner(puncut,1,1,0); 
end

rekon = p';

rekon = abs(hilbert(rekon));
% rekon = (rekon.^2);


% bmode = log10_compression2(rekon,20,-20,1); 

% mini=min(min(sO2));
% maxi=max(max(sO2));
% max_ref= max(max(sO2));
% sO2     = sO2/max_ref;                 % noramlization for 0-dB
% tl     = 10^(-30/20);               % Thershold level
% tp     = find(rekon < tl);              % Threshold positions
% A(tp)  = tl*ones(size(tp));         % Fill thershold position with threshold level
% A_th   = rekon;                         % Copy to threshold level
% rekon  = 20*log10(rekon);               % Log compression

% figure;
%Original
% imagesc((1:64)*0.245,(1:length(rekon))*(c/F),abs(bmode));colormap(hot),colorbar,x= xlabel('x(mm)'),y= ylabel('z(mm)');
% imagesc((1:128)*0.3,(1:length(rekon))*(c/F),abs(rekon));colormap(hot),colorbar,xlabel('Lateral Distance (mm)'),ylabel('Axial Depth(mm)');


%imagesc((1:64)*0.245,(1:length(log10(rekon)))*(c/F),abs(log10(rekon)));colormap(hot),colorbar,xlabel('Lateral Distance (mm)'),ylabel('Axial Depth(mm)');

title('PA Image ');
% axis equal;
% axis tight;
% set(gca,'FontWeight','bold');
% set(gca,'FontSize',14);
% set([x],'fontweight','bold','fontsize',14);
% set([y],'fontweight','bold','fontsize',14);
%     %z = title('LineProfile-Vessel');
% set(z,'fontweight','bold','fontsize',14);

rekonuncut = puncut'; 
