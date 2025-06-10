close all;
figure;

% Parameters
folder_path = '0520/';
file_name = 'deep_80delay_OBP_Laser_PA_1356.raw'; % <-- Update to your file name
full_path = fullfile(folder_path, file_name);

fid = fopen(full_path, 'rb');
if fid == -1
    error('Cannot open file: %s', full_path);
end

data = fread(fid, 'int16');
fclose(fid);

% Assume 1024 samples and 128 channels per frame
frame_size = 1024 * 128;
total_frames = floor(numel(data) / frame_size);
body = reshape(data(1:frame_size*total_frames), [1024, 128, total_frames]);

% Sweep wavelengths from 700 to 900 nm (step = 20), 11 steps
wavelengths = 700:20:900;
frames_per_wavelength = 120;

if total_frames < frames_per_wavelength * length(wavelengths)
    error('Not enough frames in the file for %d wavelengths.', length(wavelengths));
end

for idx = 1:length(wavelengths)
    
    wv = wavelengths(idx);
    start_frame = (idx-1)*frames_per_wavelength + 1;
    end_frame = start_frame + frames_per_wavelength - 1;
    % Average selected frames
    rfm = mean(body(:, :, start_frame+10:end_frame-10), 3);
    % rfm = mean(body(:, :, frame), 3);

    % Zero out noisy region
    rfm(1:200, :) = 0;

    % Flip image left-right
    rfm = fliplr(rfm);

    % Reconstruct PA image
    [img, ~] = rekon_OA_freqdom(rfm, 40, 0.315, 1.48, 0, 1, 1, 5, 8);
    caxis([0, 1e3]);

    % Display and save image
    imshow(img, []);
    title(sprintf('Wavelength %d nm', wv));
    drawnow;

    save(sprintf('0520/PA_deep_recon_%d.mat', wv), 'img');
end
