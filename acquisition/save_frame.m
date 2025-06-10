close all;
figure;

% Parameters
folder_path = '0504/';
file_name = 'ml_70delay_OBP_Laser_PA_1344.raw'; % <-- Update to your file name
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

for idx = 2:length(wavelengths)
    idx
    wv = wavelengths(idx);
    start_frame = (idx-1)*frames_per_wavelength;
    end_frame = start_frame + frames_per_wavelength;
    
    img_time = zeros(1024, 1017, (end_frame-start_frame+1));
    
    frame_num = 1;
    for frame = start_frame: end_frame
        % Average selected frames
        rfm = mean(body(:, :, frame), 3);
    
        % Zero out noisy region
        rfm(1:200, :) = 0;
    
        % Flip image left-right
        rfm = fliplr(rfm);
    
        % Reconstruct PA image
        [img, ~] = rekon_oa_freqdom(rfm, 40, 0.315, 1.48, 0, 1, 1, 5, 8);
        caxis([0, 1e3]);
    
        % Display and save image
        % imshow(img, []);
        img_time(:,:,frame_num) = img;
        frame_num = frame_num+1;
    end
    save(sprintf('0504/ml70delay_all_time_%d.mat', wv),'img_time')
    % save('0504/ml1_time_all.mat','img_time')

end
