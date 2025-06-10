clear all;
close all;

%% ==== Parameters ====
us_file = 'human/ml_70delay_OBP_Laser_B_1344.raw';    % US data file
pa_file = 'human/ml_70delay_OBP_Laser_PA_1344.raw';   % PA data file
frame_start = 150;  
frame_end = 250;  % 200 frames
pitch = 0.315;    % cm
c = 1.48;         % mm/us or cm/us depending on function

%% ==== Load US Data ====
body_us = load_raw_frame(us_file);
frame_num_us = size(body_us,3);

%% ==== Load PA Data ====
body_pa = load_raw_frame(pa_file);
frame_num_pa = size(body_pa,3);

%% ==== Check Frame Consistency ====
assert(frame_end <= min(frame_num_us, frame_num_pa), 'Frame range exceeds available data');

%% ==== Reconstruct All Frames ====
rgb_frames = {};  % store overlaid RGB images

for i = frame_start:frame_end
    disp(['Reconstructing frame ', num2str(i)]);

    % Extract US and PA RF frames
    rfm_us = double(body_us(:,:,i)) ./ 20;
    rfm_pa = double(body_pa(:,:,i)) ./ 20;

    % Clean
    rfm_us(1:200,:) = 0;
    rfm_pa(1:200,:) = 0;

    rfm_us = fliplr(rfm_us);
    rfm_pa = fliplr(rfm_pa);

    % Reconstruct images
    us_img = reconstruct_us(rfm_us, struct('F',20,'pitch',pitch,'c',c));
    [pa_img] = reconstruct_pa(rfm_pa, struct('F',40,'pitch',pitch,'c',c));

    % Normalize for overlay
    us_norm = mat2gray(us_img);
    pa_norm = mat2gray(pa_img);
    imshow(us_norm,[])
    % imshow(pa_norm,[])

    % Create overlay image: US as grayscale, PA as red channel

    rgb(:,:,1) = pa_norm(44:end,:);       % R: PA
    rgb(:,:,2) = us_norm(1:end-43,:);       % G: US
    rgb(:,:,3) = us_norm(1:end-43,:);       % B: US
    
    imshow(rgb,[]);
    rgb_frames{end+1} = rgb;
end

%% ==== Display All Frames ====
figure;
for j = 1:length(rgb_frames)
    rgb_show = rgb_frames{j};
    rgb_show(:,:,1) = rgb_show(:,:,1)*1.5;
    imshow(rgb_show);
    title(['Frame ', num2str(frame_start + j - 1)]);
    pause(0.1);  % Adjust playback speed as needed
end



