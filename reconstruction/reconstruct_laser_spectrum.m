% clear all;
close all;
figure;
for wv = 720:20:900
    %open your raw file here
    folder_path = '0503/';
    file_pattern = sprintf('ml_%d_OBP_Laser_PA_*.raw', wv);

    file_list = dir(fullfile(folder_path, file_pattern));
    file_name = fullfile(folder_path, file_list(1).name); 
    body = load_raw_frame(file_name);
    frame_num = size(body,3);
    %averaging step. Here we average 20 frames (frame rate 3 Hz), change according to your
    %requirement
    %4 KHz PRF data is already averaged 64 times in DAQ, so 20*64 = 1280, Frame
    %rate = 4000/1280 = 3.125 Hz
    % Average frames starting from 150 (adjust as needed)
    rfm = mean(body(:,:,150:frame_num), 3);
    %to avoid noise in the begining and end frames
    rfm(1:200,:) = 0;
    rfm = fliplr (rfm);
    
    % For PA reconstruction
    rekon = reconstruct_pa(rfm);
    caxis([0 1e+03]);
    
    % rekon_us_freqdom1(rfm, 20 , 0.315, 1.48 , 0,1,1,5,1,0);
    %% 
    
    imshow(rekon(1:end,:), []);
    img = rekon;
    % img =(img-min(img(:)))/(max(img(:))-min(img(:)));
    save(sprintf('0503/PA_ml_recon%d.mat',wv),'img');

end
