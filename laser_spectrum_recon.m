% clear all;
close all;
figure;
for wv = 720:20:900
    %open your raw file here
    folder_path = '0503/';
    file_pattern = sprintf('ml_%d_OBP_Laser_PA_*.raw', wv);

    file_list = dir(fullfile(folder_path, file_pattern));
    file_name = fullfile(folder_path, file_list(1).name); 
    % fid = fopen('Rf_060623_051928_OBP_B_1536.raw', 'rb');
    fid = fopen(file_name, 'rb');
    data = fread(fid, 'int16');
    fclose(fid);
    frame_num = size(data,1)/1024/128;
    body = reshape(data, [1024 128 frame_num]);
    clear data;
    %averaging step. Here we average 20 frames (frame rate 3 Hz), change according to your
    %requirement
    %4 KHz PRF data is already averaged 64 times in DAQ, so 20*64 = 1280, Frame
    %rate = 4000/1280 = 3.125 Hz
    rfm = zeros(1024,128);
    %Change 20 here for different frame rates
    for i=150:frame_num
        rfm = rfm + body(:,:,i);
    end
    %Change this 20 also for different frame rates
    rfm = rfm./20;
    %to avoid noise in the begining and end frames
    rfm(1:200,:) = 0;
    rfm = fliplr (rfm);
    
    %For PA reconstruction
    [rekon,rekonuncut] = rekon_OA_freqdom(rfm(:,:),40,.315,1.48,0,1,1,5,8); 
    caxis([0 1e+03]);
    
    % rekon_US_freqdom1(rfm, 20 , 0.315, 1.48 , 0,1,1,5,1,0);
    %% 
    
    imshow(rekon(1:end,:), []);
    img = rekon;
    % img =(img-min(img(:)))/(max(img(:))-min(img(:)));
    save(sprintf('0503/PA_ml_recon%d.mat',wv),'img');

end
