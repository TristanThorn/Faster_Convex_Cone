% Load raw frames
body = load_raw_frame('0501//ml_900_OBP_Laser_PA_674.raw');
frame_num = size(body, 3);
%averaging step. Here we average 20 frames (frame rate 3 Hz), change according to your
%requirement
%4 KHz PRF data is already averaged 64 times in DAQ, so 20*64 = 1280, Frame
%rate = 4000/1280 = 3.125 Hz
% Average all frames (change the range below for different frame rates)
rfm = mean(body(:,:,1:frame_num), 3);
%to avoid noise in the begining and end frames
rfm(1:200,:) = 0;
rfm = fliplr (rfm);
% For PA reconstruction
rekon = reconstruct_pa(rfm);
caxis([0 1e+03]);

% rekon_US_freqdom1(rfm, 20 , 0.315, 1.48 , 0,1,1,5,1,0);
%% 

imshow(rekon, [])
