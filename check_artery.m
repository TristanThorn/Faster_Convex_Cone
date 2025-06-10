clear all;
close all;

% Parameters setup
wavelengths = 720:20:900;  % Wavelength range
fps = 20;                  % Sampling rate (Hz)
total_fft = [];            % Initialize FFT accumulator

% Process each wavelength separately
for wv = wavelengths
    wv
    % Load data
    load(sprintf('human/ml1_all_time_%d.mat', wv));
    
    % Compute FFT for current wavelength
    current_fft = fft(double(img_time), [], 3);
    
    % Accumulate FFT results
    if isempty(total_fft)
        total_fft = current_fft;
    else
        total_fft = total_fft + current_fft;
    end
end

% Get dimensions from last loaded file
[x_len, y_len, t_len] = size(img_time);


%% Compute average image
[x_len, y_len, t_len] = size(img_time);

% Compute FFT along time dimension
img_fft = total_fft;
% Frequency axis
f = (0:t_len-1) * (fps / t_len);  % Frequency in Hz
half_len = floor(t_len / 2) + 1;  % Index for non-negative frequencies

% Find the index where frequency >= 0.5 Hz
start_idx = find(f >= 0.5, 1); 
f_half = f(start_idx:half_len);  % Frequency starts from 0.5 Hz
img_fft_half = img_fft(:,:,start_idx:half_len);  % Corresponding FFT data

%% Show map at heartbeat freqency 
ratio_map = zeros(x_len, y_len);

for x = 1:x_len
    x
    for y = 1:y_len
        x = 889;
        y = 607;

        x = 594;
        y = 612;


        x_range = max(1,x-3):min(x_len,x+3);
        y_range = max(1,y-10):min(y_len,y+10);
    
        % Get averaged spectrum in the region
        region_fft = img_fft_half(x_range, y_range, :);
        spectrum = mean(mean(abs(region_fft), 1), 2);  % Spatial average
        spectrum = squeeze(squeeze(spectrum));
        spectrum = (spectrum-min(spectrum(:)))/(max(spectrum(:))-min(spectrum(:)));
        % plot (spectrum);
        band_max = max(spectrum(5:7));
        band_mean = mean(spectrum(15:end));

        heart_beat_map(x,y) = band_max/band_mean;
    end
end

ratio_fig = figure(1);
imagesc(heart_beat_map);
axis image;
colorbar;
title('Intensity Ratio: (1-2Hz Max)/(0.5Hz-End Mean)');
xlabel('X position');
ylabel('Y position');

%% Check fft result
main_fig = figure(2);
xy_mean = mean(img_time, 3);  
xy_mean(1:600,:,:) = 0;  
imagesc(xy_mean);
axis image;
title({'Averaged image (720-900nm)', 'Click any point to view spectrum'});
xlabel('X position');
ylabel('Y position');
colorbar;

spectrum_fig = figure(3);
set(spectrum_fig, 'Position', [100 100 800 400]);

% Create main figure for x-y image
main_fig = figure(1);
xy_mean = mean(img_time, 3);  % Average along time dimension
imagesc(xy_mean);
axis image;  % Maintain aspect ratio
title({'Averaged image (720-900nm)', 'Click any point to view spectrum (Close window to exit)'});
xlabel('X position');
ylabel('Y position');
colorbar;

% Create figure for spectrum display
spectrum_fig = figure(2);
set(spectrum_fig, 'Position', [100 100 800 400]);  % Larger window for spectrum

% Continuous point selection until main figure is closed
while ishandle(main_fig)
    figure(main_fig);  % Bring main figure to focus
    [x_click, y_click, button] = ginput(1);
    
    % Exit if figure is closed during ginput
    if ~ishandle(main_fig)
        break;
    end
    
    % Convert click position to matrix indices
    x_idx = round(y_click);  % Row index (vertical)
    y_idx = round(x_click);  % Column index (horizontal)
    
    % Validate clicked position
    if x_idx < 1 || x_idx > x_len || y_idx < 1 || y_idx > y_len
        fprintf('Clicked position out of range! Try again.\n');
        continue;
    end
    
    % Define averaging region 
    x_range = max(1,x_idx-3):min(x_len,x_idx+3);
    y_range = max(1,y_idx-10):min(y_len,y_idx+10);
    
    % Get averaged spectrum in the region
    region_fft = img_fft_half(x_range, y_range, :);
    point_spectrum = mean(mean(abs(region_fft), 1), 2);  % Spatial average
    point_spectrum = (point_spectrum-min(point_spectrum(:)))/(max(point_spectrum)-min(point_spectrum(:)));
    
    % Update spectrum plot
    figure(spectrum_fig);
    clf;  % Clear previous spectrum
    plot(f_half, squeeze(point_spectrum));
    title({sprintf('Frequency spectrum at position (%d, %d)', x_idx, y_idx), ...
           sprintf('Averaged over %d-%dnm (f ≥ 0.5Hz)', min(wavelengths), max(wavelengths)), ...
           sprintf('Spatial average: ±3px vertical, ±5px horizontal (%dx%d area)', ...
                  length(x_range), length(y_range))});
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
    xlim([0.5 max(f_half)]);  % Set x-axis start at 0.5 Hz
    
    % Mark dominant frequency (excluding DC)
    [max_val, max_idx] = max(point_spectrum);
    hold on;
    plot(f_half(max_idx), max_val, 'ro');
    text(f_half(max_idx), max_val, sprintf(' %.2f Hz', f_half(max_idx)), ...
        'VerticalAlignment', 'bottom');
    hold off;
    
    % Set consistent axis limits for better comparison
    ylims = ylim;
    ylim([0 ylims(2)*1.1]);  % Add 10% headroom
    
    % Display frequency information in command window
    fprintf('Position (%d,%d) - Dominant frequency: %.2f Hz\n', x_idx, y_idx, f_half(max_idx));
    fprintf('Averaging region: X=%d:%d (%dpx), Y=%d:%d (%dpx)\n', ...
            min(x_range), max(x_range), length(x_range), ...
            min(y_range), max(y_range), length(y_range));
end

% Cleanup message
if ~ishandle(main_fig)
    fprintf('Main figure closed. Exiting program.\n');
end
