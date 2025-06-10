clear all;
close all;

for wv = 720
    % File reading section
    folder_path = '0507/';
    % file_pattern = sprintf('ml2_%d_OBP_Laser_PA_*.raw', wv);
    file_pattern = 'deep_phantom_150_OBP_Laser_PA_*.raw';
    %file_pattern = sprintf('check_saturation_%d_OBP_Laser_PA_*.raw', wv);
    file_list = dir(fullfile(folder_path, file_pattern));
    file_name = fullfile(folder_path, file_list(1).name); 
    fid = fopen(file_name, 'rb');
    data = fread(fid, 'int16');
    fclose(fid);
    
    % Data reshaping
    frame_num = size(data,1)/1024/128;
    body = reshape(data, [1024 128 frame_num]);
    clear data;
    
    % Initialize variables
    rfm = zeros(1024,128);
    max_p2p = -inf;  % Initialize max peak to peak
    max_i = 0;       % Frame index with max p2p
    max_j = 0;       % Channel index with max p2p
    
    % Create figure for all channels
    figure('Position', [100, 100, 1400, 900]);
    sgtitle(sprintf('All 128 Channels (Wavelength: %dnm)', wv));
    
    % Calculate subplot layout (16x8 grid for 128 channels)
    rows = 16;
    cols = 8;
    
    % Main processing loop
    for i = 150
        fprintf('Processing frame %d/%d\n', i, 500);
        
        % Temporary storage for current frame
        frame_signals = squeeze(body(:,:,i));
        
        % Plot all 128 channels in subplots
        for j = 32
            plot(frame_signals(:,j));
            subplot(rows, cols, j);
            plot(frame_signals(:,j), 'b', 'LineWidth', 0.5);
            title(sprintf('Ch %d', j), 'FontSize', 6);
            set(gca, 'XTick', [], 'YTick', []);
            box on;
            
            % Calculate peak-to-peak
            p_to_p = max(frame_signals(:,j)) - min(frame_signals(:,j));
            
            % Update max p2p tracking
            if p_to_p > max_p2p
                max_p2p = p_to_p;
                max_i = i;
                max_j = j;
            end
        end
        
        % Add overall x/y labels
        han = axes(figure(1), 'visible', 'off'); 
        han.XLabel.Visible = 'on';
        han.YLabel.Visible = 'on';
        xlabel(han, 'Sample Index');
        ylabel(han, 'Amplitude');
        
        % Pause to allow visualization (optional)
        drawnow;
    end
    
    % Highlight the channel with max p2p
    subplot(rows, cols, max_j);
    set(gca, 'XColor', 'r', 'YColor', 'r', 'LineWidth', 1.5);
    title(sprintf('Ch %d (MAX)', max_j), 'Color', 'r', 'FontSize', 6);
    
    % Display max p2p information
    fprintf('Max peak-to-peak value: %.4f at frame %d, channel %d\n', max_p2p, max_i, max_j);
    annotation('textbox', [0.4, 0.95, 0.2, 0.04], 'String', ...
        sprintf('Max p2p: %.2f (Frame %d, Ch %d)', max_p2p, max_i, max_j), ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'red');
end

%% output peak to peak
clear all;
close all;

for wv = 720:20:900
    % File reading section
    folder_path = '0507/';
    % file_pattern = sprintf('check_saturation_%d_OBP_Laser_PA_*.raw', wv);
    % file_pattern = sprintf('ml2_%d_OBP_Laser_PA_*.raw', wv);
    file_pattern = 'deep_phantom_150_OBP_Laser_PA_*.raw';
    file_list = dir(fullfile(folder_path, file_pattern));
    file_name = fullfile(folder_path, file_list(1).name)
    fid = fopen(file_name, 'rb');
    data = fread(fid, 'int16');
    fclose(fid);

    % Reshape data
    frame_num = size(data,1)/1024/128;
    body = reshape(data, [1024, 128, frame_num]);
    clear data;

    % Init max tracking
    max_p2p = -inf;
    max_i = 0;
    max_j = 0;

    % Optional plotting setup
    % figure('Position', [100, 100, 1400, 900]);
    % sgtitle(sprintf('All 128 Channels (Wavelength: %dnm)', wv));
    % rows = 16; cols = 8;

    for i = 150:650
        frame_signals = squeeze(body(200:end,:,i));

        for j = 1:128
            p_to_p = max(frame_signals(:,j)) - min(frame_signals(:,j));

            if p_to_p > max_p2p
                max_p2p = p_to_p;
                max_i = i;
                max_j = j;
            end

            % Optional plot (commented for speed)
            % subplot(rows, cols, j);
            % plot(frame_signals(:,j), 'b', 'LineWidth', 0.5);
            % title(sprintf('Ch %d', j), 'FontSize', 6);
            % set(gca, 'XTick', [], 'YTick', []); box on;
        end

        % fprintf('Processed frame %d\n', i);
    end

    fprintf('\n=== Wavelength %dnm ===\n', wv);
    fprintf('Max peak-to-peak value: %.4f\n', max_p2p);
    fprintf('Found at frame %d, channel %d\n', max_i, max_j);
end
