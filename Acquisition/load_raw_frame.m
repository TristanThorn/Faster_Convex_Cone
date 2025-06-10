function frames = load_raw_frame(filename)
%LOAD_RAW_FRAME Load RF data from a .raw file.
%   frames = LOAD_RAW_FRAME(filename) reads the binary file specified by
%   filename and returns a 3-D array of size [1024, 128, numFrames].
%
%   The helper automatically trims incomplete frames if the file size is
%   not an integer multiple of the frame size.

    fid = fopen(filename, 'rb');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    data = fread(fid, 'int16');
    fclose(fid);

    frame_size = 1024 * 128;
    total_frames = floor(numel(data) / frame_size);
    frames = reshape(data(1:frame_size*total_frames), [1024, 128, total_frames]);
end
