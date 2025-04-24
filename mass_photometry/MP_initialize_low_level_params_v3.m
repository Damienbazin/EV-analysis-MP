function params = initialize_low_level_params_v3(params, varargin)

    saveFolderName = 'eval_plots';
    get_metadata = true;

    % parse optional arguments
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'save_folder'
                saveFolderName = varargin{i+1};
            case 'get_metadata'
                get_metadata = varargin{i+1};
        end
    end

    mkdir([params.fpath, saveFolderName]);

    %% Base setup for mass photometry
    params.num_iscat_chans = length(params.iscat_chans);
    params.num_chans = params.num_iscat_chans;
    params.twocolor_iscat = (length(params.iscat_chans) > 1);
    params.localization_chan = 1;

    % Image cleanup and spot detection
    params.min_pixel_crud = 80;
    params.border = 30;
    params.min_spot_distance = 8;

    % Disable fluorescence and DEI completely
    params.num_lbl_chans = 0;
    params.num_lbl_targets = 0;
    params.num_dei_chans = 0;
    params.include_DEI = false;
    params.read_fluorescence = false;

    % File for h5 access (you must generate this before)
    params.h5(1).fid = [params.fpath, 'MP_data.h5'];
    params.h5(1).num_tps = params.num_EVinc_tps;

    % Mass photometry-specific flags
    params.analysis_mode = 'mass_photometry';
    params.num_frames = params.num_EVinc_tps;
    params.frame_time = 10 / params.num_frames; % 10 s total acquisition

    % Try to read metadata (image size, pixel size)
    if get_metadata
        try
            metadata = h5r_metadata(params.h5(1).fid);
            params.dim_x = metadata.dim_x;
            params.dim_y = metadata.dim_y;
            params.pixel_x = metadata.pixel_x;
            params.pixel_y = metadata.pixel_y;
        catch
            try
                files = dir(fullfile([params.fpath, '1'], '*.nd2'));
                [~, idx] = min([files.datenum]); % get the earliest file

                ch = read_nd_fluorescence(fullfile([params.fpath, '1'], files(idx).name));
                frame = ch(1).fov(1).raw;
                [params.dim_x, params.dim_y] = size(frame);
                clear frame
            end
        end
    end

    % Compatibility fallback
    if ~isfield(params, 'EV_localization_tp')
        params.EV_localization_tp = 2;
    end

    if isfield(params, 'RVT_thresh')
        fprintf(1, "\nWarning: RVT thresholding must be applied via 'RVT_thresh' varargin to the plot function!\n\n");
    end
end
