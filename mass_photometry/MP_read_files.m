function params = MP_read_files(params, varargin)
    %read the nd2 files and sort everything into the h5 file(s)
    if ismac
        slash = '/';
    elseif ispc
        slash = '\';
    end

    %debugging options - do not read incubation files at all, fluorescnece during incubation or DEI cycles
    include_incubation = true;
    include_fluorescence = params.read_fluorescence;
    include_dei = params.include_DEI;
    wbar = true;

    %define keywords with which files in mixed folders are identified
    flatfield_kw = 'Flatfield';     %flatfield files need to contain this
    fluo_kw = 'fluorescence';       %fluorescence files need to contain this
    %high-res iSCT files need to contain a keyword to be included and a
    %second one that preceeds the FOV number (as there is 1 iscat FOV per
    %file currently so they need to be sorted)
    iSCAT_kw = 'iSCAT';
    fov_keyword = 'fov';

    %read only specific timepoints
    read_tps = 1:params.num_EVinc_tps;

    %fluorescence frames get a pseudoflatfield becuase some channels have bleedthrough
    %and therefore a non-flat background. Some channels have really bright
    %spots and the pseudoFF then causes rings around them. We therefore
    %have the option to 'clip' the intensity of the bright spots going into
    %the psudoFF
    pseudoFF_clip_std = 50;

    %parse kwargs
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'include_incubation'
                include_incubation = varargin{i+1};
            case 'include_dei'
                include_dei = varargin{i+1};
            case 'include_fluorescence'
                include_fluorescence = varargin{i+1};
            case 'flatfield_kw'
                flatfield_kw = varargin{i+1};
            case 'iSCAT_kw'
                iSCAT_kw = varargin{i+1};
            case 'fov_keyword'
                fov_keyword = varargin{i+1};
            case 'fluo_kw'
                fluo_kw = varargin{i+1};
            case 'waitbar'
                wbar = varargin{i+1};
            case 'pseudoFF_clip_std'
                pseudoFF_clip_std = varargin{i+1};
            case 'read_tps'
                read_tps = varargin{i+1};
        end
    end

    for i = params.evaluate
        if wbar
            h = waitbar(0, ['Reading h5 file ', num2str(i)', '...']);
        end

        %params.evaluate = [1,2] for [DEI label, DEI disp]
        if include_incubation
            %loop over the EV incubation timepoints (normally 2-3
            %but can be larger for custom exps)
            
            for t = read_tps
                %first take care of all iSCAT related tings
                includes_hr_iscat = false(1, params.num_iscat_chans);  %keep track if we have the high-resolution iSCAT
                for c = 1:params.num_iscat_chans
                    if params.num_iscat_chans > 1
                        %add id to the keyword
                        f_kw = [flatfield_kw, num2str(c)];
                        i_kw = [iSCAT_kw, num2str(c)];
                    else
                        f_kw = flatfield_kw;
                        i_kw = iSCAT_kw;
                    end

                    %get all files in the folder (folder is called just '1', '2', etc)
                    fstr = dir([params.fpath, num2str(t), slash]);
                    fnames = {fstr(:).name};

                    %first figure out if there are high-res iSCAT channels
                    vld = sum(contains(fnames, i_kw)) > 0;

                    %if there is a flatfield file, read it (regardless of
                    %whether there are high_res iSCAT)
                    %if there is no flatfield, reuse from previous
                    %timepoint
                    flat_id = find(contains(fnames, f_kw));
                    if not(isempty(flat_id))
                        if size(flat_id, 2) > 1
                            ME = MException('MyComponent:incorrect_input', ['Timepoint ', num2str(t), ' multiple ', f_kw, ' files found']);
                            throw(ME);
                        end                    
                        file = [params.fpath,  num2str(t), slash, char(fnames(flat_id))];
            
                        ch = read_multidim_nd2(file);
                        %different NIS versions save FF frames in diferent
                        %dimensions
                        
                        
                        if length(ch(1).fov) > 1
                            num_flat = length(ch(1).fov);
                            ff_dim = 1;
                        elseif length(ch(1).fov(1).tp) > 1
                            num_flat = length(ch(1).fov(1).tp);
                            ff_dim = 2;
                        end
    
                        %in the first execution: get the image size
                        if not(isfield(params, 'dim_x'))
                            [params.dim_x, params.dim_y] = size(ch(1).fov(1).tp(1).raw);
                        end
    
                        %average all frames from the flatfield
                        flat = zeros(params.dim_x, params.dim_y, num_flat);
                        if ff_dim == 1
                            for k = 1:num_flat    
                                flat(:,:,k) = ch(1).fov(k).tp(1).raw;
                            end
                        elseif ff_dim == 2
                            for k = 1:num_flat    
                                flat(:,:,k) = ch(1).fov(1).tp(k).raw;
                            end
                        end
                        flatfield = median(flat, 3);
                        clear flat num_flat ch;
                    end

                    %there are high-res iSCAT files -> process them
                    if vld > 0
                        if size(flat_id, 2) < 1
                            %high_res iSCAT require good flatfield files - don't take them from other timepoints
                            ME = MException('MyComponent:incorrect_input', ['Timepoint ', num2str(t), ' - no ', f_kw, ' file found']);
                            throw(ME);           
                        end

                        %find all high-res iSCAT files in the folder and sort FOVs by creation date
                        fnames_iscat = sort_iscat_files(fnames, params.num_fovs, 'iSCAT_keyword', i_kw, 'fov_keyword', fov_keyword);
                        num_iscat = size(fnames_iscat, 2);
                        

                        %check for overexposure
                        overexposed_flatfield = find(flatfield == 2^16-1); 
          
                        % High-res iSCAT
                        % read in the files now

                        includes_hr_iscat(c) = true;
                        %go though FOVs and read the iscat data
                        for f = 1:params.num_fovs        
                            file = [params.fpath, num2str(t), slash, char(fnames_iscat(f))];
                            % f
                            [ch_raw, ch_name] = read_multidim_nd2(file);
                            if isempty(ch_raw)
                                ME = MException('MyComponent:incorrect_input', ['high-res iSCAT t=', num2str(t), ', c=', num2str(c), ', f=', num2str(f)]);
                                throw(ME);  
                            end
                
                           % Quick & dirt hardcode fix /Fabian
                             if strcmp(ch_name{1}, 'Mono')
                                if strcmp(f_kw, 'Flatfield1')
                                    ch_name = 'iSCAT_440';  
                                elseif strcmp(f_kw, 'Flatfield2')
                                    ch_name = 'iSCAT_740';
                                end
                             end


                            frame = double(ch_raw(1).fov(1).tp(1).raw);

                            %check for overexposure
                            overexposed = find(frame == 2^16-1);
                            pixel_inds = unique([overexposed_flatfield; overexposed]);
                            h5w_overexp(params.h5(i).fid, t, c, f, pixel_inds);

                            %process the frame
                            frame = frame ./ flatfield;                     %flatfield it
                            frame = frame ./ imgaussfilt(frame, 10);        %pseudo flatfield
                            frame = frame / median(frame, "all");           %normalize
    
                            %sort into correct channel
                            id = find(strcmp(ch_name, params.OC_names));
                            
                            if isempty(id)
                                fprintf(1, ['Unspecified channel ', char(ch_name),' encountered in: ', char(file), '\n']);
                            end
                            h5w(params.h5(i).fid, t, id, f, frame);    
                        end
                    else
                        %we don't have high-res iscat files
                        %instead we should have low-res iscat in the
                        %fluorescence file, we'll read those in below
                        includes_hr_iscat(c) = false;

                    end
                end
                %%
                %Fluorescence (incubation timepoints)
                if include_fluorescence
                    %find the right file
                    fstr = dir([params.fpath,  num2str(t), slash]);
                    fnames = {fstr(:).name};
                    id = find(contains(fnames, fluo_kw, 'IgnoreCase', true));

                    if size(id, 2) > 1
                        ME = MException('MyComponent:incorrect_input', ['Timepoint ', num2str(t), ' multiple Fluorescence files found']);
                        throw(ME);
                    end


                    if size(id, 2) < 1
                        %We'll allow skipping fluorescnece files in the incubation
                        fprintf(1,  ['Warning: timepoint ', num2str(t), ' no Fluorescence file found']);
                    else
                                           
                        %read it
                        file = [params.fpath,  num2str(t), slash, char(fnames(id))];
                        [ch_raw, oc_names] = read_multidim_nd2(file);
                        if isempty(ch_raw)
                                ME = MException('MyComponent:incorrect_input', ['Fluorescence file t=', num2str(t), ', c=', num2str(c), ', f=', num2str(f)]);
                                throw(ME);  
                        end

                        num_chans_in_file = size(oc_names, 2);
                                    
                        %process and sort it all
                        for c = 1:num_chans_in_file
                            %figure out what type of channel we are dealing with
                            is_iscat = contains(char(oc_names(c)), 'iscat', 'IgnoreCase', true);

                            %figure out which iSCAT channel it is (iSCAT
                            %channels are usually takes last if there are
                            %fluorescence labels present)
                            for s = 1:params.num_iscat_chans
                                iscat_OC_name = params.OC_names{params.iscat_chans(s)};
                                if contains(oc_names{c}, iscat_OC_name, 'IgnoreCase', true)
                                    chan_id = params.iscat_chans(s);
                                end
                            end

                            
                            %loop through FOVs
                            for f = 1:params.num_fovs
                                if is_iscat == false
                                    %this is a regular fluorescence channel
                                    frame = double(ch_raw(c).fov(f).tp(1).raw);     
            
                                    %Subtract background inhomogeneities from the frame. Sort of like the pseudo flatfield. 
                                    % However, the fluo channels have much greater dynamic range so we cap the
                                    % intensity values before creating the pseudo flatfield
                                    med_frame = median(frame(:));
                                    clipped_frame = frame;
                                    clipped_frame(frame > med_frame+pseudoFF_clip_std) = med_frame+pseudoFF_clip_std;
                                    clipped_frame(frame < med_frame-pseudoFF_clip_std) = med_frame-pseudoFF_clip_std;
                                    frame = frame - imgaussfilt(clipped_frame, 10);                     %pseudo flatfield
                                    frame = frame - median(frame, "all");               %normalize
    
                                    %sort into correct channel
                                    id = find(strcmp(oc_names{c}, params.OC_names));
                                    if isempty(id)
                                        fprintf(1, ['Unspecified channel encountered: ', oc_names{c}, '\n'])
                                    end
                                    h5w(params.h5(i).fid, t, id, f, frame);
    
                                elseif is_iscat == true
                                    %this is an iSCAT channel
                                    %if we already have high-res
                                    %iSCAT files, ignore this channel (legacy,
                                    %we used to record both hr and lr iSCAT).
                                    %Else read it in                            
    
                                    if includes_hr_iscat(chan_id) == false
                                        %process the frame. Use the last
                                        frame = double(ch_raw(c).fov(f).tp(1).raw);
                                        frame = frame ./ flatfield;                     %flatfield it
                                        frame = frame ./ imgaussfilt(frame, 10);        %pseudo flatfield
                                        frame = frame / median(frame, "all");           %normalize

                                        %sort into correct channel
                                        id = find(strcmp(oc_names(c), params.OC_names));
                                        h5w(params.h5(i).fid, t, chan_id, f, frame);    
                                    end
                                end
    
                                %seems like the FOV reading time is neglectable
                                %compared to flatfield processing. But that's
                                %what we got here
                                if c == num_chans_in_file
                                    progress = t*f/params.num_fovs/params.h5(i).num_tps;
                                    if wbar
                                        waitbar(progress, h, sprintf(['Reading h5 file ', num2str(i)', '... %d%%'], int32(100*progress)));
                                        drawnow
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    
        %%
        % DEI
        if include_dei
            %look if there is a flatfield file in the DEI folder
            fstr = dir([params.fpath, 'DEI', slash]);
            fnames = {fstr(:).name};
            flat_id = find(contains(fnames, flatfield_kw));
            has_flatfield = false;
            if not(isempty(flat_id))
                has_flatfield = true;
                if size(flat_id, 2) > 1
                    ME = MException('MyComponent:incorrect_input', ['Timepoint ', num2str(t), ' multiple ', f_kw, ' files found']);
                    throw(ME);
                end                    
                file = [params.fpath,  'DEI', slash, char(fnames(flat_id))];
    
                ch = read_multidim_nd2(file);
                %different NIS versions save FF frames in diferent
                %dimensions
                
                
                if length(ch(1).fov) > 1
                    num_flat = length(ch(1).fov);
                    ff_dim = 1;
                elseif length(ch(1).fov(1).tp) > 1
                    num_flat = length(ch(1).fov(1).tp);
                    ff_dim = 2;
                end

                %in the first execution: get the image size
                if not(isfield(params, 'dim_x'))
                    [params.dim_x, params.dim_y] = size(ch(1).fov(1).tp(1).raw);
                end

                %average all frames from the flatfield
                flat = zeros(params.dim_x, params.dim_y, num_flat);

                if ff_dim == 1
                    for k = 1:num_flat    
                        flat(:,:,k) = ch(1).fov(k).tp(1).raw;
                    end
                elseif ff_dim == 2
                    for k = 1:num_flat    
                        flat(:,:,k) = ch(1).fov(1).tp(k).raw;
                    end
                end
                flatfield = median(flat, 3);
                clear flat num_flat ch;
            end


            if i == 1
                fnames = params.labeling_files;
            elseif i == 2
                fnames = params.displacement_files;
            end

            for t = 1:size(fnames, 2)
                %loop through all the files (=DEI cycles)
                file = [params.fpath, char(fnames(t))];
                [ch_raw, oc_names] = read_multidim_nd2(file); 
                num_chans_in_file = size(oc_names, 2);
    
                for c = 1:num_chans_in_file
                    %figure out what type of channel we are dealing with
                    is_iscat = contains(char(oc_names(c)), 'iscat', 'IgnoreCase', true);
                    
                    %loop through FOVs
                    for f = 1:params.num_fovs
                        if is_iscat == false
                            %this is a regular fluorescence channel

                            frame = double(ch_raw(c).fov(f).tp(1).raw);     
    
                            %Subtract background inhomogeneities from the frame. Sort of like the pseudo flatfield. 
                            % However, the fluo channels have much greater dynamic range so we cap the
                            % intensity values before creating the pseudo flatfield
                            med_frame = median(frame(:));
                            std_frame = 50;
                            clipped_frame = frame;
                            clipped_frame(frame > med_frame+std_frame) = med_frame+std_frame;
                            clipped_frame(frame < med_frame-std_frame) = med_frame-std_frame;
                            frame = frame - imgaussfilt(clipped_frame, 10);                     %pseudo flatfield
                            frame = frame - median(frame, "all");               %normalize

                            %sort into correct channel
                            id = find(strcmp(oc_names(c), params.OC_names));
                            if isempty(id)
                                fprintf(1, ['Unspecified channel encountered: ', oc_names{c}, '\n'])
                            end
                            h5w(params.h5(i).fid, t+params.num_EVinc_tps, id, f, frame);

                        elseif is_iscat == true
                            frame = double(ch_raw(c).fov(f).tp(1).raw);

                            if has_flatfield
                                %process the frame
                                frame = frame ./ flatfield;                     %flatfield it
                                frame = frame ./ imgaussfilt(frame, 10);        %pseudo flatfield
                                frame = frame / median(frame, "all");           %normalize

                            else
    
                                % RoG filter instead of flatfield
                                frame = imgaussfilt(frame, 0.5) ./ imgaussfilt(frame, 3);
    
                                %frame = frame ./ imgaussfilt(frame, 10);        %pseudo flatfield
                                frame = frame / median(frame, "all");           %normalize
                            end
    
                            %figure out which chan this should go into
                            %the flatfield OC names are usually the same as
                            %for the non-averaged iSCAT, they just contain
                            %an additional description Flatfield
                            id = [];
                            for d = 1:params.num_chans
                                if contains(oc_names(c), params.OC_names(d))
                                    id = d;
                                end
                            end
                            if isempty(id)
                                fprintf(1, ['Unspecified channel encountered: ', oc_names{c}, '\n'])
                            end
                            h5w(params.h5(i).fid, t+params.num_EVinc_tps, id, f, frame);    
                        end
                        if c == num_chans_in_file
                            progress = params.num_EVinc_tps/params.h5(i).num_tps + sub2ind([params.num_fovs, size(fnames, 2)], f, t)/params.num_fovs/params.h5(i).num_tps;
                            if wbar 
                                waitbar(progress, h, sprintf(['Reading h5 file ', num2str(i)', '... %d%%'], int32(100*progress)));
                                drawnow
                            end

                        end
                    end
                end
            end
        end
        if wbar
            close(h);
        end

        clear h
    end
end
