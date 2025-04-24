% SPFI evaluation for single EV/LNP capture
clear all; clc;

params = struct();

params.fpath =  ['M:\Damien Bazin\iSCAT\20250325_ConfocaliSCAT_aperture_calibration\apt1\'];
params.num_EVinc_tps = 200;
params.num_fovs = 1;
params.valid_fovs = 1:params.num_fovs;
params.spot_detachment_tp = 200;

% Channel setup
params.OC_names = {'iSCAT440nm'};
params.channel_wavelengths = [440]*1e-9;
params.chans_to_correct = [];

% iSCAT setup
params.iscat_chans = [1];
params.iscat_thresh = 0.11;
params.RVT_thresh = 0;

% Labels and DEI - left empty
params.lbl_chans = [];
params.lbl_names = {};
params.lbl_tps = [];
params.lbl_threshs = [];

params.dei_chans = [];
params.dei_chan(1).targets = {};
params.dei_chan(1).threshs = 3 + zeros(1, length(params.dei_chan(1).targets));
params.dei_chan(2).targets = {};
params.dei_chan(2).threshs = 2 + zeros(1, length(params.dei_chan(2).targets));
params.targets_to_evaluate = 'all';

% Initialize
params = initialize_low_level_params_v3(params);
