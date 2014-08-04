clear all
close all
addpath('~/James_scripts/bruce/processing/');

Expt_name = 'M297';
% dir_name = ['~/Data/bruce/' Expt_name '/stims/'];
% dir_name = ['/media/NTlab_data2/Data/bruce/' Expt_name '/stims/'];
dir_name = ['/media/NTlab_data3/Data/bruce/' Expt_name '/stims/'];
cd(dir_name);

%% MAKE SURE THE RC FILES ARE ALL NAMED CONSISTENTLY WITH CHRONOLOGICAL ORDERING BEFORE RUNNING
file_names = dir([dir_name '*.rc*']);
n_rc_files = length(file_names);
rc_num = nan(n_rc_files,1);
not_rc = [];
for ii = 1:n_rc_files
    rc_num(ii) = str2num(file_names(ii).name(7:end));
end
[~,file_order] = sort(rc_num);
file_names = file_names(file_order);
%%
trial_trial_ids = [];
trial_file_ids = [];
trial_seed_nums = [];
trial_is_ds = [];
frame_file_ids = [];
frame_frame_nums = [];
frame_trial_nums = [];
all_left_pix = cell(n_rc_files,1);
all_right_pix = cell(n_rc_files,1);
all_file_names = cell(n_rc_files,1);
rc_has_ds = nan(n_rc_files,1);
rc_has_mtrS = nan(n_rc_files,1);
rc_n_trials = nan(n_rc_files,1);
rc_is_binoc = nan(n_rc_files,1);
for i = 1:n_rc_files
% for i = 42
    fprintf('Parsing file %d of %d. Name: %s\n',i,n_rc_files,file_names(i).name);
    
    [trial_ids,seed_ids,cur_frame_trial_nums,cur_frame_nums,left_pix,right_pix,trial_ds,rc_has_mtrS(i)] = ...
        parse_random_bar_file(file_names(i).name);
    
    %full blank frames are stored as NaNs which get parsed into -128's,
    %force them to all 0's.
    if ~isempty(left_pix)
    Lblank_rows = find(left_pix(:,1) == -128);
    Rblank_rows = find(right_pix(:,1) == -128);
    left_pix(Lblank_rows,:) = 0;
    right_pix(Rblank_rows,:) = 0;
    end

    trial_trial_ids = [trial_trial_ids; trial_ids];
    trial_seed_nums = [trial_seed_nums; seed_ids];
    trial_file_ids = [trial_file_ids; ones(length(trial_ids),1)*i];
    trial_is_ds = [trial_is_ds; trial_ds];

    frame_file_ids = [frame_file_ids; ones(length(cur_frame_nums),1)*i];
    frame_frame_nums = [frame_frame_nums; cur_frame_nums];
    frame_trial_nums = [frame_trial_nums; trial_ids(cur_frame_trial_nums)];

    rc_is_binoc(i) = any(left_pix(:) ~= right_pix(:));
    if any(trial_ds)
        rc_has_ds(i) = 1;
    elseif all(~trial_ds)
        rc_has_ds(i) = 0;
    end
    rc_n_trials(i) = length(trial_ids);
    
    all_left_pix{i} = left_pix;
    all_right_pix{i} = right_pix;
    all_file_names{i} = file_names(i).name;
    
%     pause
end

%%
save stim_data all_* trial* frame* rc_*