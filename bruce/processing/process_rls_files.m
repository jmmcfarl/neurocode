function [] = process_rls_files(dat_dir,rls_list)

if nargin < 2
    rls_list = nan;
end

cd(dat_dir);

%% MAKE SURE THE RC FILES ARE ALL NAMED CONSISTENTLY WITH CHRONOLOGICAL ORDERING BEFORE RUNNING
file_names = dir([dat_dir '/*.rc*']);
n_rc_files = length(file_names);
rc_num = nan(n_rc_files,1);
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
    rls_num = sscanf(file_names(i).name,'rls.rc%d');
    if any(isnan(rls_list)) || ismember(rls_num,rls_list)
        fprintf('Parsing file %d of %d. Name: %s\n',i,n_rc_files,file_names(i).name);
        
        [trial_ids,seed_ids,cur_frame_trial_nums,cur_frame_nums,left_pix,right_pix,trial_ds,rc_has_mtrS(i)] = ...
            parse_random_bar_file(file_names(i).name);
        
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
    end
end

%%
sname = [dat_dir '/stim_data.mat'];
save(sname,'all_*','trial*','frame*','rc_*');