clear all
close all
% addpath(genpath('C:/Code'))
% addpath(genpath('C:/WC_Germany/'))
addpath('~/James_scripts/mayank/persistent_downs/');
addpath('~/James_scripts/mayank/hsmm_state_detection/');
addpath('~/James_scripts/mayank/sven_thomas_combined/');
% load ~/Analysis/Mayank/final_pdown_analysis/compiled_data.mat
load ~/Analysis/Mayank/final_pdown_analysis/compiled_corticalMP_data.mat

force_recompute = 0; %overwrite existing files?

for dd = 1:length(data) %loop over recs in this database
    dd
        cur_dir = data(dd).dir;
        new_rawdir = map_to_new_drive_locs(cur_dir);
        new_tardir = data(dd).new_dir;
        cd(new_tardir)
    
    if ~exist('./pa_hsmm_state_seq_combined_fin_nd.mat','file') || force_recompute
        temp = find(new_rawdir == '/',1,'last');
        cur_name = new_rawdir(temp+1:end);
        all_EC_compute_state_seqs(cur_name)
    end
    is_bad(dd) = all_EC_compute_ctx_period();
end