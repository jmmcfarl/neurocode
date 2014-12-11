clear all
close all
addpath(genpath('C:/Code'))
addpath(genpath('C:/WC_Germany/'))
cd C:/WC_Germany/persistent_downs/
% load ./new_pdown_dir.mat
load ./overall_EC_dir

force_recompute = 0;

for dd = 216:length(sess_data)
% for dd = 16:length(new_pdown_dir)
% for dd = 36:length(new_pdown_dir)
%     cd(new_pdown_dir{dd})
cd(sess_data(dd).directory)
    pwd
    
    if ~exist('./all_combined_mp_uds.mat','file') || force_recompute
%     all_EC_compute_state_seqs(new_pdown_name{dd})
    all_EC_compute_state_seqs(sess_data(dd).name)
    
    is_bad(dd) = all_EC_compute_ctx_period();
    
    end
end