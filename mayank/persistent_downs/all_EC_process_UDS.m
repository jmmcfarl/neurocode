clear all
close all
addpath(genpath('C:/Code'))
addpath(genpath('C:/WC_Germany/'))

cd C:/WC_Germany/persistent_downs/
% load ./new_pdown_dir.mat
% load ./overall_EC_dir

cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat

force_recompute = 0;

for dd = 1:length(combined_dir)
% for dd = 1:length(sess_data)
% for dd = 1:length(new_pdown_dir)
    % for dd = 36:length(new_pdown_dir)
%     cd(new_pdown_dir{dd})
%         cd(sess_data(dd).directory)
        cd(combined_dir{dd})
    pwd
    
    %     if ~exist('./all_combined_mp_uds.mat','file') || force_recompute
    if ~exist('./pa_hsmm_state_seq_combined_fin_nd.mat','file') || force_recompute
%         all_EC_compute_state_seqs(new_pdown_name{dd})
%                 all_EC_compute_state_seqs(sess_data(dd).name)
temp = find(combined_dir{dd} == '\',1,'last');
cur_name = combined_dir{dd}(temp+1:end);
                all_EC_compute_state_seqs(cur_name)
    end
    is_bad(dd) = all_EC_compute_ctx_period();
end