%% HIGHRES EYE TRACKING
clear all
addpath('~/James_scripts/bruce/lfp_models/');
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
global Expt_name bar_ori

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};%NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
n_probes = 24;
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
% Expt_list = {'M296'};%NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
% ori_list = [45 nan];

for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    for ii = 1:2
        bar_ori = ori_list(Elist_cnt,ii);
        if ~isnan(bar_ori)
            k99_lfp_models_for_batch;
            clearvars -except Elist_cnt Expt_list Expt_name ori_list bar_ori ii
        end
    end
end
