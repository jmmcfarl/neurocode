clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/saccade_modulation/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
global Expt_name bar_ori
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    all_sac_trig_avgs_v4;
    clearvars -except Elist_cnt Expt_list Expt_name
end

%%
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/saccade_modulation/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
global Expt_name bar_ori
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    all_sac_trig_avgs_v4;
    clearvars -except Elist_cnt Expt_list Expt_name
end
%%

clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/saccade_modulation/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
global Expt_name bar_ori
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    all_sac_trig_avgs_v4;
    clearvars -except Elist_cnt Expt_list Expt_name
end

%%
% clear all
% addpath('~/James_scripts/bruce/eye_tracking_improvements/');
% addpath('~/James_scripts/bruce/saccade_modulation/');
% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
% global Expt_name 
% for Elist_cnt = 1:length(Expt_list)
%     Expt_name = Expt_list{Elist_cnt};
%     fit_corrected_simmodels;
%     clearvars -except Elist_cnt Expt_list Expt_name
% end
% 
%%
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/saccade_modulation/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
global Expt_name bar_ori use_MUA
bar_ori = 0;
use_MUA = false;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    full_sacmod_stimproc_v2;
    clearvars -except Elist_cnt Expt_list Expt_name
end

clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/saccade_modulation/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
global Expt_name bar_ori use_MUA
bar_ori = 90;
use_MUA = false;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    full_sacmod_stimproc_v2;
    clearvars -except Elist_cnt Expt_list Expt_name
end

clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/saccade_modulation/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
global Expt_name bar_ori use_MUA
bar_ori = 0;
use_MUA = false;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    full_sacmod_stimproc_v2;
    clearvars -except Elist_cnt Expt_list Expt_name
end