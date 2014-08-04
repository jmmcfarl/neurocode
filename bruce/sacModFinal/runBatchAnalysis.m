%% LEM EYE TRACKING 
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
global Expt_name bar_ori use_LOOXV
bar_ori = nan; 
use_LOOXV = 1;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    full_eyetracking_withblinks_v2;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori use_LOOXV
end

%% JBE EYE TRACKING ORI=0
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
% Expt_list = {'G086'};
global Expt_name bar_ori use_LOOXV
bar_ori = 0;
use_LOOXV = 1;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    full_eyetracking_withblinks_v2;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori use_LOOXV
end
%% JBE EYE TRACKING ORI=90
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
global Expt_name bar_ori use_LOOXV
bar_ori = 90;
use_LOOXV = 1;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    full_eyetracking_withblinks_v2;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori use_LOOXV
end

%% LEM TRIG AVGS
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
global Expt_name bar_ori
bar_ori = nan; 
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    trig_avg_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

%% JBE TRIG AVGS ORI=0
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
global Expt_name bar_ori
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    trig_avg_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end
%% JBE TRIG AVGS ORI=90
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
global Expt_name bar_ori
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    trig_avg_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

%% LEM FIT COR MODELS AVGS
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
global Expt_name bar_ori
bar_ori = nan; 
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    fitCorrectedModels;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

%% JBE FIT COR MODELS ORI=0
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
global Expt_name bar_ori
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    fitCorrectedModels;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end
%% JBE FIT COR MODELS ORI=90
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
global Expt_name bar_ori
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    fitCorrectedModels;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end
