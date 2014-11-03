%% LEM EYE TRACKING
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
global Expt_name bar_ori use_LOOXV
bar_ori = nan;
use_LOOXV = 1;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    full_eyetracking_withblinks_v2;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori use_LOOXV
end

Expt_name = 'M297'; bar_ori = 0;
full_eyetracking_withblinks_v2;
clearvars -except Elist_cnt Expt_list Expt_name bar_ori use_LOOXV
Expt_name = 'M297'; bar_ori = 90;
full_eyetracking_withblinks_v2;
clearvars -except Elist_cnt Expt_list Expt_name bar_ori use_LOOXV

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
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
global Expt_name bar_ori
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    trig_avg_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_name = 'M297'; bar_ori = 0;
trig_avg_analysis;
Expt_name = 'M297'; bar_ori = 90;
trig_avg_analysis;

%% JBE TRIG AVGS
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

%% GENERAL AVGS
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
global Expt_name bar_ori
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    gen_tavg_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_name = 'M297'; bar_ori = 0;
gen_tavg_analysis;
Expt_name = 'M297'; bar_ori = 90;
gen_tavg_analysis;

clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
global Expt_name bar_ori
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    gen_tavg_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
global Expt_name bar_ori
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    gen_tavg_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

%% REGRESSION SACMOD ANALYSIS
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
global Expt_name bar_ori
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sac_regression_ratemod;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_name = 'M297'; bar_ori = 0;
sac_regression_ratemod;
Expt_name = 'M297'; bar_ori = 90;
sac_regression_ratemod;

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sac_regression_ratemod;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sac_regression_ratemod;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

%% LEM FIT COR MODELS AVGS
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
global Expt_name bar_ori
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    fitCorrectedModels;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end
Expt_name = 'M297'; bar_ori = 0;
fitCorrectedModels;
Expt_name = 'M297'; bar_ori = 90;
fitCorrectedModels;

%% JBE FIT COR MODELS
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

clear all
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
global Expt_name bar_ori
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    fitCorrectedModels;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end


%% LEM sacStimMod Analysis
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
%         Expt_list = {'M294'};
global Expt_name bar_ori
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sacMod_stimProc;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_name = 'M297'; bar_ori = 0;
sacMod_stimProc;
Expt_name = 'M297'; bar_ori = 90;
sacMod_stimProc;

%% JBE sacStimMOd Analysis
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
global Expt_name bar_ori

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sacMod_stimProc;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sacMod_stimProc;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end


%% LEM sac_info_timing
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
global Expt_name bar_ori
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    info_timing_calcs;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_name = 'M297'; bar_ori = 0;
info_timing_calcs;
Expt_name = 'M297'; bar_ori = 90;
info_timing_calcs;

%% JBE sac_info_timing Analysis
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
global Expt_name bar_ori
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    info_timing_calcs;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    info_timing_calcs;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end



%% LEM LFP Analysis
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
global Expt_name bar_ori
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sac_LP_LFP_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_name = 'M297'; bar_ori = 0;
sac_LP_LFP_analysis;
Expt_name = 'M297'; bar_ori = 90;
sac_LP_LFP_analysis;


%% LEM sacStimTypeDep Analysis
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
global Expt_name bar_ori
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sacMod_typedep;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

%% JBE sacStimTypeDep Analysis
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
global Expt_name bar_ori

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sacMod_typedep;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    sacMod_typedep;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

%% Fix disparity analysis
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
addpath('~/James_scripts/bruce/sacModFinal/');
global Expt_name bar_ori

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
global Expt_name bar_ori
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    eyepos_accuracy_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_name = 'M297'; bar_ori = 0;
eyepos_accuracy_analysis;
Expt_name = 'M297'; bar_ori = 90;
eyepos_accuracy_analysis;

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    eyepos_accuracy_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end

Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    eyepos_accuracy_analysis;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori
end



%%
% clear all
% addpath('~/James_scripts/bruce/eye_tracking_improvements/');
% addpath('~/James_scripts/bruce/sacModFinal/');
% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
% %         Expt_list = {'M294'};
% global Expt_name bar_ori
% bar_ori = nan;
% for Elist_cnt = 1:length(Expt_list)
%     Expt_name = Expt_list{Elist_cnt};
%     sac_add_pre_latency;
%     clearvars -except Elist_cnt Expt_list Expt_name bar_ori
% end
% 
% Expt_name = 'M297'; bar_ori = 0;
% sac_add_pre_latency;
% Expt_name = 'M297'; bar_ori = 90;
% sac_add_pre_latency;
% 
% % JBE sacStimMOd Analysis
% clear all
% global Expt_name bar_ori
% 
% Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
% bar_ori = 0;
% for Elist_cnt = 1:length(Expt_list)
%     Expt_name = Expt_list{Elist_cnt};
%     sac_add_pre_latency;
%     clearvars -except Elist_cnt Expt_list Expt_name bar_ori
% end
% 
% Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
% bar_ori = 90;
% for Elist_cnt = 1:length(Expt_list)
%     Expt_name = Expt_list{Elist_cnt};
%     sac_add_pre_latency;
%     clearvars -except Elist_cnt Expt_list Expt_name bar_ori
% end
