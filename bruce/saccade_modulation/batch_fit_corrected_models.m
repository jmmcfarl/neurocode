clear all
addpath('~/James_scripts/bruce/saccade_modulation/');

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
% Expt_list = {'M281','M287','M289','M294'};
global Expt_name bar_ori use_MUA fit_unCor fit_rect_quad

use_MUA = true;
fit_unCor = false; %fit models to uncorrected stim?
fit_rect_quad = true; %fit quad models with split linear filter?

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    fit_corrected_simmodels_v2;
    clearvars -except ee Expt_list Expt_name bar_ori use_MUA fit_unCor fit_rect_quad
end

%%
clear all
addpath('~/James_scripts/bruce/saccade_modulation/');

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
global Expt_name bar_ori use_MUA fit_unCor fit_rect_quad

bar_ori = 0; %bar orientation to use (only for UA recs)
use_MUA = true;
fit_unCor = false; %fit models to uncorrected stim?
fit_rect_quad = true; %fit quad models with split linear filter?

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    fit_corrected_simmodels_v2;
    clearvars -except ee Expt_list Expt_name bar_ori use_MUA fit_unCor fit_rect_quad
end

%%
clear all
addpath('~/James_scripts/bruce/saccade_modulation/');

Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
global Expt_name bar_ori use_MUA fit_unCor fit_rect_quad

bar_ori = 90; %bar orientation to use (only for UA recs)
use_MUA = false;
fit_unCor = true; %fit models to uncorrected stim?
fit_rect_quad = true; %fit quad models with split linear filter?

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    fit_corrected_simmodels_v2;
    clearvars -except ee Expt_list Expt_name bar_ori use_MUA fit_unCor fit_rect_quad
end
