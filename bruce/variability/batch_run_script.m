clear all

addpath('~/James_scripts/bruce/variability/');
global Expt_name bar_ori monk_name rec_number 

Expt_list = {};
expt_oris = [];
expt_mname = {};
expt_rnum = [];

% Expt_list = cat(2,Expt_list,{'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'});
% expt_oris = cat(1,expt_oris,[80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90]);
% expt_mname = cat(2,expt_mname,repmat({'lem'},1,length(Expt_list)));
% expt_rnum = cat(1,expt_rnum,ones(length(Expt_list),2));
% 
% Expt_list = cat(2,Expt_list,{'G085','G086','G087','G088','G089','G091','G093','G095'});
% expt_oris = cat(1,expt_oris,[0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan]);
% expt_mname = cat(2,expt_mname,repmat({'jbe'},1,8));
% expt_rnum = cat(1,expt_rnum,ones(8,2));
% 
Expt_list = cat(2,Expt_list,{'M005','M309','M009','M010','M011','M012','M013','M014'});
expt_oris = cat(1,expt_oris,[50 nan; 120 nan; 0 nan; 60 nan; 160 160; 0 0; 100 nan; 40 nan]);
expt_mname = cat(2,expt_mname,{'jbe','lem','jbe','jbe','jbe','jbe','jbe','jbe'});
expt_rnum = cat(1,expt_rnum,[1 1; 1 1; 1 1; 1 1; 1 2; 1 2; 1 1';1 1]);
% 
%%
% batch_function = 'fit_corrected_models_compactData';
% batch_function = 'variability_model_compactData';
% batch_function = 'variability_rpt_anal_compact';
% batch_function = 'rpt_stimprops_sim_test';
% batch_function = 'create_processed_data';
% batch_function = 'full_eyetracking_hres_compactData_allSULOO';
% batch_function = 'fix_modelfit_ModData';
batch_function = 'full_eyetracking_compactData';
% batch_function = 'drift_grating_simulations';
% batch_function = 'add_unCorr_modfit';
% batch_function = 'fit_microsac_models_compactData';

for Elist_cnt = 1:length(Expt_list)
% for Elist_cnt = [23 24]
    Expt_name = Expt_list{Elist_cnt};
    monk_name = expt_mname{Elist_cnt}; 
    for bori_cnt = 1:2
       bar_ori = expt_oris(Elist_cnt,bori_cnt);
       rec_number = expt_rnum(Elist_cnt,bori_cnt);
       if ~isnan(bar_ori)
        fprintf('Running script %s on Expt %s ori %d rec %d\n',batch_function,Expt_name,bar_ori,rec_number);
       eval(batch_function);
       clearvars -except Elist_cnt bori_cnt Expt_list Expt_name bar_ori monk_name expt_oris expt_rnum expt_mname rec_number batch_function
       end
    end
end

