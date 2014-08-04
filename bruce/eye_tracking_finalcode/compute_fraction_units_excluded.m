clear all


expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
anal_name = 'monoc_eyecorr_hbar2.mat';

clear jbe_*
for ee = 1:length(expt_list)
    
    cur_expt_name = expt_list{ee};
    cur_expt_num = str2num(cur_expt_name(2:end));
    
    anal_dir = ['~/Analysis/bruce/' cur_expt_name '/ET_final/'];
    cd(anal_dir);
    load(anal_name,'it_LLimp','et_tr_set')
    jbe_tot_units(ee) = size(it_LLimp,2);
    jbe_used_units(ee) = length(et_tr_set);
end

%%
clear lem_*

expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
anal_name = 'monoc_eyecorr2_Cprior.mat';
for ee = 1:length(expt_list)
    
    cur_expt_name = expt_list{ee};
    cur_expt_num = str2num(cur_expt_name(2:end));
    
    anal_dir = ['~/Analysis/bruce/' cur_expt_name '/ET_final/'];
    cd(anal_dir);
    load(anal_name,'it_LLimp','et_tr_set')
    lem_tot_units(ee) = size(it_LLimp,2);
    lem_used_units(ee) = length(et_tr_set);
end

%%
lem_fract_used = lem_used_units./lem_tot_units;
jbe_fract_used = jbe_used_units./jbe_tot_units;

all_fract_used = cat(2,jbe_fract_used,lem_fract_used);