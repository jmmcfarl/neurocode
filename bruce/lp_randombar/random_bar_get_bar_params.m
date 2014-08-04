clear all
% close all

ExptNum = 239;
cd(['~/Data/bruce/M' num2str(ExptNum)]);
load ./random_bar_eyedata_ftime.mat

load(['lemM' num2str(ExptNum) 'Expts.mat']);
dt = 0.01;

if ExptNum == 235
    bar_expts(bar_expts==51) = []; %this has different set of stim positions
end
if ExptNum == 239
    bar_expts(bar_expts==40) = []; %this has different set of stim positions
end

%%
full_exptvec = [];
full_taxis = [];
full_trialvec = [];
full_bar_pos = [];
for ee = 1:length(bar_expts);
    fprintf('Processing expt %d of %d\n',ee,length(bar_expts));
    
    trial_durs{ee} = [Expts{bar_expts(ee)}.Trials(:).dur];
    n_trials(ee) = length(Expts{bar_expts(ee)}.Trials);
    all_t_axis = [];
    all_used_inds = [];
    all_expt_vec = [];
    all_trial_vec = [];
    all_bar_Op = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
                
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
                                        
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op(:)];
    end
        
    full_exptvec = [full_exptvec; ones(length(all_t_axis),1)*ee];
    full_taxis = [full_taxis; all_t_axis];
    full_trialvec = [full_trialvec; all_trial_vec];
    full_bar_pos = [full_bar_pos; all_bar_Op];    
end

%%
if ExptNum==235
    n_bars = 6;
end
if ExptNum==239
    n_bars = 6;
end
bar_is_on = find(full_bar_pos > -1000);
% bar_range = [min(full_bar_pos(bar_is_on)) max(full_bar_pos(bar_is_on))];
% bar_axis = bar_range(1):0.125:(bar_range(1)+0.125*(n_bars-1));
bar_axis = unique(full_bar_pos(bar_is_on));

bar_params.n_bars = n_bars;
% bar_params.bar_range = bar_range;
bar_params.bar_axis = bar_axis;
%%
save bar_params bar_params