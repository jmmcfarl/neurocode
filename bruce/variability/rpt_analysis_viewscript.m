clear all

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
expt_mname = repmat({'lem'},1,length(Expt_list));

Expt_list = cat(2,Expt_list,{'M005','M309','M009','M010','M011','M012'});
expt_oris = cat(1,expt_oris,[50 nan; 120 nan; 0 nan; 60 nan; 160 nan; 0 nan]);
expt_mname = cat(2,expt_mname,{'jbe','lem','jbe','jbe','jbe','jbe'});

%%
close all
base_sname = 'rpt_variability_compact';
xl = [0 0.4];

cnt = 1;
[all_spline_alpha,all_spline_alpha_base,all_mod_alpha,all_Enum,all_ntrials,all_avgrates] = deal([]);
all_mon = cell(0);
all_tp = [];
all_ep = [];
all_Fmod_alpha = [];
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    monk_name = expt_mname{Elist_cnt};
    for bori_cnt = 1:2
        bar_ori = expt_oris(Elist_cnt,bori_cnt);
        if ~isnan(bar_ori)
            fprintf('Loading %s on Expt %s ori %d\n',base_sname,Expt_name,bar_ori);
            data_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
            
            sname = [data_dir base_sname sprintf('_ori%d.mat',bar_ori)];
            load(sname);
            
            for cc = 1:length(EP_data)
                if ~isnan(EP_data(cc).ov_avg_BS)
                    fprintf('Cell %d/%d\n',cc,length(EP_data));
                    
                    %                    subplot(2,1,1);
                    %                     shadedErrorBar(EP_params.eval_xx,EP_data(cc).spline_pred_looEP,EP_data(cc).spline_boot_sd,{'color','r'});
                    %                     hold on
                    %                     xlim(xl);
                    %                     plot(EP_params.eval_xx,EP_data(cc).spline_pred_baseEP,'g')
                    %                     plot(EP_params.EP_bin_centers,EP_data(cc).var_ep_binned);
                    %                     line(xl,[0 0] + EP_data(cc).pair_psth_var,'color','k');
                    %                    subplot(2,1,2);
                    %                     shadedErrorBar(EP_params.eval_xx,EP_data(cc).spline_pred_looEP,EP_data(cc).spline_boot_sd,{'color','r'});
                    %                     hold on
                    %                     xlim([0 0.05]);
                    %                     plot(EP_params.eval_xx,EP_data(cc).spline_pred_baseEP,'g')
                    %                     plot(EP_params.EP_bin_centers,EP_data(cc).var_ep_binned);
                    %                     line(xl,[0 0] + EP_data(cc).pair_psth_var,'color','k');
                    
                    %                     fprintf('Spline:%.3f  No LOO: %.3f   Model: %.3f\n',spline_alpha,spline_alpha_base,mod_alpha);
                    %                     fprintf('PSTH FF: %.3f  Spline FF: %.3f\n',EP_data(cc).psth_FF,EP_data(cc).spline_FF);
                    %                     fprintf('Ntrials %.2f avg rate %3.f\n',sum(EP_data(cc).n_utrials),EP_data(cc).ov_avg_BS/EP_params.base_dt);
                    %                     pause
                    %                     clf
                    %
                    %                     fprintf('\n\n');
                    
                    all_tp = cat(1,all_tp,EP_data(cc).tune_props);
                    all_ep = cat(1,all_ep,EP_data(cc));
                    all_Enum = cat(1,all_Enum,Elist_cnt);
                    all_mon{cnt} = monk_name; cnt = cnt + 1;
                    
                end
            end
        end
    end
end

%%

all_ntrials = arrayfun(@(x) sum(x.n_utrials),all_ep);
all_avgrates = [all_ep(:).ov_avg_BS]'/EP_params.base_dt;
all_spline_vars = arrayfun(@(x) x.spline_pred_looEP(1),all_ep);
all_psth_vars = [all_ep(:).pair_psth_var]';
all_ball_vars = cat(2,all_ep.eps_ball_var)';
all_ball_alphas = bsxfun(@rdivide,all_psth_vars,all_ball_vars);
all_spline_alphas = all_psth_vars./all_spline_vars;

all_mod_psth_vars = arrayfun(@(x) mean(x.mod_psth_vars),all_ep);
all_mod_tot_vars = arrayfun(@(x) mean(x.mod_tot_vars),all_ep);
all_mod_alphas = all_mod_psth_vars./all_mod_tot_vars;

RF_ecc = [all_tp(:).RF_ecc]';

%%
% close all

uset = find(all_ntrials >= 50 & all_avgrates >= 5);
f1 = figure(); hold on
jbe_units = uset(strcmp(all_mon(uset),'jbe'));
lem_units = uset(strcmp(all_mon(uset),'lem'));
% plot(all_mod_alphas(lem_units),all_spline_alphas(lem_units),'o');
% plot(all_mod_alphas(jbe_units),all_spline_alphas(jbe_units),'ro');
plot(all_mod_alphas(lem_units),all_ball_alphas(lem_units,1),'o');
plot(all_mod_alphas(jbe_units),all_ball_alphas(jbe_units,1),'ro');
line([0 1],[0 1],'color','k');

% f2 = figure();
% plot(RF_ecc(uset),all_spline_alpha(uset),'o')


%%

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
expt_oris = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
expt_mname = repmat({'lem'},1,length(Expt_list));

Expt_list = cat(2,Expt_list,{'G085','G086','G087','G088','G089','G091','G093','G095'});
expt_oris = cat(1,expt_oris,[0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan]);
expt_mname = cat(2,expt_mname,repmat({'jbe'},1,8));

Expt_list = cat(2,Expt_list,{'M005','M309','M009','M010','M011','M012'});
expt_oris = cat(1,expt_oris,[50 nan; 120 nan; 0 nan; 60 nan; 160 nan; 0 nan]);
expt_mname = cat(2,expt_mname,{'jbe','lem','jbe','jbe','jbe','jbe'});


new_sname = [data_dir 'model_variability_compact' sprintf('_ori%d',bar_ori)];
temp = load(new_sname);
mod_Data = temp.EP_data(temp.targs);

for cc = 1:length(mod_Data)
    if ~isempty(mod_Data(cc).ModData.tune_props)
        all_Fmod_alpha = cat(1,all_Fmod_alpha,mod_Data(cc).alpha_funs(end));
        all_tp = cat(1,all_tp,mod_Data(cc).ModData.tune_props);
        all_Enum = cat(1,all_Enum,Elist_cnt);
        all_mon{cnt} = monk_name; cnt = cnt + 1;
    end
end
