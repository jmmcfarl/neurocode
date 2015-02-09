clear all
close all

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};%NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
n_probes = 24;
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];

enum = 9; 
ii = 1;

Expt_name = Expt_list{enum};
bar_ori = ori_list(enum,ii);

anal_dir = ['~/Analysis/bruce/' Expt_name '/lfp_models/'];
cd(anal_dir);

% sname = 'lfp_models2';
% sname = [sname sprintf('_ori%d',bar_ori)];

% load(sname);

mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];
mod_data_name = 'corrected_models2';
fprintf('Loading model fits\n');
mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];
cd(mod_data_dir)
load(mod_data_name);

%%

close all
fig_dir = '/home/james/Desktop/K99_figures/';

all_NEfilts = [];
all_NIfilts = [];
all_base_imp = [];
all_lfp_imp = [];
all_off_imp = [];
all_ESD = [];
all_ISD = [];
all_Epow = [];
all_Ipow = [];
all_Opow = [];
all_Cind = [];
all_Eind = [];
all_models = [];
all_rate_modfac = [];
for enum = 1:10
  enum
  for ii = 1:2
        
        Expt_name = Expt_list{enum};
        bar_ori = ori_list(enum,ii);
        if ~isnan(bar_ori)
            anal_dir = ['~/Analysis/bruce/' Expt_name '/lfp_models/'];
            cd(anal_dir);
            
            sname = 'lfp_models4';
            sname = [sname sprintf('_ori%d',bar_ori)];
            
            load(sname);
            
            mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];
            mod_data_name = 'corrected_models2';
            fprintf('Loading model fits\n');
            mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];
            cd(mod_data_dir)
            load(mod_data_name);            
            
            for cc = 25:length(lfp_models)
                if ~isempty(lfp_models(cc).Egain_SD)
                   cur_mod = ModData(cc).rectGQM;
                   
                   cur_rate_modfac = lfp_models(cc).gain_rate_diff_SD/lfp_models(cc).base_rate_SD;
                   
                   mod_filt_signs = [cur_mod.mods(:).sign];
                   all_Cind = cat(1,all_Cind,cc);
                   all_Eind = cat(1,all_Eind,enum);
                   all_models = cat(1,all_models,lfp_models(cc).ModData);
                   
                   all_NEfilts = cat(1,all_NEfilts,sum(mod_filt_signs==1));
                   all_NIfilts = cat(1,all_NIfilts,sum(mod_filt_signs==-1));
                                      
                   all_rate_modfac = cat(1,all_rate_modfac,cur_rate_modfac);
                   all_base_imp = cat(1,all_base_imp,lfp_models(cc).base_xvImp);
                   all_lfp_imp = cat(1,all_lfp_imp,lfp_models(cc).lfp_xvImp);
                   all_off_imp = cat(1,all_off_imp,lfp_models(cc).off_xvImp);
                   all_ESD = cat(1,all_ESD,lfp_models(cc).Egain_SD);
                   all_ISD = cat(1,all_ISD,lfp_models(cc).Igain_SD);
                   all_Epow = cat(1,all_Epow,lfp_models(cc).pow_spectra(:,1)');
                   all_Ipow = cat(1,all_Ipow,lfp_models(cc).pow_spectra(:,2)');
                   all_Opow = cat(1,all_Opow,lfp_models(cc).pow_spectra(:,3)');
                end
            end
        end
    end
end

%%
[max_Efilt_norm,max_Ifilt_norm,max_Efilt_peak,max_Ifilt_peak] = deal(nan(length(all_Eind),1));
for ii = 1:length(all_Eind)
   cur_mod = all_models(ii).rectGQM;
   cur_mod_signs = [cur_mod.mods(:).sign];
   cur_filts = [cur_mod.mods(:).filtK];
   filt_norms = sqrt(sum(cur_filts.^2));
   filt_peaks = max(abs(cur_filts));
   max_Efilt_norm(ii) = max(filt_norms(cur_mod_signs == 1));
   max_Ifilt_norm(ii) = max(filt_norms(cur_mod_signs == -1));
   max_Efilt_peak(ii) = max(filt_peaks(cur_mod_signs == 1));
   max_Ifilt_peak(ii) = max(filt_peaks(cur_mod_signs == -1));
end

C = unique([all_Eind all_Cind],'rows');
n_un = size(C,1);
to_use = true(length(all_Eind),1);
for ii = 1:n_un
   curset = find(all_Eind == C(ii,1) & all_Cind == C(ii,2));
   if length(curset) > 1
       [~,best_loc] = max(all_base_imp(curset));
       curset(best_loc) = [];
       to_use(curset) = false;
   end
end

all_rates = arrayfun(@(x) x.unit_data.avg_rate,all_models);
tot_spks = arrayfun(@(x) x.unit_data.tot_spikes,all_models);

rel_lfp_imp = (all_lfp_imp - all_base_imp)./all_base_imp;
rel_off_imp = (all_off_imp - all_base_imp)./all_base_imp;

rel_stim_def = all_base_imp./all_lfp_imp;

% use_mods = find(all_NEfilts > 1 & all_NIfilts > 1 & to_use);
% use_bmods = find(all_lfp_imp > all_off_imp & all_rates > 5 & to_use);
% use_bmods = find(max_Efilt_norm > 0.4 & max_Ifilt_norm > 0.4 & to_use);

%select usable units based on a minimum firing rate, as well as a minimum
%filter amplitude (peak amplitude of biggest E or I filter). This seems to
%be a reasonable way to weed out cells without a clear E or I filter.
min_rate = 5;
min_famp = 0.1;
use_bmods = find(all_rates >= min_rate & max_Efilt_peak >= min_famp & max_Ifilt_peak >= min_famp & to_use);
use_Emods = find(all_rates >= min_rate & max_Efilt_peak >= min_famp & to_use);
use_Imods = find(all_rates >= min_rate & max_Ifilt_peak >= min_famp & to_use);

close all
f1 = figure();
plot(all_ESD(use_bmods),all_ISD(use_bmods),'o');
xlim([0 0.55]); ylim([0 0.55]);
line([0 0.55],[0 0.55],'color','k')
xlabel('Excitatory gain SD');
ylabel('Inhibitory gain SD');

% fname = [fig_dir 'gain_SD_scatter.pdf'];
% fig_width = 4; rel_height = 1;
% figufy(f1);
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


f = lfp_models(25).spectra_f;
all_norm_Epow = bsxfun(@rdivide,all_Epow,nanmean(all_Epow,2));
all_norm_Ipow = bsxfun(@rdivide,all_Ipow,nanmean(all_Ipow,2));
all_norm_Opow = bsxfun(@rdivide,all_Opow,nanmean(all_Opow,2));

[~,Emaxpowloc] = max(all_Epow,[],2);
[~,Imaxpowloc] = max(all_Ipow,[],2);

f2 = figure();
plot(f(Emaxpowloc(use_Emods)),all_ESD(use_Emods),'o');hold on
plot(f(Imaxpowloc(use_Imods)),all_ISD(use_Imods),'ro');
% xlim([0 0.45]);
xlim([0 20]);
ylim([0 0.55]);
xlabel('Peak frequency (Hz)');
ylabel('Gain modulation (SD)');

% fname = [fig_dir 'gain_freq_scatter.pdf'];
% fig_width = 4; rel_height = 1;
% figufy(f2);
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

f3 = figure();
boxplot(rel_lfp_imp(use_bmods)+1);
box off
ylim([1 2.2]);
ylabel('Relative performance');
fname = [fig_dir 'xvLL_imp.pdf'];
fig_width = 4; rel_height = 0.8;
figufy(f3);
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);


%%
% lags = lfp_models(25).lags;
% new_lags = -0.2:0.0025:0.2;
% all_xcov = [];
% for cc = 25:length(lfp_models)
%     cur_inter_xc = interp1(lags,lfp_models(cc).xcov,new_lags,'spline');
%     all_xcov = cat(2,all_xcov,cur_inter_xc');
% end
% 
% close all
% fig_dir = '/home/james/Desktop/K99_figures/';
% f1 = figure();
% uset = [1 2 3 6];
% cmap = jet(4);
% cmap(3,:) = [0.2 0.8 0.2];
% for ii = 1:4
%     plot(new_lags,all_xcov(:,uset(ii)),'color',cmap(ii,:));
%     hold on
% end
% xlabel('Lag (s)');
% ylabel('Correlation');
% 
% fname = [fig_dir 'xcorr_examples.pdf'];
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);
% 
% %%
% clc
% close all
% for cc = 25:length(lfp_models)
%     cc
%     if ~isempty(lfp_models(cc).Egain_SD)
%     lfp_models(cc)
%     
%     lfp_ampmodels(cc)
%     
%     figure(1)
%     plot(lfp_models(cc).spectra_f,lfp_models(cc).pow_spectra(:,1:2));
%     
%     figure(2)
%     plot(lfp_models(cc).lags,lfp_models(cc).xcov);
%     
%         figure(3)
%     plot(lfp_models(cc).spectra_f,lfp_models(cc).pow_spectra(:,3),'k');
% 
%     h = NMMdisplay_model(ModData(cc).rectGQM);
%     
%     pause
%     figure(1);clf;figure(2);clf; figure(3);clf;close(h.stim_filts);
%     end
% end
% 
