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

sname = 'lfp_models';
sname = [sname sprintf('_ori%d',bar_ori)];

load(sname);

mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];
mod_data_name = 'corrected_models2';
fprintf('Loading model fits\n');
mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];
cd(mod_data_dir)
load(mod_data_name);

%%
clc
close all
for cc = 25:length(lfp_models)
    cc
    if ~isempty(lfp_models(cc).Egain_SD)
    lfp_models(cc)
    
    lfp_ampmodels(cc)
    
    figure(1)
    plot(lfp_models(cc).spectra_f,lfp_models(cc).pow_spectra(:,1:2));
    
    figure(2)
    plot(lfp_models(cc).lags,lfp_models(cc).xcov);
    
        figure(3)
    plot(lfp_models(cc).spectra_f,lfp_models(cc).pow_spectra(:,3),'k');

    h = NMMdisplay_model(ModData(cc).rectGQM);
    
    pause
    figure(1);clf;figure(2);clf; figure(3);clf;close(h.stim_filts);
    end
end

%%
lags = lfp_models(25).lags;
new_lags = -0.2:0.0025:0.2;
all_xcov = [];
for cc = 25:length(lfp_models)
    cur_inter_xc = interp1(lags,lfp_models(cc).xcov,new_lags,'spline');
    all_xcov = cat(2,all_xcov,cur_inter_xc');
end

close all
fig_dir = '/home/james/Desktop/K99_figures/';
f1 = figure();
uset = [1 2 3 6];
cmap = jet(4);
cmap(3,:) = [0.2 0.8 0.2];
for ii = 1:4
    plot(new_lags,all_xcov(:,uset(ii)),'color',cmap(ii,:));
    hold on
end
xlabel('Lag (s)');
ylabel('Correlation');

fname = [fig_dir 'xcorr_examples.pdf'];
fig_width = 4; rel_height = 0.8;
figufy(f1);
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


%%
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
for enum = 1:10
  enum
  for ii = 1:2
        
        Expt_name = Expt_list{enum};
        bar_ori = ori_list(enum,ii);
        if ~isnan(bar_ori)
            anal_dir = ['~/Analysis/bruce/' Expt_name '/lfp_models/'];
            cd(anal_dir);
            
            sname = 'lfp_models';
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
                   mod_filt_signs = [cur_mod.mods(:).sign];
                   all_NEfilts = cat(1,all_NEfilts,sum(mod_filt_signs==1));
                   all_NIfilts = cat(1,all_NIfilts,sum(mod_filt_signs==-1));
                   
                   all_base_imp = cat(1,all_base_imp,lfp_models(cc).base_xvImp);
                   all_lfp_imp = cat(1,all_lfp_imp,lfp_models(cc).lfp_xvImp);
                   all_off_imp = cat(1,all_off_imp,lfp_models(cc).off_xvImp);
                   all_ESD = cat(1,all_ESD,lfp_models(cc).Egain_SDxv);
                   all_ISD = cat(1,all_ISD,lfp_models(cc).Igain_SDxv);
                   all_Epow = cat(1,all_Epow,lfp_models(cc).pow_spectra(:,1)');
                   all_Ipow = cat(1,all_Ipow,lfp_models(cc).pow_spectra(:,2)');
                   all_Opow = cat(1,all_Opow,lfp_models(cc).pow_spectra(:,3)');
                end
            end
        end
    end
end
f = lfp_models(25).spectra_f;
rel_lfp_imp = (all_lfp_imp - all_base_imp)./all_base_imp;
rel_off_imp = (all_off_imp - all_base_imp)./all_base_imp;

all_norm_Epow = bsxfun(@rdivide,all_Epow,nanmean(all_Epow,2));
all_norm_Ipow = bsxfun(@rdivide,all_Ipow,nanmean(all_Ipow,2));
all_norm_Opow = bsxfun(@rdivide,all_Opow,nanmean(all_Opow,2));

use_mods = find(all_NEfilts > 1 & all_NIfilts > 1);

use_bmods = find(all_base_imp > 0.05 & all_base_imp < 1 & all_NEfilts > 1);

close all
f1 = figure();
plot(all_ESD(use_bmods),all_ISD(use_bmods),'o');
xlim([0 0.45]); ylim([0 0.45]);
line([0 0.45],[0 0.45],'color','k')
xlabel('Excitatory gain SD');
ylabel('Inhibitory gain SD');

fname = [fig_dir 'gain_SD_scatter.pdf'];
fig_width = 4; rel_height = 1;
figufy(f1);
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
