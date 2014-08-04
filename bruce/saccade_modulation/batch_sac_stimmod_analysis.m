%%
close all
clear all
% Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
% Expt_list = {'M281','M287','M289','M294'};
% Expt_list = {'M266','M270','M275','M277'};
% Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095','M266','M270','M275','M277'};
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095','M266','M270','M275','M277','M281','M287','M289','M294'};
% Expt_list = {'M266'};
% sname = 'sacStimProcMUA';
sname = 'sacStimProc_v2';
% rmfield_list = {'gsac_sta','gsac_sta_norm','gsac_sta_sm','gsac_sta_sm_norm','gsac_abssta',...
%     'gsac_abssta_norm','gsac_abssta_sm','gsac_abssta_sm_norm','msac_sta','msac_sta_norm','msac_sta_sm',...
%     'msac_sta_sm_norm','gsac_phaseDep_filt','gsac_phaseInd_filt','msac_phaseDep_filt','msac_phaseInd_filt'};
rmfield_list = {};

fig_dir = '/home/james/Analysis/bruce/saccade_modulation/';

all_jbe_data = [];
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
    cd(save_dir)
    load(sname)
    
    ucells = arrayfun(@(x) length(x.ModData),sacStimProc) > 0;
    cur_data = sacStimProc(ucells);
    cur_data = rmfields(cur_data,rmfield_list);
    [cur_data.expt_num] = deal(Expt_num);
    
    all_jbe_data = cat(2,all_jbe_data,cur_data);
end
%%
for ii = 1:length(all_jbe_data)
    cur_mod = all_jbe_data(ii).ModData.rectGQM;
    mod_signs = [cur_mod.mods(:).sign];
    mod_signs(end) = 0; %try excluding rectified suppressive filter
    sd = cur_mod.stim_params(1).stim_dims;
    mod_filts = reshape([cur_mod.mods(:).filtK],[sd(1) sd(2) length(mod_signs)]);
    tkerns = squeeze(std(mod_filts,[],2));
    all_Ekerns(ii,:) = mean(tkerns(:,mod_signs == 1),2);
    all_Ikerns(ii,:) = mean(tkerns(:,mod_signs == -1),2);
    
    cur_mod = all_jbe_data(ii).gsac_post_EImod;
    all_pmod_off(ii,:) = cur_mod.mods(2).filtK;
    all_pmod_Ekern(ii,:) = cur_mod.mods(3).filtK;
    all_pmod_Ikern(ii,:) = cur_mod.mods(4).filtK;

    cur_mod = all_jbe_data(ii).gsac_post_singmod;
    all_spmod_off(ii,:) = cur_mod.mods(2).filtK;
    all_spmod_gain(ii,:) = cur_mod.mods(3).filtK;
    
%     cur_mod = all_jbe_data(ii).msac_post_fullmod;
%     all_pMmod_off(ii,:) = cur_mod.mods(2).filtK;
%     all_pMmod_Ekern(ii,:) = cur_mod.mods(3).filtK;
%     all_pMmod_Ikern(ii,:) = cur_mod.mods(4).filtK;

end
flen = sd(1);

[ikern_max,ikern_loc] = max(all_Ikerns,[],2);
[ekern_max,ekern_loc] = max(all_Ekerns,[],2);
kern_peak_range = [4:12];
min_kern_peak = 0.02;
usable_ikern = find(ismember(ikern_loc,kern_peak_range) & ikern_max >= min_kern_peak);
usable_ekern = find(ismember(ekern_loc,kern_peak_range) & ekern_max >= min_kern_peak);

%%
avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_jbe_data);
rec_dur = arrayfun(@(x) x.ModData.unit_data.N_used_samps,all_jbe_data)*dt/60; %in min
clust_iso_dist = arrayfun(@(x) x.ModData.unit_data.SU_isodist,all_jbe_data);
nullLL = arrayfun(@(x) x.ModData.unit_data.nullLL,all_jbe_data);
modLL = arrayfun(@(x) x.ModData.rectGQM.LL_seq(end),all_jbe_data);
LLimp = (modLL-nullLL)/log(2);
all_ov_smod_info = [all_jbe_data(:).gsac_spost_ov_modinfo];
all_ov_submod_info = [all_jbe_data(:).gsac_sub_ov_modinfo];
all_ov_TBmod_info = [all_jbe_data(:).gsac_ov_TB_info];
% all_ov_info_rate = all_ov_info.*avg_rates;

min_rate = 5;
min_dur = 20; 
min_ss_info = 0.01;
% used_units = 1:length(avg_rates);
used_units = find(avg_rates >= min_rate & rec_dur >= min_dur & all_ov_smod_info >= min_ss_info);

usable_ikern = usable_ikern(ismember(usable_ikern,used_units));
usable_ekern = usable_ekern(ismember(usable_ekern,used_units));
usable_bothkerns = intersect(usable_ikern,usable_ekern);

% for ii = 1:length(all_jbe_data)
%     NMMdisplay_model(all_jbe_data(ii).stimMod);
%     ii
%     [avg_rates(ii) rec_dur(ii) all_ov_info(ii) all_ov_info_rate(ii) ]
%     [ismember(ii,used_units) ismember(ii,usable_ekern) ismember(ii,usable_ikern)]
%     pause
%     close 
% end


%%
all_gsac_rates = reshape([all_jbe_data(:).gsac_avg_rate],[],length(all_jbe_data))';
all_gsac_nrates = bsxfun(@rdivide,all_gsac_rates,avg_rates'*dt);

all_gsac_smod_info = reshape([all_jbe_data(:).gsac_spost_modinfo],[],length(all_jbe_data))';
all_gsac_submod_info = reshape([all_jbe_data(:).gsac_sub_modinfo],[],length(all_jbe_data))';
% all_gsac_minfo = reshape([all_jbe_data(:).gsac_postmod_info],[],length(all_jbe_data))';
% all_gsac_minfo = reshape([all_jbe_data(:).gsac_postLL_info],[],length(all_jbe_data))';
all_gsac_TBinfo = reshape([all_jbe_data(:).gsac_TB_info],[],length(all_jbe_data))';

all_gsac_smod_inforate = all_gsac_smod_info.*all_gsac_rates;
all_gsac_submod_inforate = all_gsac_submod_info.*all_gsac_rates;
all_gsac_TBinforate = all_gsac_TBinfo.*all_gsac_rates;

all_gsac_rel_smod_info = bsxfun(@rdivide,all_gsac_smod_info,all_ov_smod_info');
all_gsac_rel_submod_info = bsxfun(@rdivide,all_gsac_submod_info,all_ov_submod_info');
all_gsac_rel_TBinfo = bsxfun(@rdivide,all_gsac_TBinfo,all_ov_TBmod_info');
all_gsac_rel_smod_inforate = bsxfun(@rdivide,all_gsac_smod_inforate,all_ov_smod_info'.*avg_rates'*dt);
all_gsac_rel_submod_inforate = bsxfun(@rdivide,all_gsac_submod_inforate,all_ov_submod_info'.*avg_rates'*dt);
all_gsac_rel_TBinforate = bsxfun(@rdivide,all_gsac_TBinforate,all_ov_TBmod_info'.*avg_rates'*dt);

% all_gsac_gaink = cell2mat(arrayfun(@(x) x.gsacGainMod.gain_kernel,all_jbe_data,'uniformoutput',0))';
% all_gsac_stimk = cell2mat(arrayfun(@(x) x.gsacGainMod.stim_kernel,all_jbe_data,'uniformoutput',0))';
% all_gsac_offk = cell2mat(arrayfun(@(x) x.gsacGainMod.off_kernel,all_jbe_data,'uniformoutput',0))';

%%
ulags = find(slags*dt >= 0);
[gsac_peak_mods,gsac_peak_locs] = max(all_gsac_nrates(:,ulags),[],2);
[gsac_val_mods,gsac_val_locs] = min(all_gsac_nrates(:,ulags),[],2);

[~,gsac_peak_ord] = sort(gsac_peak_mods(used_units));
[~,gsac_val_ord] = sort(1-gsac_val_mods(used_units));

h = figure();
subplot(1,3,1)
imagesc(slags*dt,1:length(used_units),all_gsac_nrates(used_units(gsac_peak_ord),:));
caxis([0.25 1.75])
ylabel('Unit number');
xlabel('Time (s)');
title('Relative rates');

subplot(1,3,2)
imagesc(slags*dt,1:length(used_units),all_gsac_rel_TBinfo(used_units(gsac_peak_ord),:));
caxis([0 2])
ylabel('Unit number');
xlabel('Time (s)');
title('SS Info');

subplot(1,3,3)
imagesc(slags*dt,1:length(used_units),all_gsac_rel_TBinforate(used_units(gsac_peak_ord),:));
caxis([0 2])
ylabel('Unit number');
xlabel('Time (s)');
title('Info rate');

% fig_width = 8; rel_height = 0.3;
% figufy(h);
% fname = [fig_dir 'Gsac_enhance_sort.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

% [~,ord] = sort(all_ov_info_rate(used_units));
% % [~,ord] = sort(all_ov_info_rate(used_units));
% % [~,ord] = sort(clust_iso_dist(used_units));
% figure;
% imagesc(slags*dt,1:length(used_units),1+all_gsac_stimk(used_units(ord),:));

%%
gsac_peak_times = slags(ulags(gsac_peak_locs))*dt;
gsac_val_times = slags(ulags(gsac_val_locs))*dt;
good_peaks = find(~ismember(gsac_peak_times,slags(ulags([1 end]))*dt));
good_vals = find(~ismember(gsac_val_times,slags(ulags([1 end]))*dt));
good_peaks(~ismember(good_peaks,used_units)) = [];
good_vals(~ismember(good_vals,used_units)) = [];

tax = (0:(flen-1))*dt + dt/2;
ekern_temp_loc = tax(ekern_loc);

jit_amp = 0.001;
jit_vals = randn(length(gsac_peak_times),2)*jit_amp;

h = figure; 
subplot(2,1,1);hold on
plot(ekern_temp_loc(good_vals)'+jit_vals(good_vals,1),gsac_val_times(good_vals)'+jit_vals(good_vals,2),'r.','markersize',10);
xlabel('Stimulus filter delay (s)');
ylabel('Sac-suppression timing (s)');
title('Suppression');
[rV,statsV] = robustfit(ekern_temp_loc(good_vals),gsac_val_times(good_vals));
xx = linspace(0.03,0.08,50);
plot(xx,rV(1)+rV(2)*xx,'k');

subplot(2,1,2);hold on
plot(ekern_temp_loc(good_peaks)'+jit_vals(good_peaks,1),gsac_peak_times(good_peaks)'+jit_vals(good_peaks,2),'.','markersize',10);
title('Enhancement');
xlabel('Stimulus filter delay (s)');
ylabel('Sac-enhancement timing (s)');
[rP,statsP] = robustfit(ekern_temp_loc(good_peaks),gsac_peak_times(good_peaks));
xx = linspace(0.03,0.08,50);
plot(xx,rP(1)+rP(2)*xx,'k');

% fig_width = 4; rel_height = 1.6;
% figufy(h);
% fname = [fig_dir 'Gsac_filt_sacmod_timing.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%% compare pre-post gains
h = figure; 
subplot(2,1,1); hold on
h1=shadedErrorBar(slags*dt,1+mean(all_gsac_gaink(used_units,:)),std(all_gsac_gaink(used_units,:))/sqrt(length(used_units)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_gsac_stimk(used_units,:)),std(all_gsac_stimk(used_units,:))/sqrt(length(used_units)),{'color','k'});
legend([h1.mainLine h2.mainLine],{'Pre-filt','Post-filt'},'Location','southeast');
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Gains');

subplot(2,1,2)
shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');


% fig_width = 4; rel_height = 1.6;
% figufy(h);
% fname = [fig_dir 'Gsac_prepost_gains.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%% information compare
h = figure; 
subplot(2,1,1); hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_TBinfo(used_units,:)),std(all_gsac_rel_TBinfo(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_gsac_rel_TBinforate(used_units,:)),std(all_gsac_rel_TBinforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_smod_info(used_units,:)),std(all_gsac_rel_smod_info(used_units,:))/sqrt(length(used_units)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_rel_smod_inforate(used_units,:)),std(all_gsac_rel_smod_inforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
% h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_submod_info(used_units,:)),std(all_gsac_rel_submod_info(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_gsac_rel_submod_inforate(used_units,:)),std(all_gsac_rel_submod_inforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
xl = xlim(); yl = ylim();
legend([h1.mainLine h2.mainLine],{'SS Info','Info rate'},'Location','Southeast');
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative information');

subplot(2,1,2)
shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

% fig_width = 4; rel_height = 1.6;
% figufy(h);
% fname = [fig_dir 'Gsac_info.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%% Post gain 
h = figure; 
subplot(2,1,1);hold on
h1=shadedErrorBar(slags*dt,1+mean(all_spmod_off(used_units,:)),std(all_spmod_off(used_units,:))/sqrt(length(used_units)),{'color','r'});
h2=shadedErrorBar(slags*dt,1+mean(all_spmod_gain(used_units,:)),std(all_spmod_gain(used_units,:))/sqrt(length(used_units)),{'color','b'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

subplot(2,1,2);hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_minfo(used_units,:)),std(all_gsac_rel_minfo(used_units,:))/sqrt(length(used_units)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_rel_minforate(used_units,:)),std(all_gsac_rel_minforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

fig_width = 4; rel_height = 1.6;
figufy(h);
fname = [fig_dir 'Gsac_post_gain.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);

%% E-gain vs I-gain
h = figure; 
subplot(3,1,1);hold on
% h1=shadedErrorBar(slags*dt,1+mean(all_pmod_Ekern(used_units,:)),std(all_pmod_Ekern(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(slags*dt,1+mean(all_pmod_Ikern(used_units,:)),std(all_pmod_Ikern(used_units,:))/sqrt(length(used_units)),{'color','r'});
h1=shadedErrorBar(slags*dt,1+mean(all_pmod_Ekern(usable_ekern,:)),std(all_pmod_Ekern(usable_ekern,:))/sqrt(length(usable_ekern)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_pmod_Ikern(usable_ikern,:)),std(all_pmod_Ikern(usable_ikern,:))/sqrt(length(usable_ikern)),{'color','r'});
xl = xlim(); yl = ylim();
legend([h1.mainLine h2.mainLine],{'E gain','I gain'});
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Gain');

subplot(3,1,2)
shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

tax = (0:(flen-1))*dt + dt/2;
subplot(3,1,3); hold on
% h1=shadedErrorBar(tax,mean(all_Ekerns(used_units,:)),std(all_Ekerns(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(tax,mean(all_Ikerns(used_units,:)),std(all_Ikerns(used_units,:))/sqrt(length(used_units)),{'color','r'});
h1=shadedErrorBar(tax,mean(all_Ekerns(usable_ekern,:)),std(all_Ekerns(usable_ekern,:))/sqrt(length(usable_ekern)),{'color','b'});
h2=shadedErrorBar(tax,mean(all_Ikerns(usable_ikern,:)),std(all_Ikerns(usable_ikern,:))/sqrt(length(usable_ikern)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'E kern','I kern'});
xlabel('Lag (s)');
ylabel('Temporal kernel');

% fig_width = 4; rel_height = 2.2;
% figufy(h);
% fname = [fig_dir 'Gsac_EI_gains.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%% EI ratio 
avg_gsac_nrate = mean(all_gsac_nrates(used_units,:));
[~,minloc] = min(avg_gsac_nrate);
[~,maxloc] = max(avg_gsac_nrate);
all_EI_ratio = (1+all_pmod_Ekern)./(1+all_pmod_Ikern);

h = figure; 
subplot(2,1,1);hold on
% h1=shadedErrorBar(slags*dt,mean(all_EI_ratio(used_units,:)),std(all_EI_ratio(used_units,:))/sqrt(length(used_units)),{'color','b'});
h1=shadedErrorBar(slags*dt,mean(all_EI_ratio(usable_bothkerns,:)),std(all_EI_ratio(usable_bothkerns,:))/sqrt(length(usable_bothkerns)),{'color','b'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
line(slags([minloc minloc])*dt,yl,'color','r','linestyle','--');
line(slags([maxloc maxloc])*dt,yl,'color','g','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('E/I ratio');

subplot(2,1,2);hold on
% shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
shadedErrorBar(slags*dt,mean(all_gsac_nrates(usable_bothkerns,:)),std(all_gsac_nrates(usable_bothkerns,:))/sqrt(length(usable_bothkerns)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');
line(slags([minloc minloc])*dt,yl,'color','r','linestyle','--');
line(slags([maxloc maxloc])*dt,yl,'color','g','linestyle','--');

% fig_width = 4; rel_height = 1.6;
% figufy(h);
% fname = [fig_dir 'Gsac_EI_ratio.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%%

all_msac_rates = reshape([all_jbe_data(:).msac_avg_rate],[],length(all_jbe_data))';
all_msac_nrates = bsxfun(@rdivide,all_msac_rates,avg_rates'*dt);

% all_msac_minfo = reshape([all_jbe_data(:).msac_mod_info],[],length(all_jbe_data))';
all_msac_minfo = reshape([all_jbe_data(:).msac_postmod_info],[],length(all_jbe_data))';
all_msac_TBinforate = reshape([all_jbe_data(:).msac_TB_info],[],length(all_jbe_data))';
all_msac_TBinfo = all_msac_TBinforate./all_msac_rates;

all_msac_minforate = all_msac_minfo.*all_msac_rates;

all_msac_rel_minfo = bsxfun(@rdivide,all_msac_minfo,all_ov_info');
all_msac_rel_TBinfo = bsxfun(@rdivide,all_msac_TBinfo,all_ov_info');
all_msac_rel_minforate = bsxfun(@rdivide,all_msac_minforate,all_ov_info'.*avg_rates'*dt);
all_msac_rel_TBinforate = bsxfun(@rdivide,all_msac_TBinforate,all_ov_info'.*avg_rates'*dt);

% all_msac_gaink = cell2mat(arrayfun(@(x) x.msacGainMod.gain_kernel,all_jbe_data,'uniformoutput',0))';
% all_msac_stimk = cell2mat(arrayfun(@(x) x.msacGainMod.stim_kernel,all_jbe_data,'uniformoutput',0))';
% all_msac_offk = cell2mat(arrayfun(@(x) x.msacGainMod.off_kernel,all_jbe_data,'uniformoutput',0))';


%%
figure; 
subplot(2,1,1); hold on
shadedErrorBar(slags*dt,mean(all_msac_offk(used_units,:)),std(all_msac_offk(used_units,:))/sqrt(length(used_units)),{'color','r'});
shadedErrorBar(slags*dt,mean(all_msac_gaink(used_units,:)),std(all_msac_gaink(used_units,:))/sqrt(length(used_units)),{'color','b'});
shadedErrorBar(slags*dt,mean(all_msac_stimk(used_units,:)),std(all_msac_stimk(used_units,:))/sqrt(length(used_units)),{'color','k'});
xl = xlim(); yl = ylim();
line(xl,[0 0],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Filter amp');

subplot(2,1,2)
shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

%%
figure; 
subplot(2,1,1); hold on
h1=shadedErrorBar(slags*dt,mean(all_msac_rel_TBinfo(used_units,:)),std(all_msac_rel_TBinfo(used_units,:))/sqrt(length(used_units)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_msac_rel_TBinforate(used_units,:)),std(all_msac_rel_TBinforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
% h3=shadedErrorBar(slags*dt,mean(all_msac_rel_minfo(used_units,:)),std(all_msac_rel_minfo(used_units,:))/sqrt(length(used_units)),{'color','r'});
% h4=shadedErrorBar(slags*dt,mean(all_msac_rel_minforate(used_units,:)),std(all_msac_rel_minforate(used_units,:))/sqrt(length(used_units)),{'color','g'});
xl = xlim(); yl = ylim();
legend([h1.mainLine h2.mainLine],{'SS Info','Info rate'});
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative information');

subplot(2,1,2)
shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

%%
avg_msac_nrate = mean(all_msac_nrates(used_units,:));
[~,minloc] = min(avg_msac_nrate);
[~,maxloc] = max(avg_msac_nrate);
all_EI_ratio = (1+all_pMmod_Ekern)./(1+all_pMmod_Ikern);
figure; 
subplot(2,1,1);hold on
h1=shadedErrorBar(slags*dt,mean(all_EI_ratio(used_units,:)),std(all_EI_ratio(used_units,:))/sqrt(length(used_units)),{'color','b'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
line(slags([minloc minloc])*dt,yl,'color','r','linestyle','--');
line(slags([maxloc maxloc])*dt,yl,'color','g','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('E/I ratio');

subplot(2,1,2);hold on
shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');
line(slags([minloc minloc])*dt,yl,'color','r','linestyle','--');
line(slags([maxloc maxloc])*dt,yl,'color','g','linestyle','--');

%%
figure; 
subplot(2,1,1);hold on
h1=shadedErrorBar(slags*dt,1+mean(all_pmod_Ekern(used_units,:)),std(all_pmod_Ekern(used_units,:))/sqrt(length(used_units)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_pmod_Ikern(used_units,:)),std(all_pmod_Ikern(used_units,:))/sqrt(length(used_units)),{'color','r'});
xl = xlim(); yl = ylim();
legend([h1.mainLine h2.mainLine],{'E gain','I gain'});
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Gain');

subplot(2,1,2)
shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');


%%
figure;
hold on
shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','b'});
shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

% figure;
% hold on
% shadedErrorBar(slags*dt,mean(all_gsac_offk(used_units,:)),std(all_gsac_offk(used_units,:))/sqrt(length(used_units)),{'color','b'});
% shadedErrorBar(slags*dt,mean(all_msac_offk(used_units,:)),std(all_msac_offk(used_units,:))/sqrt(length(used_units)),{'color','r'});

%%
h = figure; 
subplot(3,1,1); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_TBinfo(used_units,:)),std(all_gsac_rel_TBinfo(used_units,:))/sqrt(length(used_units)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_msac_rel_TBinfo(used_units,:)),std(all_msac_rel_TBinfo(used_units,:))/sqrt(length(used_units)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Big sac','Micro sac'},'Location','Southeast');
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('SS info');

subplot(3,1,2); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_TBinforate(used_units,:)),std(all_gsac_rel_TBinforate(used_units,:))/sqrt(length(used_units)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_msac_rel_TBinforate(used_units,:)),std(all_msac_rel_TBinforate(used_units,:))/sqrt(length(used_units)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Info rate');

subplot(3,1,3); hold on
shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','b'});
shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
ylim([0.6 1.4])
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

fig_width = 4; rel_height = 2.4;
figufy(h);
fname = [fig_dir 'Gsac_Msac_sacinfo.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);


%%
% tax = (0:(flen-1))*dt + dt/2;
% figure; hold on
% shadedErrorBar(tax,mean(all_Ekerns(used_units,:)),std(all_Ekerns(used_units,:))/sqrt(length(used_units)),{'color','b'});
% shadedErrorBar(tax,mean(all_Ikerns(used_units,:)),std(all_Ikerns(used_units,:))/sqrt(length(used_units)),{'color','r'});

%%
% figure; hold on
% shadedErrorBar(slags*dt,mean(all_msac_rel_minfo(used_units,:)),std(all_msac_rel_minfo(used_units,:))/sqrt(length(used_units)),{'color','r'});
% shadedErrorBar(slags*dt,mean(all_msac_rel_minforate(used_units,:)),std(all_msac_rel_minforate(used_units,:))/sqrt(length(used_units)),{'color','g'});

