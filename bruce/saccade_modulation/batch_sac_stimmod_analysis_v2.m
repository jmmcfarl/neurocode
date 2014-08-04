%%
close all
clear all
fig_dir = '/home/james/Analysis/bruce/saccade_modulation/';

%% LOAD JBE HORIZONTAL
% sname = 'sacStimProc_sta';
sname = 'sacStimProc_v3';
% sname = 'sacStimProc_v2';
% sname = 'sacStimProc_mua';
tname = 'sac_trig_avg_data6';
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
rmfield_list = {};

all_hor_data = [];
all_hor_tdata = [];
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
    [cur_data.bar_ori] = deal(0);
    
    all_hor_data = cat(2,all_hor_data,cur_data);
    cur_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,cur_data);

    tdata = load(tname);
    tdat_SU_numbers = arrayfun(@(x) x.unit_data.SU_numbers,tdata.sua_trig_avgs);
    uset = find(ismember(tdat_SU_numbers,cur_SU_numbers));
    all_hor_tdata = cat(2,all_hor_tdata,tdata.sua_trig_avgs(uset));
end

%% LOAD JBE VERTICAL
% sname = 'sacStimProc_sta_vbars_unCor';
% sname = 'sacStimProc_sta_vbars';
sname = 'sacStimProc_v3_vbars';
% sname = 'sacStimProc_v2_vbars';
% sname = 'sacStimProc_mua_vbars';
tname = 'sac_trig_avg_data6_vbars';
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
rmfield_list = {};

all_ver_data = [];
all_ver_tdata = [];
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
    [cur_data.bar_ori] = deal(90);
    
    all_ver_data = cat(2,all_ver_data,cur_data);
    cur_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,cur_data);

    tdata = load(tname);
    tdat_SU_numbers = arrayfun(@(x) x.unit_data.SU_numbers,tdata.sua_trig_avgs);
    uset = find(ismember(tdat_SU_numbers,cur_SU_numbers));
    all_ver_tdata = cat(2,all_ver_tdata,tdata.sua_trig_avgs(uset));
end

% all_jbe_data = all_ver_data(:);
% all_jbe_tdata = all_ver_tdata(:);
% [all_jbe_data.animal] = deal('jbe');
% all_SU_data = all_jbe_data;
% all_SU_tdata = all_jbe_tdata;

%% COMBINE JBE DATA
hor_SU_number = arrayfun(@(x) x.ModData.unit_data.SU_number,all_hor_data);
hor_probe_num = arrayfun(@(x) x.ModData.unit_data.probe_number,all_hor_data);
hor_expt_num = [all_hor_data(:).expt_num];
hor_mod_info = [all_hor_data(:).gsac_spost_ov_modinfo];
hor_avg_rate = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_hor_data);
hor_mod_inforate = hor_mod_info.*hor_avg_rate;

ver_SU_number = arrayfun(@(x) x.ModData.unit_data.SU_number,all_ver_data);
ver_probe_num = arrayfun(@(x) x.ModData.unit_data.probe_number,all_ver_data);
ver_expt_num = [all_ver_data(:).expt_num];
ver_mod_info = [all_ver_data(:).gsac_spost_ov_modinfo];
ver_avg_rate = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_ver_data);
ver_mod_inforate = ver_mod_info.*ver_avg_rate;

% all_id_vecs = [hor_expt_num' hor_SU_number'; ver_expt_num' ver_SU_number'];
all_id_vecs = [hor_expt_num' hor_probe_num' hor_SU_number'; ver_expt_num' ver_probe_num' ver_SU_number'];
all_id_vecs(isnan(all_id_vecs(:,3)),3) = 0; %set NAN SU numbers to 0
all_info_rate = [hor_mod_inforate'; ver_mod_inforate'];
all_comb_data = [all_hor_data'; all_ver_data'];
all_comb_tdata = [all_hor_tdata'; all_ver_tdata'];
[C,IA,IC] = unique(all_id_vecs,'rows');

%for each array SU, use data from either the horizontal or vertical bar
%stimuli, depending on which had the better model (in terms of info rate)
all_jbe_data = [];
all_jbe_tdata = [];
all_jbe_nonpref_data = [];
all_jbe_nonpref_tdata = [];
all_jbe_pref_data = [];
all_jbe_pref_tdata = [];
for ii = 1:length(C)
%     curset = find(all_id_vecs(:,1) == C(ii,1) & all_id_vecs(:,2) == C(ii,2));
    curset = find(all_id_vecs(:,1) == C(ii,1) & all_id_vecs(:,2) == C(ii,2) & all_id_vecs(:,3) == C(ii,3));
    if length(curset) > 1
        [~,better] = max(all_info_rate(curset));
        worse = setdiff(curset,curset(better));
        curset = curset(better);
        
        all_jbe_nonpref_data = cat(1,all_jbe_nonpref_data,all_comb_data(worse));
        all_jbe_pref_data = cat(1,all_jbe_pref_data,all_comb_data(curset));
%         all_jbe_nonpref_tdata = cat(1,all_jbe_nonpref_tdata,all_comb_tdata(worse));
%         all_jbe_pref_tdata = cat(1,all_jbe_pref_tdata,all_comb_tdata(curset));
    end
    all_jbe_data = cat(1,all_jbe_data,all_comb_data(curset));
%     all_jbe_tdata = cat(1,all_jbe_tdata,all_comb_tdata(curset));
end

[all_jbe_data.animal] = deal('jbe');



%% LOAD LEM DATA
% sname = 'sacStimProc_sta';
sname = 'sacStimProc_v3';
% sname = 'sacStimProc_v2';
% sname = 'sacStimProc_mua';
tname = 'sac_trig_avg_data6';
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
% Expt_list = {'M266','M270','M275','M277','M294'};
bar_oris = [80 60 135 70 140 90 160 40];
rmfield_list = {};

all_lem_data = [];
all_lem_tdata = [];
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
    [cur_data.bar_ori] = deal(bar_oris(ee));
    
    all_lem_data = cat(1,all_lem_data,cur_data');
    cur_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,cur_data);

    tdata = load(tname);
    tdat_SU_numbers = arrayfun(@(x) x.unit_data.SU_numbers,tdata.sua_trig_avgs);
    uset = find(ismember(tdat_SU_numbers,cur_SU_numbers));
    all_lem_tdata = cat(1,all_lem_tdata,tdata.sua_trig_avgs(uset)');
end

[all_lem_data.animal] = deal('lem');

%% COMBINE JBE AND LEM DATA
tlags = tdata.trig_avg_params.lags*tdata.trig_avg_params.dt;

all_SU_data = [all_jbe_data; all_lem_data];
% all_SU_tdata = [all_jbe_tdata; all_lem_tdata];

%% SELECT USABLE CELLS
avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_SU_data);
tot_spikes = arrayfun(@(x) x.ModData.unit_data.tot_spikes,all_SU_data);
rec_dur = arrayfun(@(x) x.ModData.unit_data.N_used_samps,all_SU_data)*dt/60; %in min
expt_nums = [all_SU_data(:).expt_num];

clust_iso_dist = arrayfun(@(x) x.ModData.unit_data.SU_isodist,all_SU_data);
clust_Lratio = arrayfun(@(x) x.ModData.unit_data.SU_Lratio,all_SU_data);
clust_refract = arrayfun(@(x) x.ModData.unit_data.SU_refract,all_SU_data);
clust_dprime = arrayfun(@(x) x.ModData.unit_data.SU_dprime,all_SU_data);

xvLLimp = arrayfun(@(x) x.ModData.rectGQM.xvLLimp,all_SU_data);
LLimp = arrayfun(@(x) x.ModData.rectGQM.LLimp,all_SU_data);

min_rate = 5; %in Hz
min_rec_dur = 20; %in min
% min_rec_dur = 45; %in min
min_xvLLimp = 0.0; %in bits/spk
used_units = find(avg_rates >= min_rate & rec_dur >= min_rec_dur & xvLLimp >= min_xvLLimp);

use_SU_data = all_SU_data(used_units);
% use_SU_tdata = all_SU_tdata(used_units);
n_SUs = length(used_units);
jbe_SUs = find(strcmp('jbe',{use_SU_data(:).animal}));
lem_SUs = find(strcmp('lem',{use_SU_data(:).animal}));

lem_fov_expt_nums = [266 275];
lem_orig_expt_nums = [266 270 275 277];
fov_SUs = find([use_SU_data(:).expt_num] < 200 | ismember([use_SU_data(:).expt_num],lem_fov_expt_nums));
parafov_SUs = find([use_SU_data(:).expt_num] > 200 & ~ismember([use_SU_data(:).expt_num],lem_fov_expt_nums));
lem_fov_SUs = find(ismember([use_SU_data(:).expt_num],lem_fov_expt_nums));
lem_orig_SUs = find(ismember([use_SU_data(:).expt_num],lem_orig_expt_nums));
lem_new_SUs = find([use_SU_data(:).expt_num] > 200 & ~ismember([use_SU_data(:).expt_num],lem_orig_expt_nums));

%%
RF_eccs = arrayfun(@(x) x.ModData.tune_props.RF_ecc,use_SU_data);
RF_sigma = arrayfun(@(x) x.ModData.tune_props.RF_sigma,use_SU_data);
RF_gSF = arrayfun(@(x) x.ModData.tune_props.RF_gSF,use_SU_data);
RF_dirsel = arrayfun(@(x) x.ModData.tune_props.RF_dirsel,use_SU_data);
RF_FTF = arrayfun(@(x) x.ModData.tune_props.RF_FTF,use_SU_data);
PRM = arrayfun(@(x) x.ModData.tune_props.PRM,use_SU_data);
net_phase_polarity = arrayfun(@(x) x.ModData.tune_props.net_phase_polarity,use_SU_data);
RF_oris = [use_SU_data(:).bar_ori];


%% LAYER DEPENDENCE OF TEMPORAL KERNELS
load('/home/james/Analysis/bruce/saccade_modulation/layer_boundaries/layer_classification.mat')
boundary_enums = [boundary_class(:).Expt_num];

un_lem_expts = unique([use_SU_data(lem_SUs).expt_num]);
n_lem_expts = length(un_lem_expts);

all_gran_units = [];
all_supra_units = [];
all_infra_units = [];
for ee = 1:n_lem_expts
   cur_mua_set = lem_SUs([use_SU_data(lem_SUs).expt_num] == un_lem_expts(ee));
   cur_probe_nums = arrayfun(@(x) x.ModData.unit_data.probe_number,use_SU_data(cur_mua_set));
   cur_bound_info = find(boundary_enums == un_lem_expts(ee));
   cur_ub = boundary_class(cur_bound_info).ub;
   cur_lb = boundary_class(cur_bound_info).lb;
    
   gran_probes = (cur_ub+1):(cur_lb-1);
   supra_probes = 1:(cur_ub-1);
   infra_probes = (cur_lb+1):24;
   
   if ~isempty(gran_probes)
       all_gran_units = cat(2,all_gran_units,cur_mua_set(ismember(cur_probe_nums,gran_probes)));
   end
   if ~isempty(supra_probes)
       all_supra_units = cat(2,all_supra_units,cur_mua_set(ismember(cur_probe_nums,supra_probes)));
   end
   if ~isempty(infra_probes)
       all_infra_units = cat(2,all_infra_units,cur_mua_set(ismember(cur_probe_nums,infra_probes)));
   end
end

all_Ekerns = nan(length(use_SU_data),flen);
all_Ikerns = nan(length(use_SU_data),flen);
for ii = 1:length(use_SU_data)
    cur_mod = use_SU_data(ii).ModData.rectGQM;
    mod_signs = [cur_mod.mods(:).sign];
    sd = cur_mod.stim_params(1).stim_dims;
    mod_filts = reshape([cur_mod.mods(:).filtK],[sd(1) sd(2) length(mod_signs)]);
    tkerns = squeeze(std(mod_filts,[],2));
    if sum(mod_signs == 1) > 0
        all_Ekerns(ii,:) = mean(tkerns(:,mod_signs == 1),2);
    end
    if sum(mod_signs == -1) > 0
        all_Ikerns(ii,:) = mean(tkerns(:,mod_signs == -1),2);
    end
end

all_nEkerns = bsxfun(@rdivide,all_Ekerns,sqrt(sum(all_Ekerns.^2,2)));
all_nIkerns = bsxfun(@rdivide,all_Ikerns,sqrt(sum(all_Ikerns.^2,2)));
tax = (0:(flen-1))*dt + dt/2;

h1 = figure;
subplot(2,1,1); hold on
shadedErrorBar(tax,mean(all_nEkerns(all_gran_units,:)),std(all_nEkerns(all_gran_units,:))/sqrt(length(all_gran_units)),{'color','r'});
shadedErrorBar(tax,mean(all_nEkerns(all_infra_units,:)),std(all_nEkerns(all_infra_units,:))/sqrt(length(all_infra_units)),{'color','b'});
shadedErrorBar(tax,mean(all_nEkerns(all_supra_units,:)),std(all_nEkerns(all_supra_units,:))/sqrt(length(all_supra_units)),{'color','k'});
xlim([0 0.15]);
xlabel('Lag (s)');
ylabel('Filter temporal env');
title('E-filters');

subplot(2,1,2); hold on
shadedErrorBar(tax,mean(all_nIkerns(all_gran_units,:)),std(all_nIkerns(all_gran_units,:))/sqrt(length(all_gran_units)),{'color','r'});
shadedErrorBar(tax,mean(all_nIkerns(all_infra_units,:)),std(all_nIkerns(all_infra_units,:))/sqrt(length(all_infra_units)),{'color','b'});
shadedErrorBar(tax,mean(all_nIkerns(all_supra_units,:)),std(all_nIkerns(all_supra_units,:))/sqrt(length(all_supra_units)),{'color','k'});
xlim([0 0.15]);
xlabel('Lag (s)');
ylabel('Filter temporal env');
title('I-filters');

fig_width = 4; rel_height = 1.6;
figufy(h1);
fname = [fig_dir 'Layerdep_Tkerns.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

%%
close all

for ii = 1:n_SUs
    cur_GQM = use_SU_data(ii).ModData.rectGQM;
    stim_params = cur_GQM.stim_params.stim_dims;
    %compute best time lag over E filts
    [Xinds,Tinds] = meshgrid(1:stim_params(2),1:stim_params(1));
    cur_filts = reshape([cur_GQM.mods(1:end).filtK],[stim_params(1) stim_params(2) length(cur_GQM.mods)]);
    cur_tfilt = squeeze(mean(std(cur_filts,[],2),3));
    [~,best_lag] = max(cur_tfilt);
        
    cur_gsac_phaseDep = reshape(use_SU_data(ii).gsac_phaseDep_subfilt,[length(slags) stim_params(1) stim_params(2)]);
    cur_gsac_phaseInd = reshape(use_SU_data(ii).gsac_phaseInd_subfilt,[length(slags) stim_params(1) stim_params(2)]);
    
    cur_gsac_phaseDep = squeeze(cur_gsac_phaseDep(:,best_lag,:));
    cur_gsac_phaseInd = squeeze(cur_gsac_phaseInd(:,best_lag,:));
    
    subplot(2,2,1)
    imagesc(slags*dt,1:stim_params(2),cur_gsac_phaseDep');
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    subplot(2,2,3);
    imagesc(slags*dt,1:stim_params(2),cur_gsac_phaseInd');
    ca = caxis(); cam = max(abs(ca)); caxis([-cam cam]);
    subplot(2,2,2);
    plot(slags*dt,use_SU_data(ii).gsac_avg_rate/dt);
    subplot(2,2,4); hold on
    plot(slags*dt,use_SU_data(ii).gsac_TB_info);
    plot(slags*dt,use_SU_data(ii).gsac_sub_modinfo,'r');
     plot(slags*dt,use_SU_data(ii).gsac_spost_modinfo,'k');
   
    pause
    clf
end

%% COMPARE DIFFERENT INFO CALCS
avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,use_SU_data);
all_gsac_rates = reshape([use_SU_data(:).gsac_avg_rate],[],length(use_SU_data))';
all_gsac_nrates = bsxfun(@rdivide,all_gsac_rates,avg_rates*dt);

all_gsac_TBinfo = reshape([use_SU_data(:).gsac_TB_info],[],length(use_SU_data))';
all_gsac_TBinforate = all_gsac_TBinfo.*all_gsac_rates;
all_gsac_ov_TB_info = [use_SU_data(:).gsac_ov_TB_info];
all_gsac_NTBinfo = bsxfun(@rdivide,all_gsac_TBinfo,all_gsac_ov_TB_info');
all_gsac_NTBinforate = bsxfun(@rdivide,all_gsac_TBinforate,all_gsac_ov_TB_info'.*avg_rates*dt);

all_gsac_smodinfo = reshape([use_SU_data(:).gsac_sub_modinfo],[],length(use_SU_data))';
all_gsac_ov_smodinfo = [use_SU_data(:).gsac_sub_ov_modinfo];
all_gsac_smodinforate = all_gsac_smodinfo.*all_gsac_rates;
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');
all_gsac_Nsmodinforate = bsxfun(@rdivide,all_gsac_smodinforate,all_gsac_ov_smodinfo'.*avg_rates*dt);

all_gsac_modinfo = reshape([use_SU_data(:).gsac_spost_modinfo],[],length(use_SU_data))';
all_gsac_ov_modinfo = [use_SU_data(:).gsac_spost_ov_modinfo];
all_gsac_modinforate = all_gsac_modinfo.*all_gsac_rates;
all_gsac_Nmodinfo = bsxfun(@rdivide,all_gsac_modinfo,all_gsac_ov_modinfo');
all_gsac_Nmodinforate = bsxfun(@rdivide,all_gsac_modinforate,all_gsac_ov_modinfo'.*avg_rates*dt);

all_gsac_LLinfo = reshape([use_SU_data(:).gsac_spost_LLinfo],[],length(use_SU_data))';
all_gsac_ov_LLinfo = [use_SU_data(:).gsac_spost_ov_LLinfo];
all_gsac_NLLinfo = bsxfun(@rdivide,all_gsac_LLinfo,all_gsac_ov_LLinfo');


f1 = figure(); 
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_NTBinfo),std(all_gsac_NTBinfo)/sqrt(length(use_SU_data)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinfo),std(all_gsac_Nmodinfo)/sqrt(length(use_SU_data)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'TB','Mod-pred'});
xlabel('Time (s)');
ylabel('Relative info');

f2 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_NTBinforate),std(all_gsac_NTBinforate)/sqrt(length(use_SU_data)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinforate),std(all_gsac_Nmodinforate)/sqrt(length(use_SU_data)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'TB','Mod-pred'});
xlabel('Time (s)');
ylabel('Relative info rate');


f3 = figure(); 
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_NTBinfo),std(all_gsac_NTBinfo)/sqrt(length(use_SU_data)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinfo),std(all_gsac_Nmodinfo)/sqrt(length(use_SU_data)),{'color','r'});
h3=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo),std(all_gsac_Nsmodinfo)/sqrt(length(use_SU_data)),{'color','k'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'TB','Mod-pred','Submod-pred'});
xlabel('Time (s)');
ylabel('Relative info');


fig_width = 4; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'Gsac_infotype_compare.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

fig_width = 4; rel_height = 0.8;
figufy(f2);
fname = [fig_dir 'Gsac_inforate_type_compare.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

fig_width = 4; rel_height = 0.8;
figufy(f3);
fname = [fig_dir 'Gsac_infotype_fullcompare.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);

%% COMPARE TENT_BASIS INFO
close all
% plot_set1 = fov_SUs;
plot_set1 = lem_fov_SUs;
plot_set2 = parafov_SUs;
% plot_set1 = parafov_SUs;
% plot_set2 = fov_SUs;

avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,use_SU_data);
all_gsac_rates = reshape([use_SU_data(:).gsac_avg_rate],[],length(use_SU_data))';
all_gsac_nrates = bsxfun(@rdivide,all_gsac_rates,avg_rates*dt);

% all_gsac_trig_rates = reshape([use_SU_tdata(:).gsac_gray],[],length(use_SU_tdata))';

all_gsac_TBinfo = reshape([use_SU_data(:).gsac_TB_info],[],length(use_SU_data))';
all_gsac_TBinforate = all_gsac_TBinfo.*all_gsac_rates;

all_gsac_ov_TB_info = [use_SU_data(:).gsac_ov_TB_info];
all_gsac_NTBinfo = bsxfun(@rdivide,all_gsac_TBinfo,all_gsac_ov_TB_info');
all_gsac_NTBinforate = bsxfun(@rdivide,all_gsac_TBinforate,all_gsac_ov_TB_info'.*avg_rates*dt);

f1 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_NTBinfo(plot_set1,:)),std(all_gsac_NTBinfo(plot_set1,:))/sqrt(length(plot_set1)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set1,:)),std(all_gsac_nrates(plot_set1,:))/sqrt(length(plot_set1)),{'color','r'});
h3=shadedErrorBar(slags*dt,mean(all_gsac_NTBinfo(plot_set2,:)),std(all_gsac_NTBinfo(plot_set2,:))/sqrt(length(plot_set2)),{'color','g'});
h4=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set2,:)),std(all_gsac_nrates(plot_set2,:))/sqrt(length(plot_set2)),{'color','m'});
xl = xlim();
line(xl,[1 1],'color','k');

f2 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_NTBinforate(plot_set1,:)),std(all_gsac_NTBinforate(plot_set1,:))/sqrt(length(plot_set1)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set1,:)),std(all_gsac_nrates(plot_set1,:))/sqrt(length(plot_set1)),{'color','r'});
h3=shadedErrorBar(slags*dt,mean(all_gsac_NTBinforate(plot_set2,:)),std(all_gsac_NTBinforate(plot_set2,:))/sqrt(length(plot_set2)),{'color','g'});
h4=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set2,:)),std(all_gsac_nrates(plot_set2,:))/sqrt(length(plot_set2)),{'color','m'});
xl = xlim();
line(xl,[1 1],'color','k');

ca = [0.25 1.75];
ulags = find(slags*dt >= 0);
[gsac_peak_mods,gsac_peak_locs] = max(all_gsac_nrates(:,ulags),[],2);
[gsac_val_mods,gsac_val_locs] = min(all_gsac_nrates(:,ulags),[],2);
[~,peak_ord] = sort(gsac_peak_mods(plot_set1));
[~,val_ord] = sort(gsac_val_mods(plot_set1));
f3 = figure();
subplot(3,2,1);
imagesc(slags*dt,1:n_SUs,all_gsac_NTBinfo(plot_set1(peak_ord),:));
caxis(ca);
subplot(3,2,3);
imagesc(slags*dt,1:n_SUs,all_gsac_NTBinforate(plot_set1(peak_ord),:));
caxis(ca);
subplot(3,2,5);
imagesc(slags*dt,1:n_SUs,all_gsac_nrates(plot_set1(peak_ord),:));
caxis(ca);
subplot(3,2,2);
imagesc(slags*dt,1:n_SUs,all_gsac_NTBinfo(plot_set1(val_ord),:));
caxis(ca);
subplot(3,2,4);
imagesc(slags*dt,1:n_SUs,all_gsac_NTBinforate(plot_set1(val_ord),:));
caxis(ca);
subplot(3,2,6);
imagesc(slags*dt,1:n_SUs,all_gsac_nrates(plot_set1(val_ord),:));
caxis(ca);

%% Compare post-mod components
close all

plot_set1 = jbe_SUs;
plot_set2 = lem_SUs;

all_postgains = [];
all_postoffss = [];
for ii = 1:length(use_SU_data)
    all_postgains = cat(1,all_postgains,use_SU_data(ii).gsac_post_singmod.mods(3).filtK');
    all_postoffss = cat(1,all_postoffss,use_SU_data(ii).gsac_post_singmod.mods(2).filtK');
end
avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,use_SU_data);
all_gsac_rates = reshape([use_SU_data(:).gsac_avg_rate],[],length(use_SU_data))';
all_gsac_nrates = bsxfun(@rdivide,all_gsac_rates,avg_rates*dt);

f1 = figure; hold on
h1=shadedErrorBar(slags*dt,1+mean(all_postgains(jbe_SUs,:)),std(all_postgains(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,1+mean(all_postgains(lem_SUs,:)),std(all_postgains(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','b'});
xl = xlim();
line(xl,[1 1],'color','k');
ylim([0.65 1.2]);
xlabel('Time (s)');
ylabel('Gain');

f2 = figure; hold on
h1=shadedErrorBar(slags*dt,mean(all_postoffss(jbe_SUs,:)),std(all_postoffss(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_postoffss(lem_SUs,:)),std(all_postoffss(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','b'});
xl = xlim();
line(xl,[0 0],'color','k');
xlabel('Time (s)');
ylabel('Offset');


f3 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set1,:)),std(all_gsac_nrates(plot_set1,:))/sqrt(length(plot_set1)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set2,:)),std(all_gsac_nrates(plot_set2,:))/sqrt(length(plot_set2)),{'color','b'});
xl = xlim();
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');


% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_Nsmod_gain.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'Gsac_Nsmod_offset.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'Gsac_Navgrate.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

%% Compare post-mod components GSAC VS MSAC
close all

plot_set1 = jbe_SUs;
plot_set2 = lem_SUs;

all_postgains_gsac = [];
all_postoffss_gsac = [];
all_postgains_msac = [];
all_postoffss_msac = [];
for ii = 1:length(use_SU_data)
    all_postgains_gsac = cat(1,all_postgains_gsac,use_SU_data(ii).gsac_post_singmod.mods(3).filtK');
    all_postoffss_gsac = cat(1,all_postoffss_gsac,use_SU_data(ii).gsac_post_singmod.mods(2).filtK');
    all_postgains_msac = cat(1,all_postgains_msac,use_SU_data(ii).msac_post_singmod.mods(3).filtK');
    all_postoffss_msac = cat(1,all_postoffss_msac,use_SU_data(ii).msac_post_singmod.mods(2).filtK');
end
avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,use_SU_data);
all_gsac_rates = reshape([use_SU_data(:).gsac_avg_rate],[],length(use_SU_data))';
all_gsac_nrates = bsxfun(@rdivide,all_gsac_rates,avg_rates*dt);
all_msac_rates = reshape([use_SU_data(:).msac_avg_rate],[],length(use_SU_data))';
all_msac_nrates = bsxfun(@rdivide,all_msac_rates,avg_rates*dt);

f1 = figure; hold on
h1=shadedErrorBar(slags*dt,1+mean(all_postgains_gsac(jbe_SUs,:)),std(all_postgains_gsac(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,1+mean(all_postgains_gsac(lem_SUs,:)),std(all_postgains_gsac(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','b'});
h1=shadedErrorBar(slags*dt,1+mean(all_postgains_msac(jbe_SUs,:)),std(all_postgains_msac(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','k'});
h2=shadedErrorBar(slags*dt,1+mean(all_postgains_msac(lem_SUs,:)),std(all_postgains_msac(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','m'});
xl = xlim();
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Gain');

f2 = figure; hold on
h1=shadedErrorBar(slags*dt,mean(all_postoffss_gsac(jbe_SUs,:)),std(all_postoffss_gsac(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_postoffss_gsac(lem_SUs,:)),std(all_postoffss_gsac(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','b'});
h1=shadedErrorBar(slags*dt,mean(all_postoffss_msac(jbe_SUs,:)),std(all_postoffss_msac(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','k'});
h2=shadedErrorBar(slags*dt,mean(all_postoffss_msac(lem_SUs,:)),std(all_postoffss_msac(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','m'});
xl = xlim();
line(xl,[0 0],'color','k');
xlabel('Time (s)');
ylabel('Offset');


f3 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set1,:)),std(all_gsac_nrates(plot_set1,:))/sqrt(length(plot_set1)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set2,:)),std(all_gsac_nrates(plot_set2,:))/sqrt(length(plot_set2)),{'color','b'});
h1=shadedErrorBar(slags*dt,mean(all_msac_nrates(plot_set1,:)),std(all_msac_nrates(plot_set1,:))/sqrt(length(plot_set1)),{'color','k'});
h2=shadedErrorBar(slags*dt,mean(all_msac_nrates(plot_set2,:)),std(all_msac_nrates(plot_set2,:))/sqrt(length(plot_set2)),{'color','m'});
xl = xlim();
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

% 
% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_Msac_Nsmod_gain.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'Gsac_Msac_Nsmod_offset.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'Gsac_Msac_Navgrate.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

%% COMPARE GSAC MSAC MODEL INFO
close all

plot_set1 = jbe_SUs;
plot_set2 = lem_SUs;

avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,use_SU_data);
all_gsac_rates = reshape([use_SU_data(:).gsac_avg_rate],[],length(use_SU_data))';
all_msac_rates = reshape([use_SU_data(:).msac_avg_rate],[],length(use_SU_data))';

all_gsac_modinfo = reshape([use_SU_data(:).gsac_spost_modinfo],[],length(use_SU_data))';
all_gsac_ov_modinfo = [use_SU_data(:).gsac_spost_ov_modinfo];
all_msac_modinfo = reshape([use_SU_data(:).msac_spost_modinfo],[],length(use_SU_data))';
all_msac_ov_modinfo = [use_SU_data(:).msac_spost_ov_modinfo];

all_gsac_modinforate = all_gsac_modinfo.*all_gsac_rates;
all_msac_modinforate = all_msac_modinfo.*all_msac_rates;

all_gsac_Nmodinfo = bsxfun(@rdivide,all_gsac_modinfo,all_gsac_ov_modinfo');
all_gsac_Nmodinforate = bsxfun(@rdivide,all_gsac_modinforate,all_gsac_ov_modinfo'.*avg_rates*dt);
all_msac_Nmodinfo = bsxfun(@rdivide,all_msac_modinfo,all_msac_ov_modinfo');
all_msac_Nmodinforate = bsxfun(@rdivide,all_msac_modinforate,all_msac_ov_modinfo'.*avg_rates*dt);

f1 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinfo(plot_set1,:)),std(all_gsac_Nmodinfo(plot_set1,:))/sqrt(length(plot_set1)),{'color','r'});
h3=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinfo(plot_set2,:)),std(all_gsac_Nmodinfo(plot_set2,:))/sqrt(length(plot_set2)),{'color','b'});
h1=shadedErrorBar(slags*dt,mean(all_msac_Nmodinfo(plot_set1,:)),std(all_msac_Nmodinfo(plot_set1,:))/sqrt(length(plot_set1)),{'color','k'});
h3=shadedErrorBar(slags*dt,mean(all_msac_Nmodinfo(plot_set2,:)),std(all_msac_Nmodinfo(plot_set2,:))/sqrt(length(plot_set2)),{'color','m'});
xl = xlim();
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative single-spike info');
ylim([0.4 1.4]);

f2 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinforate(plot_set1,:)),std(all_gsac_Nmodinforate(plot_set1,:))/sqrt(length(plot_set1)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinforate(plot_set2,:)),std(all_gsac_Nmodinforate(plot_set2,:))/sqrt(length(plot_set2)),{'color','b'});
h1=shadedErrorBar(slags*dt,mean(all_msac_Nmodinforate(plot_set1,:)),std(all_msac_Nmodinforate(plot_set1,:))/sqrt(length(plot_set1)),{'color','k'});
h2=shadedErrorBar(slags*dt,mean(all_msac_Nmodinforate(plot_set2,:)),std(all_msac_Nmodinforate(plot_set2,:))/sqrt(length(plot_set2)),{'color','m'});
xl = xlim();
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative info rate');


% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_Msac_Nsmod_info.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'Gsac_Msac_Nsmod_inforate.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 

%% COMPARE MODEL INFO
close all
% % plot_set1 = fov_SUs;
plot_set1 = lem_fov_SUs;
plot_set2 = parafov_SUs;

% plot_set1 = jbe_SUs;
% plot_set2 = lem_SUs;
% plot_set1 = lem_orig_SUs;
% % plot_set2 = lem_new_SUs;
% plot_set2 = jbe_SUs;

avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,use_SU_data);
all_gsac_rates = reshape([use_SU_data(:).gsac_avg_rate],[],length(use_SU_data))';
all_gsac_nrates = bsxfun(@rdivide,all_gsac_rates,avg_rates*dt);

all_gsac_modinfo = reshape([use_SU_data(:).gsac_spost_modinfo],[],length(use_SU_data))';
all_gsac_ov_modinfo = [use_SU_data(:).gsac_spost_ov_modinfo];
% all_gsac_modinfo = reshape([use_SU_data(:).gsac_sub_modinfo],[],length(use_SU_data))';
% all_gsac_ov_modinfo = [use_SU_data(:).gsac_sub_ov_modinfo];

all_gsac_LLinfo = reshape([use_SU_data(:).gsac_spost_LLinfo],[],length(use_SU_data))';
all_gsac_ov_LLinfo = [use_SU_data(:).gsac_spost_ov_LLinfo];
all_gsac_NLLinfo = bsxfun(@rdivide,all_gsac_LLinfo,all_gsac_ov_LLinfo');

all_gsac_modinforate = all_gsac_modinfo.*all_gsac_rates;

all_gsac_Nmodinfo = bsxfun(@rdivide,all_gsac_modinfo,all_gsac_ov_modinfo');
all_gsac_Nmodinforate = bsxfun(@rdivide,all_gsac_modinforate,all_gsac_ov_modinfo'.*avg_rates*dt);

f1 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinfo(plot_set1,:)),std(all_gsac_Nmodinfo(plot_set1,:))/sqrt(length(plot_set1)),{'color','r'});
h3=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinfo(plot_set2,:)),std(all_gsac_Nmodinfo(plot_set2,:))/sqrt(length(plot_set2)),{'color','b'});
xl = xlim();
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative single-spike info');
ylim([0.4 1.4]);

f2 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinforate(plot_set1,:)),std(all_gsac_Nmodinforate(plot_set1,:))/sqrt(length(plot_set1)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinforate(plot_set2,:)),std(all_gsac_Nmodinforate(plot_set2,:))/sqrt(length(plot_set2)),{'color','b'});
xl = xlim();
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative info rate');

f3 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set1,:)),std(all_gsac_nrates(plot_set1,:))/sqrt(length(plot_set1)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_nrates(plot_set2,:)),std(all_gsac_nrates(plot_set2,:))/sqrt(length(plot_set2)),{'color','b'});
xl = xlim();
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

ca = [0.25 1.75];
ulags = find(slags*dt >= 0);
[gsac_peak_mods,gsac_peak_locs] = max(all_gsac_nrates(:,ulags),[],2);
[gsac_val_mods,gsac_val_locs] = min(all_gsac_nrates(:,ulags),[],2);
[~,peak_ord] = sort(gsac_peak_mods);
[~,val_ord] = sort(gsac_val_mods);
f4 = figure();
subplot(3,2,1);
imagesc(slags*dt,1:n_SUs,all_gsac_Nmodinfo(peak_ord,:));
caxis(ca);
subplot(3,2,3);
imagesc(slags*dt,1:n_SUs,all_gsac_Nmodinforate(peak_ord,:));
caxis(ca);
subplot(3,2,5);
imagesc(slags*dt,1:n_SUs,all_gsac_nrates(peak_ord,:));
caxis(ca);
subplot(3,2,2);
imagesc(slags*dt,1:n_SUs,all_gsac_Nmodinfo(val_ord,:));
caxis(ca);
subplot(3,2,4);
imagesc(slags*dt,1:n_SUs,all_gsac_Nmodinforate(val_ord,:));
caxis(ca);
subplot(3,2,6);
imagesc(slags*dt,1:n_SUs,all_gsac_nrates(val_ord,:));
caxis(ca);

% [~,info_ord] = sort(xvLLimp(used_units));
% [~,info_ord] = sort(clust_iso_dist(used_units));
[~,info_ord] = sort(rec_dur(used_units));
f5 = figure();
subplot(3,1,1);
imagesc(slags*dt,1:n_SUs,all_gsac_Nmodinfo(info_ord,:));
caxis(ca);
subplot(3,1,2);
imagesc(slags*dt,1:n_SUs,all_gsac_Nmodinforate(info_ord,:));
caxis(ca);
subplot(3,1,3);
imagesc(slags*dt,1:n_SUs,all_gsac_nrates(info_ord,:));
caxis(ca);


% fig_width = 4; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_Nsmod_info_longrecs.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'Gsac_Nsmod_inforate_longrecs.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'Gsac_Navgrate.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

% fig_width = 4*2; rel_height = 0.8*3/2;
% figufy(f4);
% fname = [fig_dir 'Gsac_sorted_smodinfo.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% compare pre-post gains
all_gsac_gaink = reshape(cell2mat(arrayfun(@(x) x.gsacGainMod.gain_kernel,use_SU_data,'uniformoutput',0)),[],length(use_SU_data))';
all_gsac_stimk = reshape(cell2mat(arrayfun(@(x) x.gsacGainMod.stim_kernel,use_SU_data,'uniformoutput',0)),[],length(use_SU_data))';
all_gsac_offk = reshape(cell2mat(arrayfun(@(x) x.gsacGainMod.off_kernel,use_SU_data,'uniformoutput',0)),[],length(use_SU_data))';

avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,use_SU_data);
all_gsac_rates = reshape([use_SU_data(:).gsac_avg_rate],[],length(use_SU_data))';
all_gsac_nrates = bsxfun(@rdivide,all_gsac_rates,avg_rates*dt);

all_postgains = [];
for ii = 1:length(use_SU_data)
    all_postgains = cat(1,all_postgains,use_SU_data(ii).gsac_post_singmod.mods(3).filtK');
end

h = figure; 
subplot(2,1,1); hold on
h1=shadedErrorBar(slags*dt,1+mean(all_gsac_gaink(jbe_SUs,:)),std(all_gsac_gaink(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_gsac_stimk(jbe_SUs,:)),std(all_gsac_stimk(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','k'});
h3=shadedErrorBar(slags*dt,1+mean(all_gsac_gaink(lem_SUs,:)),std(all_gsac_gaink(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','g'});
h4=shadedErrorBar(slags*dt,1+mean(all_gsac_stimk(lem_SUs,:)),std(all_gsac_stimk(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','r'});
legend([h1.mainLine h2.mainLine],{'Pre-filt','Post-filt'},'Location','southeast');
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Gains');

subplot(2,1,2); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_nrates(jbe_SUs,:)),std(all_gsac_nrates(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_nrates(lem_SUs,:)),std(all_gsac_nrates(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','b'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

f2 = figure(); hold on
h1=shadedErrorBar(slags*dt,1+mean(all_gsac_gaink),std(all_gsac_gaink)/sqrt(length(use_SU_data)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_gsac_stimk),std(all_gsac_stimk)/sqrt(length(use_SU_data)),{'color','k'});
h3=shadedErrorBar(slags*dt,1+mean(all_postgains),std(all_postgains)/sqrt(length(use_SU_data)),{'color','r'});
xlabel('Time since sac onset (s)');
ylabel('Gains');

% f3 = figure();
% h1=shadedErrorBar(slags*dt,mean(all_gsac_offk),std(all_gsac_offk)/sqrt(length(use_SU_data)),{'color','b'});


% fig_width = 4; rel_height = 1.6;
% figufy(h);
% fname = [fig_dir 'Gsac_prepost_gains.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

% fig_width = 4; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Gsac_prepost_gains_beforeafter.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% E-I kernels
flen = 15;
all_Ekerns = nan(length(use_SU_data),flen);
all_Ikerns = nan(length(use_SU_data),flen);
for ii = 1:length(use_SU_data)
    cur_mod = use_SU_data(ii).ModData.rectGQM;
    mod_signs = [cur_mod.mods(:).sign];
    sd = cur_mod.stim_params(1).stim_dims;
    mod_filts = reshape([cur_mod.mods(:).filtK],[sd(1) sd(2) length(mod_signs)]);
    tkerns = squeeze(std(mod_filts,[],2));
    if sum(mod_signs == 1) > 0
        all_Ekerns(ii,:) = mean(tkerns(:,mod_signs == 1),2);
    end
    if sum(mod_signs == -1) > 0
        all_Ikerns(ii,:) = mean(tkerns(:,mod_signs == -1),2);
    end
end

all_nEkerns = bsxfun(@rdivide,all_Ekerns,sqrt(sum(all_Ekerns.^2,2)));
all_nIkerns = bsxfun(@rdivide,all_Ikerns,sqrt(sum(all_Ikerns.^2,2)));

all_gsac_Egains = reshape([use_SU_data(:).gsac_post_Egains],[],length(use_SU_data))';
all_gsac_Igains = reshape([use_SU_data(:).gsac_post_Igains],[],length(use_SU_data))';


h = figure; 
subplot(3,1,1);hold on
h1=shadedErrorBar(slags*dt,1+mean(all_gsac_Egains(jbe_SUs,:)),std(all_gsac_Egains(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_gsac_Igains(jbe_SUs,:)),std(all_gsac_Igains(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','r'});
h3=shadedErrorBar(slags*dt,1+mean(all_gsac_Egains(lem_SUs,:)),std(all_gsac_Egains(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','k'});
h4=shadedErrorBar(slags*dt,1+mean(all_gsac_Igains(lem_SUs,:)),std(all_gsac_Igains(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','g'});
xl = xlim(); yl = ylim();
legend([h1.mainLine h2.mainLine],{'E gain','I gain'});
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Gain');

subplot(3,1,2); hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_nrates(jbe_SUs,:)),std(all_gsac_nrates(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_nrates(lem_SUs,:)),std(all_gsac_nrates(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','b'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('Relative rate');

tax = (0:(flen-1))*dt + dt/2;
subplot(3,1,3); hold on
h1=shadedErrorBar(tax,mean(all_nEkerns(jbe_SUs,:)),std(all_nEkerns(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','b'});
h2=shadedErrorBar(tax,mean(all_nIkerns(jbe_SUs,:)),std(all_nIkerns(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','r'});
h3=shadedErrorBar(tax,mean(all_nEkerns(lem_SUs,:)),std(all_nEkerns(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','k'});
h4=shadedErrorBar(tax,mean(all_nIkerns(lem_SUs,:)),std(all_nIkerns(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','g'});
legend([h1.mainLine h2.mainLine],{'E kern','I kern'});
xlabel('Lag (s)');
ylabel('Temporal kernel');

% fig_width = 4; rel_height = 2.2;
% figufy(h);
% fname = [fig_dir 'Gsac_EI_gains.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);


avg_gsac_nrate = mean(all_gsac_nrates);
[~,minloc] = min(avg_gsac_nrate);
[~,maxloc] = max(avg_gsac_nrate);
all_EI_ratio = (1+all_gsac_Egains)./(1+all_gsac_Igains);

h = figure; 
subplot(2,1,1);hold on
h1=shadedErrorBar(slags*dt,mean(all_EI_ratio(jbe_SUs,:)),std(all_EI_ratio(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_EI_ratio(lem_SUs,:)),std(all_EI_ratio(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','r'});
xl = xlim(); yl = ylim();
line(xl,[1 1],'color','k','linestyle','--');
line([0 0],yl,'color','k','linestyle','--');
line(slags([minloc minloc])*dt,yl,'color','r','linestyle','--');
line(slags([maxloc maxloc])*dt,yl,'color','g','linestyle','--');
xlabel('Time since sac onset (s)');
ylabel('E/I ratio');

subplot(2,1,2);hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_nrates(jbe_SUs,:)),std(all_gsac_nrates(jbe_SUs,:))/sqrt(length(jbe_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_nrates(lem_SUs,:)),std(all_gsac_nrates(lem_SUs,:))/sqrt(length(lem_SUs)),{'color','r'});
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


%% Compare preferred and non-preferred post-mod components
close all

all_pref_postgains = [];
all_prepf_postoffss = [];
all_npref_postgains = [];
all_nprepf_postoffss = [];
for ii = 1:length(all_jbe_pref_data)
    all_pref_postgains = cat(1,all_pref_postgains,all_jbe_pref_data(ii).gsac_post_singmod.mods(3).filtK');
    all_prepf_postoffss = cat(1,all_prepf_postoffss,all_jbe_pref_data(ii).gsac_post_singmod.mods(2).filtK');
    all_npref_postgains = cat(1,all_npref_postgains,all_jbe_nonpref_data(ii).gsac_post_singmod.mods(3).filtK');
    all_nprepf_postoffss = cat(1,all_nprepf_postoffss,all_jbe_nonpref_data(ii).gsac_post_singmod.mods(2).filtK');
end
avg_pref_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_jbe_pref_data);
avg_npref_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_jbe_nonpref_data);

all_pref_gsac_rates = reshape([all_jbe_pref_data(:).gsac_avg_rate],[],length(all_jbe_pref_data))';
all_pref_gsac_nrates = bsxfun(@rdivide,all_pref_gsac_rates,avg_pref_rates*dt);
all_npref_gsac_rates = reshape([all_jbe_nonpref_data(:).gsac_avg_rate],[],length(all_jbe_pref_data))';
all_npref_gsac_nrates = bsxfun(@rdivide,all_npref_gsac_rates,avg_npref_rates*dt);

f1 = figure; hold on
h1=shadedErrorBar(slags*dt,1+mean(all_pref_postgains),std(all_pref_postgains)/sqrt(length(all_pref_postgains)),{'color','r'});
h2=shadedErrorBar(slags*dt,1+mean(all_npref_postgains),std(all_npref_postgains)/sqrt(length(all_npref_postgains)),{'color','b'});
xl = xlim();
line(xl,[1 1],'color','k');
ylim([0.65 1.2]);
xlabel('Time (s)');
ylabel('Gain');

f2 = figure; hold on
h1=shadedErrorBar(slags*dt,mean(all_prepf_postoffss),std(all_prepf_postoffss)/sqrt(length(all_pref_postgains)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_nprepf_postoffss),std(all_nprepf_postoffss)/sqrt(length(all_npref_postgains)),{'color','b'});
xl = xlim();
line(xl,[0 0],'color','k');
xlabel('Time (s)');
ylabel('Offset');


f3 = figure(); hold on
h1=shadedErrorBar(slags*dt,mean(all_pref_gsac_nrates),std(all_pref_gsac_nrates)/sqrt(length(all_pref_postgains)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_npref_gsac_nrates),std(all_npref_gsac_nrates)/sqrt(length(all_pref_postgains)),{'color','b'});
xl = xlim();
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');


%%
% ulags = find(slags*dt >= 0);
% [gsac_peak_mods,gsac_peak_locs] = max(all_gsac_nrates(:,ulags),[],2);
% [gsac_val_mods,gsac_val_locs] = min(all_gsac_nrates(:,ulags),[],2);
% 
% [~,gsac_peak_ord] = sort(gsac_peak_mods(used_units));
% [~,gsac_val_ord] = sort(1-gsac_val_mods(used_units));
% 
% h = figure();
% subplot(1,3,1)
% imagesc(slags*dt,1:length(used_units),all_gsac_nrates(used_units(gsac_peak_ord),:));
% caxis([0.25 1.75])
% ylabel('Unit number');
% xlabel('Time (s)');
% title('Relative rates');
% 
% subplot(1,3,2)
% imagesc(slags*dt,1:length(used_units),all_gsac_rel_TBinfo(used_units(gsac_peak_ord),:));
% caxis([0 2])
% ylabel('Unit number');
% xlabel('Time (s)');
% title('SS Info');
% 
% subplot(1,3,3)
% imagesc(slags*dt,1:length(used_units),all_gsac_rel_TBinforate(used_units(gsac_peak_ord),:));
% caxis([0 2])
% ylabel('Unit number');
% xlabel('Time (s)');
% title('Info rate');
% 
% % fig_width = 8; rel_height = 0.3;
% % figufy(h);
% % fname = [fig_dir 'Gsac_enhance_sort.pdf'];
% % exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h);
% 
% % [~,ord] = sort(all_ov_info_rate(used_units));
% % % [~,ord] = sort(all_ov_info_rate(used_units));
% % % [~,ord] = sort(clust_iso_dist(used_units));
% % figure;
% % imagesc(slags*dt,1:length(used_units),1+all_gsac_stimk(used_units(ord),:));
% 
% %%
% gsac_peak_times = slags(ulags(gsac_peak_locs))*dt;
% gsac_val_times = slags(ulags(gsac_val_locs))*dt;
% good_peaks = find(~ismember(gsac_peak_times,slags(ulags([1 end]))*dt));
% good_vals = find(~ismember(gsac_val_times,slags(ulags([1 end]))*dt));
% good_peaks(~ismember(good_peaks,used_units)) = [];
% good_vals(~ismember(good_vals,used_units)) = [];
% 
% tax = (0:(flen-1))*dt + dt/2;
% ekern_temp_loc = tax(ekern_loc);
% 
% jit_amp = 0.001;
% jit_vals = randn(length(gsac_peak_times),2)*jit_amp;
% 
% h = figure; 
% subplot(2,1,1);hold on
% plot(ekern_temp_loc(good_vals)'+jit_vals(good_vals,1),gsac_val_times(good_vals)'+jit_vals(good_vals,2),'r.','markersize',10);
% xlabel('Stimulus filter delay (s)');
% ylabel('Sac-suppression timing (s)');
% title('Suppression');
% [rV,statsV] = robustfit(ekern_temp_loc(good_vals),gsac_val_times(good_vals));
% xx = linspace(0.03,0.08,50);
% plot(xx,rV(1)+rV(2)*xx,'k');
% 
% subplot(2,1,2);hold on
% plot(ekern_temp_loc(good_peaks)'+jit_vals(good_peaks,1),gsac_peak_times(good_peaks)'+jit_vals(good_peaks,2),'.','markersize',10);
% title('Enhancement');
% xlabel('Stimulus filter delay (s)');
% ylabel('Sac-enhancement timing (s)');
% [rP,statsP] = robustfit(ekern_temp_loc(good_peaks),gsac_peak_times(good_peaks));
% xx = linspace(0.03,0.08,50);
% plot(xx,rP(1)+rP(2)*xx,'k');
% 
% % fig_width = 4; rel_height = 1.6;
% % figufy(h);
% % fname = [fig_dir 'Gsac_filt_sacmod_timing.pdf'];
% % exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h);
% 

%% information compare
% h = figure; 
% subplot(2,1,1); hold on
% % h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_TBinfo(used_units,:)),std(all_gsac_rel_TBinfo(used_units,:))/sqrt(length(used_units)),{'color','b'});
% % h2=shadedErrorBar(slags*dt,mean(all_gsac_rel_TBinforate(used_units,:)),std(all_gsac_rel_TBinforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
% h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_smod_info(used_units,:)),std(all_gsac_rel_smod_info(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_gsac_rel_smod_inforate(used_units,:)),std(all_gsac_rel_smod_inforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
% % h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_submod_info(used_units,:)),std(all_gsac_rel_submod_info(used_units,:))/sqrt(length(used_units)),{'color','b'});
% % h2=shadedErrorBar(slags*dt,mean(all_gsac_rel_submod_inforate(used_units,:)),std(all_gsac_rel_submod_inforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
% xl = xlim(); yl = ylim();
% legend([h1.mainLine h2.mainLine],{'SS Info','Info rate'},'Location','Southeast');
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative information');
% 
% subplot(2,1,2)
% shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% 
% % fig_width = 4; rel_height = 1.6;
% % figufy(h);
% % fname = [fig_dir 'Gsac_info.pdf'];
% % exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h);
% 
% %% Post gain 
% h = figure; 
% subplot(2,1,1);hold on
% h1=shadedErrorBar(slags*dt,1+mean(all_spmod_off(used_units,:)),std(all_spmod_off(used_units,:))/sqrt(length(used_units)),{'color','r'});
% h2=shadedErrorBar(slags*dt,1+mean(all_spmod_gain(used_units,:)),std(all_spmod_gain(used_units,:))/sqrt(length(used_units)),{'color','b'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% 
% subplot(2,1,2);hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_minfo(used_units,:)),std(all_gsac_rel_minfo(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_gsac_rel_minforate(used_units,:)),std(all_gsac_rel_minforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% 
% fig_width = 4; rel_height = 1.6;
% figufy(h);
% fname = [fig_dir 'Gsac_post_gain.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);
% 
% %% E-gain vs I-gain
% h = figure; 
% subplot(3,1,1);hold on
% % h1=shadedErrorBar(slags*dt,1+mean(all_pmod_Ekern(used_units,:)),std(all_pmod_Ekern(used_units,:))/sqrt(length(used_units)),{'color','b'});
% % h2=shadedErrorBar(slags*dt,1+mean(all_pmod_Ikern(used_units,:)),std(all_pmod_Ikern(used_units,:))/sqrt(length(used_units)),{'color','r'});
% h1=shadedErrorBar(slags*dt,1+mean(all_pmod_Ekern(usable_ekern,:)),std(all_pmod_Ekern(usable_ekern,:))/sqrt(length(usable_ekern)),{'color','b'});
% h2=shadedErrorBar(slags*dt,1+mean(all_pmod_Ikern(usable_ikern,:)),std(all_pmod_Ikern(usable_ikern,:))/sqrt(length(usable_ikern)),{'color','r'});
% xl = xlim(); yl = ylim();
% legend([h1.mainLine h2.mainLine],{'E gain','I gain'});
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Gain');
% 
% subplot(3,1,2)
% shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% 
% tax = (0:(flen-1))*dt + dt/2;
% subplot(3,1,3); hold on
% % h1=shadedErrorBar(tax,mean(all_Ekerns(used_units,:)),std(all_Ekerns(used_units,:))/sqrt(length(used_units)),{'color','b'});
% % h2=shadedErrorBar(tax,mean(all_Ikerns(used_units,:)),std(all_Ikerns(used_units,:))/sqrt(length(used_units)),{'color','r'});
% h1=shadedErrorBar(tax,mean(all_Ekerns(usable_ekern,:)),std(all_Ekerns(usable_ekern,:))/sqrt(length(usable_ekern)),{'color','b'});
% h2=shadedErrorBar(tax,mean(all_Ikerns(usable_ikern,:)),std(all_Ikerns(usable_ikern,:))/sqrt(length(usable_ikern)),{'color','r'});
% legend([h1.mainLine h2.mainLine],{'E kern','I kern'});
% xlabel('Lag (s)');
% ylabel('Temporal kernel');
% 
% % fig_width = 4; rel_height = 2.2;
% % figufy(h);
% % fname = [fig_dir 'Gsac_EI_gains.pdf'];
% % exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h);
% 
% %% EI ratio 
% avg_gsac_nrate = mean(all_gsac_nrates(used_units,:));
% [~,minloc] = min(avg_gsac_nrate);
% [~,maxloc] = max(avg_gsac_nrate);
% all_EI_ratio = (1+all_pmod_Ekern)./(1+all_pmod_Ikern);
% 
% h = figure; 
% subplot(2,1,1);hold on
% % h1=shadedErrorBar(slags*dt,mean(all_EI_ratio(used_units,:)),std(all_EI_ratio(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h1=shadedErrorBar(slags*dt,mean(all_EI_ratio(usable_bothkerns,:)),std(all_EI_ratio(usable_bothkerns,:))/sqrt(length(usable_bothkerns)),{'color','b'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% line(slags([minloc minloc])*dt,yl,'color','r','linestyle','--');
% line(slags([maxloc maxloc])*dt,yl,'color','g','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('E/I ratio');
% 
% subplot(2,1,2);hold on
% % shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
% shadedErrorBar(slags*dt,mean(all_gsac_nrates(usable_bothkerns,:)),std(all_gsac_nrates(usable_bothkerns,:))/sqrt(length(usable_bothkerns)),{'color','r'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% line(slags([minloc minloc])*dt,yl,'color','r','linestyle','--');
% line(slags([maxloc maxloc])*dt,yl,'color','g','linestyle','--');
% 
% % fig_width = 4; rel_height = 1.6;
% % figufy(h);
% % fname = [fig_dir 'Gsac_EI_ratio.pdf'];
% % exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(h);
% 
% %%
% 
% all_msac_rates = reshape([all_jbe_data(:).msac_avg_rate],[],length(all_jbe_data))';
% all_msac_nrates = bsxfun(@rdivide,all_msac_rates,avg_rates'*dt);
% 
% % all_msac_minfo = reshape([all_jbe_data(:).msac_mod_info],[],length(all_jbe_data))';
% all_msac_minfo = reshape([all_jbe_data(:).msac_postmod_info],[],length(all_jbe_data))';
% all_msac_TBinforate = reshape([all_jbe_data(:).msac_TB_info],[],length(all_jbe_data))';
% all_msac_TBinfo = all_msac_TBinforate./all_msac_rates;
% 
% all_msac_minforate = all_msac_minfo.*all_msac_rates;
% 
% all_msac_rel_minfo = bsxfun(@rdivide,all_msac_minfo,all_ov_info');
% all_msac_rel_TBinfo = bsxfun(@rdivide,all_msac_TBinfo,all_ov_info');
% all_msac_rel_minforate = bsxfun(@rdivide,all_msac_minforate,all_ov_info'.*avg_rates'*dt);
% all_msac_rel_TBinforate = bsxfun(@rdivide,all_msac_TBinforate,all_ov_info'.*avg_rates'*dt);
% 
% % all_msac_gaink = cell2mat(arrayfun(@(x) x.msacGainMod.gain_kernel,all_jbe_data,'uniformoutput',0))';
% % all_msac_stimk = cell2mat(arrayfun(@(x) x.msacGainMod.stim_kernel,all_jbe_data,'uniformoutput',0))';
% % all_msac_offk = cell2mat(arrayfun(@(x) x.msacGainMod.off_kernel,all_jbe_data,'uniformoutput',0))';
% 
% 
% %%
% figure; 
% subplot(2,1,1); hold on
% shadedErrorBar(slags*dt,mean(all_msac_offk(used_units,:)),std(all_msac_offk(used_units,:))/sqrt(length(used_units)),{'color','r'});
% shadedErrorBar(slags*dt,mean(all_msac_gaink(used_units,:)),std(all_msac_gaink(used_units,:))/sqrt(length(used_units)),{'color','b'});
% shadedErrorBar(slags*dt,mean(all_msac_stimk(used_units,:)),std(all_msac_stimk(used_units,:))/sqrt(length(used_units)),{'color','k'});
% xl = xlim(); yl = ylim();
% line(xl,[0 0],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Filter amp');
% 
% subplot(2,1,2)
% shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% 
% %%
% figure; 
% subplot(2,1,1); hold on
% h1=shadedErrorBar(slags*dt,mean(all_msac_rel_TBinfo(used_units,:)),std(all_msac_rel_TBinfo(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_msac_rel_TBinforate(used_units,:)),std(all_msac_rel_TBinforate(used_units,:))/sqrt(length(used_units)),{'color','k'});
% % h3=shadedErrorBar(slags*dt,mean(all_msac_rel_minfo(used_units,:)),std(all_msac_rel_minfo(used_units,:))/sqrt(length(used_units)),{'color','r'});
% % h4=shadedErrorBar(slags*dt,mean(all_msac_rel_minforate(used_units,:)),std(all_msac_rel_minforate(used_units,:))/sqrt(length(used_units)),{'color','g'});
% xl = xlim(); yl = ylim();
% legend([h1.mainLine h2.mainLine],{'SS Info','Info rate'});
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative information');
% 
% subplot(2,1,2)
% shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% 
% %%
% avg_msac_nrate = mean(all_msac_nrates(used_units,:));
% [~,minloc] = min(avg_msac_nrate);
% [~,maxloc] = max(avg_msac_nrate);
% all_EI_ratio = (1+all_pMmod_Ekern)./(1+all_pMmod_Ikern);
% figure; 
% subplot(2,1,1);hold on
% h1=shadedErrorBar(slags*dt,mean(all_EI_ratio(used_units,:)),std(all_EI_ratio(used_units,:))/sqrt(length(used_units)),{'color','b'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% line(slags([minloc minloc])*dt,yl,'color','r','linestyle','--');
% line(slags([maxloc maxloc])*dt,yl,'color','g','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('E/I ratio');
% 
% subplot(2,1,2);hold on
% shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% line(slags([minloc minloc])*dt,yl,'color','r','linestyle','--');
% line(slags([maxloc maxloc])*dt,yl,'color','g','linestyle','--');
% 
% %%
% figure; 
% subplot(2,1,1);hold on
% h1=shadedErrorBar(slags*dt,1+mean(all_pmod_Ekern(used_units,:)),std(all_pmod_Ekern(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(slags*dt,1+mean(all_pmod_Ikern(used_units,:)),std(all_pmod_Ikern(used_units,:))/sqrt(length(used_units)),{'color','r'});
% xl = xlim(); yl = ylim();
% legend([h1.mainLine h2.mainLine],{'E gain','I gain'});
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Gain');
% 
% subplot(2,1,2)
% shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% 
% 
% %%
% figure;
% hold on
% shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','b'});
% shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% 
% % figure;
% % hold on
% % shadedErrorBar(slags*dt,mean(all_gsac_offk(used_units,:)),std(all_gsac_offk(used_units,:))/sqrt(length(used_units)),{'color','b'});
% % shadedErrorBar(slags*dt,mean(all_msac_offk(used_units,:)),std(all_msac_offk(used_units,:))/sqrt(length(used_units)),{'color','r'});
% 
% %%
% h = figure; 
% subplot(3,1,1); hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_TBinfo(used_units,:)),std(all_gsac_rel_TBinfo(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_msac_rel_TBinfo(used_units,:)),std(all_msac_rel_TBinfo(used_units,:))/sqrt(length(used_units)),{'color','r'});
% legend([h1.mainLine h2.mainLine],{'Big sac','Micro sac'},'Location','Southeast');
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('SS info');
% 
% subplot(3,1,2); hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_rel_TBinforate(used_units,:)),std(all_gsac_rel_TBinforate(used_units,:))/sqrt(length(used_units)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_msac_rel_TBinforate(used_units,:)),std(all_msac_rel_TBinforate(used_units,:))/sqrt(length(used_units)),{'color','r'});
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Info rate');
% 
% subplot(3,1,3); hold on
% shadedErrorBar(slags*dt,mean(all_gsac_nrates(used_units,:)),std(all_gsac_nrates(used_units,:))/sqrt(length(used_units)),{'color','b'});
% shadedErrorBar(slags*dt,mean(all_msac_nrates(used_units,:)),std(all_msac_nrates(used_units,:))/sqrt(length(used_units)),{'color','r'});
% ylim([0.6 1.4])
% xl = xlim(); yl = ylim();
% line(xl,[1 1],'color','k','linestyle','--');
% line([0 0],yl,'color','k','linestyle','--');
% xlabel('Time since sac onset (s)');
% ylabel('Relative rate');
% 
% fig_width = 4; rel_height = 2.4;
% figufy(h);
% fname = [fig_dir 'Gsac_Msac_sacinfo.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);
% 
% 
% %%
% % tax = (0:(flen-1))*dt + dt/2;
% % figure; hold on
% % shadedErrorBar(tax,mean(all_Ekerns(used_units,:)),std(all_Ekerns(used_units,:))/sqrt(length(used_units)),{'color','b'});
% % shadedErrorBar(tax,mean(all_Ikerns(used_units,:)),std(all_Ikerns(used_units,:))/sqrt(length(used_units)),{'color','r'});
% 
% %%
% % figure; hold on
% % shadedErrorBar(slags*dt,mean(all_msac_rel_minfo(used_units,:)),std(all_msac_rel_minfo(used_units,:))/sqrt(length(used_units)),{'color','r'});
% % shadedErrorBar(slags*dt,mean(all_msac_rel_minforate(used_units,:)),std(all_msac_rel_minforate(used_units,:))/sqrt(length(used_units)),{'color','g'});

