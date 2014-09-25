% close all
clear all
clc

fit_unCor = 0;
include_bursts = 0;

fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';
base_tname = 'sac_trig_avg_data';
base_sname = 'sacStimProcFinR';
base_yname = 'sacTypeDep';

if include_bursts
    base_tname = strcat(base_tname,'_withbursts');
    base_sname = strcat(base_sname,'_withbursts');
    base_yname = strcat(base_yname,'_withbursts');
end

all_SU_data = [];
all_SU_NPdata = [];

%% LOAD JBE
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
n_probes = 96;
ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
rmfield_list = {};

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    
    clear ori_data
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            %load trig avg data
            tname = strcat(sac_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            %load stimulus proc data
            sname = strcat(sac_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            if fit_unCor
                sname = strcat(sname,'_unCor');
            end
            load(sname);
            %load type-dep data (STILL NEED TO ADD THIS)
            yname = strcat(sac_dir,base_yname,sprintf('_ori%d',ori_list(ee,ii)));
            if fit_unCor
                yname = strcat(yname,'_unCor');
            end
            load(yname);
            
            su_range = (n_probes+1):length(sacStimProc);
            clear SU_data
            for jj = 1:length(su_range)
                SU_data(jj).sacStimProc = sacStimProc(su_range(jj));
                SU_data(jj).trig_avg = sua_data(jj);
                SU_data(jj).type_dep = sacTypeDep(su_range(jj));
            end
            [SU_data.expt_num] = deal(Expt_num);
            [SU_data.bar_ori] = deal(ori_list(ee,ii));
            [SU_data.animal] = deal('jbe');
            
            ori_data(ii,:) = SU_data;
        end
    end
    
    [cur_Noris,cur_Nunits] = size(ori_data);
    used = arrayfun(@(x) x.sacStimProc.used,ori_data);
    xvLLimps = nan(cur_Noris,cur_Nunits);
    xvLLimps(used) = arrayfun(@(x)x.sacStimProc.ModData.rectGQM.xvLLimp,ori_data(used));
    avg_rates = nan(cur_Noris,cur_Nunits);
    avg_rates(used) = arrayfun(@(x)x.sacStimProc.ModData.unit_data.avg_rate,ori_data(used));
    xvLLrates = xvLLimps.*avg_rates;
    
    [mvals,mlocs] = max(xvLLrates,[],1);
    nplocs = mod(mlocs,2)+1; %these are indices for the non-pref ori
    for ii = 1:cur_Nunits
        if ~isnan(xvLLrates(mlocs(ii),ii))
            all_SU_data = cat(1,all_SU_data,ori_data(mlocs(ii),ii));
        end
    end
end

%% LOAD LEM
% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
Expt_list = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};%NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
n_probes = 24;
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
rmfield_list = {};

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    
    clear ori_data
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            %load trig avg data
            tname = strcat(sac_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            %load stimulus proc data
            sname = strcat(sac_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            if fit_unCor
                sname = strcat(sname,'_unCor');
            end
            load(sname);
            %load type-dep data (STILL NEED TO ADD THIS)
            yname = strcat(sac_dir,base_yname,sprintf('_ori%d',ori_list(ee,ii)));
            if fit_unCor
                yname = strcat(yname,'_unCor');
            end
            load(yname);
            
            su_range = (n_probes+1):length(sacStimProc);
            clear SU_data
            for jj = 1:length(su_range)
                SU_data(jj).sacStimProc = sacStimProc(su_range(jj));
                SU_data(jj).trig_avg = sua_data(jj);
                SU_data(jj).type_dep = sacTypeDep(su_range(jj));
            end
            [SU_data.expt_num] = deal(Expt_num);
            [SU_data.bar_ori] = deal(ori_list(ee,ii));
            [SU_data.animal] = deal('lem');
            
            ori_data(ii,:) = SU_data;
        end
    end
    
    [cur_Noris,cur_Nunits] = size(ori_data);
    used = arrayfun(@(x) x.sacStimProc.used,ori_data);
    xvLLimps = nan(cur_Noris,cur_Nunits);
    xvLLimps(used) = arrayfun(@(x)x.sacStimProc.ModData.rectGQM.xvLLimp,ori_data(used));
    avg_rates = nan(cur_Noris,cur_Nunits);
    avg_rates(used) = arrayfun(@(x)x.sacStimProc.ModData.unit_data.avg_rate,ori_data(used));
    xvLLrates = xvLLimps.*avg_rates;
    
    [mvals,mlocs] = max(xvLLrates,[],1);
    nplocs = mod(mlocs,2)+1; %these are indices for the non-pref ori
    for ii = 1:cur_Nunits
        if ~isnan(xvLLrates(mlocs(ii),ii))
            all_SU_data = cat(1,all_SU_data,ori_data(mlocs(ii),ii));
        end
    end
end
%% time axis for trig-avg data
tlags = trig_avg_params.lags;
Tdt = trig_avg_params.dt;
tlags = tlags*Tdt;

%% SELECT USABLE CELLS
%selection criteria
min_rate = 5; % min avg rate in Hz (5)
min_Nsacs = 1e3; % min number of saccades for full analysis (1e3)
min_TA_Nsacs = 250; %min number of saccades for trig avg analysis
% min_xvLLimp = 0.05; %(0.05);

tot_Nunits = length(all_SU_data);
avg_rates = arrayfun(@(x) x.sacStimProc.ModData.unit_data.avg_rate,all_SU_data);
tot_spikes = arrayfun(@(x) x.sacStimProc.ModData.unit_data.tot_spikes,all_SU_data);
rec_dur = arrayfun(@(x) x.sacStimProc.ModData.unit_data.N_used_samps,all_SU_data)*dt/60;
mod_xvLLimps = arrayfun(@(x) x.sacStimProc.ModData.rectGQM.xvLLimp,all_SU_data);
expt_nums = [all_SU_data(:).expt_num];
expt_oris = [all_SU_data(:).bar_ori];

clust_iso_dist = arrayfun(@(x) x.sacStimProc.ModData.unit_data.SU_isodist,all_SU_data);
clust_Lratio = arrayfun(@(x) x.sacStimProc.ModData.unit_data.SU_Lratio,all_SU_data);
clust_refract = arrayfun(@(x) x.sacStimProc.ModData.unit_data.SU_refract,all_SU_data);
clust_dprime = arrayfun(@(x) x.sacStimProc.ModData.unit_data.SU_dprime,all_SU_data);
rate_stability_cv = arrayfun(@(x) x.sacStimProc.ModData.unit_data.rate_stability_cv,all_SU_data);
dprime_stability_cv = arrayfun(@(x) x.sacStimProc.ModData.unit_data.dprime_stability_cv,all_SU_data);

jbe_SUs = find(strcmp('jbe',{all_SU_data(:).animal}));
lem_SUs = find(strcmp('lem',{all_SU_data(:).animal}));

lem_fov_expt_nums = [266 275];
lem_parafov_expt_nums = [270 277 281 287 289 294 229 297];
fov_SUs = find([all_SU_data(:).expt_num] < 200 | ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
parafov_SUs = find([all_SU_data(:).expt_num] > 200 & ~ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
lem_parafov_SUs = intersect(parafov_SUs,lem_SUs);
lem_fov_SUs = intersect(fov_SUs,lem_SUs);

N_gsacs = arrayfun(@(x) x.trig_avg.N_gsacs,all_SU_data);
N_msacs = arrayfun(@(x) x.trig_avg.N_msacs,all_SU_data);
N_gsacs_gray = arrayfun(@(x) x.trig_avg.N_gsacs_gray,all_SU_data);
N_msacs_gray = arrayfun(@(x) x.trig_avg.N_msacs_gray,all_SU_data);

%% GRAY-BACKGROUND SACCADE TRIGGERED AVERAGES
cur_SUs = find(avg_rates >= min_rate & N_gsacs_gray >= min_TA_Nsacs);
all_gsac_gray = cell2mat(arrayfun(@(x) x.trig_avg.gsac_gray_avg', all_SU_data(:),'uniformoutput',0));

xl = [-0.1 0.35];

close all

%FOR SEPARATE MONKEYS
% f1 = figure(); hold on
% curset = intersect(jbe_SUs,cur_SUs);
% h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','r'});
% curset = intersect(lem_SUs,cur_SUs);
% h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','b'});

%COMBINED MONKEYS WITH SEPARATE TRACE FOR SUB_POP
f1 = figure(); hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(cur_SUs,:)),nanstd(all_gsac_gray(cur_SUs,:))/sqrt(length(cur_SUs)),{'color','b'});
curset = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs);
% h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','r'});
plot(tlags,nanmean(all_gsac_gray(curset,:)),'r','linewidth',2);

xlim(xl);
ylim([0.75 1.25]);
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

% %PRINT FIGURE
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_gray_avg.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% GSAC GRAYBACK STRENGTH VS TIMING SCATTERPLOT
cur_SUs = find(avg_rates >= min_rate & N_gsacs_gray >= min_TA_Nsacs);
all_gsac_gray = cell2mat(arrayfun(@(x) x.trig_avg.gsac_gray_avg', all_SU_data(:),'uniformoutput',0));

search_range = [0 0.35];
[gsac_Efact,gsac_exctime] = get_tavg_peaks((all_gsac_gray-1),tlags,search_range);
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_gray-1),tlags,search_range);

prc = 95;
Epeak_CI = nan(tot_Nunits,1);
Ipeak_CI = nan(tot_Nunits,1);
for ii = 1:tot_Nunits
    Epeak_CI(ii) = prctile(all_SU_data(ii).trig_avg.gsac_gray_nullP,prc);
    Ipeak_CI(ii) = prctile(all_SU_data(ii).trig_avg.gsac_gray_nullV,prc);
end
sigE = cur_SUs(gsac_Efact(cur_SUs) > Epeak_CI(cur_SUs));
sigI = cur_SUs(gsac_Ifact(cur_SUs) > Ipeak_CI(cur_SUs));
sigB = intersect(sigE,sigI);

xl = [0 1]; %mod strength axis
yl = [0 0.35]; %mod timing axis
mS = 3; %marker size
tjitter = Tdt/4; %jitter time peaks when plotting to prevent dot occlusion

close all
%PLOT MOD STRENGTH VS TIMING USING ONLY SIGNIFICANT PEAKS
f1 = figure();
hold on
plot(gsac_Efact(sigE),gsac_exctime(sigE) + randn(size(sigE))*tjitter,'bo','markersize',mS);
plot(gsac_Ifact(sigI),gsac_inhtime(sigI) + randn(size(sigI))*tjitter,'ro','markersize',mS);
set(gca,'xtick',0:0.2:1);
xlabel('Modulation strength');
ylabel('Modulation timing (s)');
xlim(xl); ylim(yl);

%create marginal distribution histograms
nbins = 20;
efact_binedges = linspace(xl(1),xl(2),nbins+1);
time_binedges = linspace(yl(1),yl(2),nbins+1);
gsac_Efact_dist = histc(gsac_Efact(sigE),efact_binedges);
gsac_Efact_dist = gsac_Efact_dist/sum(gsac_Efact_dist);
gsac_Etime_dist = histc(gsac_exctime(sigE),time_binedges);
gsac_Etime_dist = gsac_Etime_dist/sum(gsac_Etime_dist);
gsac_Ifact_dist = histc(gsac_Ifact(sigI),efact_binedges);
gsac_Ifact_dist = gsac_Ifact_dist/sum(gsac_Ifact_dist);
gsac_Itime_dist = histc(gsac_inhtime(sigI),time_binedges);
gsac_Itime_dist = gsac_Itime_dist/sum(gsac_Itime_dist);

%mod strength histogram
f2 = figure();
hold on
stairs(efact_binedges,gsac_Efact_dist);
stairs(efact_binedges,gsac_Ifact_dist,'r');
xlim(xl);

%peak timing histogram
f3 = figure();
hold on
stairs(time_binedges,gsac_Etime_dist);
stairs(time_binedges,gsac_Itime_dist,'r');
xlim(yl);

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_time_mod_scatter.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'Gsac_mod_dist.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'Gsac_time_dist.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

%% FOR SUPPLEMENTARY FIGURE LOOKING AT "REVERSE POLARITY" UNITS
%SCATTERPLOT OF MODULATION TIMING
yl = [0 0.3]; %mod timing axis
f1 = figure();
plot(gsac_inhtime(sigB) + randn(size(sigB))*tjitter,gsac_exctime(sigB) + randn(size(sigB))*tjitter,'bo','markersize',mS);
line(yl,yl,'color','k');
xlim(yl); ylim(yl);
xlabel('Suppression time (s)');
ylabel('Enhancement time (s)');

xl = [-0.1 0.35];

%PLOT SAC TRIG AVGS FOR TWO SETS OF UNITS (NORMAL AND REVERSE-POLARITY)
set1 = sigB(gsac_exctime(sigB) > gsac_inhtime(sigB));
set2 = sigB(gsac_exctime(sigB) < gsac_inhtime(sigB));
f2 = figure();
hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(set1,:)),nanstd(all_gsac_gray(set1,:))/sqrt(length(set1)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(set2,:)),nanstd(all_gsac_gray(set2,:))/sqrt(length(set2)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_timing_scatter.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'Gsac_reverse_polarity_avgs.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 

%% TB GAIN AND OFFSET 
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs);
TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod.lagX*dt;

TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates

TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));

%plot relative TB offsets
f1 = figure();
h1 = shadedErrorBar(TB_Xtick,nanmean(TB_Noffset),nanstd(TB_Noffset)/sqrt(length(cur_SUs)),{'color','b'});
xlabel('Time (s)');
ylabel('Normalized fffset');
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim([-0.1 0.3]);

%plot relative TB gains
f2 = figure();
h1 = shadedErrorBar(TB_Xtick,nanmean(TB_gains),nanstd(TB_gains)/sqrt(length(cur_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Gain');
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim([-0.1 0.3]);

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_TB_offset.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'Gsac_TB_gain.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% TB INFO AND INFO RATE
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs);

TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info
gsac_avg_rates = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_avg_rate',all_SU_data(cur_SUs),'uniformoutput',0));

TB_SSI_rate = TB_SSI.*gsac_avg_rates;
ov_SSI_rate = TB_ovinfos.*gsac_ov_rates;
TB_NSSI_rate = bsxfun(@rdivide,TB_SSI_rate,ov_SSI_rate);

%plot TB SSI and SSI rate
f1 = figure(); hold on
h1 = shadedErrorBar(slags*dt,nanmean(TB_NSSI),nanstd(TB_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2 = shadedErrorBar(slags*dt,nanmean(TB_NSSI_rate),nanstd(TB_NSSI_rate)/sqrt(length(cur_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Normalized fffset');
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim([-0.1 0.3]);

%% PLOTS COMPARING STRONGLY ENHANCED VS STRONGLY SUPPRESSED UNITS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs);
all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));

search_range = [0 0.35];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_tavg-1),tlags,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_tavg-1,tlags,search_range);

EIrat = gsac_Efact./gsac_Ifact;
EIprc = prctile(EIrat,[25 50 75]);

enh_units = find(EIrat >= EIprc(3));
sup_units = find(EIrat <= EIprc(1));
% enh_units = find(EIrat >= EIprc(2));
% sup_units = find(EIrat <= EIprc(2));

TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info

TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates
TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod.lagX*dt;

close all
%COMPARE TRIG AVG RATES
xl = [-0.1 0.3];
f1 = figure();hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_tavg(enh_units,:)),nanstd(all_gsac_tavg(enh_units,:))/sqrt(length(enh_units)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(all_gsac_tavg(sup_units,:)),nanstd(all_gsac_tavg(sup_units,:))/sqrt(length(sup_units)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');

%COMPARE SSIs
f1 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(TB_NSSI(enh_units,:)),nanstd(TB_NSSI(enh_units,:))/sqrt(length(enh_units)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(TB_NSSI(sup_units,:)),nanstd(TB_NSSI(sup_units,:))/sqrt(length(sup_units)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative SSI');

%COMPARE OFFSETS
f1 = figure();hold on
h1=shadedErrorBar(TB_Xtick,nanmean(TB_Noffset(enh_units,:)),nanstd(TB_Noffset(enh_units,:))/sqrt(length(enh_units)),{'color','b'});
h2=shadedErrorBar(TB_Xtick,nanmean(TB_Noffset(sup_units,:)),nanstd(TB_Noffset(sup_units,:))/sqrt(length(sup_units)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative Offset');

%COMPARE GAINS
f1 = figure();hold on
h1=shadedErrorBar(TB_Xtick,nanmean(TB_gains(enh_units,:)),nanstd(TB_gains(enh_units,:))/sqrt(length(enh_units)),{'color','b'});
h2=shadedErrorBar(TB_Xtick,nanmean(TB_gains(sup_units,:)),nanstd(TB_gains(sup_units,:))/sqrt(length(sup_units)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Gain');


%% COMPARE SSI, OFFSET AND GAINS FOR GO AND TB MODELS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs);
% TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
% gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
% cur_SUs = cur_SUs(TB_ovinfos >= prctile(TB_ovinfos,75));
% cur_SUs = cur_SUs(gsac_ov_rates >= median(gsac_ov_rates));

TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info
TB_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovLL = arrayfun(@(x) x.sacStimProc.gsac_TBmod.LLimp,all_SU_data(cur_SUs))*log(2); %overall TB model infos
TB_NLL = bsxfun(@rdivide,TB_LL,TB_ovLL); %normalize TB SSI by overall model info

TB_xvLL = arrayfun(@(x) max(x.sacStimProc.gsac_TBmod.xvLLs),all_SU_data(cur_SUs));

TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates
TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod.lagX*dt;

GO_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
GO_Noffset = bsxfun(@rdivide,GO_offset,gsac_ov_rates);
GO_gain = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_mod.ovInfo,all_SU_data(cur_SUs));
GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info
GO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovLL = arrayfun(@(x) x.sacStimProc.gsac_post_mod.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
GO_NLL = bsxfun(@rdivide,GO_LL,GO_ovLL); %normalize TB SSI by overall model info

GO_xvLL = arrayfun(@(x) max(x.sacStimProc.gsac_post_mod.lambda_xvLLImp(:)),all_SU_data(cur_SUs));

close all
%COMPARE SSIs
xl = [-0.1 0.3];
f1 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(TB_NSSI),nanstd(TB_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(GO_NSSI),nanstd(GO_NSSI)/sqrt(length(cur_SUs)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');

% xl = [-0.1 0.3];
% f1 = figure();hold on
% h1=shadedErrorBar(slags*dt,nanmean(TB_NLL),nanstd(TB_NLL)/sqrt(length(cur_SUs)),{'color','b'});
% h2=shadedErrorBar(slags*dt,nanmean(GO_NLL),nanstd(GO_NLL)/sqrt(length(cur_SUs)),{'color','r'});
% % plot(tlags,all_gsac_gray(set2,:),'k')
% line(xl,[1 1],'color','k');
% line([0 0],ylim(),'color','k');
% xlim(xl);
% xlabel('Time (s)');
% ylabel('Relative rate');

%COMPARE OFFSETS
xl = [-0.1 0.3];
f2 = figure();hold on
h1=shadedErrorBar(TB_Xtick,nanmean(TB_Noffset),nanstd(TB_Noffset)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(GO_Noffset),nanstd(GO_Noffset)/sqrt(length(cur_SUs)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');

%COMPARE GAINS
xl = [-0.1 0.3];
f3 = figure();hold on
h1=shadedErrorBar(TB_Xtick,nanmean(TB_gains),nanstd(TB_gains)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(GO_gain),nanstd(GO_gain)/sqrt(length(cur_SUs)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');

%% COMPARE GSACS AND MSACS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & N_msacs >= min_Nsacs);
% cur_SUs = cur_SUs(ismember(cur_SUs,jbe_SUs));

all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
GO_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
GO_Noffset = bsxfun(@rdivide,GO_offset,gsac_ov_rates);
GO_gain = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_mod.ovInfo,all_SU_data(cur_SUs));
GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info

all_msac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.msac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
msac_ov_rates = arrayfun(@(x) x.sacStimProc.msac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
mGO_offset = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_Noffset = bsxfun(@rdivide,mGO_offset,msac_ov_rates);
mGO_gain = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_ovinfos = arrayfun(@(x) x.sacStimProc.msac_post_mod.ovInfo,all_SU_data(cur_SUs));
mGO_NSSI = bsxfun(@rdivide,mGO_SSI,mGO_ovinfos); %normalize GO SSI by overall model info

search_range = [0 0.35];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_tavg-1),tlags,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_tavg-1,tlags,search_range);
[msac_Ifact,msac_inhtime] = get_tavg_peaks(-(all_msac_tavg-1),tlags,search_range);
[msac_Efact,msac_exctime] = get_tavg_peaks(all_msac_tavg-1,tlags,search_range);

% close all
% %COMPARE RATES
% xl = [-0.1 0.3];
% f3 = figure();hold on
% h1=shadedErrorBar(tlags,nanmean(all_gsac_tavg),nanstd(all_gsac_tavg)/sqrt(length(cur_SUs)),{'color','b'});
% h2=shadedErrorBar(tlags,nanmean(all_msac_tavg),nanstd(all_msac_tavg)/sqrt(length(cur_SUs)),{'color','r'});
% line(xl,[1 1],'color','k');
% line([0 0],ylim(),'color','k');
% xlim(xl);
% xlabel('Time (s)');
% ylabel('Relative rate');
% ylim([0.7 1.25])
% 
%COMPARE SSI
xl = [-0.1 0.3];
f3 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(GO_NSSI),nanstd(GO_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(mGO_NSSI),nanstd(mGO_NSSI)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative SSI');

% %COMPARE OFFSETS
% xl = [-0.1 0.3];
% f3 = figure();hold on
% h1=shadedErrorBar(slags*dt,nanmean(GO_Noffset),nanstd(GO_Noffset)/sqrt(length(cur_SUs)),{'color','b'});
% h2=shadedErrorBar(slags*dt,nanmean(mGO_Noffset),nanstd(mGO_Noffset)/sqrt(length(cur_SUs)),{'color','r'});
% line(xl,[0 0],'color','k');
% line([0 0],ylim(),'color','k');
% xlim(xl);
% xlabel('Time (s)');
% ylabel('Relative offset');
% 
% %COMPARE GAINS
% xl = [-0.1 0.3];
% f3 = figure();hold on
% h1=shadedErrorBar(slags*dt,nanmean(GO_gain),nanstd(GO_gain)/sqrt(length(cur_SUs)),{'color','b'});
% h2=shadedErrorBar(slags*dt,nanmean(mGO_gain),nanstd(mGO_gain)/sqrt(length(cur_SUs)),{'color','r'});
% line(xl,[1 1],'color','k');
% line([0 0],ylim(),'color','k');
% xlim(xl);
% xlabel('Time (s)');
% ylabel('Gain');

%%
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & N_msacs >= min_Nsacs);

all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsacPreGainMod.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsacPreGainMod.ovInfo,all_SU_data(cur_SUs));
GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info

all_msac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.msac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
msac_ov_rates = arrayfun(@(x) x.sacStimProc.msac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
mGO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.msacPreGainMod.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_ovinfos = arrayfun(@(x) x.sacStimProc.msacPreGainMod.ovInfo,all_SU_data(cur_SUs));
mGO_NSSI = bsxfun(@rdivide,mGO_SSI,mGO_ovinfos); %normalize GO SSI by overall model info

search_range = [0 0.35];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_tavg-1),tlags,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_tavg-1,tlags,search_range);
[msac_Ifact,msac_inhtime] = get_tavg_peaks(-(all_msac_tavg-1),tlags,search_range);
[msac_Efact,msac_exctime] = get_tavg_peaks(all_msac_tavg-1,tlags,search_range);

close all
%COMPARE RATES
xl = [-0.1 0.3];
f3 = figure();hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_tavg),nanstd(all_gsac_tavg)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(all_msac_tavg),nanstd(all_msac_tavg)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.7 1.25])

%COMPARE SSI
xl = [-0.1 0.3];
f3 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(GO_NSSI),nanstd(GO_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(mGO_NSSI),nanstd(mGO_NSSI)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative SSI');

