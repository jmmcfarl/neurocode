% close all
clear all
clc

fit_unCor = 0;
include_bursts = 0;

fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';
% fig_dir = '/Users/james/Analysis/bruce/FINsac_mod/figures/';
base_sname = 'sacStimProcFin_noXV';
% base_tname = 'sac_trig_avg_data4'; %this data set has all the bootsrap an
base_tname = 'sac_trig_avg_data_test';
base_yname = 'sacTypeDep_noXV';
base_iname = 'sac_info_timing_noXV3';
base_dname = 'sacStimDelay_noXV';

if include_bursts
    base_tname = strcat(base_tname,'_withbursts');
    base_sname = strcat(base_sname,'_withbursts');
    base_yname = strcat(base_yname,'_withbursts');
    base_dname = strcat(base_dname,'_withbursts');
end

all_SU_data = [];
all_SU_NPdata = [];

%% LOAD JBE
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
n_probes = 96;
ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
% ori_list = [0 90; 0 90; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan];
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
            %             if fit_unCor
            %                 yname = strcat(yname,'_unCor');
            %             end
            load(yname);
            
            iname = strcat(sac_dir,base_iname,sprintf('_ori%d',ori_list(ee,ii)));
            load(iname);
            
            dname = strcat(sac_dir,base_dname,sprintf('_ori%d',ori_list(ee,ii)));
            load(dname);
            
            su_range = (n_probes+1):length(sacStimProc);
            clear SU_data
            for jj = 1:length(su_range)
                SU_data(jj).sacStimProc = sacStimProc(su_range(jj));
                SU_data(jj).trig_avg = sua_data(jj);
                SU_data(jj).type_dep = sacTypeDep(su_range(jj));
                SU_data(jj).info_time = sacInfoTiming(su_range(jj));
                SU_data(jj).sac_delay = sacDelay(su_range(jj));
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
% Expt_list = {'M266','M270','M275','M277','M281'};%NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
n_probes = 24;
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
% rmfield_list = {};

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
            if ~ismember(Expt_num,[296 297])
                %load type-dep data (STILL NEED TO ADD THIS)
                yname = strcat(sac_dir,base_yname,sprintf('_ori%d',ori_list(ee,ii)));
                %                 if fit_unCor
                %                     yname = strcat(yname,'_unCor');
                %                 end
                load(yname);
            end
            
            iname = strcat(sac_dir,base_iname,sprintf('_ori%d',ori_list(ee,ii)));
            load(iname);
            
            dname = strcat(sac_dir,base_dname,sprintf('_ori%d',ori_list(ee,ii)));
            load(dname);
            
            su_range = (n_probes+1):length(sacStimProc);
            clear SU_data
            for jj = 1:length(su_range)
                SU_data(jj).sacStimProc = sacStimProc(su_range(jj));
                SU_data(jj).trig_avg = sua_data(jj);
                if ~ismember(Expt_num,[296 297])
                    SU_data(jj).type_dep = sacTypeDep(su_range(jj));
                else
                    SU_data(jj).type_dep = nan;
                end
                SU_data(jj).info_time = sacInfoTiming(su_range(jj));
                SU_data(jj).sac_delay = sacDelay(su_range(jj));
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
min_Nsacs = 500; % min number of saccades for full analysis (1e3)
min_TA_Nsacs = 500; %min number of saccades for trig avg analysis
min_xvLLimp = 0.0; %(0.05); %minimum cross-validated LL improvement for stim-proc models

%basic stats for each cell
tot_Nunits = length(all_SU_data);
avg_rates = arrayfun(@(x) x.sacStimProc.ModData.unit_data.avg_rate,all_SU_data);
tot_spikes = arrayfun(@(x) x.sacStimProc.ModData.unit_data.tot_spikes,all_SU_data);
rec_dur = arrayfun(@(x) x.sacStimProc.ModData.unit_data.N_used_samps,all_SU_data)*dt/60;
mod_xvLLimps = arrayfun(@(x) x.sacStimProc.ModData.rectGQM.xvLLimp,all_SU_data);
% mod_xvLLimps_unCor = arrayfun(@(x) x.sacStimProc.ModData.rectGQM_unCor.xvLLimp,all_SU_data);
expt_nums = [all_SU_data(:).expt_num];
expt_oris = [all_SU_data(:).bar_ori];

%cluster properties
clust_iso_dist = arrayfun(@(x) x.sacStimProc.ModData.unit_data.SU_isodist,all_SU_data);
clust_Lratio = arrayfun(@(x) x.sacStimProc.ModData.unit_data.SU_Lratio,all_SU_data);
clust_refract = arrayfun(@(x) x.sacStimProc.ModData.unit_data.SU_refract,all_SU_data);
clust_dprime = arrayfun(@(x) x.sacStimProc.ModData.unit_data.SU_dprime,all_SU_data);
rate_stability_cv = arrayfun(@(x) x.sacStimProc.ModData.unit_data.rate_stability_cv,all_SU_data);
dprime_stability_cv = arrayfun(@(x) x.sacStimProc.ModData.unit_data.dprime_stability_cv,all_SU_data);

RF_ecc = arrayfun(@(x) x.sacStimProc.ModData.tune_props.RF_ecc,all_SU_data);
RF_sigma = arrayfun(@(x) x.sacStimProc.ModData.tune_props.RF_sigma,all_SU_data);
RF_gSF = arrayfun(@(x) x.sacStimProc.ModData.tune_props.RF_gSF,all_SU_data);
RF_FSF = arrayfun(@(x) x.sacStimProc.ModData.tune_props.RF_FSF,all_SU_data);
RF_sX = arrayfun(@(x) x.sacStimProc.ModData.tune_props.screen_X,all_SU_data);
RF_sY = arrayfun(@(x) x.sacStimProc.ModData.tune_props.screen_Y,all_SU_data);

%cell rec props
jbe_SUs = find(strcmp('jbe',{all_SU_data(:).animal}));
lem_SUs = find(strcmp('lem',{all_SU_data(:).animal}));

lem_fov_expt_nums = [266 275];
lem_parafov_expt_nums = [270 277 281 287 289 294 229 297];
fov_SUs = find([all_SU_data(:).expt_num] < 200 | ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
parafov_SUs = find([all_SU_data(:).expt_num] > 200 & ~ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
lem_parafov_SUs = intersect(parafov_SUs,lem_SUs);
lem_fov_SUs = intersect(fov_SUs,lem_SUs);

%number of each type of saccade for each SU
N_gsacs = arrayfun(@(x) x.trig_avg.N_gsacs,all_SU_data);
N_msacs = arrayfun(@(x) x.trig_avg.N_msacs,all_SU_data);
N_gsacs_gray = arrayfun(@(x) x.trig_avg.N_gsacs_gray,all_SU_data);
N_msacs_gray = arrayfun(@(x) x.trig_avg.N_msacs_gray,all_SU_data);
N_gsacs_im = arrayfun(@(x) x.trig_avg.N_gsacs_im,all_SU_data);
N_msacs_im = arrayfun(@(x) x.trig_avg.N_msacs_im,all_SU_data);
N_simsacs = arrayfun(@(x) x.trig_avg.N_simsacs,all_SU_data);

%% GRAY-BACKGROUND SACCADE TRIGGERED AVERAGES
% cur_SUs = find(avg_rates >= min_rate & N_gsacs_gray >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
% all_gsac_gray = cell2mat(arrayfun(@(x) x.trig_avg.gsac_gray_avg', all_SU_data(:),'uniformoutput',0));
%USE pooled data across grayback and image back
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
all_gsac = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(:),'uniformoutput',0));

xl = [-0.1 0.3]; %x-axis limit

close all

%FOR SEPARATE MONKEYS
% f1 = figure(); hold on
% curset = intersect(jbe_SUs,cur_SUs);
% h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','r'});
% curset = intersect(lem_SUs,cur_SUs);
% h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','b'});

%COMBINED MONKEYS WITH SEPARATE TRACE FOR SUB_POP
f1 = figure(); hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac(cur_SUs,:)),nanstd(all_gsac(cur_SUs,:))/sqrt(length(cur_SUs)),{'color','k'});
% curset = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs);
% h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','r'});
% plot(tlags,nanmean(all_gsac_gray(curset,:)),'r','linewidth',2);
xlim(xl);
ylim([0.7 1.25]);
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

%this plot is using a broadened x-axis for comparison with the sac-in-dark
%condition
xl2 = [-0.2 0.5];
f2 = figure();
plot(tlags,nanmean(all_gsac(cur_SUs,:)),'linewidth',2);
set(gca,'yAxisLocation','right');
ylim([0.65 1.3]);
xlim(xl2);

% %PRINT FIGURE
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_gray_avg.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

% figufy(f2);
% fname = [fig_dir 'Gsac_gray_avg_compare.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% GSAC GRAYBACK STRENGTH VS TIMING SCATTERPLOT
% cur_SUs = find(avg_rates >= min_rate & N_gsacs_gray >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
% all_gsac_gray = cell2mat(arrayfun(@(x) x.trig_avg.gsac_gray_avg', all_SU_data(:),'uniformoutput',0));
% all_msac_gray = cell2mat(arrayfun(@(x) x.trig_avg.msac_gray_avg', all_SU_data(:),'uniformoutput',0));
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
all_gsac = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(:),'uniformoutput',0));
all_msac = cell2mat(arrayfun(@(x) x.trig_avg.msac_avg', all_SU_data(:),'uniformoutput',0));

tlags_up = linspace(tlags(1),tlags(end),500); %up-sampled t-axis for interpolation
all_gsac_up = nan(size(all_gsac,1),length(tlags_up));
incl = find(~isnan(all_gsac(:,1)));
all_gsac_up(incl,:) = spline(tlags,all_gsac(incl,:),tlags_up); %spline interpolation

%get the timing and magnitude of enhancement and suppression peaks
search_range = [0 0.3]; %range of time lags to search for local extrema
[gsac_Efact,gsac_exctime] = get_tavg_peaks((all_gsac-1),tlags,search_range);
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac-1),tlags,search_range);
[gsac_Efact_up,gsac_exctime_up] = get_tavg_peaks((all_gsac_up-1),tlags_up,search_range);
[gsac_Ifact_up,gsac_inhtime_up] = get_tavg_peaks(-(all_gsac_up-1),tlags_up,search_range);

prc = 95; %pull out the upper prc from the null sampling dists of peak amps
Epeak_CI = nan(tot_Nunits,1);
Ipeak_CI = nan(tot_Nunits,1);
for ii = 1:tot_Nunits
    %     Epeak_CI(ii) = prctile(all_SU_data(ii).trig_avg.gsac_gray_nullP,prc);
    %     Ipeak_CI(ii) = prctile(all_SU_data(ii).trig_avg.gsac_gray_nullV,prc);
    Epeak_CI(ii) = prctile(all_SU_data(ii).trig_avg.gsac_nullP,prc);
    Ipeak_CI(ii) = prctile(all_SU_data(ii).trig_avg.gsac_nullV,prc);
end
sigE = cur_SUs(gsac_Efact(cur_SUs) > Epeak_CI(cur_SUs)); %significant enhancement peaks
sigI = cur_SUs(gsac_Ifact(cur_SUs) > Ipeak_CI(cur_SUs)); %significant suppression peaks
sigB = intersect(sigE,sigI); %units with significant E and S peaks
normal_polarity_units = sigB(gsac_exctime(sigB) > gsac_inhtime(sigB)); %units with normal I<E timing
reverse_polarity_units = sigB(gsac_exctime(sigB) < gsac_inhtime(sigB)); %units with reverse E<I timing
nsigE = setdiff(cur_SUs,sigE); %number of non-sig E peaks
nsigI = setdiff(cur_SUs,sigI); %number of non-sig I peaks
nsigB = setdiff(cur_SUs,sigB); %number of non-both sigs

xl = [0 1]; %mod strength axis
yl = [0 0.3]; %mod timing axis
mS = 4; %marker size
mS2 = 8; %marker size

tjitter = 0; %jitter time peaks when plotting to prevent dot occlusion
% close all
%PLOT MOD STRENGTH VS TIMING USING ONLY SIGNIFICANT PEAKS
f1 = figure();
hold on
% plot(gsac_Efact(sigE),gsac_exctime(sigE) + randn(size(sigE))*tjitter,'b.','markersize',mS);
% plot(gsac_Ifact(sigI),gsac_inhtime(sigI) + randn(size(sigI))*tjitter,'r.','markersize',mS);
plot(gsac_Efact_up(sigE),gsac_exctime_up(sigE) + randn(size(sigE))*tjitter,'bo','markersize',mS,'linewidth',0.5);
plot(gsac_Ifact_up(sigI),gsac_inhtime_up(sigI) + randn(size(sigI))*tjitter,'ro','markersize',mS,'linewidth',0.5);
plot(gsac_Efact_up(nsigE),gsac_exctime_up(nsigE) ,'b.','markersize',mS2);
plot(gsac_Ifact_up(nsigI),gsac_inhtime_up(nsigI),'r.','markersize',mS2);
set(gca,'xtick',0:0.2:1);
xlabel('Modulation strength');
ylabel('Modulation timing (s)');
xlim(xl); ylim(yl);

%create marginal distribution histograms
nbins = 20;
efact_binedges = linspace(xl(1),xl(2),nbins+1);
time_binedges = linspace(yl(1),yl(2),nbins+1);
% gsac_Efact_dist = histc(gsac_Efact(sigE),efact_binedges);
gsac_Efact_dist = histc(gsac_Efact_up(sigE),efact_binedges);
gsac_Efact_dist = gsac_Efact_dist/sum(gsac_Efact_dist);
% gsac_Etime_dist = histc(gsac_exctime(sigE),time_binedges);
gsac_Etime_dist = histc(gsac_exctime_up(sigE),time_binedges);
gsac_Etime_dist = gsac_Etime_dist/sum(gsac_Etime_dist);
% gsac_Ifact_dist = histc(gsac_Ifact(sigI),efact_binedges);
gsac_Ifact_dist = histc(gsac_Ifact_up(sigI),efact_binedges);
gsac_Ifact_dist = gsac_Ifact_dist/sum(gsac_Ifact_dist);
% gsac_Itime_dist = histc(gsac_inhtime(sigI),time_binedges);
gsac_Itime_dist = histc(gsac_inhtime_up(sigI),time_binedges);
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

% % %PRINT PLOTS
% fig_width = 3.5; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Gsac_time_mod_scatter.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Gsac_mod_dist.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
%
% figufy(f3);
% fname = [fig_dir 'Gsac_time_dist.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
%
%% FOR SUPPLEMENTARY FIGURE LOOKING AT "REVERSE POLARITY" UNITS
%SCATTERPLOT OF MODULATION TIMING
yl = [0 0.3]; %mod timing axis

mS = 3;
%PLOT SAC TRIG AVGS FOR TWO SETS OF UNITS (NORMAL AND REVERSE-POLARITY)
f1 = figure(); hold on
plot(gsac_inhtime(normal_polarity_units),gsac_exctime(normal_polarity_units),'bo','markersize',mS,'linewidth',0.7);
plot(gsac_inhtime(reverse_polarity_units),gsac_exctime(reverse_polarity_units),'ro','markersize',mS,'linewidth',0.7);
plot(gsac_inhtime(nsigB),gsac_exctime(nsigB),'k.','markersize',8);
line(yl,yl,'color','k');
xlim(yl); ylim(yl);
xlabel('Suppression time (s)');
ylabel('Enhancement time (s)');

xl = [-0.1 0.3];

f2 = figure();
hold on
h2=shadedErrorBar(tlags,nanmean(all_gsac(reverse_polarity_units,:)),nanstd(all_gsac(reverse_polarity_units,:))/sqrt(length(reverse_polarity_units)),{'color','r'});
h1=shadedErrorBar(tlags,nanmean(all_gsac(normal_polarity_units,:)),nanstd(all_gsac(normal_polarity_units,:))/sqrt(length(normal_polarity_units)),{'color','b'});
% plot(tlags,all_gsac_gray(reverse_polarity_units,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.6 1.6]);

% f3 = figure();
% hold on
% h2=shadedErrorBar(tlags,nanmean(all_msac_gray(reverse_polarity_units,:)),nanstd(all_msac_gray(reverse_polarity_units,:))/sqrt(length(reverse_polarity_units)),{'color','r'});
% h1=shadedErrorBar(tlags,nanmean(all_msac_gray(normal_polarity_units,:)),nanstd(all_msac_gray(normal_polarity_units,:))/sqrt(length(normal_polarity_units)),{'color','b'});
% % plot(tlags,all_gsac_gray(reverse_polarity_units,:),'k')
% line(xl,[1 1],'color','k');
% line([0 0],ylim(),'color','k');
% xlim(xl);
% xlabel('Time (s)');
% ylabel('Relative rate');
% ylim([0.6 1.4]);

% % %PRINT PLOTS
% fig_width = 3.5; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Gsac_timing_scatter.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Gsac_reverse_polarity_avgs.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Gsac_reverse_polarity_avgs.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% TB GAIN AND OFFSET
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod{1}.lagX*dt; %time axis for TB rep
base_lags = find(TB_Xtick <= 0); %pre-sac time lags

lambda_ii = 4; %selection of reg strength parameter
xl = [-0.1 0.3]; %x-axis window for plotting

TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));

TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates
% TB_Noffset = bsxfun(@rdivide,TB_Noffset,TB_gains); %normalize by gains

%normalize gain and offset by baseline avgs
TB_gains = bsxfun(@rdivide,TB_gains,mean(TB_gains(:,base_lags),2));
TB_Noffset = bsxfun(@minus,TB_Noffset,mean(TB_Noffset(:,base_lags),2));

%plot relative TB offsets
f1 = figure();
h1 = shadedErrorBar(TB_Xtick,nanmean(TB_Noffset),nanstd(TB_Noffset)/sqrt(length(cur_SUs)),{'color','b'});
xlabel('Time (s)');
ylabel('Normalized offset');
line(xl,[0 0],'color','k');
% line([0 0],ylim(),'color','k');
xlim([-0.1 0.3]);
% ylim([-0.05 0.3]);
ylim([-0.3 0.3]);

%plot relative TB gains
f2 = figure();
h1 = shadedErrorBar(TB_Xtick,nanmean(TB_gains),nanstd(TB_gains)/sqrt(length(cur_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Gain');
line(xl,[1 1],'color','k');
% line([0 0],ylim(),'color','k');
xlim([-0.1 0.3]);
% ylim([0.5 1.15]);
ylim([0.5 1.5]);
set(gca,'YaxisLocation','right');

% % PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_TB_offset2.pdf'];
% if fit_unCor
% fname = [fig_dir 'Gsac_TB_offset_unCor.pdf'];
% end
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'Gsac_TB_gain2.pdf'];
% if fit_unCor
% fname = [fig_dir 'Gsac_TB_gain_unCor.pdf'];
% end
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% TB INFO AND INFO RATE
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);

gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates

lambda_ii = 4; %selection of regularization strength hyperparam
xl = [-0.1 0.3];
TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info
TB_NSSI = bsxfun(@rdivide,TB_NSSI,mean(TB_NSSI(:,base_lags),2)); %normalize by baseline avgs

% gsac_avg_rates = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_avg_rate',all_SU_data(cur_SUs),'uniformoutput',0));
% gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
% gsac_rel_rates = bsxfun(@rdivide,gsac_avg_rates,gsac_ov_rates);

%saccade-triggered avg rates, interpolated onto same time axis to use for
%normalization
gsac_rel_rates = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
gsac_rel_rates_interp = interp1(tlags,gsac_rel_rates',slags*dt)';

% TB_SSI_rate = TB_SSI.*gsac_avg_rates;
% ov_SSI_rate = TB_ovinfos.*gsac_ov_rates;
% TB_NSSI_rate = bsxfun(@rdivide,TB_SSI_rate,ov_SSI_rate);

%normalize SSI by sac-trig-avg rate
TB_NSSI_rate = bsxfun(@times,TB_NSSI,gsac_rel_rates_interp);
TB_NSSI_rate = bsxfun(@rdivide,TB_NSSI_rate,mean(TB_NSSI_rate(:,base_lags),2)); %normalize by baseline avg


%plot TB SSI and SSI rate
f1 = figure(); hold on
h1 = shadedErrorBar(slags*dt,nanmean(TB_NSSI),nanstd(TB_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2 = shadedErrorBar(slags*dt,nanmean(TB_NSSI_rate),nanstd(TB_NSSI_rate)/sqrt(length(cur_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Normalized fffset');
line(xl,[1 1],'color','k');
% line([0 0],ylim(),'color','k');
xlim([-0.1 0.3]);
ylim([0.4 1.1])
% ylim([0.5 1.1])

%plot TB SSI and SSI rate
f2 = figure(); hold on
% h1 = shadedErrorBar(slags*dt,nanmean(gsac_rel_rates),nanstd(gsac_rel_rates)/sqrt(length(cur_SUs)),{'color','b'});
h1 = shadedErrorBar(tlags,nanmean(gsac_rel_rates),nanstd(gsac_rel_rates)/sqrt(length(cur_SUs)),{'color','k'});
xlabel('Time (s)');
ylabel('Normalized fffset');
line(xl,[1 1],'color','k');
% line([0 0],ylim(),'color','k');
xlim([-0.1 0.3]);
ylim([0.7 1.25]);

% % % % %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% % fname = [fig_dir 'Gsac_TB_SSI.pdf'];
% fname = [fig_dir 'Gsac_TB_SSI2.pdf'];
% if fit_unCor
% fname = [fig_dir 'Gsac_TB_SSI_unCor2.pdf'];
% end
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% %
% figufy(f2);
% fname = [fig_dir 'Gsac_TBset_relrates.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% PLOTS COMPARING STRONGLY ENHANCED VS STRONGLY SUPPRESSED UNITS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);

lambda_ii = 4; %reg hyper selection
xl = [-0.1 0.3];

TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod{lambda_ii}.lagX*dt; %TB time axis
TBbase_lags = find(TB_Xtick <= 0);

%sac-trig avg firing rate
all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
search_range = [0 0.3];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_tavg-1),tlags,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_tavg-1,tlags,search_range);

%get normalized TB-model SSI (computing using standard t-axis)
TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info
TB_NSSI = bsxfun(@rdivide,TB_NSSI,mean(TB_NSSI(:,base_lags),2));

%get normalized TB-model gain and offset (computed using TB-taxis)
TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates
% TB_Noffset = bsxfun(@rdivide,TB_Noffset,TB_gains); %normalize by gain
TB_Noffset = bsxfun(@minus,TB_Noffset,mean(TB_Noffset(:,TBbase_lags),2));
TB_gains = bsxfun(@rdivide,TB_gains,mean(TB_gains(:,TBbase_lags),2));

%find timing and magnitude of TB-model SSI suppression peaks
search_range = [0 0.15]; %search range for TB SSI suppression peaks
[SSI_Ifact,SSI_inhtime] = get_tavg_peaks(-(TB_NSSI-1),slags*dt,search_range);
[Gain_Ifact,GAIN_inhtime] = get_tavg_peaks(-(TB_gains-1),TB_Xtick,search_range);
% [Off_Efact,OFF_exctime] = get_tavg_peaks((TB_Noffset),TB_Xtick,search_range);

%compute partial corr coefs between SSI and gain/offset
interp_TB_gains = interp1(TB_Xtick,TB_gains',slags*dt)';
interp_TB_offset = interp1(TB_Xtick,TB_Noffset',slags*dt)';
clear gain_SSI_corr offset_SSI_corr
for ii = 1:length(cur_SUs)
    gain_SSI_corr(ii) = partialcorr(interp_TB_gains(ii,:)',TB_NSSI(ii,:)',interp_TB_offset(ii,:)','type','spearman');
    offset_SSI_corr(ii) = partialcorr(interp_TB_offset(ii,:)',TB_NSSI(ii,:)',interp_TB_gains(ii,:)','type','spearman');
end

% EIrat = gsac_Efact./gsac_Ifact;
% EIrat = Off_Efact;
% EIrat = gsac_Ifact;
% EIrat = Gain_Ifact;
% EIrat = gsac_Efact;
% EIrat = gsac_Ifact;
% EIprc = prctile(EIrat,[25 50 75]);
EIrat = SSI_Ifact; %divide neurons into groups based on magnitude of SSI suppression
EIprc = prctile(EIrat,[33 50 67]);

% enh_units = find(gsac_Efact > prctile(gsac_Efact,67));
% sup_units = find(gsac_Ifact > prctile(gsac_Ifact,67));

% enh_units = find(EIrat >= EIprc(3));
% sup_units = find(EIrat <= EIprc(1));
%divide into two groups based on median
enh_units = find(EIrat >= EIprc(2));
sup_units = find(EIrat <= EIprc(2));


close all

%COMPARE TRIG AVG RATES
f1 = figure();hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_tavg(enh_units,:)),nanstd(all_gsac_tavg(enh_units,:))/sqrt(length(enh_units)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(all_gsac_tavg(sup_units,:)),nanstd(all_gsac_tavg(sup_units,:))/sqrt(length(sup_units)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.7 1.25]);

%COMPARE SSIs
f2 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(TB_NSSI(enh_units,:)),nanstd(TB_NSSI(enh_units,:))/sqrt(length(enh_units)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(TB_NSSI(sup_units,:)),nanstd(TB_NSSI(sup_units,:))/sqrt(length(sup_units)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative SSI');
ylim([0.35 1.1]);

%COMPARE OFFSETS
f3 = figure();hold on
h1=shadedErrorBar(TB_Xtick,nanmean(TB_Noffset(enh_units,:)),nanstd(TB_Noffset(enh_units,:))/sqrt(length(enh_units)),{'color','b'});
h2=shadedErrorBar(TB_Xtick,nanmean(TB_Noffset(sup_units,:)),nanstd(TB_Noffset(sup_units,:))/sqrt(length(sup_units)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative Offset');
ylim([-0.1 0.4]);

%COMPARE GAINS
f4 = figure();hold on
h1=shadedErrorBar(TB_Xtick,nanmean(TB_gains(enh_units,:)),nanstd(TB_gains(enh_units,:))/sqrt(length(enh_units)),{'color','b'});
h2=shadedErrorBar(TB_Xtick,nanmean(TB_gains(sup_units,:)),nanstd(TB_gains(sup_units,:))/sqrt(length(sup_units)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Gain');
ylim([0.4 1.2]);

% % %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_SW_rates2.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'Gsac_SW_SSI2.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
%
% figufy(f3);
% fname = [fig_dir 'Gsac_SW_OFFSET2.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
%
% figufy(f4);
% fname = [fig_dir 'Gsac_SW_GAIN2.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% COMPARE SSI, OFFSET AND GAINS FOR GO AND TB MODELS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod{1}.lagX*dt;
base_lags = find(slags <= 0);
TBbase_lags = find(TB_Xtick <= 0);

GO_color = [0.25 0.8 0.25]; %color to plot GO model

%regularization hyper param selection
TB_lambda_ii = 4;
GO_lambda_off = 4;
GO_lambda_gain = 3;

%get TB-model SSI, normalized
TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info
TB_NSSI = bsxfun(@rdivide,TB_NSSI,mean(TB_NSSI(:,base_lags),2));

%get TB-model log-likelihood
TB_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovLL = arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.LLimp,all_SU_data(cur_SUs)); %overall TB model infos
TB_NLL = bsxfun(@rdivide,TB_LL,TB_ovLL); %normalize TB SSI by overall model info
% TB_NLL = bsxfun(@rdivide,TB_NLL,mean(TB_NLL(:,base_lags),2));

%get TB-model offset (normalized)
TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates
TB_Noffset = bsxfun(@minus,TB_Noffset,mean(TB_Noffset(:,TBbase_lags),2));

%get TB-model gains, normalized
TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
TB_gains = bsxfun(@rdivide,TB_gains,mean(TB_gains(:,TBbase_lags),2));

%GO model offset, normalized
GO_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
GO_Noffset = bsxfun(@rdivide,GO_offset,gsac_ov_rates);
GO_Noffset = bsxfun(@minus,GO_Noffset,mean(GO_Noffset(:,base_lags),2));

%GO model gain normalized
GO_gain = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
GO_gain = bsxfun(@rdivide,GO_gain,mean(GO_gain(:,base_lags),2));

%GO model SSI, normalized
GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info
GO_NSSI = bsxfun(@rdivide,GO_NSSI,mean(GO_NSSI(:,base_lags),2));

%GO-model LL
GO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovLL = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
GO_NLL = bsxfun(@rdivide,GO_LL,GO_ovLL); %normalize TB SSI by overall model info
% GO_NLL = bsxfun(@rdivide,GO_NLL,mean(GO_NLL(:,base_lags),2));

%slight temporal smoothing of model LLs
for ii = 1:length(cur_SUs)
    GO_NLL(ii,:) = jmm_smooth_1d_cor(GO_NLL(ii,:),1);
    TB_NLL(ii,:) = jmm_smooth_1d_cor(TB_NLL(ii,:),1);
end

close all
%COMPARE SSIs
xl = [-0.1 0.3];
f1 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(TB_NSSI),nanstd(TB_NSSI)/sqrt(length(cur_SUs)),{'color','k'});
h2=shadedErrorBar(slags*dt,nanmean(GO_NSSI),nanstd(GO_NSSI)/sqrt(length(cur_SUs)),{'color',GO_color});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative SSI');
ylim([0.5 1.1])

xl = [-0.1 0.3];
f1 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(TB_NLL),nanstd(TB_NLL)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(GO_NLL),nanstd(GO_NLL)/sqrt(length(cur_SUs)),{'color','r'});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Log-likelihood');

%COMPARE OFFSETS
xl = [-0.1 0.3];
f2 = figure();hold on
h1=shadedErrorBar(TB_Xtick,nanmean(TB_Noffset),nanstd(TB_Noffset)/sqrt(length(cur_SUs)),{'color','k'});
h2=shadedErrorBar(slags*dt,nanmean(GO_Noffset),nanstd(GO_Noffset)/sqrt(length(cur_SUs)),{'color',GO_color});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Offset');
ylim([-0.05 0.3])

%COMPARE GAINS
xl = [-0.1 0.3];
f3 = figure();hold on
h1=shadedErrorBar(TB_Xtick,nanmean(TB_gains),nanstd(TB_gains)/sqrt(length(cur_SUs)),{'color','k'});
h2=shadedErrorBar(slags*dt,nanmean(GO_gain),nanstd(GO_gain)/sqrt(length(cur_SUs)),{'color',GO_color});
% plot(tlags,all_gsac_gray(set2,:),'k')
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Gain');
ylim([0.5 1.15])

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_TB_GO_SSI.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'Gsac_TB_GO_OFFSET.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
%
% figufy(f3);
% fname = [fig_dir 'Gsac_TB_GO_GAIN.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

%% COMPARE GSACS AND MSACS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & N_msacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);

%reg hyper param selection
GO_lambda_off = 4;
GO_lambda_gain = 3;

%sac trig avg rates
all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates

%GO-model sac offset, normalized
GO_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
GO_Noffset = bsxfun(@rdivide,GO_offset,gsac_ov_rates);
GO_Noffset = bsxfun(@minus,GO_Noffset,mean(GO_Noffset(:,base_lags),2));

%GO-model sac gain, normalized
GO_gain = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
GO_gain = bsxfun(@rdivide,GO_gain,mean(GO_gain(:,base_lags),2));

%GO-model sac SSI, normalized
GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info
GO_NSSI = bsxfun(@rdivide,GO_NSSI,mean(GO_NSSI(:,base_lags),2));

% %GO-model LL
% GO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
% GO_ovLL = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
% GO_NLL = bsxfun(@rdivide,GO_LL,GO_ovLL); %normalize TB SSI by overall model info
% GO_NLL = bsxfun(@rdivide,GO_NLL,mean(GO_NLL(:,base_lags),2));

%micro-sac trig avg rates
all_msac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.msac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
msac_ov_rates = arrayfun(@(x) x.sacStimProc.msac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates

%micro-sac offset, normalized
mGO_offset = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_Noffset = bsxfun(@rdivide,mGO_offset,msac_ov_rates);
mGO_Noffset = bsxfun(@minus,mGO_Noffset,mean(mGO_Noffset(:,base_lags),2));

%micro-sac gains, normalized
mGO_gain = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_gain = bsxfun(@rdivide,mGO_gain,mean(mGO_gain(:,base_lags),2));

%micro-sac SSI, normalized
mGO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_ovinfos = arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
% mGO_NSSI = bsxfun(@rdivide,mGO_SSI,mGO_ovinfos); %normalize GO SSI by overall model info
mGO_NSSI = bsxfun(@rdivide,mGO_SSI,GO_ovinfos); %normalize GO SSI by overall model info
mGO_NSSI = bsxfun(@rdivide,mGO_NSSI,mean(mGO_NSSI(:,base_lags),2));

% mGO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
% mGO_ovLL = arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
% mGO_NLL = bsxfun(@rdivide,mGO_LL,mGO_ovLL); %normalize TB SSI by overall model info
% mGO_NLL = bsxfun(@rdivide,mGO_NLL,mean(mGO_NLL(:,base_lags),2));

%spline interpolate gsac and msac trig avgs
tlags_up = linspace(tlags(1),tlags(end),500);
all_gsac_tavg_up = spline(tlags,all_gsac_tavg,tlags_up);
all_msac_tavg_up = spline(tlags,all_msac_tavg,tlags_up);

%find enh and sup peaks in gsac and msac trig avgs
search_range = [0 0.3];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_tavg_up-1),tlags_up,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_tavg_up-1,tlags_up,search_range);
[msac_Ifact,msac_inhtime] = get_tavg_peaks(-(all_msac_tavg_up-1),tlags_up,search_range);
[msac_Efact,msac_exctime] = get_tavg_peaks(all_msac_tavg_up-1,tlags_up,search_range);

%find SSI suppression peaks
search_range = [0 0.15];
[gsac_SSI_Ifact] = get_tavg_peaks(-(GO_NSSI-1),slags*dt,search_range);
[msac_SSI_Ifact] = get_tavg_peaks(-(mGO_NSSI-1),slags*dt,search_range);

close all
%COMPARE RATES
xl = [-0.1 0.3];
f1 = figure();hold on
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
f2 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(GO_NSSI),nanstd(GO_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(mGO_NSSI),nanstd(mGO_NSSI)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative SSI');
ylim([0.5 1.1]);
% ylim([0.7 1.1]);

%COMPARE OFFSETS
xl = [-0.1 0.3];
f3 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(GO_Noffset),nanstd(GO_Noffset)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(mGO_Noffset),nanstd(mGO_Noffset)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative offset');
ylim([-0.05 0.3]);

%COMPARE GAINS
xl = [-0.1 0.3];
f4 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(GO_gain),nanstd(GO_gain)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(mGO_gain),nanstd(mGO_gain)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Gain');
ylim([0.5 1.2]);

% %PRINT PLOTS
fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_MSAC_rates.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% % fname = [fig_dir 'Gsac_MSAC_SSI.pdf'];
% fname = [fig_dir 'Gsac_MSAC_SSI_compare2_unCor.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
%
% figufy(f3);
% fname = [fig_dir 'Gsac_MSAC_OFFSET.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
%
% figufy(f4);
% fname = [fig_dir 'Gsac_MSAC_GAIN.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% COMPARE GRAY AND IMAGE BACKGROUNDS
cur_min_Nsacs = 250; %require half as many saccades in each background type condition
cur_SUs = find(avg_rates >= min_rate & N_gsacs_gray >= cur_min_Nsacs & N_gsacs_im >= cur_min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);

%reg hyperparam selection
GO_lambda_off = 4;
GO_lambda_gain = 3;

%gray back trig avg rate
all_gsac_gr_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_gray_avg', all_SU_data(cur_SUs),'uniformoutput',0));

%image back trig avg rate
all_gsac_im_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_im_avg', all_SU_data(cur_SUs),'uniformoutput',0));

%gray back GO-model SSI, normalized
GO_gr_SSI = cell2mat(arrayfun(@(x) x.type_dep.gsacGR_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_gr_ovinfos = arrayfun(@(x) x.type_dep.gsacGR_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_gr_NSSI = bsxfun(@rdivide,GO_gr_SSI,GO_gr_ovinfos); %normalize GO SSI by overall model info
GO_gr_NSSI = bsxfun(@rdivide,GO_gr_NSSI,mean(GO_gr_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info

%image back GO-model SSI, normalized
GO_im_SSI = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_im_ovinfos = arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_im_NSSI = bsxfun(@rdivide,GO_im_SSI,GO_im_ovinfos); %normalize GO SSI by overall model info
GO_im_NSSI = bsxfun(@rdivide,GO_im_NSSI,mean(GO_im_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info

%find SSI suppresion peaks for each condition
search_range = [0 0.15];
[GR_SSI_Ifact,GR_SSI_inhtime] = get_tavg_peaks(-(GO_gr_NSSI-1),slags*dt,search_range);
[IM_SSI_Ifact,IM_SSI_inhtime] = get_tavg_peaks(-(GO_im_NSSI-1),slags*dt,search_range);

%spline interpolate sac-trig-avgs (upsample)
tlags_up = linspace(tlags(1),tlags(end),500);
all_gsac_gr_tavg_up = spline(tlags,all_gsac_gr_tavg,tlags_up);
all_gsac_im_tavg_up = spline(tlags,all_gsac_im_tavg,tlags_up);

%find enh and sup timing and mag for each condition
search_range = [0 0.3];
[GR_Efact,GR_exctime] = get_tavg_peaks((all_gsac_gr_tavg_up-1),tlags_up,search_range);
[IM_Efact,IM_exctime] = get_tavg_peaks((all_gsac_im_tavg_up-1),tlags_up,search_range);
[GR_Ifact,GR_inhtime] = get_tavg_peaks(-(all_gsac_gr_tavg_up-1),tlags_up,search_range);
[IM_Ifact,IM_inhtime] = get_tavg_peaks(-(all_gsac_im_tavg_up-1),tlags_up,search_range);

%compute GO-model image back offset, normalized
IM_ov_rates = arrayfun(@(x) x.type_dep.gsac_IM_ovavg_rate,all_SU_data(cur_SUs));
IM_offset = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
IM_Noffset = bsxfun(@rdivide,IM_offset,IM_ov_rates);
IM_Noffset = bsxfun(@minus,IM_Noffset,mean(IM_Noffset(:,base_lags),2));

%compute GO-model gray back offset, normalized
GR_ov_rates = arrayfun(@(x) x.type_dep.gsac_GR_ovavg_rate,all_SU_data(cur_SUs));
GR_offset = cell2mat(arrayfun(@(x) x.type_dep.gsacGR_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
GR_Noffset = bsxfun(@rdivide,GR_offset,GR_ov_rates);
GR_Noffset = bsxfun(@minus,GR_Noffset,mean(GR_Noffset(:,base_lags),2));

%compute GO-model image back gain, normalized
IM_gain = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
IM_gain = bsxfun(@rdivide,IM_gain,mean(IM_gain(:,base_lags),2));

%compute GO-model gray back gain, normalized
GR_gain = cell2mat(arrayfun(@(x) x.type_dep.gsacGR_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
GR_gain = bsxfun(@rdivide,GR_gain,mean(GR_gain(:,base_lags),2));

%COMPARE TRIG AVG RATES
xl = [-0.1 0.3];
f1 = figure();hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_gr_tavg),nanstd(all_gsac_gr_tavg)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(all_gsac_im_tavg),nanstd(all_gsac_im_tavg)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.65 1.3])

%COMPARE SSI
xl = [-0.1 0.3];
f2 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(GO_gr_NSSI),nanstd(GO_gr_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(GO_im_NSSI),nanstd(GO_im_NSSI)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative SSI');
ylim([0.5 1.2])

%COMPARE OFFSET
xl = [-0.1 0.3];
f3 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(GR_Noffset),nanstd(GR_Noffset)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(IM_Noffset),nanstd(IM_Noffset)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative offset');
ylim([-0.1 0.3])

%COMPARE GAIN
xl = [-0.1 0.3];
f4 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(GR_gain),nanstd(GR_gain)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(IM_gain),nanstd(IM_gain)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Gain');
ylim([0.5 1.15])

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_GRIM_rates.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'Gsac_GRIM_SSI.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
%
% figufy(f3);
% fname = [fig_dir 'Gsac_GRIM_OFFSET.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
%
% figufy(f4);
% fname = [fig_dir 'Gsac_GRIM_GAIN.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% COMPARE REAL (WITH IMAGE BACK) AND SIMULATED SACCADES
cur_min_Nsacs = 250; %min nsacs (or sim sacs) required for each condition
cur_SUs = find(avg_rates >= min_rate & N_gsacs_im >= cur_min_Nsacs & N_simsacs >= cur_min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);

%reg param selection
GO_lambda_off = 4;
GO_lambda_gain = 3;

%sac trig avgs for each condition
all_gsac_im_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_im_avg', all_SU_data(cur_SUs),'uniformoutput',0));
all_simsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.simsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));

%compute GO-model SSI for real sacs, normalized
GO_im_SSI = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_im_ovinfos = arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_im_NSSI = bsxfun(@rdivide,GO_im_SSI,GO_im_ovinfos); %normalize GO SSI by overall model info
GO_im_NSSI = bsxfun(@rdivide,GO_im_NSSI,mean(GO_im_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info

%compute GO_model SSI for sim sacs, normalized
GO_sim_SSI = cell2mat(arrayfun(@(x) x.type_dep.simsac_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_sim_ovinfos = arrayfun(@(x) x.type_dep.simsac_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_sim_NSSI = bsxfun(@rdivide,GO_sim_SSI,GO_sim_ovinfos); %normalize GO SSI by overall model info
GO_sim_NSSI = bsxfun(@rdivide,GO_sim_NSSI,mean(GO_sim_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info

%compute SSI suppression peaks for each condition
search_range = [0 0.15];
[IM_SSI_Ifact,IM_SSI_inhtime] = get_tavg_peaks(-(GO_im_NSSI-1),slags*dt,search_range);
[SIM_SSI_Ifact,SIM_SSI_inhtime] = get_tavg_peaks(-(GO_sim_NSSI-1),slags*dt,search_range);

%spline interpolation for trig avg rates
tlags_up = linspace(tlags(1),tlags(end),500);
all_simsac_tavg_up = spline(tlags,all_simsac_tavg,tlags_up);
all_gsac_im_tavg_up = spline(tlags,all_gsac_im_tavg,tlags_up);

%find enh and sup peak timing and mag
search_range = [0 0.3];
[IM_Ifact,IM_inhtime] = get_tavg_peaks(-(all_gsac_im_tavg_up-1),tlags_up,search_range);
[SIM_Ifact,SIM_inhtime] = get_tavg_peaks(-(all_simsac_tavg_up-1),tlags_up,search_range);
[IM_Efact,IM_exctime] = get_tavg_peaks((all_gsac_im_tavg_up-1),tlags_up,search_range);
[SIM_Efact,SIM_exctime] = get_tavg_peaks((all_simsac_tavg_up-1),tlags_up,search_range);

%compute GO-model offsets for real sacs, normalized
IM_ov_rates = arrayfun(@(x) x.type_dep.gsac_IM_ovavg_rate,all_SU_data(cur_SUs));
IM_offset = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
IM_Noffset = bsxfun(@rdivide,IM_offset,IM_ov_rates);
IM_Noffset = bsxfun(@minus,IM_Noffset,mean(IM_Noffset(:,base_lags),2));

%compute GO-model offsets for sim sacs, normalized
SIM_ov_rates = arrayfun(@(x) x.type_dep.simsac_ovavg_rate,all_SU_data(cur_SUs));
SIM_offset = cell2mat(arrayfun(@(x) x.type_dep.simsac_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
SIM_Noffset = bsxfun(@rdivide,SIM_offset,SIM_ov_rates);
SIM_Noffset = bsxfun(@minus,SIM_Noffset,mean(SIM_Noffset(:,base_lags),2));

%compute GO-model gains for real sacs, normalized
IM_gain = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
IM_gain = bsxfun(@rdivide,IM_gain,mean(IM_gain(:,base_lags),2));

%compute GO-model gains for sim sacs, normalized
SIM_gain = cell2mat(arrayfun(@(x) x.type_dep.simsac_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
SIM_gain = bsxfun(@rdivide,SIM_gain,mean(SIM_gain(:,base_lags),2));

%COMPARE trig avg rates
xl = [-0.1 0.3];
f1 = figure();hold on
h1=shadedErrorBar(tlags,nanmean(all_simsac_tavg),nanstd(all_simsac_tavg)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(all_gsac_im_tavg),nanstd(all_gsac_im_tavg)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.65 1.3])

%COMPARE SSI
xl = [-0.1 0.3];
f2 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(GO_sim_NSSI),nanstd(GO_sim_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(GO_im_NSSI),nanstd(GO_im_NSSI)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative SSI');
ylim([0.5 1.2])

%compare offsets
f3 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(IM_Noffset),nanstd(IM_Noffset)/sqrt(length(cur_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,nanmean(SIM_Noffset),nanstd(SIM_Noffset)/sqrt(length(cur_SUs)),{'color','b'});
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative Offset');
ylim([-0.1 0.3]);

%compare gains
f4 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(IM_gain),nanstd(IM_gain)/sqrt(length(cur_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,nanmean(SIM_gain),nanstd(SIM_gain)/sqrt(length(cur_SUs)),{'color','b'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Gain');
ylim([0.5 1.15]);

% % PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_SIM_rates.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'Gsac_SIM_SSI.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
%
% figufy(f3);
% fname = [fig_dir 'Gsac_SIM_OFFSET.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
%
% figufy(f4);
% fname = [fig_dir 'Gsac_SIM_GAIN.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% COMPARE PRE AND POST MODELS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= -0.025); %use slightly earlier def of backgnd time points since the upstream kernels start slightly earlier
close all

%hyperparameter selection
post_lambda_off = 4;
pre_lambda = 3; %this is just the gain kernel lambda, the offset kernel lambda is fixed equal to post-model (4)
post_lambda_gain = 3;

%pre and post gain kernels
gsac_pre_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.stim_kernel',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_post_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{post_lambda_off,post_lambda_gain}.mods(3).filtK',all_SU_data(cur_SUs),'uniformoutput',0));

%noramlize gain kernels
gsac_post_gain = bsxfun(@rdivide,gsac_post_gain,mean(gsac_post_gain(:,base_lags),2));
gsac_pre_gain = bsxfun(@rdivide,gsac_pre_gain,mean(gsac_pre_gain(:,base_lags),2));

%compute pre and post offset kernels
gsac_pre_off =cell2mat(arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.off_kernel',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_post_off =cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{post_lambda_off,post_lambda_gain}.mods(2).filtK',all_SU_data(cur_SUs),'uniformoutput',0));

%normalize offset kernels
gsac_pre_off = bsxfun(@minus,gsac_pre_off,mean(gsac_pre_off(:,base_lags),2));
gsac_post_off = bsxfun(@minus,gsac_post_off,mean(gsac_post_off(:,base_lags),2));

%compute LLs of different models
base_LLs = arrayfun(@(x) x.sacStimProc.gsac_base_LLimp,all_SU_data(cur_SUs));
off_LLs = arrayfun(@(x) x.sacStimProc.gsac_off_mod{post_lambda_off}.ovLLimp,all_SU_data(cur_SUs));
post_LLs = arrayfun(@(x) x.sacStimProc.gsac_post_mod{post_lambda_off,post_lambda_gain}.ovLLimp,all_SU_data(cur_SUs));
pre_LLs = arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.ovLLimp,all_SU_data(cur_SUs));

%COMPARE PRE AND POST GAINS
xl = [-0.1 0.2];
zxl = [-0.05 0.15];
f1 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(gsac_pre_gain),nanstd(gsac_pre_gain)/sqrt(length(cur_SUs)),{'color','k'});
h2=shadedErrorBar(slags*dt,nanmean(gsac_post_gain),nanstd(gsac_post_gain)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Gain');
ylim([0.55 1.05])

xl = [-0.1 0.3];
f2 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(gsac_pre_off),nanstd(gsac_pre_off)/sqrt(length(cur_SUs)),{'color','k'});
h2=shadedErrorBar(slags*dt,nanmean(gsac_post_off),nanstd(gsac_post_off)/sqrt(length(cur_SUs)),{'color','b'});
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Offset');

%rel LL improvement of pre model over post model (normalized by improvement
%of post-model over a sac-mod null model)
rel_LLimps = (pre_LLs - post_LLs)./(post_LLs - off_LLs);
f3 = figure();
hist(rel_LLimps,15);
h = findobj(gca,'type','patch');
set(h,'facecolor','k','edgecolor','w');
xlim([-0.9 0.9]);
% boxplot(rel_LLimps);
% hold on

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_pre_post_GAIN.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% figure(f1);
% % xlim(zxl);
% % fname = [fig_dir 'Gsac_pre_post_GAIN_zoom.pdf'];
% % exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

% % figufy(f2);
% % fname = [fig_dir 'Gsac_pre_post_OFFSET.pdf'];
% % exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f2);
% %
% figufy(f3);
% fname = [fig_dir 'Gsac_pre_post_LLhist.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

%% COMPARE GSAC AND MSAC PRE MODELS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & N_msacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= -0.025); %again, use slightly earlier definition of backgnd time points to handle upstream gain kernel timing

%reg param selection
post_lambda_off = 4;
pre_lambda = 3;
post_lambda_gain = 3;

%gsac and msac upstream gain kernels
gsac_pre_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.stim_kernel',all_SU_data(cur_SUs),'uniformoutput',0));
msac_pre_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.msacPreGainMod{pre_lambda}.stim_kernel',all_SU_data(cur_SUs),'uniformoutput',0));

%normalize upstream gain kernels
msac_pre_gain = bsxfun(@rdivide,msac_pre_gain,mean(msac_pre_gain(:,base_lags),2));
gsac_pre_gain = bsxfun(@rdivide,gsac_pre_gain,mean(gsac_pre_gain(:,base_lags),2));

%get sac model LLs
gsac_base_LLs = arrayfun(@(x) x.sacStimProc.gsac_base_LLimp,all_SU_data(cur_SUs));
gsac_off_LLs = arrayfun(@(x) x.sacStimProc.gsac_off_mod{post_lambda_off}.ovLLimp,all_SU_data(cur_SUs));
gsac_post_LLs = arrayfun(@(x) x.sacStimProc.gsac_post_mod{post_lambda_off,post_lambda_gain}.ovLLimp,all_SU_data(cur_SUs));
gsac_pre_LLs = arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.ovLLimp,all_SU_data(cur_SUs));

%get msac model LLs
msac_base_LLs = arrayfun(@(x) x.sacStimProc.msac_base_LLimp,all_SU_data(cur_SUs));
msac_off_LLs = arrayfun(@(x) x.sacStimProc.msac_off_mod{post_lambda_off}.ovLLimp,all_SU_data(cur_SUs));
msac_post_LLs = arrayfun(@(x) x.sacStimProc.msac_post_mod{post_lambda_off,post_lambda_gain}.ovLLimp,all_SU_data(cur_SUs));
msac_pre_LLs = arrayfun(@(x) x.sacStimProc.msacPreGainMod{pre_lambda}.ovLLimp,all_SU_data(cur_SUs));

%compute relative LL improvements
gsac_rel_LLimps = (gsac_pre_LLs - gsac_post_LLs)./(gsac_post_LLs - gsac_off_LLs);
msac_rel_LLimps = (msac_pre_LLs - msac_post_LLs)./(msac_post_LLs - msac_off_LLs);

%COMPARE PRE AND POST GAINS
xl = [-0.1 0.2];
xlz = [-0.05 0.1];
f1 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(gsac_pre_gain),nanstd(gsac_pre_gain)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(msac_pre_gain),nanstd(msac_pre_gain)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Gain');

% % %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_msac_pre_GAIN.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% figure(f1);
% xlim(xlz);
% fname = [fig_dir 'Gsac_msac_pre_GAIN_zoom.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% COMPARE STIM TIMING AND SAC-MOD TIMING
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);

%get the E and I-filter temporal kernels (hilbert based) for each units
%model
all_Ekerns = [];
all_Ikerns = [];
for ii = 1:length(cur_SUs)
    [~,cur_Ekern,cur_Ikern] = get_hilbert_tempkerns(all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM);
    all_Ekerns = cat(1,all_Ekerns,cur_Ekern');
    all_Ikerns = cat(1,all_Ikerns,cur_Ikern');
end

%sac trig avgs, spline interpolated for up-sampling
all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
tlags_up = linspace(tlags(1),tlags(end),500);
all_gsac_up = spline(tlags,all_gsac_tavg,tlags_up);

%find E and I sta peaks
search_range = [0 0.3];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_up-1),tlags_up,search_range);
search_range = [0 0.3];
[gsac_Efact,gsac_exctime] = get_tavg_peaks((all_gsac_up-1),tlags_up,search_range);

%only use units where rate suppression comes before rate enhancement
include = find(gsac_inhtime < gsac_exctime);
opp_pol = find(gsac_exctime <= gsac_inhtime);

% [earlier_timing,peak_earlier] = min([gsac_inhtime gsac_exctime],[],2);

%stim filter time axis
flen = 15;
lag_ax = ((1:flen)*dt - dt/2);
up_lagax = linspace(lag_ax(1),lag_ax(end),500);

%up-sample stim filter temporal kernels
all_Ekerns_up = spline(lag_ax,all_Ekerns,up_lagax);
all_Ikerns_up = spline(lag_ax,all_Ikerns,up_lagax);

%fine timing of temporal kernel peaks
search_range = [0 max(up_lagax)];
[Ekern_max,Ekern_time] = get_tavg_peaks(all_Ekerns_up,up_lagax,search_range);
[Ikern_max,Ikern_time] = get_tavg_peaks(all_Ikerns_up,up_lagax,search_range);

% cur_sigI = find(ismember(cur_SUs,sigI));
% cur_sigE = find(ismember(cur_SUs,sigE));
% cur_rP = find(ismember(cur_SUs,reverse_polarity_units));
% cur_O = false(size(cur_SUs));
% cur_O(cur_sigI(peak_earlier(cur_sigI) == 1)) = true;
% cur_O(cur_sigE(peak_earlier(cur_sigE) == 1)) = true;

%plot E-filter temp-kernel timing vs sac-suppression timing
xl = [0.02 0.08];
yl = [0.02 0.16];
mS = 3;
f1 = figure(); hold on
plot(Ekern_time(include),gsac_inhtime(include),'b.','markersize',12);
% plot(Ekern_time(opp_pol),gsac_exctime(opp_pol),'ro');
xlim(xl); ylim(yl);
% r = robustfit(Ekern_time,gsac_inhtime);
r = robustfit(Ekern_time(include),gsac_inhtime(include));
xx = linspace(xl(1),xl(2),50);
plot(xx,r(1)+r(2)*xx,'r')
xlabel('Stim latency (s)');
ylabel('Suppression timing (s)');
line([0 0.15],[0 0.15],'color','k');

SU_exptnums = [all_SU_data(cur_SUs).expt_num];
[H,ATAB,CTAB,STATS]=aoctool(Ekern_time(include),gsac_inhtime(include),SU_exptnums(include)',0.05,'stim-response','sac-sup','rec_num','off',4);


% f2 = figure(); hold on
% plot(Ekern_time,gsac_exctime,'.')
% plot(Ekern_time(cur_sigE),gsac_exctime(cur_sigE),'r.');
%
% f3 = figure(); hold on
% plot(Ikern_time,gsac_exctime,'.')
% plot(Ikern_time(cur_sigE),gsac_exctime(cur_sigE),'r.');

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Supp_Stim_timing_scatter_new.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% COMPARE TEMPKERNS AND MODEL GAIN LATENCIES
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);
noise_lags = find(slags <= 0 | slags*dt > 0.15); %time lags used to define 'background SD'
% noise_lags = find(slags <= 0); %time lags used to define 'background SD'
close all

%set regularization params
lambda_L2_ii = 1;
lambda_d2T_ii = 3;
GO_lambda_off = 4;
GO_lambda_gain = 3;

%create stim-filter taxis
flen = 15;
lag_ax = ((1:flen)*dt - dt/2);
up_lagax = linspace(lag_ax(1),lag_ax(end),500);
slags_up = linspace(slags(1)*dt,slags(end)*dt,500); %sta up-sampled taxis
jit_amp = 0.000;%jitter data points with gaussian RV to avoid point occlusion

search_range = [0 0.15]; %range of tlags to search for peaks in sac kernels
tsearch_range = [0 max(up_lagax)];
noise_SD_thresh = 2; %min peak magnitude beyond background SD to include in analysis

all_tempkerns = [];
all_gainkerns = [];
all_relweights = [];
all_modsigns = [];
cell_tg_slope = nan(length(cur_SUs),1);
for ii = 1:length(cur_SUs)
    %get stim-filter temporal kernels
    cur_tkerns = get_hilbert_tempkerns(all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM);
    %get gain kernels from full model
    cur_gainkerns = reshape([all_SU_data(cur_SUs(ii)).sacStimProc.gsac_post_Fullmod{1,lambda_d2T_ii,lambda_L2_ii}.mods(3).filtK],length(slags),[])';
    cur_relweights = all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM.rel_filt_weights; %relative weight of each stim filter
    cur_modsigns = [all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM.mods(:).sign]; %sign of each stim filter
    
    clear cur_tkerns_up
    uset = find(cur_relweights > 0); %set of filters with non-zero magnitude
    cur_tkerns_up(uset,:) = spline(lag_ax,cur_tkerns(:,uset)',up_lagax); %spline interpolate temporal kernels
    [~,cur_tk_time] = get_tavg_peaks(cur_tkerns_up,up_lagax,tsearch_range); %get peaks of temporal kernels
    cur_gkerns_up = spline(slags*dt,cur_gainkerns,slags_up); %spline interpolate gain kernels
    [cur_gkern_max,cur_gk_time] = get_tavg_peaks(-(cur_gkerns_up-1),slags_up,search_range); %get peaks of gain kernels
    cur_noise_level = std(cur_gainkerns(:,noise_lags),[],2); %compute 'background noise' of gain kernels
    uset = uset(cur_gkern_max(uset) >= cur_noise_level(uset)*noise_SD_thresh); %only take filters with sufficient peak SNR
    
    %analyze within-cell relationship between filter and gain kernel timing
    if length(uset) >= 5
        %         r = robustfit(cur_tk_time(uset),cur_gk_time(uset));
        %         cell_tg_slope(ii) = r(2);
        cell_tg_slope(ii) = corr(cur_tk_time(uset),cur_gk_time(uset),'type','spearman');
    end
    
    all_tempkerns = cat(1,all_tempkerns,cur_tkerns');
    all_gainkerns = cat(1,all_gainkerns,1+cur_gainkerns);
    all_relweights = cat(1,all_relweights,cur_relweights');
    all_modsigns = cat(1,all_modsigns,cur_modsigns');
end

%normalize gain kernels
all_gainkerns = bsxfun(@rdivide,all_gainkerns,mean(all_gainkerns(:,base_lags),2));

%spline interp all temporal kernels and find peaks
all_tkerns_up = nan(length(all_relweights),length(up_lagax));
nzero_filts = find(all_relweights > 0);
all_tkerns_up(nzero_filts,:) = spline(lag_ax,all_tempkerns(nzero_filts,:),up_lagax);
[tkern_max,tkern_time] = get_tavg_peaks(all_tkerns_up,up_lagax,tsearch_range);

%spline interpolate all gain kernels and find peaks
all_gkerns_up = spline(slags*dt,all_gainkerns,slags_up);
[gkern_max,gkern_time] = get_tavg_peaks(-(all_gkerns_up-1),slags_up,search_range);

%find gain kernels with sufficient SNR of peaks
noise_level =std(all_gainkerns(:,noise_lags),[],2);
really_poss_ukerns = find(all_relweights > 0 & ~isnan(tkern_max));
poss_ukerns = find(all_relweights > 0 & ~isnan(tkern_max) & ~isnan(gkern_max));
ukerns = find(all_relweights > 0 & gkern_max > noise_SD_thresh*noise_level & ~isnan(tkern_max)); %only use subunits with nonzero stim filters, and significant gain kernel peaks
[a,b] = corr(gkern_time(ukerns),tkern_time(ukerns),'type','spearman') %compute correlation between gain kernel suppression timing and stimulus filter temporal kernel timing across usable subunits

%usabl E and I subunits
esubs = ukerns(all_modsigns(ukerns)==1);
isubs = ukerns(all_modsigns(ukerns)==-1);

%if using jittering to prevent occlusion
if jit_amp > 0
    gkern_time = gkern_time + randn(size(gkern_time))*jit_amp;
    tkern_time = tkern_time + randn(size(tkern_time))*jit_amp;
end

mS = 8;
f1 = figure(); hold on
% plot(gkern_time(ukerns), tkern_time(ukerns),'.')
plot(gkern_time(esubs), tkern_time(esubs) ,'k.','markersize',mS)
plot(gkern_time(isubs) , tkern_time(isubs),'r.','markersize',mS)
r_ov = robustfit(gkern_time(ukerns),tkern_time(ukerns));
xx = linspace(0,0.2,50);
plot(xx,r_ov(1) + r_ov(2)*xx,'k');
xlim([0 0.15]);
xlabel('Gain suppression timing (s)');
ylabel('Temporal response timing (s)');

%separate correlations for E and I subunits
[a,b] = corr(gkern_time(esubs),tkern_time(esubs),'type','spearman')
[a,b] = corr(gkern_time(isubs),tkern_time(isubs),'type','spearman')


% xl = [-0.1 0.3];
% f2 = figure();hold on
% shadedErrorBar(slags*dt,mean(all_gainkerns),std(all_gainkerns)/sqrt(size(all_gainkerns,1)),{'color','k'});
% line(xl,[1 1],'color','k');
%
% GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
% GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
% GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info
% GO_NSSI = bsxfun(@rdivide,GO_NSSI,mean(GO_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info
%
% FULL_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_Fullmod{1,lambda_d2T_ii,lambda_L2_ii}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
% FULL_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_Fullmod{1,lambda_d2T_ii,lambda_L2_ii}.ovInfo,all_SU_data(cur_SUs));
% FULL_NSSI = bsxfun(@rdivide,FULL_SSI,FULL_ovinfos); %normalize GO SSI by overall model info
% FULL_NSSI = bsxfun(@rdivide,FULL_NSSI,mean(FULL_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info
%
% GO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
% GO_ovLL = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
% GO_NLL = bsxfun(@rdivide,GO_LL,GO_ovLL); %normalize TB SSI by overall model info
% GO_NLL = bsxfun(@rdivide,GO_NLL,mean(GO_NLL(:,base_lags),2)); %normalize GO SSI by overall model info
%
% FULL_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_Fullmod{1,lambda_d2T_ii,lambda_L2_ii}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
% FULL_ovLL = arrayfun(@(x) x.sacStimProc.gsac_post_Fullmod{1,lambda_d2T_ii,lambda_L2_ii}.ovLLimp,all_SU_data(cur_SUs));
% FULL_NLL = bsxfun(@rdivide,FULL_LL,FULL_ovLL); %normalize GO SSI by overall model info
% FULL_NLL = bsxfun(@rdivide,FULL_NLL,mean(FULL_NLL(:,base_lags),2)); %normalize GO SSI by overall model info
%
% for ii = 1:length(cur_SUs)
%     GO_NLL(ii,:) = jmm_smooth_1d_cor(GO_NLL(ii,:),1);
%     FULL_NLL(ii,:) = jmm_smooth_1d_cor(FULL_NLL(ii,:),1);
% end
%
% f3 = figure();hold on
% shadedErrorBar(slags*dt,mean(GO_NSSI),std(GO_NSSI)/sqrt(length(cur_SUs)),{'color','k'});
% shadedErrorBar(slags*dt,mean(FULL_NSSI),std(FULL_NSSI)/sqrt(size(cur_SUs,1)),{'color','r'});
% line(xl,[1 1],'color','k');
%
% f4 = figure();hold on
% shadedErrorBar(slags*dt,mean(GO_NLL),std(GO_NLL)/sqrt(length(cur_SUs)),{'color','k'});
% shadedErrorBar(slags*dt,mean(FULL_NLL),std(FULL_NLL)/sqrt(size(cur_SUs,1)),{'color','r'});

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 1;
% figufy(f1);
% fname = [fig_dir 'Full_mod_timing_scatter.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'Full_mod_avg_gkern.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
%

%% COMPARE E-I Tempkerns and Gain kerns
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);

%get up-sampled temporal filter t-axis
flen = 15;
lag_ax = ((1:flen)*dt - dt/2)*1e3;
up_lagax = linspace(lag_ax(1),lag_ax(end),500);
slags_up = linspace(slags(1)*dt,slags(end)*dt,500);

base_lags = find(slags <= 0);
noise_lags = find(slags <= 0 | slags*dt > 0.15); %times defining gain kernel background 'noise level'
noise_SD_thresh = 0; %minimum peak amplitude for gain kernels above background
search_range = [0 0.15];
tsearch_range = up_lagax([1 end]);

lambda_ii = 3; %reg parameter selection

%pull out E and I gain kernels for each unit
gsac_EIgain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_EImod{lambda_ii}.mods(3).filtK',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_Egain = gsac_EIgain(:,1:length(slags));
gsac_Igain = gsac_EIgain(:,(length(slags)+1):end);

%normalize E and I gain kernels
gsac_Egain = bsxfun(@rdivide,gsac_Egain,mean(gsac_Egain(:,base_lags),2));
gsac_Igain = bsxfun(@rdivide,gsac_Igain,mean(gsac_Igain(:,base_lags),2));

all_Ekerns = nan(length(cur_SUs),flen);
all_Ikerns = nan(length(cur_SUs),flen);
all_NIfilts = nan(length(cur_SUs),1);
all_NEfilts = nan(length(cur_SUs),1);
for ii = 1:length(cur_SUs)
    %get E and I stim filter temporal kernels
    [cur_tkerns,cur_ekern,cur_ikern] = get_hilbert_tempkerns(all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM);
    all_Ekerns(ii,:) = cur_ekern;
    all_Ikerns(ii,:) = cur_ikern;
    
    %keep count of the total number of E and I filters for each cell
    cur_modsigns = [all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM.mods(:).sign];
    all_NIfilts(ii) = sum(cur_modsigns==-1);
    all_NEfilts(ii) = sum(cur_modsigns==1);
end

%spline interpolate stimulus filter temporal kernels and find peaks
all_Ekerns_up = spline(lag_ax,all_Ekerns,up_lagax);
all_Ikerns_up = spline(lag_ax,all_Ikerns,up_lagax);
[Ekern_max,Ekern_time] = get_tavg_peaks(all_Ekerns_up,up_lagax,tsearch_range);
[Ikern_max,Ikern_time] = get_tavg_peaks(all_Ikerns_up,up_lagax,tsearch_range);

%interpolate gain kernels and find peaks
gsac_Egain_up = spline(slags*dt,gsac_Egain,slags_up);
gsac_Igain_up = spline(slags*dt,gsac_Igain,slags_up);
[Egain_max,Egain_time] = get_tavg_peaks(-(gsac_Egain_up-1),slags_up,search_range);
[Igain_max,Igain_time] = get_tavg_peaks(-(gsac_Igain_up-1),slags_up,search_range);

%find gain kernel noise levels for E and I filters
Enoise_level =std(gsac_Egain(:,noise_lags),[],2);
Inoise_level =std(gsac_Igain(:,noise_lags),[],2);

%only use E and I gain kernels where peak is sufficiently high above bckgnd
%noise level
use_E = find(Egain_max >= noise_SD_thresh*Enoise_level);
use_I = find(Igain_max >= noise_SD_thresh*Inoise_level);
use_B = intersect(use_E,use_I);


%COMPARE PRE AND POST GAINS
xl = [-0.1 0.3];
xlz = [0 0.15];
f1 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(gsac_Egain),nanstd(gsac_Egain)/sqrt(length(cur_SUs)),{'color',[0.2 0.8 0.2]});
h2=shadedErrorBar(slags*dt,nanmean(gsac_Igain),nanstd(gsac_Igain)/sqrt(length(cur_SUs)),{'color',[0.8 0.1 0.8]});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative rate');

%COMPARE PRE AND POST GAINS
f2 = figure();hold on
h1=shadedErrorBar(up_lagax,nanmean(all_Ekerns_up),nanstd(all_Ekerns_up)/sqrt(length(cur_SUs)),{'color',[0.2 0.8 0.2]});
h2=shadedErrorBar(up_lagax,nanmean(all_Ikerns_up),nanstd(all_Ikerns_up)/sqrt(length(cur_SUs)),{'color',[0.8 0.1 0.8]});
xlabel('Time (s)');

mS = 12;
f3 = figure(); hold on
plot(Ekern_time,Ikern_time,'k.','markersize',mS);
plot(Egain_time*1e3,Igain_time*1e3,'r.','markersize',mS);
line([0 150],[0 150],'color','k');
xlabel('Excitatory inputs');
ylabel('Inhibitory inputs');

% f4 = figure(); hold on
% plot(Ekern_time,Egain_time*1e3,'.','markersize',mS);
% plot(Ikern_time,Igain_time*1e3,'r.','markersize',mS);
% line([0 150],[0 150],'color','k');
% xlabel('Stim response');
% ylabel('Sac suppression');

f5 = figure(); hold on
xx = linspace(-75,75,25);
gain_diff_hist = histc(1e3*(Egain_time-Igain_time),xx);
gain_diff_hist = gain_diff_hist/sum(gain_diff_hist);
kern_diff_hist = histc((Ekern_time-Ikern_time),xx);
kern_diff_hist = kern_diff_hist/sum(kern_diff_hist);
stairs(xx,gain_diff_hist,'k');
stairs(xx,kern_diff_hist,'r');
xlim([-75 75]);

% % % %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'EI_gainkerns.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% figure(f1);
% xlim(xlz);
% fname = [fig_dir 'EI_gainkerns_zoom.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'EI_tempkerns.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
%
% figufy(f3);
% fname = [fig_dir 'EI_timing_scatter.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
%
% figufy(f5);
% fname = [fig_dir 'EI_timing_dists.pdf'];
% exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f5);

%% ANALYZE LAMINAR DEPENDENCIES WITH MUA

%get stim filter taxis
mod_dt = 0.01; flen = 15;
tax = (0:(flen-1))*mod_dt + mod_dt/2;
up_tax = linspace(tax(1),tax(end),500);
tsearch_range = up_tax([1 end]);

search_range = [0 0.3]; %search range for STAs
xl = [-0 0.25]; %xrange for plotting

mua_sm = 0.005/trig_avg_params.dt; %smoothing sigma for MUA rates

%load info about layer boundaries
load('~/Analysis/bruce/FINsac_mod/layer_boundaries/layer_classification.mat')
boundary_enums = [boundary_class(:).Expt_num];
all_lbs = [boundary_class(:).lb]; %location of lower L4 bound
all_ubs = [boundary_class(:).ub]; %location of upper L4 bound
mean_lb = mean(all_lbs); %across rec avg boundary loc
mean_ub = mean(all_ubs);
mean_diff = mean_lb-mean_ub; %avg size of L4
interp_ax = 1:1:24; %axis for interpolating depth profiles

%load in all MUA models and compute excitatory temporal kernels
base_mname = 'corrected_models2';
flen = 15;n_probes = 24;
Elist = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};
Olist = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
n_lem_expts = length(Elist);
expt_MUA_Ekerns = nan(2,n_lem_expts,n_probes,flen);
expt_MUA_stas = nan(2,n_lem_expts,n_probes,length(tlags));
for ii = 1:n_lem_expts
    Expt_name = Elist{ii};
    mod_dir = ['~/Analysis/bruce/' Expt_name '/models/'];
    sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    for kk = 1:2
        cur_ori = Olist(ii,kk);
        if ~isnan(cur_ori)
            sname = strcat(mod_dir,base_mname,sprintf('_ori%d',cur_ori));
            load(sname); %load in models
            tname = strcat(sac_dir,base_tname,sprintf('_ori%d',cur_ori));
            load(tname); %load in trig avg data
            expt_MUA_stas(kk,ii,:,:) = mua_data.gsac_avg'; %store MUA sac-trig avgs
            
            %cycle over MUA probes and compute the Ekern temporal profile
            %for each
            for jj = 1:n_probes
                cur_mod = ModData(jj).rectGQM;
                [avg_tkerns,avg_Ekerns] = get_hilbert_tempkerns(cur_mod);
                expt_MUA_Ekerns(kk,ii,jj,:) = avg_Ekerns;
            end
        end
    end
end
%average across multiple stim oris for each expt
expt_MUA_Ekerns = squeeze(nanmean(expt_MUA_Ekerns));
expt_MUA_stas = squeeze(nanmean(expt_MUA_stas));

%slight smoothing of MUA trig avg rates
for ii = 1:size(expt_MUA_stas,1)
    for jj = 1:size(expt_MUA_stas,2)
        expt_MUA_stas(ii,jj,:) = jmm_smooth_1d_cor(squeeze(expt_MUA_stas(ii,jj,:)),mua_sm);
    end
end

SU_probenums = arrayfun(@(x) x.sacStimProc.ModData.unit_data.probe_number,all_SU_data);
SU_exptnums = [all_SU_data(:).expt_num];

%loop over expts and assign probes to layer classes
all_SU_lclass = nan(length(all_SU_data),1); %layer classififcation for each SU
all_MU_lclass = nan(n_lem_expts,n_probes); %layer classification for each MU
all_MU_recID = nan(n_lem_expts,n_probes); %recording ID for each MU
% all_interp_mua_sta = nan(n_lem_expts,length(interp_ax),length(tlags));
% all_interp_mua_Ekern = nan(n_lem_expts,length(interp_ax),length(up_tax));
all_mua_Ekerns = [];
for ee = 1:n_lem_expts
    %boundary info for current expt
    cur_enum = str2num(Elist{ee}(2:end));
    cur_bound_info = find(boundary_enums == cur_enum,1);
    cur_ub = boundary_class(cur_bound_info).ub;
    cur_lb = boundary_class(cur_bound_info).lb;
    
    %data needed for depth profile interpolation
    cur_diff = cur_lb - cur_ub;
    cur_slope = cur_diff/mean_diff;
    cur_off = mean_ub - cur_ub*cur_slope;
    cur_ax = (1:24)*cur_slope + cur_off;
    
    %     cur_Ekern = squeeze(expt_MUA_Ekerns(ee,:,:));
    %     cur_Ekern = spline(tax,cur_Ekern,up_tax);
    %     all_interp_mua_sta(ee,:,:) = interp1(cur_ax,squeeze(expt_MUA_stas(ee,:,:)),interp_ax);
    %     all_interp_mua_Ekern(ee,:,:) = interp1(cur_ax,cur_Ekern,interp_ax);
    
    %probe layer ids
    gran_probes = (cur_ub+1):(cur_lb-1);
    supra_probes = 1:(cur_ub-1);
    infra_probes = (cur_lb+1):24;
    
    %get layer classification for each MU probe
    cur_pclass = nan(24,1);
    cur_pclass(supra_probes) = 1; cur_pclass(gran_probes) = 2; cur_pclass(infra_probes) = 3;
    
    all_MU_lclass(ee,:) = cur_pclass;
    all_MU_recID(ee,:) = ee;
    
    %find SUs from this expt and classify their layer
    cur_SU_set = find(SU_exptnums == cur_enum);
    cur_SU_probenums = SU_probenums(cur_SU_set);
    cur_SU_lclass = nan(length(cur_SU_set),1);
    cur_SU_lclass(ismember(cur_SU_probenums,supra_probes)) = 1;
    cur_SU_lclass(ismember(cur_SU_probenums,gran_probes)) = 2;
    cur_SU_lclass(ismember(cur_SU_probenums,infra_probes)) = 3;
    all_SU_lclass(cur_SU_set) = cur_SU_lclass;
end
all_MU_lclass = reshape(all_MU_lclass,[],1);
all_MU_recID = reshape(all_MU_recID,[],1);

expt_MUA_Ekerns = reshape(expt_MUA_Ekerns,[],flen);
expt_MUA_stas = reshape(expt_MUA_stas,[],length(tlags));

%spline interpolation for temporal kernels
expt_MUA_Ekerns = spline(tax,expt_MUA_Ekerns,up_tax);

supra = find(all_MU_lclass == 1);
gran = find(all_MU_lclass == 2);
infra = find(all_MU_lclass == 3);
f1 = figure();hold on
h1 = shadedErrorBar(up_tax,mean(expt_MUA_Ekerns(supra,:)),std(expt_MUA_Ekerns(supra,:))/sqrt(length(supra)),{'color','b'});
h2 = shadedErrorBar(up_tax,mean(expt_MUA_Ekerns(gran,:)),std(expt_MUA_Ekerns(gran,:))/sqrt(length(gran)),{'color','r'});
h3 = shadedErrorBar(up_tax,mean(expt_MUA_Ekerns(infra,:)),std(expt_MUA_Ekerns(infra,:))/sqrt(length(infra)),{'color','k'});
xlim([0 0.15]);

[MUA_Tpeak,MUA_Ttime] = get_tavg_peaks(expt_MUA_Ekerns,up_tax,tsearch_range);

%spline interpolate MUA stas
tlags_up = linspace(tlags(1),tlags(end),500);
expt_MUA_stas_up = spline(tlags,expt_MUA_stas,tlags_up);
% [MUA_Efact,MUA_Etime] = get_tavg_peaks((expt_MUA_stas-1),tlags,search_range);
% [MUA_Ifact,MUA_Itime] = get_tavg_peaks(-(expt_MUA_stas-1),tlags,search_range);
[MUA_Efact,MUA_Etime] = get_tavg_peaks((expt_MUA_stas_up-1),tlags_up,search_range);
[MUA_Ifact,MUA_Itime] = get_tavg_peaks(-(expt_MUA_stas_up-1),tlags_up,search_range);

% f2 = figure();hold on
% h1 = shadedErrorBar(tlags,mean(expt_MUA_stas(supra,:)),std(expt_MUA_stas(supra,:))/sqrt(length(supra)),{'color','b'});
% h2 = shadedErrorBar(tlags,mean(expt_MUA_stas(gran,:)),std(expt_MUA_stas(gran,:))/sqrt(length(gran)),{'color','r'});
% h3 = shadedErrorBar(tlags,mean(expt_MUA_stas(infra,:)),std(expt_MUA_stas(infra,:))/sqrt(length(infra)),{'color','k'});
% xlim(xl);
% ylim([0.7 1.3]);
% line(xl,[1 1],'color','k');
%
% f3 = figure();
% subplot(1,3,1);
% bounds = prctile(tempdata,[10 90]);
% tempdata(tempdata < bounds(1)) = bounds(1); tempdata(tempdata > bounds(2)) = bounds(2);
% boxplot_capped(MUA_Itime*1e3,all_MU_lclass,[10 90]);
% % ylim([40 95])
% subplot(1,3,2);
% boxplot_capped(MUA_Etime*1e3,all_MU_lclass,[10 90]);
% subplot(1,3,3);
% boxplot_capped(MUA_Ttime*1e3,all_MU_lclass,[10 90]);
%
%test whether there are layer dependent differences in suppression timing.
%Use a two-way anova with random effects on recording ID.
group = {all_MU_lclass all_MU_recID};
[p,table,stats] = anovan(MUA_Itime,group,'random',2);
[p,table,stats] = anovan(MUA_Ttime,group,'random',2);
% [p,table,stats] = anovan(MUA_Ifact,group,'random',2);
[p,table,stats] = anovan(MUA_Etime,group,'random',2);
% [p,table,stats] = anovan(MUA_Efact,group,'random',2);

% all_Ekerns = nan(length(all_SU_data),flen);
% for ii = 1:length(all_SU_data)
%     [cur_tkerns,cur_ekern,cur_ikern] = get_hilbert_tempkerns(all_SU_data(ii).sacStimProc.ModData.rectGQM);
%     all_Ekerns(ii,:) = cur_ekern;
% end
% all_Ekerns_up = spline(tax,all_Ekerns,up_tax);
%
% cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_TA_Nsacs);
% supra_SUs = cur_SUs(all_SU_lclass(cur_SUs) == 1);
% gran_SUs = cur_SUs(all_SU_lclass(cur_SUs) == 2);
% infra_SUs = cur_SUs(all_SU_lclass(cur_SUs) == 3);
%
% f3 = figure();hold on
% h1 = shadedErrorBar(up_tax,mean(all_Ekerns_up(supra_SUs,:)),std(all_Ekerns_up(supra_SUs,:))/sqrt(length(supra_SUs)),{'color','b'});
% h2 = shadedErrorBar(up_tax,mean(all_Ekerns_up(gran_SUs,:)),std(all_Ekerns_up(gran_SUs,:))/sqrt(length(gran_SUs)),{'color','r'});
% h3 = shadedErrorBar(up_tax,mean(all_Ekerns_up(infra_SUs,:)),std(all_Ekerns_up(infra_SUs,:))/sqrt(length(infra_SUs)),{'color','k'});
%
% SU_stas = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(:),'uniformoutput',0));
% search_range = [0 0.3];
% [SU_Efact,SU_Etime] = get_tavg_peaks((SU_stas-1),tlags,search_range);
% [SU_Ifact,SU_Itime] = get_tavg_peaks(-(SU_stas-1),tlags,search_range);
%
% f4 = figure();hold on
% h1 = shadedErrorBar(tlags,mean(SU_stas(supra_SUs,:)),std(SU_stas(supra_SUs,:))/sqrt(length(supra_SUs)),{'color','b'});
% h2 = shadedErrorBar(tlags,mean(SU_stas(gran_SUs,:)),std(SU_stas(gran_SUs,:))/sqrt(length(gran_SUs)),{'color','r'});
% h3 = shadedErrorBar(tlags,mean(SU_stas(infra_SUs,:)),std(SU_stas(infra_SUs,:))/sqrt(length(infra_SUs)),{'color','k'});

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'layer_tkerns.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'layer_sacrates.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

% figufy(f3);
% fname = [fig_dir 'layer_boxplots.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

%% COMPARE GO MOD AND SUBSPACE MOD
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);

%regularization parameter selection
GO_lambda_gain = 3;
GO_lambda_off = 4;
sub_lambda = 3;

%get GO-model SSI, normalized
GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info
GO_NSSI = bsxfun(@rdivide,GO_NSSI,mean(GO_NSSI(:,base_lags),2)); %normalize by pre-sac avg

%get GO-model LL
GO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovLL = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
GO_NLL = bsxfun(@rdivide,GO_LL,GO_ovLL); %normalize TB SSI by overall model info

%get subspace model SSI, normalized
sub_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_subMod{sub_lambda}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
sub_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_subMod{sub_lambda}.ovInfo,all_SU_data(cur_SUs));
sub_NSSI = bsxfun(@rdivide,sub_SSI,sub_ovinfos); %normalize GO SSI by overall model info
sub_NSSI = bsxfun(@rdivide,sub_NSSI,mean(sub_NSSI(:,base_lags),2)); %normalize by pre-sac avg

%get subspace model LL
sub_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_subMod{sub_lambda}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
sub_ovLL = arrayfun(@(x) x.sacStimProc.gsac_subMod{sub_lambda}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
sub_NLL = bsxfun(@rdivide,sub_LL,sub_ovLL); %normalize TB SSI by overall model info

%slight temporal smoothing of LLs
for ii = 1:length(cur_SUs)
    GO_NLL(ii,:) = jmm_smooth_1d_cor(GO_NLL(ii,:),1);
    sub_NLL(ii,:) = jmm_smooth_1d_cor(sub_NLL(ii,:),1);
end

xl = [-0.1 0.3];
f1 = figure();hold on
h1 = shadedErrorBar(slags*dt,mean(GO_NSSI),std(GO_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2 = shadedErrorBar(slags*dt,mean(sub_NSSI),std(sub_NSSI)/sqrt(length(cur_SUs)),{'color','r'});
ylim([0.5 1.15]);
line(xl,[1 1],'color','k');

% f4 = figure();hold on
% h1 = shadedErrorBar(slags*dt,mean(GO_NLL),std(GO_NLL)/sqrt(length(cur_SUs)),{'color','b'});
% h2 = shadedErrorBar(slags*dt,mean(sub_NLL),std(sub_NLL)/sqrt(length(cur_SUs)),{'color','r'});

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'GO_SUB_SSI.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% MUA GSAC GRAY TRIG AVG
mua_sm = 0.0075/trig_avg_params.dt; %smoothing sigma for MUA rates

%load in all MUA models and compute excitatory temporal kernels
n_probes = 24;
Elist = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};
Olist = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
n_lem_expts = length(Elist);
expt_MUA_stas = nan(2,n_lem_expts,n_probes,length(tlags));
for ii = 1:n_lem_expts
    Expt_name = Elist{ii};
    mod_dir = ['~/Analysis/bruce/' Expt_name '/models/'];
    sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    for kk = 1:2
        cur_ori = Olist(ii,kk);
        if ~isnan(cur_ori)
            tname = strcat(sac_dir,base_tname,sprintf('_ori%d',cur_ori));
            load(tname);
            expt_MUA_stas(kk,ii,:,:) = mua_data.gsac_avg';
        end
    end
end
%average across multiple stim oris for each expt
expt_MUA_stas = squeeze(nanmean(expt_MUA_stas));
expt_MUA_stas_lem = reshape(expt_MUA_stas,[],length(tlags));

%load in all MUA models and compute excitatory temporal kernels
n_probes = 96;
Elist = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};
Olist = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
Elist = {'G085','G086','G087','G088','G089','G091','G093','G095'};
Olist = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
n_lem_expts = length(Elist);
expt_MUA_stas = nan(2,n_lem_expts,n_probes,length(tlags));
for ii = 1:n_lem_expts
    Expt_name = Elist{ii};
    mod_dir = ['~/Analysis/bruce/' Expt_name '/models/'];
    sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    for kk = 1:2
        cur_ori = Olist(ii,kk);
        if ~isnan(cur_ori)
            tname = strcat(sac_dir,base_tname,sprintf('_ori%d',cur_ori));
            load(tname);
            expt_MUA_stas(kk,ii,:,:) = mua_data.gsac_avg';
        end
    end
end
%average across multiple stim oris for each expt
expt_MUA_stas = squeeze(nanmean(expt_MUA_stas));
expt_MUA_stas_jbe = reshape(expt_MUA_stas,[],length(tlags));


%slight smoothing of MUA trig avg rates
for ii = 1:size(expt_MUA_stas_jbe,1)
    expt_MUA_stas_jbe(ii,:) = jmm_smooth_1d_cor(expt_MUA_stas_jbe(ii,:),mua_sm);
end
for ii = 1:size(expt_MUA_stas_lem,1)
    expt_MUA_stas_lem(ii,:) = jmm_smooth_1d_cor(expt_MUA_stas_lem(ii,:),mua_sm);
end

all_MUA_stas = cat(1,expt_MUA_stas_lem,expt_MUA_stas_jbe);
avg_MUA_stas = 0.5*mean(expt_MUA_stas_lem) + 0.5*mean(expt_MUA_stas_jbe);

f1 = figure(); hold on
% plot(tlags,mean(all_MUA_stas),'g','linewidth',2);
% plot(tlags,mean(expt_MUA_stas_lem),'r','linewidth',2);
% plot(tlags,mean(expt_MUA_stas_jbe),'b','linewidth',2);
plot(tlags,avg_MUA_stas,'k','linewidth',2);
ylim([0.75 1.25]);
xlim([-0.25 0.55]);

%PRINT PLOTS
fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'MUA_GSAC_AVG.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% PARALLEL VS ORTHOGANOL MSACS
%number of parallel and orthoganol msacs (subsampled so they are matched)
N_msacs_Par_sub = arrayfun(@(x) x.trig_avg.N_msacs_Par_sub,all_SU_data);
N_msacs_Orth_sub = arrayfun(@(x) x.trig_avg.N_msacs_Orth_sub,all_SU_data);

cur_SUs = find(avg_rates >= min_rate & N_msacs >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
% cur_SUs = find(avg_rates >= min_rate & N_msacs_Par_sub >= min_TA_Nsacs & N_msacs_Orth_sub >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
% all_msac_par = cell2mat(arrayfun(@(x) x.trig_avg.msac_Par_avg', all_SU_data(cur_SUs),'uniformoutput',0));
% all_msac_orth = cell2mat(arrayfun(@(x) x.trig_avg.msac_Orth_avg', all_SU_data(cur_SUs),'uniformoutput',0));

%par and orth msac trig avgs, using matched sample sizes
all_msac_par = cell2mat(arrayfun(@(x) x.trig_avg.msac_Par_sub_avg', all_SU_data(cur_SUs),'uniformoutput',0));
all_msac_orth = cell2mat(arrayfun(@(x) x.trig_avg.msac_Orth_sub_avg', all_SU_data(cur_SUs),'uniformoutput',0));

%spline interpolate
tlags_up = linspace(tlags(1),tlags(end),500);
all_msac_par_up = spline(tlags,all_msac_par,tlags_up);
all_msac_orth_up = spline(tlags,all_msac_orth,tlags_up);

%find suppression peaks
search_range = [0 0.3];
[par_Ifact,par_inhtime] = get_tavg_peaks(-(all_msac_par_up-1),tlags_up,search_range);
[orth_Ifact,orth_inhtime] = get_tavg_peaks(-(all_msac_orth_up-1),tlags_up,search_range);

%find excitation peaks
search_range = [0 0.3];
[par_Efact,par_exctime] = get_tavg_peaks((all_msac_par_up-1),tlags_up,search_range);
[orth_Efact,orth_exctime] = get_tavg_peaks((all_msac_orth_up-1),tlags_up,search_range);

xl = [-0.1 0.3];

close all
%COMBINED MONKEYS WITH SEPARATE TRACE FOR SUB_POP
f1 = figure(); hold on
h1=shadedErrorBar(tlags,nanmean(all_msac_par),nanstd(all_msac_par)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(all_msac_orth),nanstd(all_msac_orth)/sqrt(length(cur_SUs)),{'color','k'});
xlim(xl);
ylim([0.7 1.25]);
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

% % %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'MUA_par_orth3.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
%% PARALLEL VS ORTHOGANOL GSACS

cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);

%par and orth msac trig avgs, using matched sample sizes
all_gsac_par = cell2mat(arrayfun(@(x) x.trig_avg.gsac_small_ocomp', all_SU_data(cur_SUs),'uniformoutput',0));
all_gsac_orth = cell2mat(arrayfun(@(x) x.trig_avg.gsac_large_ocomp', all_SU_data(cur_SUs),'uniformoutput',0));

%spline interpolate
tlags_up = linspace(tlags(1),tlags(end),500);
all_gsac_par_up = spline(tlags,all_gsac_par,tlags_up);
all_gsac_orth_up = spline(tlags,all_gsac_orth,tlags_up);

%find suppression peaks
search_range = [0 0.3];
[par_Ifact,par_inhtime] = get_tavg_peaks(-(all_gsac_par_up-1),tlags_up,search_range);
[orth_Ifact,orth_inhtime] = get_tavg_peaks(-(all_gsac_orth_up-1),tlags_up,search_range);

%find excitation peaks
search_range = [0 0.3];
[par_Efact,par_exctime] = get_tavg_peaks((all_gsac_par_up-1),tlags_up,search_range);
[orth_Efact,orth_exctime] = get_tavg_peaks((all_gsac_orth_up-1),tlags_up,search_range);

xl = [-0.1 0.3];

close all
%COMBINED MONKEYS WITH SEPARATE TRACE FOR SUB_POP
f1 = figure(); hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_par),nanstd(all_gsac_par)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(all_gsac_orth),nanstd(all_gsac_orth)/sqrt(length(cur_SUs)),{'color','k'});
xlim(xl);
ylim([0.7 1.25]);
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

% % %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'GSAC_OCOMP.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


%% STA ANALYSIS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs);

clear z_D
for ii = 1:length(cur_SUs)
    ov_STA{ii} = all_SU_data(cur_SUs(ii)).sacStimProc.ov_phaseDep_sta;
    ov_STA_len(ii) = sqrt(sum(ov_STA{ii}(:).^2))/numel(ov_STA{ii});
    ov_STA_max(ii) = max(abs(ov_STA{ii}(:)));
    avg_spks(ii) = mean(all_SU_data(cur_SUs(ii)).sacStimProc.gsac_nspks);
    
    space_profile = mean(abs(hilbert(ov_STA{ii})));
    space_profile = space_profile/sum(space_profile);
    [~,opt_lag] = max(std(ov_STA{ii},[],2));
    [XX,TT] = meshgrid(1:size(ov_STA{ii},2),1:size(ov_STA{ii},1));
    uset = find(TT(:) == opt_lag);
    %    a_STA = all_SU_data(cur_SUs(ii)).sacStimProc.gsac_phaseDep_sta;
    %    z_STA = (a_STA - all_SU_data(cur_SUs(ii)).sacStimProc.gsac_phaseDep_sta_nullMean);
    %    z_STA = z_STA./all_SU_data(cur_SUs(ii)).sacStimProc.gsac_phaseDep_sta_nullStd;
    z_STA = (all_SU_data(cur_SUs(ii)).sacStimProc.gsac_phaseDep_sta_nullStd);
    
    z_STA = reshape(z_STA,length(slags),[]);
    z_D(ii,:) = sqrt(sum(bsxfun(@times,z_STA(:,uset).^2,space_profile),2))/length(uset);
end

beg_slags = find(slags <= 0);
z_D_avg = mean(z_D(:,beg_slags),2);
z_Dn = bsxfun(@rdivide,z_D,z_D_avg);

for ii = 1:length(cur_SUs)
    z_Dn(ii,:) = jmm_smooth_1d_cor(z_Dn(ii,:),1);
end

min_nspks = 250;

prc = prctile(ov_STA_max(avg_spks >= min_nspks),[25 50 75]);

xl = [-0.1 0.3];
f1 = figure();hold on
cur = find(avg_spks >= min_nspks);
shadedErrorBar(slags*dt,mean(z_Dn(cur,:)),std(z_Dn(cur,:))/sqrt(length(cur)),{'color','k'});
cur = find(avg_spks >= min_nspks & ov_STA_max >= prc(2));
shadedErrorBar(slags*dt,mean(z_Dn(cur,:)),std(z_Dn(cur,:))/sqrt(length(cur)),{'color','b'});
% cur = find(avg_spks >= min_nspks & ov_STA_max >= prc(3));
% shadedErrorBar(slags*dt,mean(z_Dn(cur,:)),std(z_Dn(cur,:))/sqrt(length(cur)),{'color','r'});

% line(xl,[1 1],'color','k');
% line([0 0],ylim(),'color','k');

%% GENERATE ILLUSTRATION OF SPATIOTEMPORAL STIM
close all
npix = 12;
nT = 50;
probNGray = 1-0.12;

randset = rand(nT,npix) > probNGray;
randX = zeros(nT,npix);
randX(randset) = 1;
randset = rand(nT,npix) > 0.5;
randX(randset) = -randX(randset);

f1 = figure();
imagesc(randX); colormap(gray);
set(gca,'Xtick',[],'Ytick',[]);
%PRINT PLOTS
fig_width = 2.5; rel_height = 3.5;
figufy(f1);
fname = [fig_dir 'examp_ST_stim.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% INFO TIMING ANALYSIS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
lambda_gain = 3;
info_lags = sacInfoTiming(end).lag_axis;

gsac_info_before = cell2mat(arrayfun(@(x) x.info_time.gsac_info_before(lambda_gain,:), all_SU_data(cur_SUs),'uniformoutput',0));
gsac_info_after = cell2mat(arrayfun(@(x) x.info_time.gsac_info_after(lambda_gain,:), all_SU_data(cur_SUs),'uniformoutput',0));
gsac_info_during = cell2mat(arrayfun(@(x) x.info_time.gsac_info_during(lambda_gain,:), all_SU_data(cur_SUs),'uniformoutput',0));
gsac_info_Bbefore = cell2mat(arrayfun(@(x) x.info_time.gsac_base_info_before, all_SU_data(cur_SUs),'uniformoutput',0));
gsac_info_Bafter = cell2mat(arrayfun(@(x) x.info_time.gsac_base_info_after, all_SU_data(cur_SUs),'uniformoutput',0));
gsac_info_Bduring = cell2mat(arrayfun(@(x) x.info_time.gsac_base_info_during, all_SU_data(cur_SUs),'uniformoutput',0));
gsac_info_unshuff = cell2mat(arrayfun(@(x) x.info_time.gsac_base_info_unshuff, all_SU_data(cur_SUs),'uniformoutput',0));

gsac_avg_rate = cell2mat(arrayfun(@(x) x.info_time.gsac_avg_rate', all_SU_data(cur_SUs),'uniformoutput',0));

base_lags = find(info_lags < -0.05 & info_lags > -0.15);
% baseline_info = mean(gsac_info_Bbefore(:,1:25),2);
% baseline_info = mean(gsac_info_before(:,base_lags),2);
% baseline_info = mean(gsac_info_unshuff,2);
baseline_info = arrayfun(@(x) x.sacStimProc.ModData.rectGQM.LLimp,all_SU_data(cur_SUs));
ov_avgrate = arrayfun(@(x) x.sacStimProc.ModData.unit_data.avg_rate,all_SU_data(cur_SUs))*dt;
baseline_info = baseline_info.*ov_avgrate;

gsac_info_before = bsxfun(@times,gsac_info_before,gsac_avg_rate);
gsac_info_after = bsxfun(@times,gsac_info_after,gsac_avg_rate);
gsac_info_during = bsxfun(@times,gsac_info_during,gsac_avg_rate);

gsac_info_Bbefore = bsxfun(@times,gsac_info_Bbefore,ov_avgrate);
gsac_info_Bafter = bsxfun(@times,gsac_info_Bafter,ov_avgrate);
gsac_info_Bduring = bsxfun(@times,gsac_info_Bduring,ov_avgrate);

norm_info_before = bsxfun(@rdivide,gsac_info_before,baseline_info);
norm_info_after = bsxfun(@rdivide,gsac_info_after,baseline_info);
norm_info_during = bsxfun(@rdivide,gsac_info_during,baseline_info);

norm_Binfo_before = bsxfun(@rdivide,gsac_info_Bbefore,baseline_info);
norm_Binfo_after = bsxfun(@rdivide,gsac_info_Bafter,baseline_info);
norm_Binfo_during = bsxfun(@rdivide,gsac_info_Bduring,baseline_info);


sm_sig = 0.75;
if sm_sig > 0
    for ii = 1:length(cur_SUs)
        norm_info_before(ii,:) = jmm_smooth_1d_cor(norm_info_before(ii,:),sm_sig);
        norm_info_after(ii,:) = jmm_smooth_1d_cor(norm_info_after(ii,:),sm_sig);
        norm_info_during(ii,:) = jmm_smooth_1d_cor(norm_info_during(ii,:),sm_sig);
        norm_Binfo_before(ii,:) = jmm_smooth_1d_cor(norm_Binfo_before(ii,:),sm_sig);
        norm_Binfo_after(ii,:) = jmm_smooth_1d_cor(norm_Binfo_after(ii,:),sm_sig);
        norm_Binfo_during(ii,:) = jmm_smooth_1d_cor(norm_Binfo_during(ii,:),sm_sig);
    end
end


info_mod_during = norm_info_during - norm_Binfo_during;
info_mod_before = norm_info_before - norm_Binfo_before;
info_mod_after = norm_info_after - norm_Binfo_after;

% for ii = 1:length(cur_SUs)
%     info_mod_during(ii,:) = jmm_smooth_1d_cor(info_mod_during(ii,:),0.5);
% end

f1 = figure(); hold on
shadedErrorBar(info_lags,mean(norm_info_before),std(norm_info_before)/sqrt(length(cur_SUs)),{'color','r'});
shadedErrorBar(info_lags,mean(norm_info_after),std(norm_info_after)/sqrt(length(cur_SUs)),{'color','b'});
shadedErrorBar(info_lags,mean(norm_info_during),std(norm_info_during)/sqrt(length(cur_SUs)),{'color','k'});
plot(info_lags,mean(norm_Binfo_before),'r--','linewidth',1);
plot(info_lags,mean(norm_Binfo_after),'b--','linewidth',1);
plot(info_lags,mean(norm_Binfo_during),'k--','linewidth',1);
line([-0.3 0.3],[1 1],'color','k')
xlim([-0.1 0.2]);
ylim([0 1.2])
line([0 0],[0 1.2],'color','k')
xlabel('Time since fixation onset (s)');
ylabel('Relative stim info');


% f2 = figure(); hold on
% shadedErrorBar(info_lags,mean(info_mod_before),std(info_mod_before)/sqrt(length(cur_SUs)),{'color','r'});
% shadedErrorBar(info_lags,mean(info_mod_after),std(info_mod_after)/sqrt(length(cur_SUs)),{'color','b'});
% shadedErrorBar(info_lags,mean(info_mod_during),std(info_mod_during)/sqrt(length(cur_SUs)),{'color','k'});
% xlim([-0.1 0.2]);
% xlabel('Time since fixation onset (s)');

% exCell = find([all_SU_data(cur_SUs).expt_num] == 93,1,'last');
% f2 = figure(); hold on
% plot(info_lags,norm_info_before(exCell,:),'r','linewidth',2);
% plot(info_lags,norm_info_after(exCell,:),'b','linewidth',2);
% plot(info_lags,norm_info_during(exCell,:),'k','linewidth',2);
% plot(info_lags,norm_Binfo_before(exCell,:),'r--','linewidth',1);
% plot(info_lags,norm_Binfo_after(exCell,:),'b--','linewidth',1);
% plot(info_lags,norm_Binfo_during(exCell,:),'k--','linewidth',1);
% xlim([-0.2 0.2]);
% ylim([0 1.2])
% line([0 0],[0 1.2],'color','k')
% xlabel('Time since fixation onset (s)');
% ylabel('Relative stim info');

% PRINT PLOTS
fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'gsac_info_timing4.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);
%
% figufy(f2);
% fname = [fig_dir 'gsac_info_supp.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% INFO TIMING ANALYSIS MSACS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
lambda_gain = 3;
info_lags = sacInfoTiming(end).lag_axis;

msac_info_before = cell2mat(arrayfun(@(x) x.info_time.msac_info_before(lambda_gain,:), all_SU_data(cur_SUs),'uniformoutput',0));
msac_info_after = cell2mat(arrayfun(@(x) x.info_time.msac_info_after(lambda_gain,:), all_SU_data(cur_SUs),'uniformoutput',0));
msac_info_during = cell2mat(arrayfun(@(x) x.info_time.msac_info_during(lambda_gain,:), all_SU_data(cur_SUs),'uniformoutput',0));
msac_info_Bbefore = cell2mat(arrayfun(@(x) x.info_time.msac_base_info_before, all_SU_data(cur_SUs),'uniformoutput',0));
msac_info_Bafter = cell2mat(arrayfun(@(x) x.info_time.msac_base_info_after, all_SU_data(cur_SUs),'uniformoutput',0));
msac_info_Bduring = cell2mat(arrayfun(@(x) x.info_time.msac_base_info_during, all_SU_data(cur_SUs),'uniformoutput',0));
msac_info_unshuff = cell2mat(arrayfun(@(x) x.info_time.msac_base_info_unshuff, all_SU_data(cur_SUs),'uniformoutput',0));

msac_avg_rate = cell2mat(arrayfun(@(x) x.info_time.msac_avg_rate', all_SU_data(cur_SUs),'uniformoutput',0));

base_lags = find(info_lags < -0.05 & info_lags > -0.15);
% baseline_info = mean(msac_info_Bbefore(:,1:25),2);
% baseline_info = mean(msac_info_before(:,base_lags),2);
baseline_info = arrayfun(@(x) x.sacStimProc.ModData.rectGQM.LLimp,all_SU_data(cur_SUs));
ov_avgrate = arrayfun(@(x) x.sacStimProc.ModData.unit_data.avg_rate,all_SU_data(cur_SUs))*dt;
baseline_info = baseline_info.*ov_avgrate;

msac_info_before = bsxfun(@times,msac_info_before,msac_avg_rate);
msac_info_after = bsxfun(@times,msac_info_after,msac_avg_rate);
msac_info_during = bsxfun(@times,msac_info_during,msac_avg_rate);

msac_info_Bbefore = bsxfun(@times,msac_info_Bbefore,ov_avgrate);
msac_info_Bafter = bsxfun(@times,msac_info_Bafter,ov_avgrate);
msac_info_Bduring = bsxfun(@times,msac_info_Bduring,ov_avgrate);


norm_info_before = bsxfun(@rdivide,msac_info_before,baseline_info);
norm_info_after = bsxfun(@rdivide,msac_info_after,baseline_info);
norm_info_during = bsxfun(@rdivide,msac_info_during,baseline_info);

norm_Binfo_before = bsxfun(@rdivide,msac_info_Bbefore,baseline_info);
norm_Binfo_after = bsxfun(@rdivide,msac_info_Bafter,baseline_info);
norm_Binfo_during = bsxfun(@rdivide,msac_info_Bduring,baseline_info);

sm_sig = 0.75;
if sm_sig > 0
    for ii = 1:length(cur_SUs)
        norm_info_before(ii,:) = jmm_smooth_1d_cor(norm_info_before(ii,:),sm_sig);
        norm_info_after(ii,:) = jmm_smooth_1d_cor(norm_info_after(ii,:),sm_sig);
        norm_info_during(ii,:) = jmm_smooth_1d_cor(norm_info_during(ii,:),sm_sig);
        norm_Binfo_before(ii,:) = jmm_smooth_1d_cor(norm_Binfo_before(ii,:),sm_sig);
        norm_Binfo_after(ii,:) = jmm_smooth_1d_cor(norm_Binfo_after(ii,:),sm_sig);
        norm_Binfo_during(ii,:) = jmm_smooth_1d_cor(norm_Binfo_during(ii,:),sm_sig);
    end
end


info_mod_during = norm_info_during - norm_Binfo_during;
info_mod_before = norm_info_before - norm_Binfo_before;
info_mod_after = norm_info_after - norm_Binfo_after;

f1 = figure(); hold on
shadedErrorBar(info_lags,mean(norm_info_before),std(norm_info_before)/sqrt(length(cur_SUs)),{'color','r'});
shadedErrorBar(info_lags,mean(norm_info_after),std(norm_info_after)/sqrt(length(cur_SUs)),{'color','b'});
shadedErrorBar(info_lags,mean(norm_info_during),std(norm_info_during)/sqrt(length(cur_SUs)),{'color','k'});
plot(info_lags,mean(norm_Binfo_before),'r--','linewidth',1);
plot(info_lags,mean(norm_Binfo_after),'b--','linewidth',1);
plot(info_lags,mean(norm_Binfo_during),'k--','linewidth',1);
line([-0.3 0.3],[1 1],'color','k')
xlim([-0.1 0.2]);
ylim([0 1.2])
line([0 0],[0 1.2],'color','k')
xlabel('Time since fixation onset (s)');
ylabel('Relative stim info');

% f2 = figure(); hold on
% shadedErrorBar(info_lags,mean(info_mod_before),std(info_mod_before)/sqrt(length(cur_SUs)),{'color','r'});
% shadedErrorBar(info_lags,mean(info_mod_after),std(info_mod_after)/sqrt(length(cur_SUs)),{'color','b'});
% shadedErrorBar(info_lags,mean(info_mod_during),std(info_mod_during)/sqrt(length(cur_SUs)),{'color','k'});
% xlim([-0.2 0.2]);
% xlabel('Time since fixation onset (s)');

% exCell = find([all_SU_data(cur_SUs).expt_num] == 93,1,'last');
% f2 = figure(); hold on
% plot(info_lags,norm_info_before(exCell,:),'r','linewidth',2);
% plot(info_lags,norm_info_after(exCell,:),'b','linewidth',2);
% plot(info_lags,norm_info_during(exCell,:),'k','linewidth',2);
% plot(info_lags,norm_Binfo_before(exCell,:),'r--','linewidth',1);
% plot(info_lags,norm_Binfo_after(exCell,:),'b--','linewidth',1);
% plot(info_lags,norm_Binfo_during(exCell,:),'k--','linewidth',1);
% xlim([-0.2 0.2]);o
% ylim([0 1.2])
% line([0 0],[0 1.2],'color','k')
% xlabel('Time since fixation onset (s)');
% ylabel('Relative stim info');
%
% % PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'msac_info_timing4.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'msac_info_supp.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% PLOTS COMPARING SACCADE MODULATION AND RF PROPS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);

lambda_ii = 4; %reg hyper selection
xl = [-0.1 0.3];

TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod{lambda_ii}.lagX*dt; %TB time axis
TBbase_lags = find(TB_Xtick <= 0);

%sac-trig avg firing rate
all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
search_range = [0 0.3];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_tavg-1),tlags,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_tavg-1,tlags,search_range);

%get normalized TB-model SSI (computing using standard t-axis)
TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info
TB_NSSI = bsxfun(@rdivide,TB_NSSI,mean(TB_NSSI(:,base_lags),2));

%get normalized TB-model gain and offset (computed using TB-taxis)
TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates
% TB_Noffset = bsxfun(@rdivide,TB_Noffset,TB_gains); %normalize by gain
TB_Noffset = bsxfun(@minus,TB_Noffset,mean(TB_Noffset(:,TBbase_lags),2));
TB_gains = bsxfun(@rdivide,TB_gains,mean(TB_gains(:,TBbase_lags),2));

%find timing and magnitude of TB-model SSI suppression peaks
search_range = [0 0.15]; %search range for TB SSI suppression peaks
[SSI_Ifact,SSI_inhtime] = get_tavg_peaks(-(TB_NSSI-1),slags*dt,search_range);
[Gain_Ifact,GAIN_inhtime] = get_tavg_peaks(-(TB_gains-1),TB_Xtick,search_range);

[Off_Efact,Off_Etime] = get_tavg_peaks(TB_Noffset,TB_Xtick,[0 0.3]);

max_RF_size = 1; 
include = find(2*RF_sigma(cur_SUs) <= max_RF_size);

cur_jbe = include(ismember(cur_SUs(include),jbe_SUs));
cur_lem = include(ismember(cur_SUs(include),lem_SUs));

rf_ticks = [0.1:0.1:1]; %make ticks at these RF size locations in log space

%plot RF size vs rate suppression strength
f1 = figure();
plot(log10(RF_sigma(cur_SUs(cur_jbe))*2),gsac_Ifact(cur_jbe),'o');
hold on
plot(log10(RF_sigma(cur_SUs(cur_lem))*2),gsac_Ifact(cur_lem),'ro');
set(gca,'xtick',log10(rf_ticks),'xticklabel',rf_ticks)
xlabel('RF size (deg)');
ylabel('Rate suppression strength');

[a,b] = corr(RF_sigma(cur_SUs(include)),gsac_Ifact(include),'type','spearman');
[aj,bj] = corr(RF_sigma(cur_SUs(cur_jbe)),gsac_Ifact(cur_jbe),'type','spearman');
[al,bl] = corr(RF_sigma(cur_SUs(cur_lem)),gsac_Ifact(cur_lem),'type','spearman');
fprintf('overall p=%.4f. JBE p=%.4f. LEM p=%.4f\n',b,bj,bl);

%plot RF size vs enhancement strength
f2 = figure();
plot(log10(RF_sigma(cur_SUs(cur_jbe))*2),gsac_Efact(cur_jbe),'o');
hold on
plot(log10(RF_sigma(cur_SUs(cur_lem))*2),gsac_Efact(cur_lem),'ro');
set(gca,'xtick',log10(rf_ticks),'xticklabel',rf_ticks)
xlabel('RF size (deg)');
ylabel('Rate enhancement strength');

[a,b] = corr(RF_sigma(cur_SUs(include)),gsac_Efact(include),'type','spearman');
[aj,bj] = corr(RF_sigma(cur_SUs(cur_jbe)),gsac_Efact(cur_jbe),'type','spearman');
[al,bl] = corr(RF_sigma(cur_SUs(cur_lem)),gsac_Efact(cur_lem),'type','spearman');
fprintf('overall p=%.4f. JBE p=%.4f. LEM p=%.4f\n',b,bj,bl);

% plot RF size vs SSI suppression strength
f3 = figure();
plot(log10(RF_sigma(cur_SUs(cur_jbe))*2),SSI_Ifact(cur_jbe),'o');
hold on
plot(log10(RF_sigma(cur_SUs(cur_lem))*2),SSI_Ifact(cur_lem),'ro');
set(gca,'xtick',log10(rf_ticks),'xticklabel',rf_ticks)
xlabel('RF size (deg)');
ylabel('SSI suppression strength');

%fit animal specific robust linear fits
r_JBE = robustfit(log10(RF_sigma(cur_SUs(cur_jbe))*2),SSI_Ifact(cur_jbe));
r_LEM = robustfit(log10(RF_sigma(cur_SUs(cur_lem))*2),SSI_Ifact(cur_lem));
xeval = log10(logspace(log10(0.1),log10(0.5),100));
% plot(xeval,xeval*r_JBE(2)+r_JBE(1),'b');
% plot(xeval,xeval*r_LEM(2)+r_LEM(1),'r');

[a,b] = corr(RF_sigma(cur_SUs(include)),SSI_Ifact(include),'type','spearman');
[aj,bj] = corr(RF_sigma(cur_SUs(cur_jbe)),SSI_Ifact(cur_jbe),'type','spearman');
[al,bl] = corr(RF_sigma(cur_SUs(cur_lem)),SSI_Ifact(cur_lem),'type','spearman');
fprintf('overall p=%.4e. JBE p=%.4g. LEM p=%.4g\n',b,bj,bl);

SU_exptnums = [all_SU_data(cur_SUs).expt_num];
[H,ATAB,CTAB,STATS]=aoctool(log10(RF_sigma(cur_SUs(include))),SSI_Ifact(include),SU_exptnums(include)',0.05,'RF Sigma','SSI-sup','rec_num','off',5);

include = 1:length(cur_SUs);
cur_jbe = include(ismember(cur_SUs(include),jbe_SUs));
cur_lem = include(ismember(cur_SUs(include),lem_SUs));

ecc_ticks = [0.5:0.5:5]; %make ticks at these RF size locations in log space
%now plot SSI suppression vs eccentricity
f4 = figure();
plot(log10(RF_ecc(cur_SUs(cur_jbe))),SSI_Ifact(cur_jbe),'o');
hold on
plot(log10(RF_ecc(cur_SUs(cur_lem))),SSI_Ifact(cur_lem),'ro');
xlabel('RF eccentricity (deg)');
ylabel('SSI suppression strength');
set(gca,'xtick',log10(ecc_ticks),'xticklabel',ecc_ticks)
%fit animal specific robust linear fits
r_JBE = robustfit(log10(RF_ecc(cur_SUs(cur_jbe))),SSI_Ifact(cur_jbe));
r_LEM = robustfit(log10(RF_ecc(cur_SUs(cur_lem))),SSI_Ifact(cur_lem));
xeval = log10(logspace(log10(0.3),log10(5),100));
% plot(xeval,xeval*r_JBE(2)+r_JBE(1),'b');
% plot(xeval,xeval*r_LEM(2)+r_LEM(1),'r');

% % PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'RFsigma_vs_ratesup.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'RFsigma_vs_rateenh.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'RFsigma_vs_SSIsup.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
% 
% figufy(f4);
% fname = [fig_dir 'RFecc_vs_SSIsup.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% SACCADE GAIN TIMING ANALYSIS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0); %use slightly earlier def of backgnd time points since the upstream kernels start slightly earlier

all_SDs = cell2mat(arrayfun(@(x) x.sac_delay.single_gSD',all_SU_data(cur_SUs),'uniformoutput',0));
all_SDs = bsxfun(@rdivide,all_SDs,max(all_SDs,[],2));
min_SD_frac = 0.2;

flen = 15;
gsac_gain = nan(length(cur_SUs),flen,length(slags));
% gsac_gain_sig = nan(length(cur_SUs),flen,length(slags));
% n_use_cells = nan(flen,1);
for ff = 1:flen;
    cur_gains = 1+cell2mat(arrayfun(@(x) x.sac_delay.gain_filts(ff,:),all_SU_data(cur_SUs),'uniformoutput',0));
    cur_gains = bsxfun(@rdivide,cur_gains,mean(cur_gains(:,base_lags),2)); %normalize by pre-sac gain strength
    gsac_gain(:,ff,:) = cur_gains;
    
%     cur_use_cells = all_SDs(:,ff) >= min_SD_frac;
%     cur_gains(~cur_use_cells,:) = nan;
%     n_use_cells(ff) = sum(cur_use_cells);
%     gsac_gain_sig(:,ff,:) = cur_gains;
end

%population avgs
avg_gains = squeeze(nanmean(gsac_gain));
sem_gains = squeeze(nanstd(gsac_gain))/sqrt(length(cur_SUs));

%interpolate onto finer temporal grid
slags_up = linspace(slags(1),slags(end),200);
avg_gains_up = spline(slags,avg_gains,slags_up);
[avg_gain_peakamps,avg_gain_peaklocs] = min(avg_gains_up,[],2);

yr = [0.6 1.1];
ulag_range = 3:10;
cmap = jet(length(ulag_range));
f1 = figure(); hold on
for ii = 1:length(ulag_range)
    plot(slags_up*dt,avg_gains_up(ulag_range(ii),:),'color',cmap(ii,:),'linewidth',1);
    line([0 0]+slags_up(avg_gain_peaklocs(ulag_range(ii)))*dt,yr,'color',cmap(ii,:));
    %     shadedErrorBar(slags*dt,avg_gains(ulag_range(ii),:),sem_gains(ulag_range(ii),:),{'color',cmap(ii,:)});
    plot(slags_up(avg_gain_peaklocs(ulag_range(ii)))*dt,avg_gain_peakamps(ulag_range(ii)),'o','color',cmap(ii,:),'linewidth',2)
end
xlim([-0.05 0.2])
ylim(yr);
line([0 0],yr,'color','k','linestyle','--');
line([-0.1 0.3],[1 1],'color','k','linestyle','--');
xlabel('Time (s)');
ylabel('Gain');

%hyperparameter selection
post_lambda_off = 4;
post_lambda_gain = 3;

%pre and post gain kernels
gsac_post_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{post_lambda_off,post_lambda_gain}.mods(3).filtK',all_SU_data(cur_SUs),'uniformoutput',0));
%noramlize gain kernels
gsac_post_gain = bsxfun(@rdivide,gsac_post_gain,mean(gsac_post_gain(:,base_lags),2));
avg_post_gain = mean(gsac_post_gain);

post_gain_up = spline(slags,avg_post_gain,slags_up);
plot(slags_up*dt,post_gain_up,'k--','linewidth',2);


search_range = [0. 0.2];
resh_gsac_gain = reshape(permute(gsac_gain,[3 1 2]),length(slags),[]);
%spline interpolate all gain kernels and find peaks
all_gkerns_up = spline(slags',resh_gsac_gain',slags_up');
all_gkerns_up = reshape(all_gkerns_up',length(slags_up),length(cur_SUs),flen);

[gkern_amps,gkern_times] = deal(nan(length(cur_SUs),flen));
for ff = 1:flen
    [cur_gkern_max,cur_gkern_time] = get_tavg_peaks(-(squeeze(all_gkerns_up(:,:,ff))'-1),slags_up*dt,search_range);
    gkern_amps(:,ff) = cur_gkern_max;
    gkern_times(:,ff) = cur_gkern_time;
end
uset = all_SDs > min_SD_frac;
gkern_times(~uset) = nan;
noise_SDs = squeeze(std(gsac_gain(:,:,base_lags),[],3));
min_SNR = 0;
bad_peaks = gkern_amps <= min_SNR*noise_SDs;
gkern_times(bad_peaks) = nan;

tax = (0:(flen-1))*mod_dt + mod_dt/2;

min_pts = 5;
[cell_corr ,cell_slope,cell_offset] = deal(nan(length(cur_SUs),1));
for ii = 1:length(cur_SUs)
    curset = find(~isnan(gkern_times(ii,:)));
    if length(curset) >= min_pts
        cell_corr(ii) = corr((1:length(curset))',gkern_times(ii,curset)','type','spearman');
        temp = robustfit(tax(curset),gkern_times(ii,curset)');
%         temp = regress(gkern_times(ii,curset)',[ones(length(curset),1) tax(curset)']);
        cell_slope(ii) = temp(2); cell_offset(ii) = temp(1);
    end
end

all_preds = bsxfun(@times,cell_slope,tax) + nanmean(cell_offset);

min_used_units = 10;
n_used_units = sum(~isnan(gkern_times));
gkern_times(:,n_used_units < min_used_units) = nan;
all_preds(:,n_used_units < min_used_units) = nan;
f2 = figure(); hold on
% shadedErrorBar(1:flen,nanmedian(gkern_times),iqr(gkern_times),{'color','r'});
% shadedErrorBar(tax,nanmean(all_preds),nanstd(all_preds)./sqrt(n_used_units),{'color','r'});
shadedErrorBar(tax,nanmean(gkern_times),nanstd(gkern_times)./sqrt(n_used_units),{'color','b'});
% plot(1:flen,gkern_times,'k.')
line([0 0.14],[0 0.14]+0.04,'color','k');
xlim([0 0.12]);
xlabel('Stimulus latency (s)');
ylabel('Suppression timing (s)');

xe = linspace(-2.5,2.5,50);
nn = histc(cell_slope,xe);
f3 = figure();
stairs(xe,nn,'b','linewidth',2);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line([1 1],yl,'color','r','linestyle','--');
axis tight
xlabel('Slope');
ylabel('Counts');

% for ex_cell = 1:84
ex_cell = 32;
f4 = figure();
imagesc(slags*dt,tax,squeeze(gsac_gain(ex_cell,:,:)))
caxis([0.2 1.8])
ylim([0 0.12]);
xlabel('Time since saccade onset (s)');
ylabel('Stimulus latency (s)');
xlim([-0.05 0.2]);
colorbar
% pause
% close
% end
% % PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'gsac_stimlatency_gainmod.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'gsac_stimlatency_gainmod_avg.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% fig_width = 3.5; rel_height = 0.8;
% figufy(f3);
% fname = [fig_dir 'gsac_stimlatency_gainmod_slopedist.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

fig_width = 3.5; rel_height = 0.8;
figufy(f4);
fname = [fig_dir sprintf('gsac_stimlatency_excell%d.pdf',ex_cell)];
exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f4);

%%
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= -0.025); %use slightly earlier def of backgnd time points since the upstream kernels start slightly earlier

%hyperparameter selection
post_lambda_off = 4;
pre_lambda = 3; %this is just the gain kernel lambda, the offset kernel lambda is fixed equal to post-model (4)

%pre and post gain kernels
gsac_pre_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.stim_kernel',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_pre_gain = bsxfun(@rdivide,gsac_pre_gain,mean(gsac_pre_gain(:,base_lags),2));

all_Tkerns = nan(length(cur_SUs),flen);
all_Ekerns = nan(length(cur_SUs),flen);
all_Ikerns = nan(length(cur_SUs),flen);
for ii = 1:length(cur_SUs)
    %get E and I stim filter temporal kernels
    [cur_tkerns,cur_ekern,cur_ikern] = get_hilbert_tempkerns(all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM);
    all_Ekerns(ii,:) = cur_ekern;
    all_Ikerns(ii,:) = cur_ikern;
    all_Tkerns(ii,:) = nanmean(cur_tkerns,2);
end


% avg_Tkerns = fliplr(nanmean(all_Tkerns));
% avg_Tkerns = (nanmean(all_Tkerns));

for ii = 1:84
examp = ii;
base_kern = all_Tkerns(examp,:); base_kern = base_kern - min(base_kern);base_kern = base_kern/max(base_kern); 
avg_Tkerns = fliplr(all_Tkerns(examp,:));

% avg_pre_gain = nanmean(gsac_pre_gain);
avg_pre_gain = gsac_pre_gain(examp,:);

avg_Tkerns = [avg_Tkerns - min(avg_Tkerns) zeros(1,flen-1)];
% avg_Tkerns = [zeros(1,flen-1) avg_Tkerns - min(avg_Tkerns)];
avg_Tkerns = avg_Tkerns/max(avg_Tkerns);

effective_gains = ones(length(slags),2*flen-1);
for ss = 1:length(slags)
   curset = (ss-flen+1):(ss+flen-1); 
   cur_uset = curset >= 1 & curset <= length(slags); 
   effective_gains(ss,cur_uset) = avg_pre_gain(curset(cur_uset));
end
effective_gains = bsxfun(@times,effective_gains,avg_Tkerns);
effective_gains = fliplr(effective_gains(:,1:flen));

[~,ploc] = max(base_kern); ploc = ploc - 1;
[~,gloc] = min(avg_pre_gain);
f1 = figure();
plot(base_kern,'k');
hold on
plot(effective_gains(ploc+gloc+2,:),'b');
plot(effective_gains(ploc+gloc,:),'r');
plot(effective_gains(ploc+gloc-2,:),'g');

ii
pause
close all
end