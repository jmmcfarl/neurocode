% close all
clear all
clc

fit_unCor = 0;
include_bursts = 0;

fig_dir = '/Users/james/Analysis/bruce/FINsac_mod/figures/';
base_sname = 'sacStimProcFin_noXV';
% base_tname = 'sac_trig_avg_data3';
base_tname = 'sac_trig_avg_data3test';
base_yname = 'sacTypeDep_noXV';

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
% ori_list = [0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan];
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
% Expt_list = {'M266','M270','M275','M277','M281'};%NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
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
            if ~ismember(Expt_num,[296 297])
                %load type-dep data (STILL NEED TO ADD THIS)
                yname = strcat(sac_dir,base_yname,sprintf('_ori%d',ori_list(ee,ii)));
%                 if fit_unCor
%                     yname = strcat(yname,'_unCor');
%                 end
                load(yname);
            end
            
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
min_xvLLimp = 0.0; %(0.05);

tot_Nunits = length(all_SU_data);
avg_rates = arrayfun(@(x) x.sacStimProc.ModData.unit_data.avg_rate,all_SU_data);
tot_spikes = arrayfun(@(x) x.sacStimProc.ModData.unit_data.tot_spikes,all_SU_data);
rec_dur = arrayfun(@(x) x.sacStimProc.ModData.unit_data.N_used_samps,all_SU_data)*dt/60;
mod_xvLLimps = arrayfun(@(x) x.sacStimProc.ModData.rectGQM.xvLLimp,all_SU_data);
mod_xvLLimps_unCor = arrayfun(@(x) x.sacStimProc.ModData.rectGQM_unCor.xvLLimp,all_SU_data);
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
N_gsacs_im = arrayfun(@(x) x.trig_avg.N_gsacs_im,all_SU_data);
N_msacs_im = arrayfun(@(x) x.trig_avg.N_msacs_im,all_SU_data);
N_simsacs = arrayfun(@(x) x.trig_avg.N_simsacs,all_SU_data);

%% GRAY-BACKGROUND SACCADE TRIGGERED AVERAGES
% cur_SUs = find(avg_rates >= min_rate & N_gsacs_gray >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
% all_gsac_gray = cell2mat(arrayfun(@(x) x.trig_avg.gsac_gray_avg', all_SU_data(:),'uniformoutput',0));
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
all_gsac_gray = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(:),'uniformoutput',0));


xl = [-0.1 0.4];

close all

%FOR SEPARATE MONKEYS
% f1 = figure(); hold on
% curset = intersect(jbe_SUs,cur_SUs);
% h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','r'});
% curset = intersect(lem_SUs,cur_SUs);
% h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','b'});

%COMBINED MONKEYS WITH SEPARATE TRACE FOR SUB_POP
f1 = figure(); hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(cur_SUs,:)),nanstd(all_gsac_gray(cur_SUs,:))/sqrt(length(cur_SUs)),{'color','k'});
curset = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs);
% h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','r'});
% plot(tlags,nanmean(all_gsac_gray(curset,:)),'r','linewidth',2);
xlim(xl);
ylim([0.7 1.25]);
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

% xl2 = [-0.25 0.55];
% f2 = figure();
% plot(tlags,nanmean(all_gsac_gray(cur_SUs,:)),'linewidth',2);
% set(gca,'yAxisLocation','right');
% ylim([0.75 1.25]);
% xlim(xl2);

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
all_gsac_gray = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(:),'uniformoutput',0));
all_msac_gray = cell2mat(arrayfun(@(x) x.trig_avg.msac_avg', all_SU_data(:),'uniformoutput',0));

tlags_up = linspace(tlags(1),tlags(end),500);
all_gsac_up = nan(size(all_gsac_gray,1),length(tlags_up));
incl = find(~isnan(all_gsac_gray(:,1)));
all_gsac_up(incl,:) = spline(tlags,all_gsac_gray(incl,:),tlags_up);

search_range = [0 0.3];
[gsac_Efact,gsac_exctime] = get_tavg_peaks((all_gsac_gray-1),tlags,search_range);
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_gray-1),tlags,search_range);
[gsac_Efact_up,gsac_exctime_up] = get_tavg_peaks((all_gsac_up-1),tlags_up,search_range);
[gsac_Ifact_up,gsac_inhtime_up] = get_tavg_peaks(-(all_gsac_up-1),tlags_up,search_range);

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
normal_polarity_units = sigB(gsac_exctime(sigB) > gsac_inhtime(sigB));
reverse_polarity_units = sigB(gsac_exctime(sigB) < gsac_inhtime(sigB));
nsigE = setdiff(cur_SUs,sigE);
nsigI = setdiff(cur_SUs,sigI);
nsigB = setdiff(cur_SUs,sigB);

xl = [0 1]; %mod strength axis
yl = [0 0.3]; %mod timing axis
mS = 4; %marker size
mS2 = 8; %marker size

% tjitter = Tdt/4; %jitter time peaks when plotting to prevent dot occlusion
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
% 
% %PRINT PLOTS
fig_width = 3.5; rel_height = 1;
figufy(f1);
fname = [fig_dir 'Gsac_time_mod_scatter.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

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

xl = [-0.1 0.4];

f2 = figure();
hold on
h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(reverse_polarity_units,:)),nanstd(all_gsac_gray(reverse_polarity_units,:))/sqrt(length(reverse_polarity_units)),{'color','r'});
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(normal_polarity_units,:)),nanstd(all_gsac_gray(normal_polarity_units,:))/sqrt(length(normal_polarity_units)),{'color','b'});
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

% %PRINT PLOTS
fig_width = 3.5; rel_height = 1;
figufy(f1);
fname = [fig_dir 'Gsac_timing_scatter.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

fig_width = 3.5; rel_height = 0.8;
figufy(f2);
fname = [fig_dir 'Gsac_reverse_polarity_avgs.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Gsac_reverse_polarity_avgs.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% TB GAIN AND OFFSET 
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod{1}.lagX*dt;
base_lags = find(TB_Xtick <= 0);

lambda_ii = 4;
xl = [-0.1 0.3];

TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));

TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates
% TB_Noffset = bsxfun(@rdivide,TB_Noffset,TB_gains); %normalize by gains

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

lambda_ii = 4;
xl = [-0.1 0.3];
TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info
TB_NSSI = bsxfun(@rdivide,TB_NSSI,mean(TB_NSSI(:,base_lags),2));

% gsac_avg_rates = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_avg_rate',all_SU_data(cur_SUs),'uniformoutput',0));
% gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
% gsac_rel_rates = bsxfun(@rdivide,gsac_avg_rates,gsac_ov_rates);

gsac_rel_rates = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
gsac_rel_rates_interp = interp1(tlags,gsac_rel_rates',slags*dt)';

% TB_SSI_rate = TB_SSI.*gsac_avg_rates;
% ov_SSI_rate = TB_ovinfos.*gsac_ov_rates;
% TB_NSSI_rate = bsxfun(@rdivide,TB_SSI_rate,ov_SSI_rate);

TB_NSSI_rate = bsxfun(@times,TB_NSSI,gsac_rel_rates_interp);
TB_NSSI_rate = bsxfun(@rdivide,TB_NSSI_rate,mean(TB_NSSI_rate(:,base_lags),2));


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

% % % %PRINT PLOTS
fig_width = 3.5; rel_height = 0.8;
figufy(f1);
% fname = [fig_dir 'Gsac_TB_SSI.pdf'];
fname = [fig_dir 'Gsac_TB_SSI2.pdf'];
if fit_unCor
fname = [fig_dir 'Gsac_TB_SSI_unCor2.pdf'];
end
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);
% % 
% figufy(f2);
% fname = [fig_dir 'Gsac_TBset_relrates.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% PLOTS COMPARING STRONGLY ENHANCED VS STRONGLY SUPPRESSED UNITS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);

base_lags = find(slags <= 0);

lambda_ii = 4;
xl = [-0.1 0.3];

TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod{lambda_ii}.lagX*dt;
TBbase_lags = find(TB_Xtick <= 0);


all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
search_range = [0 0.3];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_tavg-1),tlags,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_tavg-1,tlags,search_range);

TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info
TB_NSSI = bsxfun(@rdivide,TB_NSSI,mean(TB_NSSI(:,base_lags),2));

TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));

TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{lambda_ii}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates
% TB_Noffset = bsxfun(@rdivide,TB_Noffset,TB_gains); %normalize by gain

TB_Noffset = bsxfun(@minus,TB_Noffset,mean(TB_Noffset(:,TBbase_lags),2));
TB_gains = bsxfun(@rdivide,TB_gains,mean(TB_gains(:,TBbase_lags),2));

search_range = [0 0.15];
[SSI_Ifact,SSI_inhtime] = get_tavg_peaks(-(TB_NSSI-1),slags*dt,search_range);
[Gain_Ifact,GAIN_inhtime] = get_tavg_peaks(-(TB_gains-1),TB_Xtick,search_range);
[Off_Efact,OFF_exctime] = get_tavg_peaks((TB_Noffset),TB_Xtick,search_range);

interp_TB_gains = interp1(TB_Xtick,TB_gains',slags*dt)';
interp_TB_offset = interp1(TB_Xtick,TB_Noffset',slags*dt)';
clear gain_SSI_corr offset_SSI_corr
for ii = 1:length(cur_SUs)
    gain_SSI_corr(ii) = partialcorr(interp_TB_gains(ii,:)',TB_NSSI(ii,:)',interp_TB_offset(ii,:)','type','spearman');
    offset_SSI_corr(ii) = partialcorr(interp_TB_offset(ii,:)',TB_NSSI(ii,:)',interp_TB_gains(ii,:)','type','spearman');
end

% EIrat = gsac_Efact./gsac_Ifact;
% EIrat = Off_Efact;
EIrat = SSI_Ifact;
% EIrat = gsac_Ifact;
% EIrat = Gain_Ifact;
% EIrat = gsac_Efact;
% EIrat = gsac_Ifact;
% EIprc = prctile(EIrat,[25 50 75]);
EIprc = prctile(EIrat,[33 50 67]);

% enh_units = find(gsac_Efact > prctile(gsac_Efact,67));
% sup_units = find(gsac_Ifact > prctile(gsac_Ifact,67));

% enh_units = find(EIrat >= EIprc(3));
% sup_units = find(EIrat <= EIprc(1));
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

% %PRINT PLOTS
fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'Gsac_SW_rates2.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'Gsac_SW_SSI2.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

figufy(f3);
fname = [fig_dir 'Gsac_SW_OFFSET2.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);

figufy(f4);
fname = [fig_dir 'Gsac_SW_GAIN2.pdf'];
exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f4);

%% COMPARE SSI, OFFSET AND GAINS FOR GO AND TB MODELS
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
TB_Xtick = all_SU_data(1).sacStimProc.gsac_TBmod{1}.lagX*dt;
base_lags = find(slags <= 0);
TBbase_lags = find(TB_Xtick <= 0);

GO_color = [0.25 0.8 0.25];

TB_lambda_ii = 4;
GO_lambda_off = 4;
GO_lambda_gain = 3;

TB_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.ovInfo,all_SU_data(cur_SUs)); %overall TB model infos
TB_NSSI = bsxfun(@rdivide,TB_SSI,TB_ovinfos); %normalize TB SSI by overall model info
TB_NSSI = bsxfun(@rdivide,TB_NSSI,mean(TB_NSSI(:,base_lags),2));

TB_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
TB_ovLL = arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.LLimp,all_SU_data(cur_SUs)); %overall TB model infos
TB_NLL = bsxfun(@rdivide,TB_LL,TB_ovLL); %normalize TB SSI by overall model info
TB_NLL = bsxfun(@rdivide,TB_NLL,mean(TB_NLL(:,base_lags),2));

TB_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
TB_Noffset = bsxfun(@rdivide,TB_offset,gsac_ov_rates); %normalize offset by overall avg rates
TB_Noffset = bsxfun(@minus,TB_Noffset,mean(TB_Noffset(:,TBbase_lags),2));

TB_gains = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_TBmod{TB_lambda_ii}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
TB_gains = bsxfun(@rdivide,TB_gains,mean(TB_gains(:,TBbase_lags),2));


GO_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
GO_Noffset = bsxfun(@rdivide,GO_offset,gsac_ov_rates);
GO_Noffset = bsxfun(@minus,GO_Noffset,mean(GO_Noffset(:,base_lags),2));

GO_gain = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
GO_gain = bsxfun(@rdivide,GO_gain,mean(GO_gain(:,base_lags),2));

GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info
GO_NSSI = bsxfun(@rdivide,GO_NSSI,mean(GO_NSSI(:,base_lags),2));

GO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovLL = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
GO_NLL = bsxfun(@rdivide,GO_LL,GO_ovLL); %normalize TB SSI by overall model info
GO_NLL = bsxfun(@rdivide,GO_NLL,mean(GO_NLL(:,base_lags),2));

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

% xl = [-0.1 0.3];
% f1 = figure();hold on
% h1=shadedErrorBar(slags*dt,nanmean(TB_NLL),nanstd(TB_NLL)/sqrt(length(cur_SUs)),{'color','b'});
% h2=shadedErrorBar(slags*dt,nanmean(GO_NLL),nanstd(GO_NLL)/sqrt(length(cur_SUs)),{'color','r'});
% % plot(tlags,all_gsac_gray(set2,:),'k')
% line(xl,[1 1],'color','k');
% line([0 0],ylim(),'color','k');
% xlim(xl);
% xlabel('Time (s)');
% ylabel('Log-likelihood');

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

GO_lambda_off = 4;
GO_lambda_gain = 3;

all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
gsac_ov_rates = arrayfun(@(x) x.sacStimProc.gsac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates
GO_offset = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
GO_Noffset = bsxfun(@rdivide,GO_offset,gsac_ov_rates);
GO_Noffset = bsxfun(@minus,GO_Noffset,mean(GO_Noffset(:,base_lags),2));

GO_gain = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
GO_gain = bsxfun(@rdivide,GO_gain,mean(GO_gain(:,base_lags),2));

GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info
GO_NSSI = bsxfun(@rdivide,GO_NSSI,mean(GO_NSSI(:,base_lags),2));

GO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovLL = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
GO_NLL = bsxfun(@rdivide,GO_LL,GO_ovLL); %normalize TB SSI by overall model info
GO_NLL = bsxfun(@rdivide,GO_NLL,mean(GO_NLL(:,base_lags),2));

all_msac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.msac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
msac_ov_rates = arrayfun(@(x) x.sacStimProc.msac_ovavg_rate,all_SU_data(cur_SUs)); %overall average rates

mGO_offset = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_Noffset = bsxfun(@rdivide,mGO_offset,msac_ov_rates);
mGO_Noffset = bsxfun(@minus,mGO_Noffset,mean(mGO_Noffset(:,base_lags),2));

mGO_gain = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_gain = bsxfun(@rdivide,mGO_gain,mean(mGO_gain(:,base_lags),2));

mGO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_ovinfos = arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
mGO_NSSI = bsxfun(@rdivide,mGO_SSI,mGO_ovinfos); %normalize GO SSI by overall model info
% mGO_NSSI = bsxfun(@rdivide,mGO_NSSI,mean(mGO_NSSI(:,base_lags),2));

mGO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
mGO_ovLL = arrayfun(@(x) x.sacStimProc.msac_post_mod{GO_lambda_off,GO_lambda_gain}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
mGO_NLL = bsxfun(@rdivide,mGO_LL,mGO_ovLL); %normalize TB SSI by overall model info
mGO_NLL = bsxfun(@rdivide,mGO_NLL,mean(mGO_NLL(:,base_lags),2));

tlags_up = linspace(tlags(1),tlags(end),500);
all_gsac_tavg_up = spline(tlags,all_gsac_tavg,tlags_up);
all_msac_tavg_up = spline(tlags,all_msac_tavg,tlags_up);

search_range = [0 0.3];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_tavg_up-1),tlags_up,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_tavg_up-1,tlags_up,search_range);
[msac_Ifact,msac_inhtime] = get_tavg_peaks(-(all_msac_tavg_up-1),tlags_up,search_range);
[msac_Efact,msac_exctime] = get_tavg_peaks(all_msac_tavg_up-1,tlags_up,search_range);

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
% h1=shadedErrorBar(slags*dt,nanmean(GO_NSSI),nanstd(GO_NSSI)/sqrt(length(cur_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,nanmean(mGO_NSSI),nanstd(mGO_NSSI)/sqrt(length(cur_SUs)),{'color','r'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative SSI');
% ylim([0.5 1.1]);
ylim([0.7 1.1]);

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
figufy(f2);
% fname = [fig_dir 'Gsac_MSAC_SSI.pdf'];
fname = [fig_dir 'Gsac_MSAC_SSI_compare2_unCor.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);
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
cur_min_Nsacs = 250;
cur_SUs = find(avg_rates >= min_rate & N_gsacs_gray >= cur_min_Nsacs & N_gsacs_im >= cur_min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);

GO_lambda_off = 4;
GO_lambda_gain = 3;

all_gsac_gr_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_gray_avg', all_SU_data(cur_SUs),'uniformoutput',0));
GO_gr_SSI = cell2mat(arrayfun(@(x) x.type_dep.gsacGR_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_gr_ovinfos = arrayfun(@(x) x.type_dep.gsacGR_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_gr_NSSI = bsxfun(@rdivide,GO_gr_SSI,GO_gr_ovinfos); %normalize GO SSI by overall model info
GO_gr_NSSI = bsxfun(@rdivide,GO_gr_NSSI,mean(GO_gr_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info

all_gsac_im_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_im_avg', all_SU_data(cur_SUs),'uniformoutput',0));
GO_im_SSI = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_im_ovinfos = arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_im_NSSI = bsxfun(@rdivide,GO_im_SSI,GO_im_ovinfos); %normalize GO SSI by overall model info
GO_im_NSSI = bsxfun(@rdivide,GO_im_NSSI,mean(GO_im_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info

search_range = [0 0.15];
[GR_SSI_Ifact,GR_SSI_inhtime] = get_tavg_peaks(-(GO_gr_NSSI-1),slags*dt,search_range);
[IM_SSI_Ifact,IM_SSI_inhtime] = get_tavg_peaks(-(GO_im_NSSI-1),slags*dt,search_range);

tlags_up = linspace(tlags(1),tlags(end),500);
all_gsac_gr_tavg_up = spline(tlags,all_gsac_gr_tavg,tlags_up);
all_gsac_im_tavg_up = spline(tlags,all_gsac_im_tavg,tlags_up);
search_range = [0 0.3];
[GR_Efact,GR_exctime] = get_tavg_peaks((all_gsac_gr_tavg_up-1),tlags_up,search_range);
[IM_Efact,IM_exctime] = get_tavg_peaks((all_gsac_im_tavg_up-1),tlags_up,search_range);
[GR_Ifact,GR_inhtime] = get_tavg_peaks(-(all_gsac_gr_tavg_up-1),tlags_up,search_range);
[IM_Ifact,IM_inhtime] = get_tavg_peaks(-(all_gsac_im_tavg_up-1),tlags_up,search_range);

IM_ov_rates = arrayfun(@(x) x.type_dep.gsac_IM_ovavg_rate,all_SU_data(cur_SUs));
IM_offset = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
IM_Noffset = bsxfun(@rdivide,IM_offset,IM_ov_rates);
IM_Noffset = bsxfun(@minus,IM_Noffset,mean(IM_Noffset(:,base_lags),2));

GR_ov_rates = arrayfun(@(x) x.type_dep.gsac_GR_ovavg_rate,all_SU_data(cur_SUs));
GR_offset = cell2mat(arrayfun(@(x) x.type_dep.gsacGR_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
GR_Noffset = bsxfun(@rdivide,GR_offset,GR_ov_rates);
GR_Noffset = bsxfun(@minus,GR_Noffset,mean(GR_Noffset(:,base_lags),2));

IM_gain = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
IM_gain = bsxfun(@rdivide,IM_gain,mean(IM_gain(:,base_lags),2));

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
ylim([0.6 1.3])

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
ylim([0.5 1.1])

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
ylim([-0.05 0.3])

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

%% COMPARE REAL AND SIMULATED SACCADES
cur_min_Nsacs = 250;
cur_SUs = find(avg_rates >= min_rate & N_gsacs_im >= cur_min_Nsacs & N_simsacs >= cur_min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);

GO_lambda_off = 4;
GO_lambda_gain = 3;

all_gsac_im_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_im_avg', all_SU_data(cur_SUs),'uniformoutput',0));
GO_im_SSI = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_im_ovinfos = arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_im_NSSI = bsxfun(@rdivide,GO_im_SSI,GO_im_ovinfos); %normalize GO SSI by overall model info
GO_im_NSSI = bsxfun(@rdivide,GO_im_NSSI,mean(GO_im_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info

all_simsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.simsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));
GO_sim_SSI = cell2mat(arrayfun(@(x) x.type_dep.simsac_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_sim_ovinfos = arrayfun(@(x) x.type_dep.simsac_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_sim_NSSI = bsxfun(@rdivide,GO_sim_SSI,GO_sim_ovinfos); %normalize GO SSI by overall model info
GO_sim_NSSI = bsxfun(@rdivide,GO_sim_NSSI,mean(GO_sim_NSSI(:,base_lags),2)); %normalize GO SSI by overall model info

search_range = [0 0.15];
[IM_SSI_Ifact,IM_SSI_inhtime] = get_tavg_peaks(-(GO_im_NSSI-1),slags*dt,search_range);
[SIM_SSI_Ifact,SIM_SSI_inhtime] = get_tavg_peaks(-(GO_sim_NSSI-1),slags*dt,search_range);

tlags_up = linspace(tlags(1),tlags(end),500);
all_simsac_tavg_up = spline(tlags,all_simsac_tavg,tlags_up);
all_gsac_im_tavg_up = spline(tlags,all_gsac_im_tavg,tlags_up);
search_range = [0 0.3];
[IM_Ifact,IM_inhtime] = get_tavg_peaks(-(all_gsac_im_tavg_up-1),tlags_up,search_range);
[SIM_Ifact,SIM_inhtime] = get_tavg_peaks(-(all_simsac_tavg_up-1),tlags_up,search_range);
[IM_Efact,IM_exctime] = get_tavg_peaks((all_gsac_im_tavg_up-1),tlags_up,search_range);
[SIM_Efact,SIM_exctime] = get_tavg_peaks((all_simsac_tavg_up-1),tlags_up,search_range);

IM_ov_rates = arrayfun(@(x) x.type_dep.gsac_IM_ovavg_rate,all_SU_data(cur_SUs));
IM_offset = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
IM_Noffset = bsxfun(@rdivide,IM_offset,IM_ov_rates);
IM_Noffset = bsxfun(@minus,IM_Noffset,mean(IM_Noffset(:,base_lags),2));

SIM_ov_rates = arrayfun(@(x) x.type_dep.simsac_ovavg_rate,all_SU_data(cur_SUs));
SIM_offset = cell2mat(arrayfun(@(x) x.type_dep.simsac_mod{GO_lambda_off,GO_lambda_gain}.sac_offset',all_SU_data(cur_SUs),'uniformoutput',0));
SIM_Noffset = bsxfun(@rdivide,SIM_offset,SIM_ov_rates);
SIM_Noffset = bsxfun(@minus,SIM_Noffset,mean(SIM_Noffset(:,base_lags),2));

IM_gain = cell2mat(arrayfun(@(x) x.type_dep.gsacIM_mod{GO_lambda_off,GO_lambda_gain}.sac_gain',all_SU_data(cur_SUs),'uniformoutput',0));
IM_gain = bsxfun(@rdivide,IM_gain,mean(IM_gain(:,base_lags),2));

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
ylim([0.6 1.3])

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
ylim([0.5 1.25])

%compare offsets
f3 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(IM_Noffset),nanstd(IM_Noffset)/sqrt(length(cur_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,nanmean(SIM_Noffset),nanstd(SIM_Noffset)/sqrt(length(cur_SUs)),{'color','b'});
line(xl,[0 0],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative Offset');
ylim([-0.1 0.25]);

%compare gains
f4 = figure();hold on
h1=shadedErrorBar(slags*dt,nanmean(IM_gain),nanstd(IM_gain)/sqrt(length(cur_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,nanmean(SIM_gain),nanstd(SIM_gain)/sqrt(length(cur_SUs)),{'color','b'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Gain');
ylim([0.5 1.2]);

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
base_lags = find(slags <= -0.025);
close all

post_lambda_off = 4;
pre_lambda =3;
post_lambda_gain = 3;

gsac_pre_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.stim_kernel',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_post_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{post_lambda_off,post_lambda_gain}.mods(3).filtK',all_SU_data(cur_SUs),'uniformoutput',0));

gsac_post_gain = bsxfun(@rdivide,gsac_post_gain,mean(gsac_post_gain(:,base_lags),2));
gsac_pre_gain = bsxfun(@rdivide,gsac_pre_gain,mean(gsac_pre_gain(:,base_lags),2));

gsac_pre_off =cell2mat(arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.off_kernel',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_post_off =cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{post_lambda_off,post_lambda_gain}.mods(2).filtK',all_SU_data(cur_SUs),'uniformoutput',0));

gsac_pre_off = bsxfun(@minus,gsac_pre_off,mean(gsac_pre_off(:,base_lags),2));
gsac_post_off = bsxfun(@minus,gsac_post_off,mean(gsac_post_off(:,base_lags),2));

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
base_lags = find(slags <= -0.025);

post_lambda_off = 4;
pre_lambda = 3;
post_lambda_gain = 3;

gsac_pre_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.stim_kernel',all_SU_data(cur_SUs),'uniformoutput',0));
msac_pre_gain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.msacPreGainMod{pre_lambda}.stim_kernel',all_SU_data(cur_SUs),'uniformoutput',0));

msac_pre_gain = bsxfun(@rdivide,msac_pre_gain,mean(msac_pre_gain(:,base_lags),2));
gsac_pre_gain = bsxfun(@rdivide,gsac_pre_gain,mean(gsac_pre_gain(:,base_lags),2));

gsac_base_LLs = arrayfun(@(x) x.sacStimProc.gsac_base_LLimp,all_SU_data(cur_SUs));
gsac_off_LLs = arrayfun(@(x) x.sacStimProc.gsac_off_mod{post_lambda_off}.ovLLimp,all_SU_data(cur_SUs));
gsac_post_LLs = arrayfun(@(x) x.sacStimProc.gsac_post_mod{post_lambda_off,post_lambda_gain}.ovLLimp,all_SU_data(cur_SUs));
gsac_pre_LLs = arrayfun(@(x) x.sacStimProc.gsacPreGainMod{pre_lambda}.ovLLimp,all_SU_data(cur_SUs));

msac_base_LLs = arrayfun(@(x) x.sacStimProc.msac_base_LLimp,all_SU_data(cur_SUs));
msac_off_LLs = arrayfun(@(x) x.sacStimProc.msac_off_mod{post_lambda_off}.ovLLimp,all_SU_data(cur_SUs));
msac_post_LLs = arrayfun(@(x) x.sacStimProc.msac_post_mod{post_lambda_off,post_lambda_gain}.ovLLimp,all_SU_data(cur_SUs));
msac_pre_LLs = arrayfun(@(x) x.sacStimProc.msacPreGainMod{pre_lambda}.ovLLimp,all_SU_data(cur_SUs));

gsac_rel_LLimps = (gsac_pre_LLs - gsac_post_LLs)./(gsac_pre_LLs - gsac_off_LLs);
msac_rel_LLimps = (msac_pre_LLs - msac_post_LLs)./(msac_pre_LLs - msac_off_LLs);

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

all_Ekerns = [];
all_Ikerns = [];
for ii = 1:length(cur_SUs)
    [~,cur_Ekern,cur_Ikern] = get_hilbert_tempkerns(all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM);
    all_Ekerns = cat(1,all_Ekerns,cur_Ekern');
    all_Ikerns = cat(1,all_Ikerns,cur_Ikern');
end

all_gsac_tavg = cell2mat(arrayfun(@(x) x.trig_avg.gsac_avg', all_SU_data(cur_SUs),'uniformoutput',0));

tlags_up = linspace(tlags(1),tlags(end),500);
all_gsac_up = spline(tlags,all_gsac_tavg,tlags_up);

search_range = [0 0.3];
[gsac_Ifact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_up-1),tlags_up,search_range);
search_range = [0 0.3];
[gsac_Efact,gsac_exctime] = get_tavg_peaks((all_gsac_up-1),tlags_up,search_range);

include = find(gsac_inhtime < gsac_exctime);


[earlier_timing,peak_earlier] = min([gsac_inhtime gsac_exctime],[],2);
flen = 15;
lag_ax = ((1:flen)*dt - dt/2);
up_lagax = linspace(lag_ax(1),lag_ax(end),500);
all_Ekerns_up = spline(lag_ax,all_Ekerns,up_lagax);
all_Ikerns_up = spline(lag_ax,all_Ikerns,up_lagax);
search_range = [0 max(up_lagax)];
[Ekern_max,Ekern_time] = get_tavg_peaks(all_Ekerns_up,up_lagax,search_range);
[Ikern_max,Ikern_time] = get_tavg_peaks(all_Ikerns_up,up_lagax,search_range);

% cur_sigI = find(ismember(cur_SUs,sigI));
% cur_sigE = find(ismember(cur_SUs,sigE));
% cur_rP = find(ismember(cur_SUs,reverse_polarity_units));
% cur_O = false(size(cur_SUs));
% cur_O(cur_sigI(peak_earlier(cur_sigI) == 1)) = true;
% cur_O(cur_sigE(peak_earlier(cur_sigE) == 1)) = true;

xl = [0.02 0.08];
yl = [0.02 0.16];

mS = 3;
f1 = figure(); hold on
% plot(Ekern_time,gsac_inhtime,'b.','markersize',12)
% plot(Ekern_time(cur_sigI),gsac_inhtime(cur_sigI),'r.');
plot(Ekern_time(include),gsac_inhtime(include),'b.','markersize',12);
% plot(Ekern_time,earlier_timing,'b.','markersize',12)
% plot(Ekern_time(cur_O),earlier_timing(cur_O),'r.','markersize',12)
% plot(Ekern_time(cur_sigI),gsac_inhtime(cur_sigI),'r.');
% plot(Ekern_time(cur_rP),gsac_inhtime(cur_rP),'ro','markersize',12);
xlim(xl); ylim(yl);
% line(xl,xl,'color','k')
r = robustfit(Ekern_time,gsac_inhtime);
xx = linspace(xl(1),xl(2),50);
plot(xx,r(1)+r(2)*xx,'r')
xlabel('Stim latency (s)');
ylabel('Suppression timing (s)');

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
% fname = [fig_dir 'Supp_Stim_timing_scatter.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% COMPARE TEMPKERNS ADN MODEL GAIN LATENCIES
cur_SUs = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs & mod_xvLLimps > min_xvLLimp);
base_lags = find(slags <= 0);
noise_lags = find(slags <= 0 | slags*dt > 0.15);
close all

lambda_L2_ii = 2;
lambda_d2T_ii = 3;
GO_lambda_off = 4;
GO_lambda_gain = 3;

flen = 15;
lag_ax = ((1:flen)*dt - dt/2);
up_lagax = linspace(lag_ax(1),lag_ax(end),500);
slags_up = linspace(slags(1)*dt,slags(end)*dt,500);
jit_amp = 0.000;%jitter data points with gaussian RV to avoid point occlusion

search_range = [0 0.15];
tsearch_range = [0 max(up_lagax)];
noise_SD_thresh = 0;

all_tempkerns = [];
all_gainkerns = [];
all_relweights = [];
all_modsigns = [];
cell_tg_slope = nan(length(cur_SUs),1);
for ii = 1:length(cur_SUs)
    %get stim-filter temporal kernels
    cur_tkerns = get_hilbert_tempkerns(all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM);
    %get gain kerels from full model
    cur_gainkerns = reshape([all_SU_data(cur_SUs(ii)).sacStimProc.gsac_post_Fullmod{1,lambda_d2T_ii,lambda_L2_ii}.mods(3).filtK],length(slags),[])';
    cur_relweights = all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM.rel_filt_weights;
    cur_modsigns = [all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM.mods(:).sign];
    
    clear cur_tkerns_up
    uset = find(cur_relweights > 0);
    cur_tkerns_up(uset,:) = spline(lag_ax,cur_tkerns(:,uset)',up_lagax); %spline interpoltae temporal kernels
    [~,cur_tk_time] = get_tavg_peaks(cur_tkerns_up,up_lagax,tsearch_range); %get peaks of temporal kernels
    cur_gkerns_up = spline(slags*dt,cur_gainkerns,slags_up); %spline interpolate gain kernels
    [cur_gkern_max,cur_gk_time] = get_tavg_peaks(-(cur_gkerns_up-1),slags_up,search_range); %get peaks of gain kernels
    cur_noise_level = std(cur_gainkerns(:,noise_lags),[],2); %compute 'background noise' of gain kernels
    uset = uset(cur_gkern_max(uset) >= cur_noise_level(uset)*noise_SD_thresh); %only take filters with sufficient peak SNR
    %fit within-cell robust linear reg between temporal kernel and gain
    %kernel peak timing
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
ukerns = find(all_relweights > 0 & gkern_max > noise_SD_thresh*noise_level & ~isnan(tkern_max));
[a,b] = corr(gkern_time(ukerns),tkern_time(ukerns),'type','spearman')

%E and I subunits
esubs = ukerns(all_modsigns(ukerns)==1);
isubs = ukerns(all_modsigns(ukerns)==-1);

gkern_time = gkern_time + randn(size(gkern_time))*jit_amp;
tkern_time = tkern_time + randn(size(tkern_time))*jit_amp;

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

[a,b] = corr(gkern_time(esubs),tkern_time(esubs),'type','spearman')
[a,b] = corr(gkern_time(isubs),tkern_time(isubs),'type','spearman')


xl = [-0.1 0.3];
f2 = figure();hold on
shadedErrorBar(slags*dt,mean(all_gainkerns),std(all_gainkerns)/sqrt(size(all_gainkerns,1)),{'color','k'});
line(xl,[1 1],'color','k');
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
% fig_width = 3.5; rel_height = 0.8;
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
flen = 15;
lag_ax = ((1:flen)*dt - dt/2)*1e3;
up_lagax = linspace(lag_ax(1),lag_ax(end),500);
slags_up = linspace(slags(1)*dt,slags(end)*dt,500);

base_lags = find(slags <= 0);
noise_lags = find(slags <= 0 | slags*dt > 0.15);
noise_SD_thresh = 0;
search_range = [0 0.15];
tsearch_range = up_lagax([1 end]);

lambda_ii = 3;

gsac_EIgain = 1+cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_EImod{lambda_ii}.mods(3).filtK',all_SU_data(cur_SUs),'uniformoutput',0));
gsac_Egain = gsac_EIgain(:,1:length(slags));
gsac_Igain = gsac_EIgain(:,(length(slags)+1):end);

gsac_Egain = bsxfun(@rdivide,gsac_Egain,mean(gsac_Egain(:,base_lags),2));
gsac_Igain = bsxfun(@rdivide,gsac_Igain,mean(gsac_Igain(:,base_lags),2));

flen = 15;
all_Ekerns = nan(length(cur_SUs),flen);
all_Ikerns = nan(length(cur_SUs),flen);
all_NIfilts = nan(length(cur_SUs),1);
all_NEfilts = nan(length(cur_SUs),1);
for ii = 1:length(cur_SUs)
    [cur_tkerns,cur_ekern,cur_ikern] = get_hilbert_tempkerns(all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM);
    all_Ekerns(ii,:) = cur_ekern;
    all_Ikerns(ii,:) = cur_ikern;
    cur_modsigns = [all_SU_data(cur_SUs(ii)).sacStimProc.ModData.rectGQM.mods(:).sign];
    all_NIfilts(ii) = sum(cur_modsigns==-1);
    all_NEfilts(ii) = sum(cur_modsigns==1);
end

all_Ekerns_up = spline(lag_ax,all_Ekerns,up_lagax);
all_Ikerns_up = spline(lag_ax,all_Ikerns,up_lagax);
[Ekern_max,Ekern_time] = get_tavg_peaks(all_Ekerns_up,up_lagax,tsearch_range);
[Ikern_max,Ikern_time] = get_tavg_peaks(all_Ikerns_up,up_lagax,tsearch_range);

gsac_Egain_up = spline(slags*dt,gsac_Egain,slags_up);
gsac_Igain_up = spline(slags*dt,gsac_Igain,slags_up);
[Egain_max,Egain_time] = get_tavg_peaks(-(gsac_Egain_up-1),slags_up,search_range);
[Igain_max,Igain_time] = get_tavg_peaks(-(gsac_Igain_up-1),slags_up,search_range);
Enoise_level =std(gsac_Egain(:,noise_lags),[],2);
Inoise_level =std(gsac_Igain(:,noise_lags),[],2);
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

% % %PRINT PLOTS
fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'EI_gainkerns.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
figure(f1);
xlim(xlz);
fname = [fig_dir 'EI_gainkerns_zoom.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'EI_tempkerns.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

figufy(f3);
fname = [fig_dir 'EI_timing_scatter.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);

figufy(f5);
fname = [fig_dir 'EI_timing_dists.pdf'];
exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f5);

%% ANALYZE LAMINAR DEPENDENCIES WITH MUA
mod_dt = 0.01;
flen = 15;
tax = (0:(flen-1))*mod_dt + mod_dt/2;
up_tax = linspace(tax(1),tax(end),500);
tsearch_range = up_tax([1 end]);
search_range = [0 0.3];
xl = [-0 0.25];

mua_sm = 0.005/trig_avg_params.dt; %smoothing sigma for MUA rates

%load info about layer boundaries
load('~/Analysis/bruce/FINsac_mod/layer_boundaries/layer_classification.mat')
boundary_enums = [boundary_class(:).Expt_num];
all_lbs = [boundary_class(:).lb];
all_ubs = [boundary_class(:).ub];
mean_lb = mean(all_lbs);
mean_ub = mean(all_ubs);
mean_diff = mean_lb-mean_ub;
interp_ax = 1:1:24;

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
            load(sname);
            tname = strcat(sac_dir,base_tname,sprintf('_ori%d',cur_ori));
            load(tname);
            expt_MUA_stas(kk,ii,:,:) = mua_data.gsac_avg';
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
all_SU_lclass = nan(length(all_SU_data),1);
all_MU_lclass = nan(n_lem_expts,n_probes);
all_MU_recID = nan(n_lem_expts,n_probes);
% all_interp_mua_sta = nan(n_lem_expts,length(interp_ax),length(tlags));
% all_interp_mua_Ekern = nan(n_lem_expts,length(interp_ax),length(up_tax));
all_mua_Ekerns = [];
for ee = 1:n_lem_expts
    %boundary info for current expt
    cur_enum = str2num(Elist{ee}(2:end));
    cur_bound_info = find(boundary_enums == cur_enum,1);
    cur_ub = boundary_class(cur_bound_info).ub;
    cur_lb = boundary_class(cur_bound_info).lb;
    
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
% [p,table,stats] = anovan(MUA_Ttime,group,'random',2);
% [p,table,stats] = anovan(MUA_Ifact,group,'random',2);
% [p,table,stats] = anovan(MUA_Etime,group,'random',2);
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

GO_lambda_gain = 3;
GO_lambda_off = 4;
sub_lambda = 3;

GO_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovInfo,all_SU_data(cur_SUs));
GO_NSSI = bsxfun(@rdivide,GO_SSI,GO_ovinfos); %normalize GO SSI by overall model info
GO_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
GO_ovLL = arrayfun(@(x) x.sacStimProc.gsac_post_mod{GO_lambda_off,GO_lambda_gain}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
GO_NLL = bsxfun(@rdivide,GO_LL,GO_ovLL); %normalize TB SSI by overall model info

sub_SSI = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_subMod{sub_lambda}.sac_modinfo',all_SU_data(cur_SUs),'uniformoutput',0));
sub_ovinfos = arrayfun(@(x) x.sacStimProc.gsac_subMod{sub_lambda}.ovInfo,all_SU_data(cur_SUs));
sub_NSSI = bsxfun(@rdivide,sub_SSI,sub_ovinfos); %normalize GO SSI by overall model info
sub_LL = cell2mat(arrayfun(@(x) x.sacStimProc.gsac_subMod{sub_lambda}.sac_LLimp',all_SU_data(cur_SUs),'uniformoutput',0));
sub_ovLL = arrayfun(@(x) x.sacStimProc.gsac_subMod{sub_lambda}.ovLLimp,all_SU_data(cur_SUs)); %overall TB model infos
sub_NLL = bsxfun(@rdivide,sub_LL,sub_ovLL); %normalize TB SSI by overall model info

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
N_msacs_Par_sub = arrayfun(@(x) x.trig_avg.N_msacs_Par_sub,all_SU_data);
N_msacs_Orth_sub = arrayfun(@(x) x.trig_avg.N_msacs_Orth_sub,all_SU_data);

cur_SUs = find(avg_rates >= min_rate & N_msacs >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
% cur_SUs = find(avg_rates >= min_rate & N_msacs_Par_sub >= min_TA_Nsacs & N_msacs_Orth_sub >= min_TA_Nsacs & mod_xvLLimps > min_xvLLimp);
% all_msac_par = cell2mat(arrayfun(@(x) x.trig_avg.msac_Par_avg', all_SU_data(cur_SUs),'uniformoutput',0));
% all_msac_orth = cell2mat(arrayfun(@(x) x.trig_avg.msac_Orth_avg', all_SU_data(cur_SUs),'uniformoutput',0));
all_msac_par = cell2mat(arrayfun(@(x) x.trig_avg.msac_Par_sub_avg', all_SU_data(cur_SUs),'uniformoutput',0));
all_msac_orth = cell2mat(arrayfun(@(x) x.trig_avg.msac_Orth_sub_avg', all_SU_data(cur_SUs),'uniformoutput',0));

tlags_up = linspace(tlags(1),tlags(end),500);
all_msac_par_up = spline(tlags,all_msac_par,tlags_up);
all_msac_orth_up = spline(tlags,all_msac_orth,tlags_up);

search_range = [0 0.3];
[par_Ifact,par_inhtime] = get_tavg_peaks(-(all_msac_par_up-1),tlags_up,search_range);
[orth_Ifact,orth_inhtime] = get_tavg_peaks(-(all_msac_orth_up-1),tlags_up,search_range);

search_range = [0 0.3];
[par_Efact,par_exctime] = get_tavg_peaks((all_msac_par_up-1),tlags_up,search_range);
[orth_Efact,orth_exctime] = get_tavg_peaks((all_msac_orth_up-1),tlags_up,search_range);

xl = [-0.1 0.4];

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
% fname = [fig_dir 'MUA_par_orth2.pdf'];
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
