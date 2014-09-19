% close all
clear all
clc

fit_unCor = 0;

fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';
base_tname = 'sac_trig_avg_data';
base_sname = 'sacStimProcFin';

all_SU_data = [];
all_SU_NPdata = [];
all_SU_tdata = [];
all_SU_timedata = [];

%% LOAD JBE
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
n_probes = 96;
ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
rmfield_list = {};

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            sname = strcat(sac_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            temp = load(sname);
            tname = strcat(sac_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            if fit_unCor
                tname = strcat(tname,'_unCor');
            end
            load(tname);            
                       
            %SU numbers for trig avg data
            tavg_SU_numbers = [temp.sua_data(:).SU_numbers];
            
            %find units where we have computed sacStimMod data
            ucells = arrayfun(@(x) length(x.gsac_avg_rate),sacStimProc) > 0;
            sua_data = sacStimProc(ucells);
            SM_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,sua_data);
            SM_SU_xvLLimp = arrayfun(@(x) x.ModData.rectGQM.xvLLimp,sua_data);
            SM_SU_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,sua_data);
            
            %find which SUs we have sacMod data for
            [lia,locb] = ismember(tavg_SU_numbers,SM_SU_numbers);
            lia_inds = find(lia);
            base_xvLLimps = ones(size(tavg_SU_numbers))*-Inf;
            base_xvLLimps(lia) = SM_SU_xvLLimp(locb(lia)).*SM_SU_rates(locb(lia));
            for jj = 1:length(sua_data);
                sua_data(jj).xvLLimp = base_xvLLimps(lia_inds(jj));
                sua_data(jj).N_gsacs = temp.sua_data(lia_inds(jj)).N_gsacs;
                sua_data(jj).N_msacs = temp.sua_data(lia_inds(jj)).N_msacs;
            end
            
            [sua_data.expt_num] = deal(Expt_num);
            [sua_data.bar_ori] = deal(ori_list(ee,ii));
            [sua_data.animal] = deal('jbe');
            
            ori_SU_nums(ii,:) = tavg_SU_numbers;
            ori_xvLLimps(ii,:) = base_xvLLimps;
            ori_sua_data{ii}(lia) = sua_data;
            ori_tavg_data{ii} = temp.sua_data;
            clear ModData
            
        end
    end
    %if there was only one ori, set the xvLLimps to -Inf for a placeholder
    if size(ori_xvLLimps,1) == 1
        ori_xvLLimps = [ori_xvLLimps; ones(size(ori_xvLLimps))*-Inf];
    end
    [mvals,mlocs] = max(ori_xvLLimps,[],1);
    nplocs = mod(mlocs,2)+1; %these are indices for the NP ori
    for ii = 1:length(mvals)
        if mvals(ii) > 0
            if ori_xvLLimps(nplocs(ii),ii) > 0 %if the NP ori had usable data
                all_SU_NPdata = cat(1,all_SU_NPdata,ori_sua_data{nplocs(ii)}(ii));
                has_NP = true;
            else
                has_NP = false;
            end
            cur_struct = ori_sua_data{mlocs(ii)}(ii);
            cur_struct.has_NP = has_NP;
            all_SU_data = cat(1,all_SU_data,cur_struct);
            all_SU_tdata = cat(1,all_SU_tdata,ori_tavg_data{mlocs(ii)}(ii));
        end
    end
    
    clear ori_SU_nums ori_xvLLimps ori_sua_data ori_tavg_data
end

clear ori_*
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
    
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            sname = strcat(sac_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            temp = load(sname);
            tname = strcat(sac_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            if fit_unCor
                tname = strcat(tname,'_unCor');
            end
            load(tname);
            
            tavg_SU_numbers = [temp.sua_data(:).SU_numbers];
            
            ucells = arrayfun(@(x) length(x.gsac_avg_rate),sacStimProc) > 0;
            sua_data = sacStimProc(ucells);
            SM_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,sua_data);
            SM_SU_xvLLimp = arrayfun(@(x) x.ModData.rectGQM.xvLLimp,sua_data);
            SM_SU_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,sua_data);
            
            [lia,locb] = ismember(tavg_SU_numbers,SM_SU_numbers);
            lia_inds = find(lia);
            base_xvLLimps = ones(size(tavg_SU_numbers))*-Inf;
            base_xvLLimps(lia) = SM_SU_xvLLimp(locb(lia)).*SM_SU_rates(locb(lia));
            for jj = 1:length(sua_data);
                sua_data(jj).xvLLimp = base_xvLLimps(lia_inds(jj));
                sua_data(jj).N_gsacs = temp.sua_data(lia_inds(jj)).N_gsacs;
                sua_data(jj).N_msacs = temp.sua_data(lia_inds(jj)).N_msacs;
            end
            
            [sua_data.expt_num] = deal(Expt_num);
            [sua_data.bar_ori] = deal(ori_list(ee,ii));
            [sua_data.animal] = deal('lem');
            
            ori_SU_nums(ii,:) = tavg_SU_numbers;
            ori_xvLLimps(ii,:) = base_xvLLimps;
            ori_sua_data{ii}(lia) = sua_data;
            ori_tavg_data{ii} = temp.sua_data;
            clear ModData
            
        end
    end
    %if there was only one ori, set the xvLLimps to -Inf for a placeholder
    if size(ori_xvLLimps,1) == 1
        ori_xvLLimps = [ori_xvLLimps; ones(size(ori_xvLLimps))*-Inf];
    end
    [mvals,mlocs] = max(ori_xvLLimps,[],1);
    nplocs = mod(mlocs,2)+1; %these are indices for the NP ori
    for ii = 1:length(mvals)
        if mvals(ii) > 0
            if ori_xvLLimps(nplocs(ii),ii) > 0 %if the NP ori had usable data
                all_SU_NPdata = cat(1,all_SU_NPdata,ori_sua_data{nplocs(ii)}(ii));
                has_NP = true;
            else
                has_NP = false;
            end
            cur_struct = ori_sua_data{mlocs(ii)}(ii);
            cur_struct.has_NP = has_NP;
            all_SU_data = cat(1,all_SU_data,cur_struct);
            all_SU_tdata = cat(1,all_SU_tdata,ori_tavg_data{mlocs(ii)}(ii));
        end
    end
    
    clear ori_SU_nums ori_xvLLimps ori_sua_data ori_tavg_data 
end

%% time axis for trig-avg data
tlags = temp.trig_avg_params.lags;
Tdt = temp.trig_avg_params.dt;
tlags = tlags*Tdt;

%% SELECT USABLE CELLS
%selection criteria
min_rate = 5; %in Hz (5)
min_Nsacs = 1e3; % (1e3)
min_TA_Nsacs = 250;
% min_xvLLimp = 0.05; %(0.05);

tot_Nunits = length(all_SU_data);
avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_SU_data)';
tot_spikes = arrayfun(@(x) x.ModData.unit_data.tot_spikes,all_SU_data)';
rec_dur = arrayfun(@(x) x.ModData.unit_data.N_used_samps,all_SU_data)'*dt/60;
expt_nums = [all_SU_data(:).expt_num];
expt_oris = [all_SU_data(:).bar_ori];
% xvLLimps = [all_SU_data(:).xvLLimp];
% mod_LLimps = [all_SU_data(:).gsac_spost_ov_modinfo];
mod_LLimps = [all_SU_data(:).xvLLimp];

clust_iso_dist = arrayfun(@(x) x.ModData.unit_data.SU_isodist,all_SU_data)';
clust_Lratio = arrayfun(@(x) x.ModData.unit_data.SU_Lratio,all_SU_data)';
clust_refract = arrayfun(@(x) x.ModData.unit_data.SU_refract,all_SU_data)';
clust_dprime = arrayfun(@(x) x.ModData.unit_data.SU_dprime,all_SU_data)';
rate_stability_cv = arrayfun(@(x) x.ModData.unit_data.rate_stability_cv,all_SU_data)';
dprime_stability_cv = arrayfun(@(x) x.ModData.unit_data.dprime_stability_cv,all_SU_data)';

lem_fov_SUs = intersect(fov_SUs,lem_SUs);
jbe_SUs = find(strcmp('jbe',{all_SU_data(:).animal}));
lem_SUs = find(strcmp('lem',{all_SU_data(:).animal}));

lem_fov_expt_nums = [266 275];
lem_parafov_expt_nums = [270 277 281 287 289 294 229 297];
fov_SUs = find([all_SU_data(:).expt_num] < 200 | ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
parafov_SUs = find([all_SU_data(:).expt_num] > 200 & ~ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
lem_parafov_SUs = intersect(parafov_SUs,lem_SUs);

N_gsacs = [all_SU_data(:).N_gsacs];
N_msacs = [all_SU_data(:).N_msacs];
N_gsacs_gray = [all_SU_tdata(:).N_gsacs_gray];
N_msacs_gray = [all_SU_tdata(:).N_msacs_gray];

%% GRAY-BACKGROUND SACCADE TRIGGERED AVERAGES
cur_SUs = find(avg_rates >= min_rate & N_gsacs_gray >= min_TA_Nsacs);
all_gsac_gray = reshape([all_SU_tdata(:).gsac_gray_avg],[],tot_Nunits)';

xl = [-0.15 0.4];

close all
f1 = figure(); hold on
curset = intersect(jbe_SUs,cur_SUs);
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','r'});
curset = intersect(lem_SUs,cur_SUs);
h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curset,:)),nanstd(all_gsac_gray(curset,:))/sqrt(length(curset)),{'color','b'});

curset = find(avg_rates >= min_rate & N_gsacs >= min_Nsacs);
plot(tlags,nanmean(all_gsac_gray(curset,:)),'k','linewidth',2);

xlim(xl);
ylim([0.7 1.3]);
legend([h1.mainLine h2.mainLine],{'JBE','LEM'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
title('Gsac TA Grayback');

%% GSAC GRAYBACK STRENGTH VS TIMING SCATTERPLOT
cur_SUs = find(avg_rates >= min_rate & N_gsacs_gray >= min_TA_Nsacs);

all_gsac_gray = reshape([all_SU_tdata(:).gsac_gray_avg],[],tot_Nunits)';
%compute z-score of trig-avgs using bootstrap SE estimates
all_gsac_gray_SE = reshape([all_SU_tdata(:).gsac_gray_std],[],tot_Nunits)';
all_gsac_gray_Z = (all_gsac_gray - 1)./all_gsac_gray_SE;

search_range = [0 0.35];

[gsac_Efact,gsac_exctime] = get_tavg_peaks((all_gsac_gray-1),tlags,search_range);
[gsac_Sfact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_gray-1),tlags,search_range);
[gsac_exc_Z] = get_tavg_peaks(all_gsac_gray_Z,tlags,search_range);
[gsac_inh_Z] = get_tavg_peaks(-all_gsac_gray_Z,tlags,search_range);

%USE CIS HERE RATHER THAN ZCORES
% Zthresh = 3; %minimum z-score value for peak to be counted as 'significant'
% sig_E = cur_SUs(gsac_exc_Z(cur_SUs) >=  Zthresh);
% sig_S = cur_SUs(gsac_inh_Z(cur_SUs) >=  Zthresh);
% sig_both = intersect(sig_E,sig_S);

xl = [0 1];
yl = [0 0.3];
mS = 3; %marker size

close all
f1 = figure(); 
hold on
plot(gsac_Efact(sig_E),gsac_exctime(sig_E),'bo','markersize',mS);
plot(gsac_Sfact(sig_S),gsac_inhtime(sig_S),'ro','markersize',mS);
set(gca,'xtick',0:0.2:1);

xlabel('Modulation strength');
ylabel('Modulation timing (s)');
xlim(xl); ylim(yl);

f2 = figure();
plot(gsac_inhtime(sig_both),gsac_exctime(sig_both),'bo','markersize',mS);
line(yl,yl,'color','k');
xlim(yl); ylim(yl);

f3 = figure();
plot(gsac_Sfact(sig_both),gsac_Efact(sig_both),'bo','markersize',mS);
line(xl,xl,'color','k');
xlim(xl); ylim(xl);

nbins = 20;
efact_binedges = linspace(xl(1),xl(2),nbins+1);
time_binedges = linspace(yl(1),yl(2),nbins+1);
gsac_Efact_dist = histc(gsac_Efact(sig_E),efact_binedges);
gsac_Efact_dist = gsac_Efact_dist/sum(gsac_Efact_dist);
gsac_Etime_dist = histc(gsac_exctime(sig_E),time_binedges);
gsac_Etime_dist = gsac_Etime_dist/sum(gsac_Etime_dist);
gsac_Sfact_dist = histc(gsac_Sfact(sig_S),efact_binedges);
gsac_Sfact_dist = gsac_Sfact_dist/sum(gsac_Sfact_dist);
gsac_Stime_dist = histc(gsac_inhtime(sig_S),time_binedges);
gsac_Stime_dist = gsac_Stime_dist/sum(gsac_Stime_dist);

f4 = figure();
hold on
stairs(efact_binedges,gsac_Efact_dist);
stairs(efact_binedges,gsac_Sfact_dist,'r');
xlim(xl);
f5 = figure();
hold on
stairs(time_binedges,gsac_Etime_dist);
stairs(time_binedges,gsac_Stime_dist,'r');
xlim(yl);

set1 = sig_both(gsac_exctime(sig_both) > gsac_inhtime(sig_both));
set2 = sig_both(gsac_exctime(sig_both) < gsac_inhtime(sig_both));
f3 = figure();
hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(set1,:)),nanstd(all_gsac_gray(set1,:))/sqrt(length(set1)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(set2,:)),nanstd(all_gsac_gray(set2,:))/sqrt(length(set2)),{'color','r'});