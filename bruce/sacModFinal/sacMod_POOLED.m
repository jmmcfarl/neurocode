%%
% close all
clear all
clc

fit_unCor = 0;

fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';
base_tname = 'sac_trig_avg_data';
base_sname = 'sacStimProcFin';
% base_sname = 'sacStimProcTest';
% base_sname = 'sacStimProc2';
% base_ename = 'sacStimProc';
base_timename = 'sac_info_timing';

all_SU_data = [];
all_SU_NPdata = [];
all_SU_tdata = [];
all_SU_timedata = [];

%% LOAD JBE
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
n_probes = 96;
% ori_list = [0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan];
ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
rmfield_list = {};

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    osac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
        sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
%     sac_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod/'];
    tavg_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            sname = strcat(tavg_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            temp = load(sname);
            tname = strcat(sac_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)),sprintf('_ori%d',ori_list(ee,ii)));
            if fit_unCor
                tname = strcat(tname,'_unCor');
            end
            load(tname);
            timename = strcat(osac_dir,base_timename,sprintf('_ori%d',ori_list(ee,ii)));
            load(timename);
            
            
%             if isfield(sacStimProc,'msac_ovavg_rate')
%                 sacStimProc = rmfield(sacStimProc,'msac_ovavg_rate');
%             end

%             tname = strcat(osac_dir,base_ename,sprintf('_ori%d',ori_list(ee,ii)));
%             temp2 = load(tname);
%             for bbb = 1:length(sacStimProc)
%                 sacStimProc(bbb).gsacGainMod = temp2.sacStimProc(bbb).gsacGainMod;
%             end
            
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
%             base_xvLLimps(lia) = SM_SU_xvLLimp(locb(lia));
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
            ori_time_data{ii}(lia) = sacInfoTiming(ucells);
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
            all_SU_timedata = cat(1,all_SU_timedata,ori_time_data{mlocs(ii)}(ii));
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
    osac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
        sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
%     sac_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod/'];
    tavg_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            sname = strcat(tavg_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            temp = load(sname);
            tname = strcat(sac_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)),sprintf('_ori%d',ori_list(ee,ii)));
             if fit_unCor
                tname = strcat(tname,'_unCor');
            end
           load(tname);
            timename = strcat(osac_dir,base_timename,sprintf('_ori%d',ori_list(ee,ii)));
            load(timename);
            
            
%             if isfield(sacStimProc,'msac_ovavg_rate')
%             sacStimProc = rmfield(sacStimProc,'msac_ovavg_rate');
%             end
%             
%                         tname = strcat(osac_dir,base_ename,sprintf('_ori%d',ori_list(ee,ii)));
%             temp2 = load(tname);
%             for bbb = 1:length(sacStimProc)
%                 sacStimProc(bbb).gsacGainMod = temp2.sacStimProc(bbb).gsacGainMod;
%             end

            tavg_SU_numbers = [temp.sua_data(:).SU_numbers];
            
            ucells = arrayfun(@(x) length(x.gsac_avg_rate),sacStimProc) > 0;
            sua_data = sacStimProc(ucells);
            SM_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,sua_data);
            SM_SU_xvLLimp = arrayfun(@(x) x.ModData.rectGQM.xvLLimp,sua_data);
            SM_SU_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,sua_data);
            
            [lia,locb] = ismember(tavg_SU_numbers,SM_SU_numbers);
            lia_inds = find(lia);
            base_xvLLimps = ones(size(tavg_SU_numbers))*-Inf;
%             base_xvLLimps(lia) = SM_SU_xvLLimp(locb(lia));
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
            ori_time_data{ii}(lia) = sacInfoTiming(ucells);
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
            all_SU_timedata = cat(1,all_SU_timedata,ori_time_data{mlocs(ii)}(ii));
        end
    end
    
    clear ori_SU_nums ori_xvLLimps ori_sua_data ori_tavg_data 
end

%% time axis for trig-avg data
tlags = temp.trig_avg_params.lags;
Tdt = temp.trig_avg_params.dt;
tlags = tlags*Tdt;

%% SELECT USABLE CELLS
avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_SU_data);
tot_spikes = arrayfun(@(x) x.ModData.unit_data.tot_spikes,all_SU_data);
rec_dur = arrayfun(@(x) x.ModData.unit_data.N_used_samps,all_SU_data)*dt/60;
expt_nums = [all_SU_data(:).expt_num];
expt_oris = [all_SU_data(:).bar_ori];
xvLLimps = [all_SU_data(:).xvLLimp];
mod_LLimps = [all_SU_data(:).gsac_spost_ov_modinfo];

clust_iso_dist = arrayfun(@(x) x.ModData.unit_data.SU_isodist,all_SU_data);
clust_Lratio = arrayfun(@(x) x.ModData.unit_data.SU_Lratio,all_SU_data);
clust_refract = arrayfun(@(x) x.ModData.unit_data.SU_refract,all_SU_data);
clust_dprime = arrayfun(@(x) x.ModData.unit_data.SU_dprime,all_SU_data);
rate_stability_cv = arrayfun(@(x) x.ModData.unit_data.rate_stability_cv,all_SU_data);
dprime_stability_cv = arrayfun(@(x) x.ModData.unit_data.dprime_stability_cv,all_SU_data);

min_rate = 5; %in Hz (5)
min_Nsacs = 1e3; % (1000)
min_xvLLimp = 0.05; %(0.05);
jbe_SUs = find(strcmp('jbe',{all_SU_data(:).animal}));
lem_SUs = find(strcmp('lem',{all_SU_data(:).animal}));

lem_fov_expt_nums = [266 275];
lem_parafov_expt_nums = [270 277 281 287 289 294 229 297];
fov_SUs = find([all_SU_data(:).expt_num] < 200 | ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
parafov_SUs = find([all_SU_data(:).expt_num] > 200 & ~ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
lem_fov_SUs = intersect(fov_SUs,lem_SUs);
lem_parafov_SUs = intersect(parafov_SUs,lem_SUs);

N_gsacs = [all_SU_data(:).N_gsacs];
N_msacs = [all_SU_data(:).N_msacs];

% use_gsac_SUs = find(avg_rates' >= min_rate & N_gsacs >= min_Nsacs & xvLLimps > min_xvLLimp);
use_gsac_SUs = find(avg_rates' >= min_rate & N_gsacs >= min_Nsacs & mod_LLimps > min_xvLLimp);
use_msac_SUs = find(avg_rates' >= min_rate & N_msacs >= min_Nsacs & mod_LLimps > min_xvLLimp);
use_jbe_SUs = intersect(use_gsac_SUs,jbe_SUs);
use_lem_SUs = intersect(use_gsac_SUs,lem_SUs);

all_gsac_rates = reshape([all_SU_data(:).gsac_avg_rate]/dt,[],length(all_SU_data))';
all_gsac_nrates = bsxfun(@rdivide,all_gsac_rates,avg_rates);

all_gsac_Trates = reshape([all_SU_tdata(:).gsac_avg],[],length(all_SU_tdata))';
all_msac_Trates = reshape([all_SU_tdata(:).msac_avg],[],length(all_SU_tdata))';

%% COMPARE GAIN/OFFSET FILTERS AND AVG RATES
close all
all_postgains = [];
all_postoffs = [];
for ii = 1:length(all_SU_data)
    all_postgains = cat(1,all_postgains,all_SU_data(ii).gsac_post_singmod.mods(3).filtK');
    all_postoffs = cat(1,all_postoffs,all_SU_data(ii).gsac_post_singmod.mods(2).filtK');
end

all_gsac_spost_offset = reshape([all_SU_data(:).gsac_spost_offset],[],length(all_SU_data))';
all_gsac_spost_offset = bsxfun(@rdivide,all_gsac_spost_offset,[all_SU_data(:).gsac_ovavg_rate]');
all_gsac_spost_gain = reshape([all_SU_data(:).gsac_spost_gain],[],length(all_SU_data))';


all_gsac_TB_offset = reshape([all_SU_data(:).gsac_TB_offset],[],length(all_SU_data))';
all_gsac_TB_offset = bsxfun(@rdivide,all_gsac_TB_offset,[all_SU_data(:).gsac_ovavg_rate]');
all_gsac_TB_gain = reshape([all_SU_data(:).gsac_TB_gain],[],length(all_SU_data))';
TB_slags = all_SU_data(1).gsac_TB_lagX;

xl = [-0.1 0.3];

yl = [0.7 1.3];
f1 = figure();
hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_nrates(use_gsac_SUs,:)),std(all_gsac_nrates(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','k'});
h2=shadedErrorBar(tlags,mean(all_gsac_Trates(use_gsac_SUs,:)),std(all_gsac_Trates(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','k'});
xlabel('Time (s)');
ylabel('Relative gain');
xlim(xl);
line(xl,[1 1],'color','k');
ylim(yl);
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

f2 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_postgains(use_gsac_SUs,:)),std(all_postgains(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_postoffs(use_gsac_SUs,:)),std(all_postoffs(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Relative gain');
line(xl,[0 0],'color','k');
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Gain');

f3 = figure();
hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_spost_offset(use_gsac_SUs,:)),std(all_gsac_spost_offset(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
h1=shadedErrorBar(TB_slags*dt,mean(all_gsac_TB_offset(use_gsac_SUs,:)),std(all_gsac_TB_offset(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Offset');
line(xl,[0 0],'color','k');
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Gain');

f4 = figure();
hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_spost_gain(use_gsac_SUs,:)),std(all_gsac_spost_gain(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
h1=shadedErrorBar(TB_slags*dt,mean(all_gsac_TB_gain(use_gsac_SUs,:)),std(all_gsac_TB_gain(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Relative gain');
line(xl,[1 1],'color','k');
xlim(xl);
yl = ylim();
line([1 1],yl,'color','k');
xlabel('Time (s)');
ylabel('Gain');

% f3 = figure();
% hold on
% h2=shadedErrorBar(slags*dt,mean(all_postoffs(use_gsac_SUs,:)),std(all_postoffs(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
% xlabel('Time (s)');
% ylabel('Offset');
% line(xl,[0 0],'color','k');
% xlim(xl);
% yl = ylim();
% line([0 0],yl,'color','k');
% xlabel('Time (s)');
% ylabel('Offset');

%
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_avg_rates.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'Gsac_postgains.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

% figufy(f3);
% fname = [fig_dir 'Gsac_respoffs.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

% figufy(f3);
% fname = [fig_dir 'Gsac_respoffs.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);

% figufy(f4);
% fname = [fig_dir 'Gsac_respgains.pdf'];
% exportfig(f4,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f4);

%% COMPARE GAINS AND RATES FOR STRONGLY ENHANCED VS STRONGLY SUPPRESSED UNITS
close all
search_range = [0 0.35];
[gsac_Sfact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_Trates-1),tlags,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_Trates-1,tlags,search_range);

EIrat = gsac_Efact./gsac_Sfact;
EIprc = prctile(EIrat(use_gsac_SUs),[25 50 75]);

% stronger_E = use_gsac_SUs(gsac_Efact(use_gsac_SUs) > gsac_Sfact(use_gsac_SUs));
% stronger_I = use_gsac_SUs(gsac_Sfact(use_gsac_SUs) > gsac_Efact(use_gsac_SUs));
stronger_E = use_gsac_SUs(EIrat(use_gsac_SUs) >= EIprc(3));
stronger_I = use_gsac_SUs(EIrat(use_gsac_SUs) <= EIprc(1));
% stronger_E = use_gsac_SUs(EIrat(use_gsac_SUs) >= EIprc(2));
% stronger_I = use_gsac_SUs(EIrat(use_gsac_SUs) < EIprc(2));



tot_mod = max([gsac_Sfact gsac_Efact],[],2);
mod_thresh = prctile(tot_mod,25);
weak_set = use_gsac_SUs(tot_mod(use_gsac_SUs) <= mod_thresh);

all_gsac_smodinfo = reshape([all_SU_data(:).gsac_spost_modinfo],[],length(all_SU_data))';
all_gsac_ov_smodinfo = [all_SU_data(:).gsac_spost_ov_modinfo];
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');
all_gsac_smodinforate = all_gsac_smodinfo.*all_gsac_rates*dt;
all_gsac_Nsmodinforate = bsxfun(@rdivide,all_gsac_smodinforate,all_gsac_ov_smodinfo'.*avg_rates*dt);

search_range = [0 0.2];
[gsac_info_Sfact,gsac_info_inhtime] = get_tavg_peaks(-(all_gsac_Nsmodinfo-1),slags*dt,search_range);

xl = [-0.1 0.3];

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,1+mean(all_postgains(stronger_E,:)),std(all_postgains(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_postgains(stronger_I,:)),std(all_postgains(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Gain');

f2 = figure();
hold on
h1=shadedErrorBar(tlags,mean(all_gsac_Trates(stronger_E,:)),std(all_gsac_Trates(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(tlags,mean(all_gsac_Trates(stronger_I,:)),std(all_gsac_Trates(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

f3 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(stronger_E,:)),std(all_gsac_Nsmodinfo(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(stronger_I,:)),std(all_gsac_Nsmodinfo(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative info');

f4 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_postoffs(stronger_E,:)),std(all_postoffs(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_postoffs(stronger_I,:)),std(all_postoffs(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[0 0],'color','k');
xlabel('Time (s)');
ylabel('Offset');

f5 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_spost_offset(stronger_E,:)),std(all_gsac_spost_offset(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_spost_offset(stronger_I,:)),std(all_gsac_spost_offset(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[0 0],'color','k');
xlabel('Time (s)');
ylabel('Offset');

f6 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_spost_gain(stronger_E,:)),std(all_gsac_spost_gain(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_spost_gain(stronger_I,:)),std(all_gsac_spost_gain(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[0 0],'color','k');
xlabel('Time (s)');
ylabel('Offset');

% f5 = figure();
% hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinforate(stronger_E,:)),std(all_gsac_Nsmodinforate(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinforate(stronger_I,:)),std(all_gsac_Nsmodinforate(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
% xlim(xl);
% line(xl,[1 1],'color','k');
% xlabel('Time (s)');
% ylabel('Relative info');

fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_ESdep_gains.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


% figufy(f2);
% fname = [fig_dir 'Gsac_ESdep_rates.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'Gsac_ESdep_info.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);


%% COMPARE SS INFO AND INFO RATES (POST MODEL BASED)

all_gsac_smodinfo = reshape([all_SU_data(:).gsac_spost_modinfo],[],length(all_SU_data))';
all_gsac_ov_smodinfo = [all_SU_data(:).gsac_spost_ov_modinfo];
all_gsac_smodinforate = all_gsac_smodinfo.*all_gsac_rates*dt;
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');
all_gsac_Nsmodinforate = bsxfun(@rdivide,all_gsac_smodinforate,all_gsac_ov_smodinfo'.*avg_rates*dt);

all_gsac_TBinfo = reshape([all_SU_data(:).gsac_TB_info],[],length(all_SU_data))';
all_gsac_TBrates = reshape([all_SU_data(:).gsac_TB_avg_rate],[],length(all_SU_data))';
all_gsac_TBinforate = all_gsac_TBinfo.*all_gsac_TBrates;
all_gsac_ov_TB_info = [all_SU_data(:).gsac_ov_TB_info];
all_gsac_NTBinfo = bsxfun(@rdivide,all_gsac_TBinfo,all_gsac_ov_TB_info');
all_gsac_NTBinforate = bsxfun(@rdivide,all_gsac_TBinforate,all_gsac_ov_TB_info'.*avg_rates*dt);
TB_slags = all_SU_data(1).gsac_TB_lagX;

% search_range = [0 0.2];
% TB_info_Sfact = get_tavg_peaks(-(all_gsac_NTBinfo-1),TB_slags*dt,search_range);
% TB_gain_Sfact = get_tavg_peaks(-(all_gsac_TB_gain-1),TB_slags*dt,search_range);
% search_range = [0 0.3];
% TB_offset_Efact = get_tavg_peaks(all_gsac_TB_offset,TB_slags*dt,search_range);

xl = [-0.1 0.3];

f1 = figure();
hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(use_gsac_SUs,:)),std(all_gsac_Nsmodinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinforate(use_gsac_SUs,:)),std(all_gsac_Nsmodinforate(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
h1=shadedErrorBar(TB_slags*dt,mean(all_gsac_NTBinfo(use_gsac_SUs,:)),std(all_gsac_NTBinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
h2=shadedErrorBar(TB_slags*dt,mean(all_gsac_NTBinforate(use_gsac_SUs,:)),std(all_gsac_NTBinforate(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','k'});
% legend([h1.mainLine h2.mainLine],{'TB','Mod-pred'});
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Relative information');

% f2 = figure();
% hold on
% % h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(use_gsac_SUs,:)),std(all_gsac_Nsmodinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
% h1=shadedErrorBar(TB_slags*dt,mean(all_gsac_TBinfo(use_gsac_SUs,:)),std(all_gsac_TBinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
% % legend([h1.mainLine h2.mainLine],{'TB','Mod-pred'});
% xlabel('Time (s)');
% ylabel('Relative info');
% xlim(xl);
% yl = ylim();
% line([0 0],yl,'color','k');
% xlabel('Time (s)');
% ylabel('Relative information');

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_SSinfo_inforate_uncor.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% COMPARE POST G/O MODEL AND TB MODEL

all_gsac_TBinfo = reshape([all_SU_data(:).gsac_TB_info],[],length(all_SU_data))';
all_gsac_TBrates = reshape([all_SU_data(:).gsac_TB_avg_rate],[],length(all_SU_data))';
all_gsac_TBinforate = all_gsac_TBinfo.*all_gsac_TBrates;
all_gsac_ov_TB_info = [all_SU_data(:).gsac_ov_TB_info];
all_gsac_NTBinfo = bsxfun(@rdivide,all_gsac_TBinfo,all_gsac_ov_TB_info');
all_gsac_NTBinforate = bsxfun(@rdivide,all_gsac_TBinforate,all_gsac_ov_TB_info'.*avg_rates*dt);
TB_slags = all_SU_data(1).gsac_TB_lagX;

all_gsac_smodinfo = reshape([all_SU_data(:).gsac_spost_modinfo],[],length(all_SU_data))';
all_gsac_ov_smodinfo = [all_SU_data(:).gsac_spost_ov_modinfo];
all_gsac_smodinforate = all_gsac_smodinfo.*all_gsac_rates;
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');
all_gsac_Nsmodinforate = bsxfun(@rdivide,all_gsac_smodinforate,all_gsac_ov_smodinfo'.*avg_rates*dt);

xl = [-0.1 0.3];

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(use_gsac_SUs,:)),std(all_gsac_Nsmodinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
h2=shadedErrorBar(TB_slags*dt,mean(all_gsac_NTBinfo(use_gsac_SUs,:)),std(all_gsac_NTBinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
% legend([h1.mainLine h2.mainLine],{'TB','Mod-pred'});
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative information');

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_smod_TB_infocompare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% PLOT SPIKE-WEIGHTED GEN SIGNALS


all_gsac_spk_gmean = reshape([all_SU_data(:).gsac_spkCondI],[],length(all_SU_data))';
ov_gmean = [all_SU_data(:).gsac_ovspkCondI];

% all_gsac_spk_gmean = bsxfun(@rdivide,all_gsac_spk_gmean,ov_gmean');
all_gsac_spk_gmean = bsxfun(@minus,all_gsac_spk_gmean,ov_gmean');

TB_slags = all_SU_data(1).gsac_TB_lagX;

xl = [-0.1 0.3];

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_spk_gmean(use_gsac_SUs,:)),std(all_gsac_spk_gmean(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
xlabel('Time (s)');
ylabel('Relative gain change');
line(xl,[0 0],'color','k');
xlim(xl);
xlabel('Time (s)');

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_spk_weight_g.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


%% COMPARE POST G/O MODEL AND SUBSPACE MODEL

all_gsac_submodinfo = reshape([all_SU_data(:).gsac_sub_LLinfo],[],length(all_SU_data))';
all_gsac_ov_submodinfo = [all_SU_data(:).gsac_sub_ov_LLinfo];
% all_gsac_ov_submodinfo = mean(all_gsac_submodinfo,2)';
all_gsac_submodinforate = all_gsac_submodinfo.*all_gsac_rates;
all_gsac_Nsubmodinfo = bsxfun(@rdivide,all_gsac_submodinfo,all_gsac_ov_submodinfo');
all_gsac_Nsubmodinforate = bsxfun(@rdivide,all_gsac_submodinforate,all_gsac_ov_submodinfo'.*avg_rates*dt);

all_gsac_smodinfo = reshape([all_SU_data(:).gsac_spost_LLinfo],[],length(all_SU_data))';
all_gsac_ov_smodinfo = [all_SU_data(:).gsac_spost_ov_LLinfo];
all_gsac_smodinforate = all_gsac_smodinfo.*all_gsac_rates;
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');
all_gsac_Nsmodinforate = bsxfun(@rdivide,all_gsac_smodinforate,all_gsac_ov_smodinfo'.*avg_rates*dt);


all_gsac_submodinfo2 = reshape([all_SU_data(:).gsac_sub_modinfo],[],length(all_SU_data))';
all_gsac_ov_submodinfo2 = [all_SU_data(:).gsac_sub_ov_modinfo];
% all_gsac_submodinforate = all_gsac_submodinfo.*all_gsac_rates;
% all_gsac_Nsubmodinfo = bsxfun(@rdivide,all_gsac_submodinfo,all_gsac_ov_submodinfo');
% all_gsac_Nsubmodinforate = bsxfun(@rdivide,all_gsac_submodinforate,all_gsac_ov_submodinfo'.*avg_rates*dt);
% 
all_gsac_smodinfo2 = reshape([all_SU_data(:).gsac_spost_modinfo],[],length(all_SU_data))';
all_gsac_ov_smodinfo2 = [all_SU_data(:).gsac_spost_ov_modinfo];
% all_gsac_smodinforate = all_gsac_smodinfo.*all_gsac_rates;
% all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');
% all_gsac_Nsmodinforate = bsxfun(@rdivide,all_gsac_smodinforate,all_gsac_ov_smodinfo'.*avg_rates*dt);
% 

xl = [-0.1 0.3];

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(use_gsac_SUs,:)),std(all_gsac_Nsmodinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_Nsubmodinfo(use_gsac_SUs,:)),std(all_gsac_Nsubmodinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
% legend([h1.mainLine h2.mainLine],{'TB','Mod-pred'});
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);
xlabel('Time (s)');
ylabel('Relative information');

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_smod_submod_infocompare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%% COMPARE COMPONENTS OF PRE/POST MODEL
% all_gsac_gaink = reshape(cell2mat(arrayfun(@(x) x.gsacGainMod.gain_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';
% all_gsac_stimk = reshape(cell2mat(arrayfun(@(x) x.gsacGainMod.stim_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';
% all_gsac_offk = reshape(cell2mat(arrayfun(@(x) x.gsacGainMod.off_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';
% all_gsac_gainoff = arrayfun(@(x) x.gsacGainMod.gain_offset,all_SU_data)-1;

all_gsac_stimk = reshape(cell2mat(arrayfun(@(x) x.gsacPreGainMod.stim_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';
all_gsac_offk = reshape(cell2mat(arrayfun(@(x) x.gsacPreGainMod.off_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';
all_msac_stimk = reshape(cell2mat(arrayfun(@(x) x.msacPreGainMod.stim_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';
all_msac_offk = reshape(cell2mat(arrayfun(@(x) x.msacPreGainMod.off_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';

% all_gsac_gaink = bsxfun(@minus,all_gsac_gaink,all_gsac_gainoff);

all_gsac_postgains = [];
all_gsac_postoffs = [];
all_msac_postgains = [];
all_msac_postoffs = [];
for ii = 1:length(all_SU_data)
    all_gsac_postgains = cat(1,all_gsac_postgains,all_SU_data(ii).gsac_post_singmod.mods(3).filtK');
    all_gsac_postoffs = cat(1,all_gsac_postoffs,all_SU_data(ii).gsac_post_singmod.mods(2).filtK');
    all_msac_postgains = cat(1,all_msac_postgains,all_SU_data(ii).msac_post_singmod.mods(3).filtK');
    all_msac_postoffs = cat(1,all_msac_postoffs,all_SU_data(ii).msac_post_singmod.mods(2).filtK');
end

xl = [-0.1 0.3];
yl = [0.6 1.2];

uu = use_gsac_SUs;
% uu = use_jbe_SUs;
% uu = use_lem_SUs;
f1 = figure();
hold on
% h1=shadedErrorBar(slags*dt,1+mean(all_gsac_gaink(uu,:)),std(all_gsac_gaink(uu,:))/sqrt(length(uu)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_gsac_stimk(uu,:)),std(all_gsac_stimk(uu,:))/sqrt(length(uu)),{'color','k'});
h2=shadedErrorBar(slags*dt,1+mean(all_msac_stimk(uu,:)),std(all_msac_stimk(uu,:))/sqrt(length(uu)),{'color','b'});
h3=plot(slags*dt,1+mean(all_gsac_postgains(uu,:)),'m','linewidth',2);
h3=plot(slags*dt,1+mean(all_msac_postgains(uu,:)),'r','linewidth',2);
% h3=plot(slags*dt,1+mean(all_gsac_postoffs(uuuub','linewidth',2);
% h4=plot(slags*dt,1+mean(all_gsac_offk(use_gsac_SUs,:)),'r','linewidth',2);
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);
ylim(yl);
line([0 0],yl,'color','k');

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_PREPOST_modcompare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
%% COMPARE DIFFERENT INFO CALCS

all_gsac_TBinfo = reshape([all_SU_data(:).gsac_TB_info],[],length(all_SU_data))';
all_gsac_TBrates = reshape([all_SU_data(:).gsac_TB_avg_rate],[],length(all_SU_data))';
all_gsac_TBinforate = all_gsac_TBinfo.*all_gsac_TBrates;
all_gsac_ov_TB_info = [all_SU_data(:).gsac_ov_TB_info];
all_gsac_NTBinfo = bsxfun(@rdivide,all_gsac_TBinfo,all_gsac_ov_TB_info');
all_gsac_NTBinforate = bsxfun(@rdivide,all_gsac_TBinforate,all_gsac_ov_TB_info'.*avg_rates*dt);

all_gsac_smodinfo = reshape([all_SU_data(:).gsac_spost_modinfo],[],length(all_SU_data))';
all_gsac_ov_smodinfo = [all_SU_data(:).gsac_spost_ov_modinfo];
all_gsac_smodinforate = all_gsac_smodinfo.*all_gsac_rates;
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');
all_gsac_Nsmodinforate = bsxfun(@rdivide,all_gsac_smodinforate,all_gsac_ov_smodinfo'.*avg_rates*dt);

all_gsac_modinfo = reshape([all_SU_data(:).gsac_modinfo],[],length(all_SU_data))';
all_gsac_ov_modinfo = [all_SU_data(:).gsac_ov_modinfo];
all_gsac_modinforate = all_gsac_modinfo.*all_gsac_rates;
all_gsac_Nmodinfo = bsxfun(@rdivide,all_gsac_modinfo,all_gsac_ov_modinfo');
all_gsac_Nmodinforate = bsxfun(@rdivide,all_gsac_modinforate,all_gsac_ov_modinfo'.*avg_rates*dt);

xl = [-0.1 0.3];

TB_slags = all_SU_data(1).gsac_TB_lagX;

f1 = figure();
hold on
h1=shadedErrorBar(TB_slags*dt,mean(all_gsac_NTBinfo(use_gsac_SUs,:)),std(all_gsac_NTBinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_Nmodinfo(use_gsac_SUs,:)),std(all_gsac_Nmodinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
h3=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(use_gsac_SUs,:)),std(all_gsac_Nsmodinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','k'});
legend([h1.mainLine h2.mainLine h3.mainLine],{'TB','Mod-pred','Submod-pred'});
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);

%% COMPARE MOD-based and LL-based INFO Calcs

all_gsac_smodinfo = reshape([all_SU_data(:).gsac_spost_modinfo],[],length(all_SU_data))';
all_gsac_ov_smodinfo = [all_SU_data(:).gsac_spost_ov_modinfo];
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');

all_gsac_smodLL = reshape([all_SU_data(:).gsac_spost_LLinfo],[],length(all_SU_data))';
all_gsac_ov_smodLL = [all_SU_data(:).gsac_spost_ov_LLinfo];
all_gsac_NsmodLL = bsxfun(@rdivide,all_gsac_smodLL,all_gsac_ov_smodLL');

all_gsac_subLL = reshape([all_SU_data(:).gsac_sub_LLinfo],[],length(all_SU_data))';
all_gsac_ov_subLL = [all_SU_data(:).gsac_sub_ov_LLinfo];
all_gsac_NsubLL = bsxfun(@rdivide,all_gsac_subLL,all_gsac_ov_subLL');

all_gsac_TBLL = reshape([all_SU_data(:).gsac_TB_LLinfo],[],length(all_SU_data))';
all_gsac_ov_TBLL = [all_SU_data(:).gsac_TB_ov_LLinfo];
all_gsac_NTBLL = bsxfun(@rdivide,all_gsac_TBLL,all_gsac_ov_TBLL');

xl = [-0.1 0.3];

TB_slags = all_SU_data(1).gsac_TB_lagX;

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(use_gsac_SUs,:)),std(all_gsac_Nsmodinfo(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_NsmodLL(use_gsac_SUs,:)),std(all_gsac_NsmodLL(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
h3=shadedErrorBar(slags*dt,mean(all_gsac_NsubLL(use_gsac_SUs,:)),std(all_gsac_NsubLL(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','k'});
h4=shadedErrorBar(slags*dt,mean(all_gsac_NTBLL(use_gsac_SUs,:)),std(all_gsac_NTBLL(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','g'});
% legend([h1.mainLine h2.mainLine h3.mainLine],{'TB','Mod-pred','Submod-pred'});
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);


%% PARSE TIMING OF EXC AND INH TUNED COMPONENTS
flen = 15;
all_Ekerns = nan(length(all_SU_data),flen);
all_Ikerns = nan(length(all_SU_data),flen);
[n_Efilts,n_Ifilts] = deal(nan(length(all_SU_data),1));
for ii = 1:length(all_SU_data)
    cur_mod = all_SU_data(ii).ModData.rectGQM;
    mod_signs = [cur_mod.mods(:).sign];
    sd = cur_mod.stim_params(1).stim_dims;
    mod_filts = reshape([cur_mod.mods(:).filtK],[sd(1) sd(2) length(mod_signs)]);
    
    qfilts = find(strcmp({cur_mod.mods(:).NLtype},'quad'));
    
    tkerns = squeeze(std(mod_filts,[],2));
    tkerns = bsxfun(@rdivide,tkerns,sqrt(sum(tkerns.^2)));
    tkerns = bsxfun(@times,tkerns,all_SU_data(ii).ModData.rectGQM.rel_filt_weights);
    tkerns(:,all_SU_data(ii).ModData.rectGQM.rel_filt_weights==0) = 0;
    %     tkerns(:,qfilts) = tkerns(:,qfilts).^2;
    
    
    if sum(mod_signs == 1) > 0
        all_Ekerns(ii,:) = mean(tkerns(:,mod_signs == 1),2);
    end
    if sum(mod_signs == -1) > 0
        all_Ikerns(ii,:) = mean(tkerns(:,mod_signs == -1),2);
    end
    %      cur_mod = all_SU_data(ii).ModData.bestGQM;
    %     mod_signs = [cur_mod.mods(:).sign];
    n_Efilts(ii) = sum(mod_signs == 1);
    n_Ifilts(ii) = sum(mod_signs == -1);
end
net_Estrength = mean(all_Ekerns,2);
net_Istrength = mean(all_Ikerns,2);

all_EI_xc = reshape([all_SU_data(:).EI_xc],[],length(all_SU_data))';
EI_lags = all_SU_data(1).EI_xc_lags;

all_nEkerns = bsxfun(@rdivide,all_Ekerns,sqrt(sum(all_Ekerns.^2,2)));
all_nIkerns = bsxfun(@rdivide,all_Ikerns,sqrt(sum(all_Ikerns.^2,2)));

all_postgains = [];
all_postoffs = [];
for ii = 1:length(all_SU_data)
    all_postgains = cat(1,all_postgains,all_SU_data(ii).gsac_post_singmod.mods(3).filtK');
    all_postoffs = cat(1,all_postoffs,all_SU_data(ii).gsac_post_singmod.mods(2).filtK');
end

tax = (0:(flen-1))*dt + dt/2;

min_NIfilts = 1;

all_gsac_Egain = reshape([all_SU_data(:).gsac_post_Egains],[],length(all_SU_data))';
all_gsac_Igain = reshape([all_SU_data(:).gsac_post_Igains],[],length(all_SU_data))';
use_gsac_Ikerns = use_gsac_SUs(n_Ifilts(use_gsac_SUs) >= min_NIfilts);

xr = [0 0.2];
poss_lagrange = find(slags*dt > xr(1) & slags*dt < xr(2));
[Egain_minloc,Igain_minloc] = deal(nan(size(all_gsac_Egain,1),1));
for ii = 1:size(all_gsac_Egain,1)
    [temp,temploc] = findpeaks(-all_gsac_Egain(ii,poss_lagrange),'sortstr','descend');
    if ~isempty(temploc)
        Egain_minloc(ii) = temploc(1);
    end
    [temp,temploc] = findpeaks(-all_gsac_Igain(ii,poss_lagrange),'sortstr','descend');
    if ~isempty(temploc)
        Igain_minloc(ii) = temploc(1);
    end
end

all_EI_ratio = (1+all_gsac_Egain)./(1+all_gsac_Igain);
% all_EI_ratio = (1+all_gsac_Egain) - (1+all_gsac_Igain)+1;
all_EI_mean = 0.5*(1+all_gsac_Egain)+0.5*(1+all_gsac_Igain);
all_EI_diff = (all_gsac_Egain) - (all_gsac_Igain);
all_EI_ratio = all_EI_diff;

% all_gsac_offk = reshape(cell2mat(arrayfun(@(x) x.gsacGainMod.off_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';

all_gsac_smodinfo = reshape([all_SU_data(:).gsac_spost_modinfo],[],length(all_SU_data))';
all_gsac_ov_smodinfo = [all_SU_data(:).gsac_spost_ov_modinfo];
all_gsac_smodinforate = all_gsac_smodinfo.*all_gsac_rates*dt;
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');

close all

xl = [-0.1 0.3];
% xl = [-0.05 0.2];

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,1+mean(all_gsac_Egain(use_gsac_Ikerns,:)),std(all_gsac_Egain(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_gsac_Igain(use_gsac_Ikerns,:)),std(all_gsac_Igain(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','r'});
h1=plot(slags*dt,mean(all_gsac_Nsmodinfo(use_gsac_Ikerns,:)),'k','linewidth',2);
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);

tax = (0:(flen-1))*dt + dt/2;
up_tax = linspace(tax(1),tax(end),500);
all_nEkerns_up = spline(tax,all_nEkerns,up_tax);
all_nIkerns_up = spline(tax,all_nIkerns,up_tax);

f2 = figure();
hold on
% h1=shadedErrorBar(up_tax,mean(all_nEkerns_up(use_gsac_SUs,:)),std(all_nEkerns_up(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
% h2=shadedErrorBar(up_tax,mean(all_nIkerns_up(use_gsac_Ikerns,:)),std(all_nIkerns_up(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','r'});
h1=shadedErrorBar(tax,mean(all_nEkerns(use_gsac_Ikerns,:)),std(all_nEkerns(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','b'});
h2=shadedErrorBar(tax,mean(all_nIkerns(use_gsac_Ikerns,:)),std(all_nIkerns(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','r'});
xlim([0 tax(end)])

f3 = figure();
subplot(2,1,1)
h1=shadedErrorBar(slags*dt,mean(all_EI_ratio(use_gsac_Ikerns,:)),std(all_EI_ratio(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','b'});
% ylim([0.8 1.2]);
ylim([-0.2 0.2]);
% line(xl,[1 1],'color','k');
line(xl,[0 0],'color','k');
xlim(xl);

subplot(2,1,2)
hold on
% h2=shadedErrorBar(slags*dt,mean(all_EI_mean(use_gsac_Ikerns,:)),std(all_EI_mean(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','r'});
h2=plot(slags*dt,1+mean(all_postgains(use_gsac_Ikerns,:)),'r');
h2=plot(slags*dt,mean(all_gsac_Nsmodinfo(use_gsac_Ikerns,:)),'k');
% h2=shadedErrorBar(slags*dt,1+mean(all_postgains(use_gsac_Ikerns,:)),std(all_postgains(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','r'});
% h2=shadedErrorBar(slags*dt,1+mean(all_postoffs(use_gsac_Ikerns,:)),std(all_postoffs(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','k'});
ylim([0.6 1.4]);
set(gca,'YAxisLocation','right');
line(xl,[1 1],'color','k');
xlim(xl);


% f4 = figure();
% plot3(slags*dt,mean(all_gsac_Egain(use_gsac_Ikerns,:))+1,mean(all_gsac_Igain(use_gsac_Ikerns,:))+1);
% 
% f5 = figure();
% h1=shadedErrorBar(EI_lags*dt,-mean(all_EI_xc(use_gsac_Ikerns,:)),std(all_EI_xc(use_gsac_Ikerns,:))/sqrt(length(use_gsac_Ikerns)),{'color','b'});

% sz_rate = 0;
% sz_offset =50;
% avg_Egains = 1 + all_gsac_Egain(use_gsac_Ikerns,:);
% avg_Igains = 1 + all_gsac_Igain(use_gsac_Ikerns,:);
% net_gain = sqrt(avg_Egains.^2 + avg_Igains.^2);
% % avg_Egains = avg_Egains./net_gain;
% % avg_Igains = avg_Igains./net_gain;
% 
% ulags = find(slags*dt >= -0.1 & slags*dt <= 0.3);
% cmap = jet(length(ulags));
% f6 = figure();hold on
% % plot3(slags*dt,mean(all_gsac_Egain(use_gsac_Ikerns,:))+1,mean(all_gsac_Igain(use_gsac_Ikerns,:))+1,'k','linewidth',1);
% % scatter3(slags(ulags)*dt,avg_Egains(ulags),avg_Igains(ulags),sz_rate*(1:length(ulags))+sz_offset,cmap,'linewidth',2);
% plot3(slags*dt,avg_Igains,avg_Egains,'k','linewidth',1);
% scatter3(slags(ulags)*dt,avg_Igains(ulags),avg_Egains(ulags),sz_rate*(1:length(ulags))+sz_offset,cmap,'linewidth',2);
% grid
% % ylim([0.6 1.1]); zlim([0.6 1.1]);
% 
% cur_Egains = 1 + all_gsac_Egain(use_gsac_Ikerns,:);
% cur_Igains = 1 + all_gsac_Igain(use_gsac_Ikerns,:);
% cur_EIangle = atan2(cur_Egains,cur_Igains);



% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_EIgains.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);
% figure(f1);
% xlim([0 0.15])
% fname = [fig_dir 'Gsac_EIgains_zoom.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%
% figufy(f2);
% fname = [fig_dir 'EI_tempkerns.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

fig_width = 3.5; rel_height = 1.6;
figufy(f3);
fname = [fig_dir 'Gsac_EI_ratio_new.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);

%% ANALYZE LAMINAR DEPENDENCIES
load('/home/james/Analysis/bruce/FINsac_mod/layer_boundaries/layer_classification.mat')
boundary_enums = [boundary_class(:).Expt_num];

cur_min_rate = 5; %in Hz (5)
cur_min_Nsacs = 250; % (100)
cur_use_gsac_SUs = find(avg_rates' >= cur_min_rate & N_gsacs >= cur_min_Nsacs);

flen = 15;
tax = (0:(flen-1))*dt + dt/2;
up_tax = linspace(tax(1),tax(end),500);
all_nEkerns_up = spline(tax,all_nEkerns,up_tax);
[~,emloc] = max(all_nEkerns_up,[],2);
emtime = up_tax(emloc);

SU_expt_nums = [all_SU_data(:).expt_num];
un_lem_expts = unique(SU_expt_nums(lem_SUs));
n_lem_expts = length(un_lem_expts);

all_gran_SUs = [];
all_infra_SUs = [];
all_supra_SUs = [];
for ee = 1:n_lem_expts
    cur_bound_info = find(boundary_enums == un_lem_expts(ee),1);
    cur_ub = boundary_class(cur_bound_info).ub;
    cur_lb = boundary_class(cur_bound_info).lb;
    
    gran_probes = (cur_ub+1):(cur_lb-1);
    supra_probes = 1:(cur_ub-1);
    infra_probes = (cur_lb+1):24;
    
    cur_SU_set = find(SU_expt_nums == un_lem_expts(ee));
    cur_SU_probenums = arrayfun(@(x) x.ModData.unit_data.probe_number,all_SU_data(cur_SU_set));
    all_gran_SUs = cat(2,all_gran_SUs,cur_SU_set(ismember(cur_SU_probenums,gran_probes)));
    all_infra_SUs = cat(2,all_infra_SUs,cur_SU_set(ismember(cur_SU_probenums,infra_probes)));
    all_supra_SUs = cat(2,all_supra_SUs,cur_SU_set(ismember(cur_SU_probenums,supra_probes)));
end

all_gran_SUs = all_gran_SUs(ismember(all_gran_SUs,cur_use_gsac_SUs));
all_supra_SUs = all_supra_SUs(ismember(all_supra_SUs,cur_use_gsac_SUs));
all_infra_SUs = all_infra_SUs(ismember(all_infra_SUs,cur_use_gsac_SUs));

f1 = figure();
hold on
shadedErrorBar(up_tax,nanmean(all_nEkerns_up(all_gran_SUs,:)),nanstd(all_nEkerns_up(all_gran_SUs,:))/sqrt(length(all_gran_SUs)),{'color','k'});
shadedErrorBar(up_tax,nanmean(all_nEkerns_up(all_supra_SUs,:)),nanstd(all_nEkerns_up(all_supra_SUs,:))/sqrt(length(all_supra_SUs)),{'color','b'});
shadedErrorBar(up_tax,nanmean(all_nEkerns_up(all_infra_SUs,:)),nanstd(all_nEkerns_up(all_infra_SUs,:))/sqrt(length(all_infra_SUs)),{'color','r'});
xlim([0 0.12])

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Lamdep_stimkerns.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


%% SCATTERPLOT OF STIM TIMING VS SACMOD TIMING
cur_min_rate = 5; %in Hz (5)
cur_min_Nsacs = 250; % (100)
cur_use_gsac_SUs = find(avg_rates' >= cur_min_rate & N_gsacs >= cur_min_Nsacs);


flen = 15;
tax = (0:(flen-1))*dt + dt/2;
up_tax = linspace(tax(1),tax(end),500);
all_nEkerns_up = spline(tax,all_nEkerns,up_tax);
[~,emloc] = max(all_nEkerns_up,[],2);
emtime = up_tax(emloc);

up_tlags = linspace(tlags(1),tlags(end),5000);
all_gsac_Trates_up = spline(tlags,all_gsac_Trates,up_tlags);

search_range = [0 0.2];
% [gsac_Sfact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_Trates-1),tlags,search_range);
[gsac_Sfact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_Trates_up-1),up_tlags,search_range);
search_range = [0 0.3];
% [gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_Trates-1,tlags,search_range);
[gsac_Efact,gsac_exctime] = get_tavg_peaks(all_gsac_Trates_up-1,up_tlags,search_range);

reverse_polarity = find(gsac_exctime < gsac_inhtime);
uu = cur_use_gsac_SUs(~ismember(cur_use_gsac_SUs,reverse_polarity));

emtime = emtime*1e3;
gsac_inhtime = gsac_inhtime*1e3;

xl = [0.02 0.1]*1e3;
yl = [0.02 0.1]*1e3;

f1 = figure();
hold on
% plot(emtime(cur_use_gsac_SUs),gsac_inhtime(cur_use_gsac_SUs),'k.')
plot(emtime(uu),gsac_inhtime(uu),'r.')
% plot(emtime(all_gran_SUs),gsac_inhtime(all_gran_SUs),'k.')
% plot(emtime(all_supra_SUs),gsac_inhtime(all_supra_SUs),'b.')
% plot(emtime(all_infra_SUs),gsac_inhtime(all_infra_SUs),'r.')
% r = robustfit(emtime(cur_use_gsac_SUs),gsac_inhtime(cur_use_gsac_SUs));
r = robustfit(emtime(uu),gsac_inhtime(uu));
% r = regress(gsac_inhtime(uu),[ones(length(uu),1) emtime(uu)']);
xlim(xl); ylim(yl);
xx = linspace(xl(1),xl(2),100);
plot(xx,r(1)+r(2)*xx,'r');
line(xl,yl,'color','k');
xlabel('Time lag (s)');
ylabel('Suppression timing (s)');

lem_uset = uu(ismember(uu,lem_SUs));
group_id = nan(length(all_SU_data));
group_id(all_supra_SUs) = 1;
group_id(all_gran_SUs) = 2;
group_id(all_infra_SUs) = 3;
aoctool(emtime(lem_uset),gsac_inhtime(lem_uset),group_id(lem_uset),0.05,[],[],[],'on',4);

% xl = [0.02 0.12];
% yl = [0.1 0.3];
%
% f2 = figure();
% hold on
% plot(emtime(cur_use_gsac_SUs),gsac_exctime(cur_use_gsac_SUs),'k.')
% r = robustfit(emtime(cur_use_gsac_SUs),gsac_exctime(cur_use_gsac_SUs));
% xlim(xl); ylim(yl);
% xx = linspace(xl(1),xl(2),100);
% plot(xx,r(1)+r(2)*xx,'r');
% line(xl,yl,'color','k');
% xlabel('Time lag (s)');
% ylabel('Suppression timing (s)');


% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Stim_sac_timing_compare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
%%
PRM = arrayfun(@(x) x.ModData.tune_props.PRM,all_SU_data);
RF_FTF = arrayfun(@(x) x.ModData.tune_props.RF_FTF,all_SU_data);
RF_FSF = arrayfun(@(x) x.ModData.tune_props.RF_FSF,all_SU_data);
net_phase_polarity = arrayfun(@(x) x.ModData.tune_props.net_phase_polarity,all_SU_data);
RF_ecc = arrayfun(@(x) x.ModData.tune_props.RF_ecc,all_SU_data);
RF_dirsel = arrayfun(@(x) x.ModData.tune_props.RF_dirsel,all_SU_data);
RF_gSF = arrayfun(@(x) x.ModData.tune_props.RF_gSF,all_SU_data);


%% INFORMATION TIMING ANALYSIS
info_tax = all_SU_timedata(1).lag_axis;
uu = use_gsac_SUs;
uu = use_jbe_SUs;
% uu = use_lem_SUs;
all_info_before = reshape([all_SU_timedata(uu).info_before],[],length(uu))';
all_info_after = reshape([all_SU_timedata(uu).info_after],[],length(uu))';
all_info_during = reshape([all_SU_timedata(uu).info_during],[],length(uu))';
all_Binfo_before = reshape([all_SU_timedata(uu).base_info_before],[],length(uu))';
all_Binfo_after = reshape([all_SU_timedata(uu).base_info_after],[],length(uu))';
all_Binfo_during = reshape([all_SU_timedata(uu).base_info_during],[],length(uu))';

baseline_info = mean(all_Binfo_before(:,1:20),2);

norm_info_before = bsxfun(@rdivide,all_info_before,baseline_info);
norm_info_after = bsxfun(@rdivide,all_info_after,baseline_info);
norm_info_during = bsxfun(@rdivide,all_info_during,baseline_info);

norm_Binfo_before = bsxfun(@rdivide,all_Binfo_before,baseline_info);
norm_Binfo_after = bsxfun(@rdivide,all_Binfo_after,baseline_info);
norm_Binfo_during = bsxfun(@rdivide,all_Binfo_during,baseline_info);

f1 = figure(); hold on
shadedErrorBar(info_tax,mean(norm_info_before),std(norm_info_before)/sqrt(length(uu)),{'color','r'});
shadedErrorBar(info_tax,mean(norm_info_after),std(norm_info_after)/sqrt(length(uu)),{'color','b'});
shadedErrorBar(info_tax,mean(norm_info_during),std(norm_info_during)/sqrt(length(uu)),{'color','k'});

% f2 = figure(); hold on
% shadedErrorBar(info_tax,mean(norm_Binfo_before(uu,:)),std(norm_Binfo_before(uu,:))/sqrt(length(uu)),{'color','r'});
% shadedErrorBar(info_tax,mean(norm_Binfo_after(uu,:)),std(norm_Binfo_after(uu,:))/sqrt(length(uu)),{'color','b'});
% shadedErrorBar(info_tax,mean(norm_Binfo_during(uu,:)),std(norm_Binfo_during(uu,:))/sqrt(length(uu)),{'color','k'});
plot(info_tax,mean(norm_Binfo_before),'r--','linewidth',1);
plot(info_tax,mean(norm_Binfo_after),'b--','linewidth',1);
plot(info_tax,mean(norm_Binfo_during),'k--','linewidth',1);
xlim([-0.2 0.2]);
ylim([0 1.2])
line([0 0],[0 1.2],'color','k')
xlabel('Time since fixation onset (s)');
ylabel('Relative stim info');

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Stim_info_timing.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


%% COMPARE GSAC AND MSAC GAIN/OFFSET FILTERS AND AVG RATES

all_gsac_postgains = [];
all_gsac_postoffs = [];
all_msac_postgains = [];
all_msac_postoffs = [];
for ii = use_gsac_SUs
    all_gsac_postgains = cat(1,all_gsac_postgains,all_SU_data(ii).gsac_post_singmod.mods(3).filtK');
    all_gsac_postoffs = cat(1,all_gsac_postoffs,all_SU_data(ii).gsac_post_singmod.mods(2).filtK');
end
for ii = use_msac_SUs
    all_msac_postgains = cat(1,all_msac_postgains,all_SU_data(ii).msac_post_singmod.mods(3).filtK');
    all_msac_postoffs = cat(1,all_msac_postoffs,all_SU_data(ii).msac_post_singmod.mods(2).filtK');
end

all_gsac_spost_offset = reshape([all_SU_data(use_gsac_SUs).gsac_spost_offset],[],length(use_gsac_SUs))';
all_gsac_spost_offset = bsxfun(@rdivide,all_gsac_spost_offset,[all_SU_data(use_gsac_SUs).gsac_ovavg_rate]');
all_gsac_spost_gain = reshape([all_SU_data(use_gsac_SUs).gsac_spost_gain],[],length(use_gsac_SUs))';
all_msac_spost_offset = reshape([all_SU_data(use_msac_SUs).msac_spost_offset],[],length(use_msac_SUs))';
all_msac_spost_offset = bsxfun(@rdivide,all_msac_spost_offset,[all_SU_data(use_msac_SUs).gsac_ovavg_rate]');
all_msac_spost_gain = reshape([all_SU_data(use_msac_SUs).msac_spost_gain],[],length(use_msac_SUs))';

search_range = [0 0.2];
[gsac_Sfact,gsac_inhtime] = get_tavg_peaks(-(all_gsac_Trates-1),tlags,search_range);
[msac_Sfact,msac_inhtime] = get_tavg_peaks(-(all_msac_Trates-1),tlags,search_range);
[gsac_gain_Sfact,gsac_gain_inhtime] = get_tavg_peaks(-(all_gsac_postgains),slags*dt,search_range);
[msac_gain_Sfact,msac_gain_inhtime] = get_tavg_peaks(-(all_msac_postgains),slags*dt,search_range);
search_range = [0 0.35];
[gsac_Efact,gsac_exctime] = get_tavg_peaks((all_gsac_Trates-1),tlags,search_range);
[msac_Efact,msac_exctime] = get_tavg_peaks((all_msac_Trates-1),tlags,search_range);

use_both_SUs = intersect(use_gsac_SUs,use_msac_SUs);
gboth_SUs = find(ismember(use_gsac_SUs,use_both_SUs));
mboth_SUs = find(ismember(use_msac_SUs,use_both_SUs));

xl = [-0.1 0.3];

yl = [0.7 1.3];
f1 = figure();
hold on
h1=shadedErrorBar(tlags,mean(all_gsac_Trates(gboth_SUs,:)),std(all_gsac_Trates(gboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','b'});
h2=shadedErrorBar(tlags,mean(all_msac_Trates(mboth_SUs,:)),std(all_msac_Trates(mboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Relative gain');
xlim(xl);
line(xl,[1 1],'color','k');
ylim(yl);
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

f2 = figure();
hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_postgains(gboth_SUs,:)),std(all_gsac_postgains(gboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_msac_postgains(mboth_SUs,:)),std(all_msac_postgains(mboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','r'});
h1=shadedErrorBar(slags*dt,mean(all_gsac_spost_gain(gboth_SUs,:)),std(all_gsac_spost_gain(gboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_msac_spost_gain(mboth_SUs,:)),std(all_msac_spost_gain(mboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Relative gain');
line(xl,[1 1],'color','k');
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Gain');


f3 = figure();
hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_postoffs(gboth_SUs,:)),std(all_gsac_postoffs(use_both_SUs,:))/sqrt(length(use_both_SUs)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_msac_postoffs(mboth_SUs,:)),std(all_msac_postoffs(mboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','r'});
h1=shadedErrorBar(slags*dt,mean(all_gsac_spost_offset(gboth_SUs,:)),std(all_gsac_spost_offset(gboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_msac_spost_offset(mboth_SUs,:)),std(all_msac_spost_offset(mboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Relative gain');
line(xl,[0 0],'color','k');
xlim(xl);
yl = ylim();
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Gain');

% figufy(f1);
% fig_width = 3.5; rel_height = 0.8;
% fname = [fig_dir 'Gsac_msac_rate.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'Gsac_msac_gains.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'Gsac_msac_offset.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
% 
%% COMPARE GSAC AND MSAC SS INFO AND INFO RATES (POST MODEL BASED)

all_gsac_smodinfo = reshape([all_SU_data(use_gsac_SUs).gsac_spost_modinfo],[],length(use_gsac_SUs))';
all_gsac_ov_smodinfo = [all_SU_data(use_gsac_SUs).gsac_spost_ov_modinfo];
all_gsac_smodinforate = all_gsac_smodinfo.*all_gsac_rates(use_gsac_SUs,:)*dt;
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');
all_gsac_Nsmodinforate = bsxfun(@rdivide,all_gsac_smodinforate,all_gsac_ov_smodinfo'.*avg_rates(use_gsac_SUs)*dt);


all_msac_rates = reshape([all_SU_data(use_msac_SUs).msac_avg_rate],[],length(use_msac_SUs))'/dt;
all_msac_smodinfo = reshape([all_SU_data(use_msac_SUs).msac_spost_modinfo],[],length(use_msac_SUs))';
all_msac_ov_smodinfo = [all_SU_data(use_msac_SUs).msac_spost_ov_modinfo];
all_msac_smodinforate = all_msac_smodinfo.*all_msac_rates*dt;
all_msac_Nsmodinfo = bsxfun(@rdivide,all_msac_smodinfo,all_msac_ov_smodinfo');
all_msac_Nsmodinforate = bsxfun(@rdivide,all_msac_smodinforate,all_msac_ov_smodinfo'.*avg_rates(use_msac_SUs)*dt);

use_both_SUs = intersect(use_gsac_SUs,use_msac_SUs);
use_both_jbe = intersect(use_both_SUs,jbe_SUs);
use_both_lem = intersect(use_both_SUs,lem_SUs);
use_msac_jbe = intersect(use_msac_SUs,jbe_SUs);
use_msac_lem = intersect(use_msac_SUs,lem_SUs);
gboth_SUs = find(ismember(use_gsac_SUs,use_both_SUs));
mboth_SUs = find(ismember(use_msac_SUs,use_both_SUs));

all_msac_rates = reshape([all_SU_data(use_msac_SUs).msac_avg_rate]/dt,[],length(use_msac_SUs))';
all_msac_nrates = bsxfun(@rdivide,all_msac_rates,avg_rates(use_msac_SUs));

search_range = [0 0.2];
[gsac_info_Sfact] = get_tavg_peaks(-(all_gsac_Nsmodinfo-1),slags*dt,search_range);
[msac_info_Sfact] = get_tavg_peaks(-(all_msac_Nsmodinfo-1),slags*dt,search_range);
[gsac_Sfact,gsac_Stime] = get_tavg_peaks(-(all_gsac_nrates-1),slags*dt,search_range);
[msac_Sfact,msac_Stime] = get_tavg_peaks(-(all_msac_nrates-1),slags*dt,search_range);

search_range = [0 0.35];
[gsac_Efact,gsac_Etime] = get_tavg_peaks((all_gsac_nrates-1),slags*dt,search_range);
[msac_Efact,msac_Etime] = get_tavg_peaks((all_msac_nrates-1),slags*dt,search_range);


xl = [-0.1 0.3];
ms = 3;

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(gboth_SUs,:)),std(all_gsac_Nsmodinfo(gboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_msac_Nsmodinfo(mboth_SUs,:)),std(all_msac_Nsmodinfo(mboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','r'});
% legend([h1.mainLine h2.mainLine],{'TB','Mod-pred'});
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);
ylim([0.5 1.1])
yl = ylim();
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Relative information');

f2 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinforate(gboth_SUs,:)),std(all_gsac_Nsmodinforate(gboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_msac_Nsmodinforate(mboth_SUs,:)),std(all_msac_Nsmodinforate(mboth_SUs,:))/sqrt(length(use_both_SUs)),{'color','r'});
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);
% ylim([0.5 1.1])
yl = ylim();
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Relative information');

% figufy(f1);
% fig_width = 3.5; rel_height = 0.8;
% fname = [fig_dir 'Gsac_msac_info.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fig_width = 3.5; rel_height = 0.8;
% fname = [fig_dir 'Gsac_msac_inforate.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% RF WIDTH VS SUPPRESSION STRENGTH
close all
rf_sigma = arrayfun(@(x) x.ModData.tune_props.RF_sigma,all_SU_data);
rf_ecc = arrayfun(@(x) x.ModData.tune_props.RF_ecc,all_SU_data);
TB_slags = all_SU_data(1).gsac_TB_lagX;

% all_gsac_smodinfo = reshape([all_SU_data(:).gsac_spost_modinfo],[],length(all_SU_data))';
% all_gsac_ov_smodinfo = [all_SU_data(:).gsac_spost_ov_modinfo];
% all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');
all_gsac_smodinfo = reshape([all_SU_data(:).gsac_TB_info],[],length(all_SU_data))';
all_gsac_ov_smodinfo = [all_SU_data(:).gsac_ov_TB_info];
all_gsac_Nsmodinfo = bsxfun(@rdivide,all_gsac_smodinfo,all_gsac_ov_smodinfo');

all_gsac_TB_offset = reshape([all_SU_data(:).gsac_TB_offset],[],length(all_SU_data))';
all_gsac_TB_offset = bsxfun(@rdivide,all_gsac_TB_offset,[all_SU_data(:).gsac_ovavg_rate]');
all_gsac_TB_gain = reshape([all_SU_data(:).gsac_TB_gain],[],length(all_SU_data))';

search_range = [0 0.2];
[gsac_info_Sfact] = get_tavg_peaks(-(all_gsac_Nsmodinfo-1),TB_slags*dt,search_range);
[gsac_gain_Sfact] = get_tavg_peaks(-(all_gsac_TB_gain-1),TB_slags*dt,search_range);
[gsac_Sfact] = get_tavg_peaks(-(all_gsac_nrates-1),slags*dt,search_range);
search_range = [0 0.35];
[gsac_Efact] = get_tavg_peaks((all_gsac_nrates-1),slags*dt,search_range);
[gsac_offset_Efact] = get_tavg_peaks((all_gsac_TB_offset),TB_slags*dt,search_range);


ms = 3;

cur_exclude = find(rf_sigma>0.5);
cur_use_jbe = use_jbe_SUs(~ismember(use_jbe_SUs,cur_exclude));
cur_use_lem = use_lem_SUs(~ismember(use_lem_SUs,cur_exclude));
cur_use = use_gsac_SUs(~ismember(use_gsac_SUs,cur_exclude));


f1 = figure();
hold on
plot(rf_sigma(cur_use_jbe)*2,gsac_info_Sfact(cur_use_jbe),'o','markersize',ms);
plot(rf_sigma(cur_use_lem)*2,gsac_info_Sfact(cur_use_lem),'ro','markersize',ms);
xlabel('RF size (deg)');
ylabel('Information suppression');
xlim([0.08 0.5])
set(gca,'xscale','log');

ff1 = figure();
hold on
plot(rf_ecc(cur_use_jbe),gsac_info_Sfact(cur_use_jbe),'o','markersize',ms);
plot(rf_ecc(cur_use_lem),gsac_info_Sfact(cur_use_lem),'ro','markersize',ms);
xlabel('RF eccentricity (deg)');
ylabel('Information suppression');
xlim([0 5])

ff2 = figure();
hold on
plot(rf_ecc(cur_use_jbe),gsac_offset_Efact(cur_use_jbe),'o','markersize',ms);
plot(rf_ecc(cur_use_lem),gsac_offset_Efact(cur_use_lem),'ro','markersize',ms);
xlabel('RF eccentricity (deg)');
ylabel('Information suppression');
xlim([0 5])

f2 = figure();
hold on
plot(rf_sigma(cur_use_jbe)*2,gsac_Sfact(cur_use_jbe),'o','markersize',ms);
plot(rf_sigma(cur_use_lem)*2,gsac_Sfact(cur_use_lem),'ro','markersize',ms);
xlabel('RF size (deg)');
ylabel('Firing rate suppression');
xlim([0.08 0.5])
set(gca,'xscale','log');

f3 = figure();
hold on
plot(rf_ecc(cur_use_jbe),gsac_info_Sfact(cur_use_jbe),'o','markersize',ms);
plot(rf_ecc(cur_use_lem),gsac_info_Sfact(cur_use_lem),'ro','markersize',ms);
xlabel('RF eccentricity (deg)');
ylabel('Information suppression');

% f4 = figure();
% hold on
% plot(rf_ecc(cur_use_jbe),gsac_Sfact(cur_use_jbe),'o','markersize',ms);
% plot(rf_ecc(cur_use_lem),gsac_Sfact(cur_use_lem),'ro','markersize',ms);
% xlabel('RF eccentricity (deg)');
% ylabel('Firing rate suppression');

f5 = figure();hold on
plot(rf_ecc(cur_use_jbe),rf_sigma(cur_use_jbe)*2,'o','markersize',ms);
plot(rf_ecc(cur_use_lem),rf_sigma(cur_use_lem)*2,'ro','markersize',ms);
ylim([0.08 0.5])
xlabel('Eccentricity (deg)');
ylabel('RF size (deg)');
set(gca,'yscale','log');

% all_msac_rates = reshape([all_SU_data(:).msac_avg_rate]/dt,[],length(all_SU_data))';
% all_msac_nrates = bsxfun(@rdivide,all_msac_rates,avg_rates);
% all_msac_smodinfo = reshape([all_SU_data(:).msac_spost_modinfo],[],length(all_SU_data))';
% all_msac_ov_smodinfo = [all_SU_data(:).msac_spost_ov_modinfo];
% all_msac_Nsmodinfo = bsxfun(@rdivide,all_msac_smodinfo,all_msac_ov_smodinfo');
% 
% search_range = [0 0.2];
% [msac_info_Sfact] = get_tavg_peaks(-(all_msac_Nsmodinfo-1),slags*dt,search_range);
% [msac_Sfact] = get_tavg_peaks(-(all_msac_nrates-1),slags*dt,search_range);
% [msac_Efact] = get_tavg_peaks((all_msac_nrates-1),slags*dt,search_range);


rf_sigma_prc = prctile(rf_sigma(cur_use),[25 50 75]);
% set1 = cur_use(rf_sigma(cur_use) <= rf_sigma_prc(1));
% set2 = cur_use(rf_sigma(cur_use) >= rf_sigma_prc(end));
set1 = cur_use(rf_sigma(cur_use) <= rf_sigma_prc(2));
set2 = cur_use(rf_sigma(cur_use) > rf_sigma_prc(2));

xl = [-0.1 0.3];
f6 = figure;hold on
h1=shadedErrorBar(tlags,mean(all_gsac_Trates(set1,:)),std(all_gsac_Trates(set1,:))/sqrt(length(set1)),{'color','b'});
h2=shadedErrorBar(tlags,mean(all_gsac_Trates(set2,:)),std(all_gsac_Trates(set2,:))/sqrt(length(set2)),{'color','r'});
xlim(xl);

xl = [-0.1 0.3];
TB_slags = all_SU_data(1).gsac_TB_lagX;
f7 = figure;hold on
h1=shadedErrorBar(TB_slags*dt,mean(all_gsac_Nsmodinfo(set1,:)),std(all_gsac_Nsmodinfo(set1,:))/sqrt(length(set1)),{'color','b'});
h2=shadedErrorBar(TB_slags*dt,mean(all_gsac_Nsmodinfo(set2,:)),std(all_gsac_Nsmodinfo(set2,:))/sqrt(length(set2)),{'color','r'});
xlim(xl);


fig_width = 3.5; rel_height = 1;

% figufy(f1);
% fname = [fig_dir 'rfwidth_info.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'rfwidth_rate.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'rfecc_rate.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
% 
% figufy(f5);
% fname = [fig_dir 'rfecc_width.pdf'];
% exportfig(f5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f5);
% 
% fig_width = 3.5; rel_height = 1;
% figufy(f6);
% fname = [fig_dir 'rfwidthdep_rate.pdf'];
% exportfig(f6,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f6);
% 
% figufy(f7);
% fname = [fig_dir 'rfwidthdep_info.pdf'];
% exportfig(f7,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f7);
% 
%% analyze upstream changes in temporal integration

clear max_* weight_* *_mlocs base_*

flen = 15;
fpast = 0:-1:(-flen+1);
base_kern = ones(flen,1);

for cc = 1:length(all_SU_data)
    
    stim_mod_signs = [all_SU_data(cc).ModData.rectGQM.mods(:).sign];
    cur_filts = reshape([all_SU_data(cc).ModData.rectGQM.mods(:).filtK],flen,[],length(stim_mod_signs));
    filt_tkerns = squeeze(std(cur_filts,[],2));
    filt_tkerns = bsxfun(@rdivide,filt_tkerns,sqrt(sum(filt_tkerns.^2)));
    filt_tkerns = bsxfun(@times,filt_tkerns,all_SU_data(cc).ModData.rectGQM.rel_filt_weights);
    filt_tkerns(:,all_SU_data(cc).rel_filt_weights==0) = 0;
    
    avg_tkern = squeeze(mean(filt_tkerns,2));
    avg_ekern = squeeze(mean(filt_tkerns(:,stim_mod_signs==1),2));
    avg_ikern = squeeze(mean(filt_tkerns(:,stim_mod_signs==-1),2));
    
    sacGainMod = all_SU_data(cc).gsacGainMod;
    base_filt = sacGainMod.stim_kernel;
    sacdep_tkerns = nan(length(slags),flen);
    for ii = 1:length(slags)
        cur_kern = base_kern;
        [lia,lib] = ismember(fpast+slags(ii),slags);
        cur_kern(lia) = 1+base_filt(lib(lia));
        sacdep_tkerns(ii,:) = cur_kern;
    end
    
    tlags = 0:(flen-1);
    up_tlags = 0:0.1:(flen-1);
    spline_sac_tkerns = spline(tlags,bsxfun(@times,sacdep_tkerns,avg_tkern'),up_tlags);
    spline_sac_ekerns = spline(tlags,bsxfun(@times,sacdep_tkerns,avg_ekern'),up_tlags);
    spline_sac_ikerns = spline(tlags,bsxfun(@times,sacdep_tkerns,avg_ikern'),up_tlags);
    
    weight_avg_tlag(cc,:) = sum(bsxfun(@times,spline_sac_tkerns,up_tlags),2)./sum(spline_sac_tkerns,2);
    weight_avg_elag(cc,:) = sum(bsxfun(@times,spline_sac_ekerns,up_tlags),2)./sum(spline_sac_ekerns,2);
    weight_avg_ilag(cc,:) = sum(bsxfun(@times,spline_sac_ikerns,up_tlags),2)./sum(spline_sac_ikerns,2);
    [mvals,tkern_mlocs] = max(spline_sac_tkerns,[],2);
    [mvals,ekern_mlocs] = max(spline_sac_ekerns,[],2);
    [mvals,ikern_mlocs] = max(spline_sac_ikerns,[],2);
    
    max_tkerns(cc,:) = up_tlags(tkern_mlocs);
    max_ekerns(cc,:) = up_tlags(ekern_mlocs);
    max_ikerns(cc,:) = up_tlags(ikern_mlocs);
    
    base_tkern_weight(cc) = sum(tlags.*avg_tkern')/sum(avg_tkern);
    base_ekern_weight(cc) = sum(tlags.*avg_ekern')/sum(avg_ekern);
    base_ikern_weight(cc) = sum(tlags.*avg_ikern')/sum(avg_ikern);
    
    [~,temp] = max(spline(tlags,avg_tkern',up_tlags));
    base_tkern_mloc(cc) = up_tlags(temp);
    [~,temp] = max(spline(tlags,avg_ekern',up_tlags));
    base_ekern_mloc(cc) = up_tlags(temp);
    [~,temp] = max(spline(tlags,avg_ikern',up_tlags));
    base_ikern_mloc(cc) = up_tlags(temp);
    
end

weight_avg_tlag = bsxfun(@minus,weight_avg_tlag,base_tkern_weight')*dt*1e3;
weight_avg_elag = bsxfun(@minus,weight_avg_elag,base_ekern_weight')*dt*1e3;
weight_avg_ilag = bsxfun(@minus,weight_avg_ilag,base_ikern_weight')*dt*1e3;

max_tkerns = bsxfun(@minus,max_tkerns,base_tkern_mloc')*dt*1e3;
max_ekerns = bsxfun(@minus,max_ekerns,base_ekern_mloc')*dt*1e3;
max_ikerns = bsxfun(@minus,max_ikerns,base_ikern_mloc')*dt*1e3;

f1 = figure();
plot(slags*dt,weight_avg_tlag(use_gsac_SUs,:),'k')
hold on
h1=shadedErrorBar(slags*dt,mean(weight_avg_tlag(use_gsac_SUs,:)),std(weight_avg_tlag(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
h1=shadedErrorBar(slags*dt,mean(weight_avg_elag(use_gsac_SUs,:)),std(weight_avg_elag(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
h1=shadedErrorBar(slags*dt,mean(weight_avg_ilag(use_gsac_SUs,:)),std(weight_avg_ilag(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','g'});
xlabel('Time since saccade onset (s)');
ylabel('Change in stimulus response latency');

% f1 = figure();
% hold on
% h1=shadedErrorBar(slags*dt,mean(max_tkerns(use_gsac_SUs,:)),std(max_tkerns(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','b'});
% h1=shadedErrorBar(slags*dt,median(max_ekerns(use_gsac_SUs,:)),std(max_ekerns(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','r'});
% h1=shadedErrorBar(slags*dt,mean(max_ikerns(use_gsac_SUs,:)),std(max_ikerns(use_gsac_SUs,:))/sqrt(length(use_gsac_SUs)),{'color','g'});

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_tempdelay_changes.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);



%% analyze any mixing of filter amps in subspace models
all_inds = [];
all_filts = [];
all_qfilts1 = [];
all_qfilts2 = [];
all_angs = [];
unit_avg_angs = [];
unit_other_angs = [];
unit_nmods = [];
max_dims = 8;
base_filts = nan(length(slags),max_dims,max_dims);
for cc = 1:length(all_SU_data)
    cur_submod = all_SU_data(cc).gsac_submod;
    cur_mod_signs = [cur_submod.mods(2:end).sign];
    if ismember(cc,use_gsac_SUs)
        cur_filts = reshape([cur_submod.mods(2:end).filtK],length(slags),[],length(cur_mod_signs));
        filts_to_add = base_filts;
        filts_to_add(:,1:length(cur_mod_signs),1:length(cur_mod_signs)) = cur_filts;
        all_filts = cat(4,all_filts,filts_to_add);
        all_inds = cat(1,all_inds,cc);
        
        cur_nmods = length(cur_mod_signs);
        filt_mags = sqrt(sum(cur_filts.^2,3));
        filt_ang = nan(cur_nmods,length(slags));
        otherfilt_ang = nan(cur_nmods,length(slags));
        for jj = 1:cur_nmods
            temp = bsxfun(@rdivide,squeeze(cur_filts(:,:,jj)),filt_mags(:,jj));
            filt_ang(jj,:) = temp(:,jj);
            others = setdiff(1:cur_nmods,jj);
            otherfilt_ang(jj,:) = mean(temp(:,others),2);
        end
        all_angs = cat(1,all_angs,filt_ang);
        unit_avg_angs = cat(1,unit_avg_angs,nanmean(filt_ang));
        unit_other_angs = cat(1,unit_other_angs,nanmean(otherfilt_ang));
        unit_nmods = cat(1,unit_nmods,cur_nmods);
    end
end

avg_subfilts = squeeze(nanmean(all_filts,4));
% avg_subfilts = avg_subfilts(:,[1 end 2:end-1],:);
avg_subfilts = permute(avg_subfilts,[2 3 1]);

for ii = 1:max_dims
    avg_subfilts(:,ii,:) = circshift(squeeze(avg_subfilts(:,ii,:)),-(ii-1));
end
% close all
% f1 = figure();
% for ii = 1:max_dims
% subplot(4,2,ii)
% imagesc(slags*dt,1:8,squeeze(avg_subfilts(:,ii,:))')
% caxis([-1.2 1.2])
% end
% 
avg_nulllev = mean(1./sqrt(unit_nmods));
f2 = figure();
hold on
shadedErrorBar(slags*dt,nanmean(unit_avg_angs),nanstd(unit_avg_angs)/size(unit_avg_angs,1),{'color','r'});
shadedErrorBar(slags*dt,nanmean(unit_other_angs),nanstd(unit_other_angs)/size(unit_other_angs,1));
ylim([0 1]);
xl = xlim();
% line(xl,[avg_nulllev avg_nulllev]);

% 
% fig_width = 5; rel_height = 1.2;
% figufy(f1);
% fname = [fig_dir 'Gsac_subspace_coefs.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
fig_width = 3.5; rel_height = 0.8;
figufy(f2);
fname = [fig_dir 'Gsac_subspace_angles.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

%% Microsac PRE/POST MODEL
all_msac_gaink = reshape(cell2mat(arrayfun(@(x) x.msacGainMod.gain_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';
all_msac_stimk = reshape(cell2mat(arrayfun(@(x) x.msacGainMod.stim_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';
all_msac_offk = reshape(cell2mat(arrayfun(@(x) x.msacGainMod.off_kernel,all_SU_data,'uniformoutput',0)),[],length(all_SU_data))';
all_msac_gainoff = arrayfun(@(x) x.msacGainMod.gain_offset,all_SU_data)-1;

% all_msac_gaink = bsxfun(@minus,all_msac_gaink,all_msac_gainoff);

all_postgains = [];
all_postoffs = [];
for ii = 1:length(all_SU_data)
    all_postgains = cat(1,all_postgains,all_SU_data(ii).msac_post_singmod.mods(3).filtK');
    all_postoffs = cat(1,all_postoffs,all_SU_data(ii).msac_post_singmod.mods(2).filtK');
end

xl = [-0.1 0.3];
yl = [0.6 1.2];

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,1+mean(all_msac_gaink(use_msac_SUs,:)),std(all_msac_gaink(use_msac_SUs,:))/sqrt(length(use_msac_SUs)),{'color','b'});
h2=shadedErrorBar(slags*dt,1+mean(all_msac_stimk(use_msac_SUs,:)),std(all_msac_stimk(use_msac_SUs,:))/sqrt(length(use_msac_SUs)),{'color','k'});
h3=plot(slags*dt,1+mean(all_postgains(use_msac_SUs,:)),'m','linewidth',2);
h4=plot(slags*dt,1+mean(all_msac_offk(use_msac_SUs,:)),'r','linewidth',2);
xlabel('Time (s)');
ylabel('Relative info');
line(xl,[1 1],'color','k');
xlim(xl);
ylim(yl);
line([0 0],yl,'color','k');

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_PREPOST_modcompare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
%


%% compare stim timing and gain supp for full models
flen = 15;
tax = (0:(flen-1))*dt + dt/2;
search_range = [0 0.2];
all_stim_kerns = [];
all_fpost_gains = [];
all_rel_weights = [];
all_slopes = [];
for ii = 1:length(all_SU_data)
    if ismember(ii,use_gsac_SUs)
        cur_mod = all_SU_data(ii).ModData.rectGQM;
        mod_signs = [cur_mod.mods(:).sign];
        sd = cur_mod.stim_params(1).stim_dims;
        mod_filts = reshape([cur_mod.mods(:).filtK],[sd(1) sd(2) length(mod_signs)]);
        
        rel_weights = [all_SU_data(ii).rel_filt_weights];
        
%         tkerns = squeeze(std(mod_filts,[],2));
        tkerns = nan(flen,length(mod_signs));
        for jj = 1:length(mod_signs)
           tkerns(:,jj) = squeeze(mean(abs(hilbert(squeeze(mod_filts(:,:,jj)))),2)); 
        end

        tkerns = bsxfun(@rdivide,tkerns,sqrt(sum(tkerns.^2)));
        tkerns = bsxfun(@times,tkerns,rel_weights);
        
        
        tkerns(:,rel_weights == 0) = [];
        
        all_stim_kerns = cat(1,all_stim_kerns,tkerns');
        
        
        fpost_gains = reshape([all_SU_data(ii).gsac_post_Fmod.mods(3).filtK],length(slags),length(mod_signs));
        fpost_gains(:,rel_weights == 0) = [];
        all_fpost_gains = cat(1,all_fpost_gains,fpost_gains');
        
        mod_signs(rel_weights == 0) = [];
        rel_weights(rel_weights==0) = [];
        all_rel_weights = cat(1,all_rel_weights,rel_weights');
        if length(mod_signs) >= 5
            [~,tmaxloc] = max(tkerns);
            [cur_gainmax,cur_gainloc] = get_tavg_peaks(-(fpost_gains'-1),slags*dt,search_range);
            B = regress(cur_gainloc,[tax(tmaxloc)' ones(length(mod_signs),1)]);
            all_slopes = cat(1,all_slopes,B(1));
        end
    end
end


up_tax = linspace(tax(1),tax(end),500);
all_tkerns_up = spline(tax,all_stim_kerns,up_tax);
[tkern_max,tkern_maxloc] = max(all_tkerns_up,[],2);

up_tax = linspace(slags(1)*dt,slags(end)*dt,1e3);
all_fpost_gains_up = spline(slags*dt,all_fpost_gains,up_tax);
[fpost_gain_max,fpost_gain_loc] = get_tavg_peaks(-(all_fpost_gains_up-1),up_tax*dt,search_range);

fpost_minamp = 1.25;
tkern_amp = 0.05;

uset = find(fpost_gain_max > fpost_minamp & all_rel_weights > tkern_amp);

%%
clust_prop = clust_iso_dist;
clust_prc = prctile(clust_prop(use_gsac_SUs),[25 50 75]);
set1 = use_gsac_SUs(clust_prop(use_gsac_SUs) >= clust_prc(3));
set2 = use_gsac_SUs(clust_prop(use_gsac_SUs) <= clust_prc(1));

all_gsac_TBmodinfo = reshape([all_SU_data(:).gsac_TB_info],[],length(all_SU_data))';
all_gsac_ov_TBmodinfo = [all_SU_data(:).gsac_ov_TB_info];
all_gsac_NTBmodinfo = bsxfun(@rdivide,all_gsac_TBmodinfo,all_gsac_ov_TBmodinfo');
TB_lag = [all_SU_data(1).gsac_TB_lagX];

search_range = [0 0.2];
gsac_Sfac = get_tavg_peaks(-(all_gsac_Trates-1),tlags,search_range);
gsac_info_Sfac = get_tavg_peaks(-(all_gsac_NTBmodinfo-1),TB_lag*dt,search_range);
search_range = [0 0.3];
gsac_Efac = get_tavg_peaks(all_gsac_Trates,tlags,search_range);

xl = [-0.1 0.3];

f1 = figure();
hold on
h1=shadedErrorBar(tlags,mean(all_gsac_Trates(set1,:)),std(all_gsac_Trates(set1,:))/sqrt(length(set1)),{'color','b'});
h2=shadedErrorBar(tlags,mean(all_gsac_Trates(set2,:)),std(all_gsac_Trates(set2,:))/sqrt(length(set2)),{'color','r'});
xlim(xl);
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

f3 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(stronger_E,:)),std(all_gsac_Nsmodinfo(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinfo(stronger_I,:)),std(all_gsac_Nsmodinfo(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative info');

f4 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_postoffs(stronger_E,:)),std(all_postoffs(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_postoffs(stronger_I,:)),std(all_postoffs(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[0 0],'color','k');
xlabel('Time (s)');
ylabel('Offset');

f5 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_spost_offset(stronger_E,:)),std(all_gsac_spost_offset(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_spost_offset(stronger_I,:)),std(all_gsac_spost_offset(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[0 0],'color','k');
xlabel('Time (s)');
ylabel('Offset');

f6 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsac_spost_gain(stronger_E,:)),std(all_gsac_spost_gain(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsac_spost_gain(stronger_I,:)),std(all_gsac_spost_gain(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
xlim(xl);
line(xl,[0 0],'color','k');
xlabel('Time (s)');
ylabel('Offset');

% f5 = figure();
% hold on
% h1=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinforate(stronger_E,:)),std(all_gsac_Nsmodinforate(stronger_E,:))/sqrt(length(stronger_E)),{'color','b'});
% h2=shadedErrorBar(slags*dt,mean(all_gsac_Nsmodinforate(stronger_I,:)),std(all_gsac_Nsmodinforate(stronger_I,:))/sqrt(length(stronger_I)),{'color','r'});
% xlim(xl);
% line(xl,[1 1],'color','k');
% xlabel('Time (s)');
% ylabel('Relative info');

fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_ESdep_gains.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);


% figufy(f2);
% fname = [fig_dir 'Gsac_ESdep_rates.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
% figufy(f3);
% fname = [fig_dir 'Gsac_ESdep_info.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);
