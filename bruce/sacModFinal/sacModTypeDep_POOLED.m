close all
clear all
clc
fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';
base_tname = 'sac_trig_avg_data';
base_sname = 'sacStimTypeDep';

all_SU_data = [];
all_SU_NPdata = [];
all_SU_tdata = [];
all_SU_timedata = [];

%% LOAD JBE
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
n_probes = 96;
% ori_list = [0 90; 0 90; 0 90; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan];
ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
rmfield_list = {};

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
        sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    tavg_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            sname = strcat(tavg_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            temp = load(sname);
            tname = strcat(sac_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            
            %SU numbers for trig avg data
            tavg_SU_numbers = [temp.sua_data(:).SU_numbers];
            
            %find units where we have computed sacStimMod data
            ucells = arrayfun(@(x) length(x.gsac_GR_ovavg_rate),sacStimProc) > 0;
            sua_data = sacStimProc(ucells);
            
            if ~isfield(sua_data,'simsac_avg_rate')
                [sua_data.simsac_avg_rate] = deal(nan);
                [sua_data.simsac_post_LLinfo] = deal(nan);
                [sua_data.simsac_post_gain] = deal(nan);
                 [sua_data.simsac_post_mod] = deal(nan);
                [sua_data.simsac_post_modinfo] = deal(nan);
                [sua_data.simsac_post_offset] = deal(nan);
                [sua_data.simsac_post_ov_LLinfo] = deal(nan);
                [sua_data.simsac_post_ov_modinfo] = deal(nan);
            end
            
            SM_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,sua_data);
            SM_SU_xvLLimp = arrayfun(@(x) x.ModData.rectGQM.xvLLimp,sua_data);
            
            %find which SUs we have sacMod data for
            [lia,locb] = ismember(tavg_SU_numbers,SM_SU_numbers);
            lia_inds = find(lia);
            base_xvLLimps = ones(size(tavg_SU_numbers))*-Inf;
            base_xvLLimps(lia) = SM_SU_xvLLimp(locb(lia));
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
Expt_list = {'M266','M270','M275','M277','M281','M287','M294'};%NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
n_probes = 24;
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
rmfield_list = {};

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    osac_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod/'];
        sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
%     sac_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod/'];
    tavg_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            sname = strcat(tavg_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            temp = load(sname);
            tname = strcat(sac_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            
            tavg_SU_numbers = [temp.sua_data(:).SU_numbers];
            
            ucells = arrayfun(@(x) length(x.gsac_GR_ovavg_rate),sacStimProc) > 0;
            sua_data = sacStimProc(ucells);
 
                        if ~isfield(sua_data,'simsac_avg_rate')
                [sua_data.simsac_avg_rate] = deal(nan);
                [sua_data.simsac_post_LLinfo] = deal(nan);
                [sua_data.simsac_post_gain] = deal(nan);
                 [sua_data.simsac_post_mod] = deal(nan);
                [sua_data.simsac_post_modinfo] = deal(nan);
                [sua_data.simsac_post_offset] = deal(nan);
                [sua_data.simsac_post_ov_LLinfo] = deal(nan);
                [sua_data.simsac_post_ov_modinfo] = deal(nan);
            end

            SM_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,sua_data);
            SM_SU_xvLLimp = arrayfun(@(x) x.ModData.rectGQM.xvLLimp,sua_data);
            
            [lia,locb] = ismember(tavg_SU_numbers,SM_SU_numbers);
            lia_inds = find(lia);
            base_xvLLimps = ones(size(tavg_SU_numbers))*-Inf;
            base_xvLLimps(lia) = SM_SU_xvLLimp(locb(lia));
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
avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_SU_data);
tot_spikes = arrayfun(@(x) x.ModData.unit_data.tot_spikes,all_SU_data);
rec_dur = arrayfun(@(x) x.ModData.unit_data.N_used_samps,all_SU_data)*dt/60;
expt_nums = [all_SU_data(:).expt_num];
expt_oris = [all_SU_data(:).bar_ori];
xvLLimps = [all_SU_data(:).xvLLimp];

clust_iso_dist = arrayfun(@(x) x.ModData.unit_data.SU_isodist,all_SU_data);
clust_Lratio = arrayfun(@(x) x.ModData.unit_data.SU_Lratio,all_SU_data);
clust_refract = arrayfun(@(x) x.ModData.unit_data.SU_refract,all_SU_data);
clust_dprime = arrayfun(@(x) x.ModData.unit_data.SU_dprime,all_SU_data);
rate_stability_cv = arrayfun(@(x) x.ModData.unit_data.rate_stability_cv,all_SU_data);
dprime_stability_cv = arrayfun(@(x) x.ModData.unit_data.dprime_stability_cv,all_SU_data);

min_rate = 5; %in Hz (5)
min_Nsacs = 500; % (1000)
min_xvLLimp = 0.05; %(0.05);
jbe_SUs = find(strcmp('jbe',{all_SU_data(:).animal}));
lem_SUs = find(strcmp('lem',{all_SU_data(:).animal}));

lem_fov_expt_nums = [266 275];
lem_parafov_expt_nums = [270 277 281 287 289 294 229 297];
fov_SUs = find([all_SU_data(:).expt_num] < 200 | ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
parafov_SUs = find([all_SU_data(:).expt_num] > 200 & ~ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
lem_fov_SUs = intersect(fov_SUs,lem_SUs);
lem_parafov_SUs = intersect(parafov_SUs,lem_SUs);

N_simsacs = [all_SU_data(:).N_simsacs];
N_gsacs = [all_SU_data(:).N_gsacs];
N_gsacs_GR = [all_SU_data(:).N_gsacs_GR];
N_gsacs_IM = [all_SU_data(:).N_gsacs_IM];

N_msacs = [all_SU_data(:).N_msacs];
N_smMsacs = [all_SU_data(:).N_smMsacs];
N_bigMsacs = [all_SU_data(:).N_bigMsacs];

% use_gsac_SUs = find(avg_rates' >= min_rate & N_gsacs >= min_Nsacs & xvLLimps > min_xvLLimp);
use_gsac_SUs = find(avg_rates' >= min_rate & N_gsacs >= min_Nsacs & xvLLimps > min_xvLLimp);
use_simsac_SUs = find(avg_rates' >= min_rate & N_simsacs >= min_Nsacs & xvLLimps > min_xvLLimp);
use_msac_SUs = find(avg_rates' >= min_rate & N_msacs >= min_Nsacs & xvLLimps > min_xvLLimp);
use_jbe_SUs = intersect(use_gsac_SUs,jbe_SUs);
use_lem_SUs = intersect(use_gsac_SUs,lem_SUs);


%%
use_both = use_gsac_SUs(N_gsacs_GR(use_gsac_SUs) >= min_Nsacs & N_gsacs_IM(use_gsac_SUs) >= min_Nsacs);

all_gsacGR_gain = reshape([all_SU_data(use_both).gsacGR_post_gain],[],length(use_both))';
all_gsacIM_gain = reshape([all_SU_data(use_both).gsacIM_post_gain],[],length(use_both))';
all_gsacGR_offset = reshape([all_SU_data(use_both).gsacGR_post_offset],[],length(use_both))';
all_gsacIM_offset = reshape([all_SU_data(use_both).gsacIM_post_offset],[],length(use_both))';

all_gsacGR_rate = reshape([all_SU_data(use_both).gsacGR_avg_rate],[],length(use_both))';
all_gsacIM_rate = reshape([all_SU_data(use_both).gsacIM_avg_rate],[],length(use_both))';
all_gsacGR_rate = bsxfun(@rdivide,all_gsacGR_rate,[all_SU_data(use_both).gsac_GR_ovavg_rate]');    
all_gsacIM_rate = bsxfun(@rdivide,all_gsacIM_rate,[all_SU_data(use_both).gsac_IM_ovavg_rate]');    

all_gsacGR_info = reshape([all_SU_data(use_both).gsacGR_post_modinfo],[],length(use_both))';
all_gsacIM_info = reshape([all_SU_data(use_both).gsacIM_post_modinfo],[],length(use_both))';
all_gsacGR_info = bsxfun(@rdivide,all_gsacGR_info,[all_SU_data(use_both).gsacGR_post_ov_modinfo]');    
all_gsacIM_info = bsxfun(@rdivide,all_gsacIM_info,[all_SU_data(use_both).gsacIM_post_ov_modinfo]');    

search_range = [0 0.2];
[gsacGR_gain_Sfact,gsacGR_gain_inhtime] = get_tavg_peaks(-(all_gsacGR_gain),slags*dt,search_range);
[gsacIM_gain_Sfact,gsacIM_gain_inhtime] = get_tavg_peaks(-(all_gsacIM_gain),slags*dt,search_range);

[gsacGR_offset_Sfact] = get_tavg_peaks(-(all_gsacGR_offset),slags*dt,search_range);
[gsacIM_offset_Sfact] = get_tavg_peaks(-(all_gsacIM_offset),slags*dt,search_range);

[gsacGR_Sfact] = get_tavg_peaks(-(all_gsacGR_rate-1),slags*dt,search_range);
[gsacIM_Sfact] = get_tavg_peaks(-(all_gsacIM_rate-1),slags*dt,search_range);

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsacGR_gain)+1,std(all_gsacGR_gain)/sqrt(length(use_both)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsacIM_gain)+1,std(all_gsacIM_gain)/sqrt(length(use_both)),{'color','b'});
xlabel('Time (s)');
ylabel('Gain');

f2 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsacGR_offset),std(all_gsacGR_offset)/sqrt(length(use_both)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsacIM_offset),std(all_gsacIM_offset)/sqrt(length(use_both)),{'color','b'});
xlabel('Time (s)');
ylabel('Gain');

f3 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_gsacGR_info),std(all_gsacGR_info)/sqrt(length(use_both)),{'color','r'});
h2=shadedErrorBar(slags*dt,mean(all_gsacIM_info),std(all_gsacIM_info)/sqrt(length(use_both)),{'color','b'});
xlabel('Time (s)');
ylabel('Info');

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'sacGain_GR_IM_compare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'sacOffset_GR_IM_compare.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%%
use_both = intersect(use_simsac_SUs,use_gsac_SUs);
use_both = use_both(N_gsacs_IM(use_both) >= min_Nsacs);

all_simsac_gain = reshape([all_SU_data(use_both).simsac_post_gain],[],length(use_both))';
all_gsacIM_gain = reshape([all_SU_data(use_both).gsacIM_post_gain],[],length(use_both))';

all_simsac_rate = reshape([all_SU_data(use_both).simsac_avg_rate],[],length(use_both))';
all_gsacIM_rate = reshape([all_SU_data(use_both).gsacIM_avg_rate],[],length(use_both))';
all_simsac_rate = bsxfun(@rdivide,all_simsac_rate,[all_SU_data(use_both).simsac_ovavg_rate]');    
all_gsacIM_rate = bsxfun(@rdivide,all_gsacIM_rate,[all_SU_data(use_both).gsac_IM_ovavg_rate]');    

% all_simsac_info = reshape([all_SU_data(use_both).simsac_post_modinfo],[],length(use_both))';
% all_gsacIM_info = reshape([all_SU_data(use_both).gsacIM_post_modinfo],[],length(use_both))';
% all_simsac_info = bsxfun(@rdivide,all_simsac_info,[all_SU_data(use_both).simsac_post_ov_modinfo]');    
% all_gsacIM_info = bsxfun(@rdivide,all_gsacIM_info,[all_SU_data(use_both).gsacIM_post_ov_modinfo]');    
all_simsac_info = reshape([all_SU_data(use_both).simsac_post_LLinfo],[],length(use_both))';
all_gsacIM_info = reshape([all_SU_data(use_both).gsacIM_post_LLinfo],[],length(use_both))';
all_simsac_info = bsxfun(@rdivide,all_simsac_info,[all_SU_data(use_both).simsac_post_ov_LLinfo]');    
all_gsacIM_info = bsxfun(@rdivide,all_gsacIM_info,[all_SU_data(use_both).gsacIM_post_ov_LLinfo]');    

all_simsac_off = reshape([all_SU_data(use_both).simsac_post_offset],[],length(use_both))';
all_gsacIM_off = reshape([all_SU_data(use_both).gsacIM_post_offset],[],length(use_both))';

search_range = [0 0.2];
[simsac_gain_Sfact] = get_tavg_peaks(-(all_simsac_gain),slags*dt,search_range);
[gsacIM_gain_Sfact] = get_tavg_peaks(-(all_gsacIM_gain),slags*dt,search_range);

[simsac_Sfact] = get_tavg_peaks(-(all_simsac_rate-1),slags*dt,search_range);
[gsacIM_Sfact] = get_tavg_peaks(-(all_gsacIM_rate-1),slags*dt,search_range);

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_simsac_gain),std(all_simsac_gain)/sqrt(length(use_both)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsacIM_gain),std(all_gsacIM_gain)/sqrt(length(use_both)),{'color','r'});


f2 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_simsac_off),std(all_simsac_off)/sqrt(length(use_both)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsacIM_off),std(all_gsacIM_off)/sqrt(length(use_both)),{'color','r'});

f3 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_simsac_info),std(all_simsac_info)/sqrt(length(use_both)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_gsacIM_info),std(all_gsacIM_info)/sqrt(length(use_both)),{'color','r'});

% f4 = figure();
% hold on
% h1=shadedErrorBar(slags*dt,mean(all_simsac_rate),std(all_simsac_rate)/sqrt(length(use_both)),{'color','k'});
% h2=shadedErrorBar(slags*dt,mean(all_gsacIM_rate),std(all_gsacIM_rate)/sqrt(length(use_both)),{'color','r'});

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'sacGain_simreal_compare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'sacOffset_simreal_compare.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
%%
close all
use_both = use_msac_SUs(N_smMsacs(use_msac_SUs) >= min_Nsacs & N_bigMsacs(use_msac_SUs) >= min_Nsacs);
% use_both = use_both(ismember(use_both,use_lem_SUs));
% use_both = use_both(ismember(use_both,use_jbe_SUs));

all_smMsac_gain = reshape([all_SU_data(use_both).smMsac_post_gain],[],length(use_both))';
all_bigMsac_gain = reshape([all_SU_data(use_both).bigMsac_post_gain],[],length(use_both))';

all_smMsac_rate = reshape([all_SU_data(use_both).smMsac_avg_rate],[],length(use_both))';
all_bigMsac_rate = reshape([all_SU_data(use_both).bigMsac_avg_rate],[],length(use_both))';
all_smMsac_rate = bsxfun(@rdivide,all_smMsac_rate,[all_SU_data(use_both).smMsac_ovavg_rate]');    
all_bigMsac_rate = bsxfun(@rdivide,all_bigMsac_rate,[all_SU_data(use_both).bigMsac_ovavg_rate]');    


all_smMsac_info = reshape([all_SU_data(use_both).smMsac_post_modinfo],[],length(use_both))';
all_bigMsac_info = reshape([all_SU_data(use_both).bigMsac_post_modinfo],[],length(use_both))';
all_smMsac_info = bsxfun(@rdivide,all_smMsac_info,[all_SU_data(use_both).smMsac_post_ov_modinfo]');    
all_bigMsac_info = bsxfun(@rdivide,all_bigMsac_info,[all_SU_data(use_both).bigMsac_post_ov_modinfo]');    


all_smMsac_off = reshape([all_SU_data(use_both).smMsac_post_offset],[],length(use_both))';
all_bigMsac_off = reshape([all_SU_data(use_both).bigMsac_post_offset],[],length(use_both))';

search_range = [0 0.2];
[smMsac_gain_Sfact] = get_tavg_peaks(-(all_smMsac_gain),slags*dt,search_range);
[bigMsac_gain_Sfact] = get_tavg_peaks(-(all_bigMsac_gain),slags*dt,search_range);
[smMsac_info_Sfact] = get_tavg_peaks(-(all_smMsac_info-1),slags*dt,search_range);
[bigMsac_info_Sfact] = get_tavg_peaks(-(all_bigMsac_info-1),slags*dt,search_range);

[smMsac_Sfact] = get_tavg_peaks(-(all_smMsac_rate-1),slags*dt,search_range);
[bigMsac_Sfact] = get_tavg_peaks(-(all_bigMsac_rate-1),slags*dt,search_range);
[smMsac_Efact] = get_tavg_peaks((all_smMsac_rate-1),slags*dt,search_range);
[bigMsac_Efact] = get_tavg_peaks((all_bigMsac_rate-1),slags*dt,search_range);

f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_smMsac_gain),std(all_smMsac_gain)/sqrt(length(use_both)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_bigMsac_gain),std(all_bigMsac_gain)/sqrt(length(use_both)),{'color','r'});


f1 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_smMsac_info),std(all_smMsac_info)/sqrt(length(use_both)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_bigMsac_info),std(all_bigMsac_info)/sqrt(length(use_both)),{'color','r'});

f2 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_smMsac_off),std(all_smMsac_off)/sqrt(length(use_both)),{'color','b'});
h2=shadedErrorBar(slags*dt,mean(all_bigMsac_off),std(all_bigMsac_off)/sqrt(length(use_both)),{'color','r'});

f3 = figure();
hold on
h1=shadedErrorBar(slags*dt,mean(all_smMsac_rate),std(all_smMsac_rate)/sqrt(length(use_both)),{'color','k'});
h2=shadedErrorBar(slags*dt,mean(all_bigMsac_rate),std(all_bigMsac_rate)/sqrt(length(use_both)),{'color','r'});

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'sacGain_simreal_compare.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'sacOffset_simreal_compare.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);
% 
rf_sigma = arrayfun(@(x) x.ModData.tune_props.RF_sigma,all_SU_data);
rf_ecc = arrayfun(@(x) x.ModData.tune_props.RF_ecc,all_SU_data);
