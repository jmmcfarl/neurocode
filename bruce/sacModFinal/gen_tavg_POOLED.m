%%
% close all
clear all
fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';

all_data = [];

base_name = 'gen_trig_avg_data';
include_bursts = 0;
if include_bursts
    base_name = [base_name '_withburst'];
end

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
            tname = strcat(sac_dir,base_name,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            
            gen_data.animal = 'jbe';
            gen_data.ori = ori_list(ee,ii);
            gen_data.exptnum = Expt_num;
            all_data = cat(1,all_data,gen_data);
        
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
    mod_dir = ['~/Analysis/bruce/' Expt_name '/models/'];
    
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            name = strcat(sac_dir,base_name,sprintf('_ori%d',ori_list(ee,ii)));
            load(name);
            
            gen_data.animal = 'lem';
            gen_data.ori = ori_list(ee,ii);
            gen_data.exptnum = Expt_num;
            all_data = cat(1,all_data,gen_data);

        end
    end
end

%%
N_msacs = [all_data(:).N_msacs];
Tot_time = [all_data(:).Tot_time];
dt = trig_avg_params.dt;

msac_rates = N_msacs./Tot_time/dt;
expt_nums = [all_data(:).exptnum];
unique_expts = unique(expt_nums);

expt_msac_rates = nan(length(unique_expts),1);
expt_animal = cell(length(unique_expts),1);
for ii = 1:length(unique_expts)
    eset = find(expt_nums == unique_expts(ii));
    expt_msac_rates(ii) = mean(msac_rates(eset));
    expt_animal{ii} = all_data(eset(1)).animal;
end

%%
all_gsac_tavg = [all_data(:).gsac_rawtavg_eyespeed]';
all_msac_tavg = [all_data(:).msac_rawtavg_eyespeed]';

dur_bin_cents = all_data(1).dur_bin_cents;
dur_dx = median(diff(dur_bin_cents));
dur_sm = 0.005/dur_dx;
all_gsac_durdist = reshape([all_data(:).gsac_dur_dist],length(dur_bin_cents),[])';
all_msac_durdist = reshape([all_data(:).msac_dur_dist],length(dur_bin_cents),[])';

all_gsac_durdist = bsxfun(@rdivide,all_gsac_durdist,sum(all_gsac_durdist,2));
all_msac_durdist = bsxfun(@rdivide,all_msac_durdist,sum(all_msac_durdist,2));

eye_ax = all_data(1).raw_eye_lags;
expt_gsac_eyespeed = nan(length(unique_expts),length(eye_ax));
expt_msac_eyespeed = nan(length(unique_expts),length(eye_ax));
expt_gsac_durdist = nan(length(unique_expts),length(dur_bin_cents));
expt_msac_durdist = nan(length(unique_expts),length(dur_bin_cents));
for ii = 1:length(unique_expts)
    eset = find(expt_nums == unique_expts(ii));
    expt_gsac_eyespeed(ii,:) = mean(all_gsac_tavg(eset,:),1);
    expt_msac_eyespeed(ii,:) = mean(all_msac_tavg(eset,:),1);
    expt_gsac_durdist(ii,:) = mean(all_gsac_durdist(eset,:),1);
    expt_msac_durdist(ii,:) = mean(all_msac_durdist(eset,:),1);
    if dur_sm > 0
       expt_gsac_durdist(ii,:) = jmm_smooth_1d_cor(expt_gsac_durdist(ii,:),dur_sm); 
       expt_msac_durdist(ii,:) = jmm_smooth_1d_cor(expt_msac_durdist(ii,:),dur_sm); 
    end
end

% f1= figure; hold on
% h1 = shadedErrorBar(eye_ax,mean(expt_gsac_eyespeed),std(expt_gsac_eyespeed)/sqrt(length(unique_expts)));
% h2 = shadedErrorBar(eye_ax,mean(expt_msac_eyespeed),std(expt_msac_eyespeed)/sqrt(length(unique_expts)),{'color','r'});

f2 = figure; hold on
h1 = shadedErrorBar(dur_bin_cents,mean(expt_gsac_durdist),std(expt_gsac_durdist)/sqrt(length(unique_expts)),{'color','b'});
h2 = shadedErrorBar(dur_bin_cents,mean(expt_msac_durdist),std(expt_msac_durdist)/sqrt(length(unique_expts)),{'color','r'});
xlabel('Duration (s)');
ylabel('Relative frequency');

%PRINT PLOTS
fig_width = 3.5; rel_height = 0.8;
figufy(f2);
fname = [fig_dir 'Gsac_msac_dur_dists.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);


%%
all_gsac_tavg = [all_data(:).gsac_tavg_eyespeed]';
all_msac_tavg = [all_data(:).msac_tavg_eyespeed]';
eye_ax = all_data(1).eye_lags;
expt_gsac_eyespeed = nan(length(unique_expts),length(eye_ax));
expt_msac_eyespeed = nan(length(unique_expts),length(eye_ax));
for ii = 1:length(unique_expts)
    eset = find(expt_nums == unique_expts(ii));
    expt_gsac_eyespeed(ii,:) = mean(all_gsac_tavg(eset,:),1);
    expt_msac_eyespeed(ii,:) = mean(all_msac_tavg(eset,:),1);
end

figure; hold on
h1 = shadedErrorBar(eye_ax,mean(expt_gsac_eyespeed),std(expt_gsac_eyespeed)/sqrt(length(unique_expts)));
h2 = shadedErrorBar(eye_ax,mean(expt_msac_eyespeed),std(expt_msac_eyespeed)/sqrt(length(unique_expts)),{'color','r'});

