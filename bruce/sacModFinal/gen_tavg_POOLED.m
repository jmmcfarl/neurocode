%%
% close all
clear all
fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';

all_data = [];
base_name = 'gen_trig_avg_data';
base_ename = 'eye_accuracy';
include_bursts = 0;
if include_bursts
    base_name = [base_name '_withburst'];
    base_ename = [base_ename '_withburst'];
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
             tname = strcat(sac_dir,base_ename,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
           
            gen_data.fix_data = fix_disp_data;
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
             tname = strcat(sac_dir,base_ename,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            
            gen_data.fix_data = fix_disp_data;
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

msac_dur_avgs = [all_data(:).msac_dur_avg];
gsac_dur_avgs = [all_data(:).gsac_dur_avg];
gsac_new_dur_avgs = arrayfun(@(x) mean(x.new_gsac_durs),all_data);
msac_new_dur_avgs = arrayfun(@(x) mean(x.new_msac_durs),all_data);

msac_rates = N_msacs./Tot_time/dt;
expt_nums = [all_data(:).exptnum];
unique_expts = unique(expt_nums);

%get avg values within each expt (across multiple bar oris if applicable)
expt_msac_rates = nan(length(unique_expts),1);
expt_gsac_durs = nan(length(unique_expts),1);
expt_msac_durs = nan(length(unique_expts),1);
expt_gsac_new_durs = nan(length(unique_expts),1);
expt_msac_new_durs = nan(length(unique_expts),1);
expt_animal = cell(length(unique_expts),1);
for ii = 1:length(unique_expts)
    eset = find(expt_nums == unique_expts(ii));
    expt_msac_rates(ii) = mean(msac_rates(eset));
    expt_msac_durs(ii) = mean(msac_dur_avgs(eset));
    expt_gsac_durs(ii) = mean(gsac_dur_avgs(eset));
    expt_msac_new_durs(ii) = mean(msac_new_dur_avgs(eset));
    expt_gsac_new_durs(ii) = mean(gsac_new_dur_avgs(eset));
    expt_animal{ii} = all_data(eset(1)).animal;
end

%% ORIGINAL CALCULATION OF SAC DURATION DISTRIBUTIONS
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

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Gsac_msac_dur_dists.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% COMPUTE SAC DURATION DISTRIBUTION USING NEW DURATION DEFINITION

dur_bin_cents = all_data(1).dur_bin_cents;
dur_dx = median(diff(dur_bin_cents));
dur_bin_edges = [dur_bin_cents - dur_dx/2 (dur_bin_cents(end) + dur_dx/2)];
dur_sm = 0.0025/dur_dx;

[all_gsac_durdist,all_msac_durdist] = deal(nan(length(all_data),length(dur_bin_cents)));
for ii = 1:length(all_data)
    temp = histc(all_data(ii).new_gsac_durs,dur_bin_edges);
    all_gsac_durdist(ii,:) = temp(1:end-1)'/sum(temp);
    temp = histc(all_data(ii).new_msac_durs,dur_bin_edges);
    all_msac_durdist(ii,:) = temp(1:end-1)/sum(temp);
end

expt_gsac_durdist = nan(length(unique_expts),length(dur_bin_cents));
expt_msac_durdist = nan(length(unique_expts),length(dur_bin_cents));
for ii = 1:length(unique_expts)
    eset = find(expt_nums == unique_expts(ii));
    expt_gsac_durdist(ii,:) = mean(all_gsac_durdist(eset,:),1);
    expt_msac_durdist(ii,:) = mean(all_msac_durdist(eset,:),1);
    if dur_sm > 0
       expt_gsac_durdist(ii,:) = jmm_smooth_1d_cor(expt_gsac_durdist(ii,:),dur_sm); 
       expt_msac_durdist(ii,:) = jmm_smooth_1d_cor(expt_msac_durdist(ii,:),dur_sm); 
    end
end

f2 = figure; hold on
h1 = shadedErrorBar(dur_bin_cents,mean(expt_gsac_durdist),std(expt_gsac_durdist)/sqrt(length(unique_expts)),{'color','b'});
h2 = shadedErrorBar(dur_bin_cents,mean(expt_msac_durdist),std(expt_msac_durdist)/sqrt(length(unique_expts)),{'color','r'});
xlabel('Duration (s)');
ylabel('Relative frequency');

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Gsac_msac_dur_dists_new.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%%
% all_gsac_tavg = [all_data(:).gsac_tavg_eyespeed]';
% all_msac_tavg = [all_data(:).msac_tavg_eyespeed]';
% eye_ax = all_data(1).eye_lags;
% expt_gsac_eyespeed = nan(length(unique_expts),length(eye_ax));
% expt_msac_eyespeed = nan(length(unique_expts),length(eye_ax));
% for ii = 1:length(unique_expts)
%     eset = find(expt_nums == unique_expts(ii));
%     expt_gsac_eyespeed(ii,:) = mean(all_gsac_tavg(eset,:),1);
%     expt_msac_eyespeed(ii,:) = mean(all_msac_tavg(eset,:),1);
% end
% 
% figure; hold on
% h1 = shadedErrorBar(eye_ax,mean(expt_gsac_eyespeed),std(expt_gsac_eyespeed)/sqrt(length(unique_expts)));
% h2 = shadedErrorBar(eye_ax,mean(expt_msac_eyespeed),std(expt_msac_eyespeed)/sqrt(length(unique_expts)),{'color','r'});

%% DISPARITY ANALYSIS
et_dt = 0.01;
all_gsac_fixdisp = cell2mat(arrayfun(@(x) x.fix_data.gsac_fix_disp(:,1)',all_data,'uniformoutput',0));
all_msac_fixdisp = cell2mat(arrayfun(@(x) x.fix_data.msac_fix_disp(:,1)',all_data,'uniformoutput',0));

expt_gsac_fixdisp = nan(length(unique_expts),length(slags));
expt_msac_fixdisp = nan(length(unique_expts),length(slags));
for ii = 1:length(unique_expts)
    eset = find(expt_nums == unique_expts(ii));
    expt_gsac_fixdisp(ii,:) = mean(all_gsac_fixdisp(eset,:),1);
    expt_msac_fixdisp(ii,:) = mean(all_msac_fixdisp(eset,:),1);
end

f1 = figure();hold on
h1 = shadedErrorBar(slags*et_dt,mean(expt_gsac_fixdisp),std(expt_gsac_fixdisp)/sqrt(length(unique_expts)),{'color','k'});
h2 = shadedErrorBar(slags*et_dt,mean(expt_msac_fixdisp),std(expt_msac_fixdisp)/sqrt(length(unique_expts)),{'color','r'});
xlabel('Time (s)');
ylabel('Fixation disparity (deg)');

% %PRINT PLOTS
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Gsac_msac_fix_disp.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%
eye_tax = all_data(1).raw_eye_lags;
gsac_rawtavg_eyespeed = cell2mat(arrayfun(@(x) x.gsac_rawtavg_eyespeed',all_data,'uniformoutput',0));
msac_rawtavg_eyespeed = cell2mat(arrayfun(@(x) x.msac_rawtavg_eyespeed',all_data,'uniformoutput',0));

gsac_avg_orth_frac = arrayfun(@(x) mean(abs(x.gsac_delta_Y_frac)),all_data);
[avg_inaccurate_orth,avg_accurate_orth] = deal(nan(length(all_data),1));
for ii = 1:length(all_data)
    cur_orth_frac = abs(all_data(ii).gsac_delta_Y_frac);
    med_orth_frac = median(cur_orth_frac);
    avg_inaccurate_orth(ii) = mean(cur_orth_frac(cur_orth_frac > med_orth_frac));
    avg_accurate_orth(ii) = mean(cur_orth_frac(cur_orth_frac < med_orth_frac));
end

gsac_avg_orthspeed = bsxfun(@times,gsac_rawtavg_eyespeed,gsac_avg_orth_frac);
gsac_inac_orthspeed = bsxfun(@times,gsac_rawtavg_eyespeed,avg_inaccurate_orth);
gsac_ac_orthspeed = bsxfun(@times,gsac_rawtavg_eyespeed,avg_accurate_orth);

[expt_gsac_eyespeed,expt_msac_eyespeed,expt_gsac_avgorthspeed,expt_gsac_inaccorthspeed,expt_gsac_accorthspeed] = deal(nan(length(unique_expts),length(eye_tax)));
[expt_avg_orth,expt_inacc_orth,expt_acc_orth] = deal(nan(length(unique_expts),1));
for ii = 1:length(unique_expts)
    eset = find(expt_nums == unique_expts(ii));
    expt_msac_eyespeed(ii,:) = mean(msac_rawtavg_eyespeed(eset,:),1);
    expt_gsac_eyespeed(ii,:) = mean(gsac_rawtavg_eyespeed(eset,:),1);
    expt_gsac_avgorthspeed(ii,:) = mean(gsac_avg_orthspeed(eset,:),1);
    expt_gsac_inaccorthspeed(ii,:) = mean(gsac_inac_orthspeed(eset,:),1);
    expt_gsac_accorthspeed(ii,:) = mean(gsac_ac_orthspeed(eset,:),1);
    expt_avg_orth(ii) = mean(gsac_avg_orth_frac(eset));
    expt_inacc_orth(ii) = mean(avg_inaccurate_orth(eset));
    expt_acc_orth(ii) = mean(avg_accurate_orth(eset));
end

f1 = figure(); hold on
plot(eye_tax,mean(expt_gsac_avgorthspeed));
plot(eye_tax,mean(expt_gsac_inaccorthspeed),'r')
plot(eye_tax,mean(expt_gsac_accorthspeed),'k');

f2 = figure(); 
subplot(2,1,1);hold on
plot(eye_tax,mean(expt_gsac_eyespeed));
hold on
plot(eye_tax,mean(expt_msac_eyespeed),'r')
xlabel('Time (s)'); ylabel('Speed (deg/sec)');
xlim([-0.02 0.08]);
line([-0.02 0.08],[3 3],'color','k','linestyle','--');
line([-0.02 0.08],[5 5],'color','r','linestyle','--');
line([-0.02 0.08],[10 10],'color','b','linestyle','--');

subplot(2,1,2);hold on
plot(eye_tax,mean(expt_gsac_eyespeed));
hold on
plot(eye_tax,mean(expt_msac_eyespeed),'r')
xlabel('Time (s)'); ylabel('Speed (deg/sec)');
set(gca,'yscale','log');
xlim([-0.02 0.08]);
line([-0.02 0.08],[3 3],'color','k','linestyle','--');
line([-0.02 0.08],[5 5],'color','r','linestyle','--');
line([-0.02 0.08],[10 10],'color','b','linestyle','--');

%PRINT PLOTS
fig_width = 3.5; rel_height = 1.5;
figufy(f2);
fname = [fig_dir 'Avg_speed_thresh_comparison.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

%%
orth_trajects.tax = eye_tax;
orth_trajects.avg_orth_speed = mean(expt_gsac_avgorthspeed);
orth_trajects.avg_inac_speed = mean(expt_gsac_inaccorthspeed);
orth_trajects.avg_acc_speed = mean(expt_gsac_accorthspeed);
dname = '~/Analysis/bruce/FINsac_mod/orth_eyetrajectories';
save(dname,'orth_trajects');
