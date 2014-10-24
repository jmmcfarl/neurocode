%%
% close all
clear all
fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';

all_SU_data = [];
all_SU_NPdata = [];
all_MU_data = [];
all_XC_data = [];

base_sname = 'corrected_models2';
% base_tname = 'sac_glm_data';
base_tname = 'sac_glm_data_withbursts';

%% LOAD JBE

Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
n_probes = 96;
ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
rmfield_list = {};

for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
    mod_dir = ['~/Analysis/bruce/' Expt_name '/models/'];
    
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            sname = strcat(mod_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            load(sname);
            tname = strcat(sac_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            
            xcorr_data.animal = 'jbe';
            all_XC_data = cat(1,all_XC_data,xcorr_data);
            
            ucells = arrayfun(@(x) length(x.unit_data),ModData) > 0;
            ModData = ModData(ucells);
            Mod_SU_numbers = arrayfun(@(x) x.unit_data.SU_number,ModData);
            Mod_SU_xvLLimp = arrayfun(@(x) x.rectGQM.xvLLimp,ModData);
            
            tavg_SU_numbers = [sua_data(:).SU_numbers];
            [lia,locb] = ismember(tavg_SU_numbers,Mod_SU_numbers);
            
            base_xvLLimps = ones(size(tavg_SU_numbers))*-Inf;
            base_xvLLimps(lia) = Mod_SU_xvLLimp(locb(lia));
            base_Mods = ModData;
            base_Mods(lia) = ModData(locb(lia));
            for jj = 1:length(sua_data); 
                sua_data(jj).xvLLimp = base_xvLLimps(jj); 
                if jj <= length(base_Mods)
                sua_data(jj).ModData = base_Mods(jj); 
                end
            end;
            
            [sua_data.expt_num] = deal(Expt_num);
            [sua_data.bar_ori] = deal(ori_list(ee,ii));
            [sua_data.animal] = deal('jbe');
            
            ori_SU_nums(ii,:) = tavg_SU_numbers;
            ori_xvLLimps(ii,:) = base_xvLLimps;
            ori_sua_data{ii} = sua_data;
            clear ModData
            
            [mua_data.expt_num] = deal(ones(1,n_probes)*Expt_num);
            [mua_data.bar_ori] = deal(ones(1,n_probes)*ori_list(ee,ii));
            [mua_data.animal] = deal(repmat({'jbe'},1,n_probes));
            all_MU_data = cat(1,all_MU_data,mua_data);
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
        end
    end
    
    clear ori_SU_nums ori_xvLLimps ori_sua_data
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
            sname = strcat(mod_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            load(sname);
            tname = strcat(sac_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            
            xcorr_data.animal = 'lem';
            all_XC_data = cat(1,all_XC_data,xcorr_data);

            ucells = arrayfun(@(x) length(x.unit_data),ModData) > 0;
            ModData = ModData(ucells);
            Mod_SU_numbers = arrayfun(@(x) x.unit_data.SU_number,ModData);
            Mod_SU_xvLLimp = arrayfun(@(x) x.rectGQM.xvLLimp,ModData);
            
            tavg_SU_numbers = [sua_data(:).SU_numbers];
            [lia,locb] = ismember(tavg_SU_numbers,Mod_SU_numbers);
            
            base_xvLLimps = ones(size(tavg_SU_numbers))*-Inf;
            base_xvLLimps(lia) = Mod_SU_xvLLimp(locb(lia));
            base_Mods = ModData;
            base_Mods(lia) = ModData(locb(lia));
            for jj = 1:length(sua_data); 
                sua_data(jj).xvLLimp = base_xvLLimps(jj); 
                if jj <= length(base_Mods)
                sua_data(jj).ModData = base_Mods(jj); 
                end
            end;
            
            [sua_data.expt_num] = deal(Expt_num);
            [sua_data.bar_ori] = deal(ori_list(ee,ii));
            [sua_data.animal] = deal('lem');
            
            ori_SU_nums(ii,:) = tavg_SU_numbers;
            ori_xvLLimps(ii,:) = base_xvLLimps;
            ori_sua_data{ii} = sua_data;
            clear ModData
            
            [mua_data.expt_num] = deal(ones(1,n_probes)*Expt_num);
            [mua_data.bar_ori] = deal(ones(1,n_probes)*ori_list(ee,ii));
            [mua_data.animal] = deal(repmat({'lem'},1,n_probes));
            all_MU_data = cat(1,all_MU_data,mua_data);
        end
    end
    %if there was only one ori, set the xvLLimps to -Inf for a placeholder
    if size(ori_xvLLimps,1) == 1
        ori_xvLLimps = [ori_xvLLimps; ones(size(ori_xvLLimps))*-Inf];
    end
    [mvals,mlocs] = max(ori_xvLLimps,[],1);    
    nplocs = mod(mlocs,2)+1; %these are indices for the NP ori
    for ii = 1:length(mvals)
            if ori_xvLLimps(nplocs(ii),ii) > 0 %if the NP ori had usable data
               all_SU_NPdata = cat(1,all_SU_NPdata,ori_sua_data{nplocs(ii)}(ii)); 
               has_NP = true;
            else
                has_NP = false;
            end
            cur_struct = ori_sua_data{mlocs(ii)}(ii);
            cur_struct.has_NP = has_NP;
            all_SU_data = cat(1,all_SU_data,cur_struct); 
    end
    
    clear ori_SU_nums ori_xvLLimps ori_sua_data
end

%%
tlags = trig_avg_params.lags;
dt = trig_avg_params.dt;
tlags = tlags*dt;

tot_spikes = [all_SU_data(:).tot_nspikes];
bad = find(tot_spikes == 0);
all_SU_data(bad) = [];
%% SELECT USABLE CELLS
avg_rates = [all_SU_data(:).avg_rates]/dt;
tot_spikes = [all_SU_data(:).tot_nspikes];
rec_dur = [all_SU_data(:).N_used_samps]*dt/60; %in min
expt_nums = [all_SU_data(:).expt_num];
expt_oris = [all_SU_data(:).bar_ori];

clust_iso_dist = [all_SU_data(:).SU_isodist];
clust_Lratio = [all_SU_data(:).SU_Lratio];
clust_refract = [all_SU_data(:).SU_refract];
clust_dprime = [all_SU_data(:).SU_dprime];
rate_stability_cv = [all_SU_data(:).rate_stability_cv];
dprime_stability_cv = [all_SU_data(:).dprime_stability_cv];

sm_sigma = 0;
min_rate = 5; %in Hz (5)
min_Nsacs = 250; % (100)
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
N_gray_gsacs = [all_SU_data(:).N_gsacs_gray];
N_gray_msacs = [all_SU_data(:).N_msacs_gray];
N_im_gsacs = [all_SU_data(:).N_gsacs_im];
N_simsacs = [all_SU_data(:).N_simsacs];
N_simmsacs = [all_SU_data(:).N_simmsacs]; %simulated microsacs
N_blanks = [all_SU_data(:).N_blanks];

%%
xc_lags = all_XC_data(1).lags;
xc_dt = median(diff(xc_lags));
xc_sm = 0.02/xc_dt;
msac_gsac_xcorr = reshape([all_XC_data(:).msac_gsac_xcorr],[],length(all_XC_data))';
msac_simsac_xcorr = reshape([all_XC_data(:).msac_simsac_xcorr],[],length(all_XC_data))';
msac_msac_xcorr = reshape([all_XC_data(:).msac_msac_xcorr],[],length(all_XC_data))';

for ii = 1:length(all_XC_data)
    msac_gsac_xcorr(ii,:) = jmm_smooth_1d_cor(msac_gsac_xcorr(ii,:),xc_sm);
    msac_simsac_xcorr(ii,:) = jmm_smooth_1d_cor(msac_simsac_xcorr(ii,:),xc_sm);
    msac_msac_xcorr(ii,:) = jmm_smooth_1d_cor(msac_msac_xcorr(ii,:),xc_sm);
end

f1 = figure();
hold on
h1=shadedErrorBar(xc_lags,nanmean(msac_gsac_xcorr)/xc_dt,nanstd(msac_gsac_xcorr)/sqrt(length(all_XC_data))/xc_dt,{'color','r'});
h2=shadedErrorBar(xc_lags,nanmean(msac_simsac_xcorr)/xc_dt,nanstd(msac_simsac_xcorr)/sqrt(length(all_XC_data))/xc_dt,{'color','b'});

% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Msac_xcorrs.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f1);



% jj = find(strcmp({all_XC_data(:).animal},'jbe'));
% ll = find(strcmp({all_XC_data(:).animal},'lem'));
% f2 = figure();
% hold on
% h3=shadedErrorBar(xc_lags,nanmean(msac_msac_xcorr(jj,:))/xc_dt,nanstd(msac_msac_xcorr(jj,:))/sqrt(length(jj))/xc_dt,{'color','r'});
% h3=shadedErrorBar(xc_lags,nanmean(msac_msac_xcorr(ll,:))/xc_dt,nanstd(msac_msac_xcorr(ll,:))/sqrt(length(ll))/xc_dt,{'color','k'});

% fig_width = 3.5; rel_height = 0.8;
% figufy(f2);
% fname = [fig_dir 'Msac_acorr_compare.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% % close(f2);


%% GSACS GLM VS TRIGAVG

gsac_used_SUs = find(N_gsacs >= min_Nsacs & avg_rates >= min_rate);

% all_gsac = reshape([all_SU_data(:).gsac_avg],[],length(all_SU_data))';
all_gsac = reshape([all_SU_data(gsac_used_SUs).tavg_gsac_rate],[],length(gsac_used_SUs))';
all_gsac = bsxfun(@rdivide,all_gsac,[all_SU_data(gsac_used_SUs).glm_gsac_avgrate]');

% if sm_sigma > 0
%     for ii = 1:size(all_gsac,1)
%         all_gsac(ii,:) = jmm_smooth_1d_cor(all_gsac(ii,:),sm_sigma);
%     end
% end
reg_gsac = reshape([all_SU_data(gsac_used_SUs).glm_gsac_rate],[],length(gsac_used_SUs))';
nreg_gsac = bsxfun(@rdivide,reg_gsac,[all_SU_data(gsac_used_SUs).glm_gsac_avgrate]');
reg_lags = all_SU_data(1).glm_lags*all_SU_data(1).glm_dt;

% close all
xl = [-0.15 0.4];

f1 = figure(); hold on
h1=shadedErrorBar(reg_lags,nanmean(all_gsac),nanstd(all_gsac)/sqrt(length(gsac_used_SUs)),{'color','r'});
h2=shadedErrorBar(reg_lags,nanmean(reg_gsac),nanstd(reg_gsac)/sqrt(length(gsac_used_SUs)),{'color','b'});
xlim(xl);
ylim([0.7 1.25]);
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'Gsac_glm_rates.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% MSACS GLM VS TRIGAVG
msac_used_SUs = find(N_msacs >= min_Nsacs & avg_rates >= min_rate);

all_msac = reshape([all_SU_data(msac_used_SUs).tavg_msac_rate],[],length(msac_used_SUs))';
all_msac = bsxfun(@rdivide,all_msac,[all_SU_data(msac_used_SUs).glm_msac_avgrate]');

% if sm_sigma > 0
%     for ii = 1:size(all_msac,1)
%         all_msac(ii,:) = jmm_smooth_1d_cor(all_msac(ii,:),sm_sigma);
%     end
% end

reg_msac = reshape([all_SU_data(msac_used_SUs).glm_msac_rate],[],length(msac_used_SUs))';
nreg_msac = bsxfun(@rdivide,reg_msac,[all_SU_data(msac_used_SUs).glm_msac_avgrate]');
reg_lags = all_SU_data(1).glm_lags*all_SU_data(1).glm_dt;

% close all
xl = [-0.15 0.4];

f1 = figure(); hold on
h1=shadedErrorBar(reg_lags,nanmean(all_msac),nanstd(all_msac)/sqrt(length(msac_used_SUs)),{'color','r'});
h2=shadedErrorBar(reg_lags,nanmean(reg_msac),nanstd(reg_msac)/sqrt(length(msac_used_SUs)),{'color','b'});
xlim(xl);
ylim([0.7 1.25]);
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'Msac_glm_rates.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% SIMSACS GLM VS TRIGAVG
simsac_used_SUs = find(N_simsacs >= min_Nsacs & avg_rates >= min_rate);

all_simsac = reshape([all_SU_data(simsac_used_SUs).tavg_simsac_rate],[],length(simsac_used_SUs))';
all_simsac = bsxfun(@rdivide,all_simsac,[all_SU_data(simsac_used_SUs).glm_simsac_avgrate]');

% if sm_sigma > 0
%     for ii = 1:size(all_simsac,1)
%         all_simsac(ii,:) = jmm_smooth_1d_cor(all_simsac(ii,:),sm_sigma);
%     end
% end

reg_simsac = reshape([all_SU_data(simsac_used_SUs).glm_simsac_rate],[],length(simsac_used_SUs))';
nreg_simsac = bsxfun(@rdivide,reg_simsac,[all_SU_data(simsac_used_SUs).glm_simsac_avgrate]');
reg_lags = all_SU_data(1).glm_lags*all_SU_data(1).glm_dt;

% close all
xl = [-0.15 0.4];

f1 = figure(); hold on
h1=shadedErrorBar(reg_lags,nanmean(all_simsac),nanstd(all_simsac)/sqrt(length(simsac_used_SUs)),{'color','r'});
h2=shadedErrorBar(reg_lags,nanmean(reg_simsac),nanstd(reg_simsac)/sqrt(length(simsac_used_SUs)),{'color','b'});
xlim(xl);
ylim([0.85 1.15]);
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');


fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'Simsac_glm_rates.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% MSACS GLM VS TRIGAVG
msac_used_SUs = find(N_msacs >= min_Nsacs & avg_rates >= min_rate);
msac_used_SUs = msac_used_SUs(ismember(msac_used_SUs,jbe_SUs));

all_msac = reshape([all_SU_data(msac_used_SUs).tavg_msac_rate],[],length(msac_used_SUs))';
all_msac = bsxfun(@rdivide,all_msac,[all_SU_data(msac_used_SUs).glm_msac_avgrate]');

% if sm_sigma > 0
%     for ii = 1:size(all_msac,1)
%         all_msac(ii,:) = jmm_smooth_1d_cor(all_msac(ii,:),sm_sigma);
%     end
% end

reg_msac = reshape([all_SU_data(msac_used_SUs).glm_msac_rate],[],length(msac_used_SUs))';
nreg_msac = bsxfun(@rdivide,reg_msac,[all_SU_data(msac_used_SUs).glm_msac_avgrate]');
reg_lags = all_SU_data(1).glm_lags*all_SU_data(1).glm_dt;

% close all
xl = [-0.15 0.4];

f1 = figure(); hold on
% h2=shadedErrorBar(reg_lags,nanmean(nreg_msac),nanstd(nreg_msac)/sqrt(length(msac_used_SUs)),{'color','b'});
h2=shadedErrorBar(reg_lags,nanmean(nreg_msac),nanstd(nreg_msac)/sqrt(length(msac_used_SUs)),{'color','r'});
xlim(xl);
ylim([0.7 1.25]);
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'Msac_GLM_burstcompare.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);
