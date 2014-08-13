%%
close all
clear all
fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';

all_SU_data = [];
all_SU_NPdata = [];
all_MU_data = [];

%% LOAD JBE
base_sname = 'corrected_models';
base_tname = 'sac_trig_avg_data';
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
            
            ucells = arrayfun(@(x) length(x.unit_data),ModData) > 0;
            ModData = ModData(ucells);
            Mod_SU_numbers = arrayfun(@(x) x.unit_data.SU_number,ModData);
            Mod_SU_xvLLimp = arrayfun(@(x) x.rectGQM.xvLLimp,ModData);
            
            tavg_SU_numbers = [sua_data(:).SU_numbers];
            [lia,locb] = ismember(tavg_SU_numbers,Mod_SU_numbers);
            
            base_xvLLimps = ones(size(tavg_SU_numbers))*-Inf;
            base_xvLLimps(lia) = Mod_SU_xvLLimp(locb(lia));
            for jj = 1:length(sua_data); sua_data(jj).xvLLimp = base_xvLLimps(jj); end;
            
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
base_sname = 'corrected_models';
base_tname = 'sac_trig_avg_data';
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
            
            ucells = arrayfun(@(x) length(x.unit_data),ModData) > 0;
            ModData = ModData(ucells);
            Mod_SU_numbers = arrayfun(@(x) x.unit_data.SU_number,ModData);
            Mod_SU_xvLLimp = arrayfun(@(x) x.rectGQM.xvLLimp,ModData);
            
            tavg_SU_numbers = [sua_data(:).SU_numbers];
            [lia,locb] = ismember(tavg_SU_numbers,Mod_SU_numbers);
            
            base_xvLLimps = ones(size(tavg_SU_numbers))*-Inf;
            base_xvLLimps(lia) = Mod_SU_xvLLimp(locb(lia));
            for jj = 1:length(sua_data); sua_data(jj).xvLLimp = base_xvLLimps(jj); end;
            
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

N_gray_gsacs = [all_SU_data(:).N_gsacs_gray];
N_gray_msacs = [all_SU_data(:).N_msacs_gray];
N_im_gsacs = [all_SU_data(:).N_gsacs_im];
N_simsacs = [all_SU_data(:).N_simsacs];
N_simmsacs = [all_SU_data(:).N_simmsacs]; %simulated microsacs
N_blanks = [all_SU_data(:).N_blanks];

mua_sm_sigma = 0.005/dt;
min_MUA_rate = 25;
MU_expt_nums = [all_MU_data(:).expt_num];
MU_avg_rates = [all_MU_data(:).avg_rates]/dt;
MU_bar_oris = [all_MU_data(:).bar_ori];
jbe_MUs = find(strcmp([all_MU_data(:).animal],'jbe') & MU_avg_rates >= min_MUA_rate);
lem_MUs = find(strcmp([all_MU_data(:).animal],'lem') & MU_avg_rates >= min_MUA_rate);
fov_MUs = find(MU_expt_nums < 200 | ismember(MU_expt_nums,lem_fov_expt_nums));
parafov_MUs = find(MU_expt_nums > 200 & ~ismember(MU_expt_nums,lem_fov_expt_nums));
lem_fov_MUs = intersect(fov_MUs,lem_MUs);
lem_parafov_MUs = intersect(parafov_MUs,lem_MUs);

%% GSACS AND MSACS MADE ON GRAY BACKGROUNDS

gsac_used_SUs = find(N_gray_gsacs >= min_Nsacs & avg_rates >= min_rate);
msac_used_SUs = find(N_gray_msacs >= min_Nsacs & avg_rates >= min_rate);

all_gsac_gray = reshape([all_SU_data(:).gsac_gray_avg],[],length(all_SU_data))';
all_msac_gray = reshape([all_SU_data(:).msac_gray_avg],[],length(all_SU_data))';

if sm_sigma > 0
    for ii = 1:size(all_gsac_gray,1)
        all_gsac_gray(ii,:) = jmm_smooth_1d_cor(all_gsac_gray(ii,:),sm_sigma);
        all_msac_gray(ii,:) = jmm_smooth_1d_cor(all_msac_gray(ii,:),sm_sigma);
    end
end

% close all
xl = [-0.15 0.4];

f1 = figure(); hold on
curSUs = intersect(jbe_SUs,gsac_used_SUs);
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','r'});
curSUs = intersect(lem_SUs,gsac_used_SUs);
h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','b'});
xlim(xl);
ylim([0.7 1.3]);
legend([h1.mainLine h2.mainLine],{'JBE','LEM'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
title('Gsac TA Grayback');

f2 = figure(); hold on
curSUs = intersect(jbe_SUs,msac_used_SUs);
h1=shadedErrorBar(tlags,nanmean(all_msac_gray(curSUs,:)),nanstd(all_msac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','r'});
curSUs = intersect(lem_SUs,msac_used_SUs);
h2=shadedErrorBar(tlags,nanmean(all_msac_gray(curSUs,:)),nanstd(all_msac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','b'});
curSUs = intersect(jbe_SUs,gsac_used_SUs);
h3=plot(tlags,nanmean(all_gsac_gray(curSUs,:)),'m','linewidth',2);
curSUs = intersect(lem_SUs,gsac_used_SUs);
h4=plot(tlags,nanmean(all_gsac_gray(curSUs,:)),'k','linewidth',2);
xlim(xl);
ylim([0.7 1.3]);
% legend([h1.mainLine h2.mainLine],{'JBE','LEM'});
line(xl,[1 1],'color','k');
line([0 0],ylim(),'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
title('Msac TA Grayback');


fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'SUA_Gsac_TA_Gback.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
figufy(f2);
fname = [fig_dir 'SUA_Msac_TA_Gback.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

%% SCATTERPLOT OF GSAC/MSAC EXC/SUP magnitude and timing
xr = [0 0.35];
poss_lagrange = find(tlags > xr(1) & tlags < xr(2));
% [gsac_exc,gsac_excloc] = max(all_gsac_gray(:,poss_lagrange),[],2);
% [gsac_inh,gsac_inhloc] = min(all_gsac_gray(:,poss_lagrange),[],2);
% [msac_exc,msac_excloc] = max(all_msac_gray(:,poss_lagrange),[],2);
% [msac_inh,msac_inhloc] = min(all_msac_gray(:,poss_lagrange),[],2);

[gsac_exc,gsac_inh,gsac_excloc,gsac_inhloc] = deal(nan(size(all_gsac_gray,1),1));
for ii = 1:size(all_gsac_gray,1)
    if ~isnan(all_gsac_gray(ii,poss_lagrange))
        [temp,temploc] = findpeaks(all_gsac_gray(ii,poss_lagrange),'sortstr','descend');
        gsac_exc(ii) = temp(1); gsac_excloc(ii) = temploc(1);
        [temp,temploc] = findpeaks(-all_gsac_gray(ii,poss_lagrange),'sortstr','descend');
        gsac_inh(ii) = -temp(1); gsac_inhloc(ii) = temploc(1);
    else
        gsac_exc(ii) = nan; gsac_excloc(ii) = 1;
        gsac_inh(ii) = nan; gsac_inhloc(ii) = 1;
    end
end

[msac_exc,msac_inh,msac_excloc,msac_inhloc] = deal(nan(size(all_msac_gray,1),1));
for ii = 1:size(all_msac_gray,1)
    if ~isnan(all_msac_gray(ii,poss_lagrange))
        [temp,temploc] = findpeaks(all_msac_gray(ii,poss_lagrange),'sortstr','descend');
        msac_exc(ii) = temp(1); msac_excloc(ii) = temploc(1);
        [temp,temploc] = findpeaks(-all_msac_gray(ii,poss_lagrange),'sortstr','descend');
        msac_inh(ii) = -temp(1); msac_inhloc(ii) = temploc(1);
    else
        msac_exc(ii) = nan; msac_excloc(ii) = 1;
        msac_inh(ii) = nan; msac_inhloc(ii) = 1;
    end
end

gsac_Efact = gsac_exc - 1;
gsac_Sfact = 1-gsac_inh;
msac_Efact = msac_exc - 1;
msac_Sfact = 1 - msac_inh;
gsac_exctime = tlags(poss_lagrange(gsac_excloc));
msac_exctime = tlags(poss_lagrange(msac_excloc));
gsac_inhtime = tlags(poss_lagrange(gsac_inhloc));
msac_inhtime = tlags(poss_lagrange(msac_inhloc));

mS = 3; %marker size

%plot E/S modulation strengths
f1 = figure(); 
subplot(2,1,1);hold on
curSUs = intersect(jbe_SUs,gsac_used_SUs);
plot(gsac_Efact(curSUs),gsac_Sfact(curSUs),'ro','markersize',mS);
curSUs = intersect(lem_SUs,gsac_used_SUs);
plot(gsac_Efact(curSUs),gsac_Sfact(curSUs),'bo','markersize',mS);
% legend('JBE','LEM');
line([0 1],[0 1],'color','k');
set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);
xlabel('Enhancement strength');
ylabel('Suppression strength');
title('Guided Sacs');
subplot(2,1,2);hold on
curSUs = intersect(jbe_SUs,msac_used_SUs);
plot(msac_Efact(curSUs),msac_Sfact(curSUs),'ro','markersize',mS);
curSUs = intersect(lem_SUs,msac_used_SUs);
plot(msac_Efact(curSUs),msac_Sfact(curSUs),'bo','markersize',mS);
line([0 1],[0 1],'color','k');
set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);
xlabel('Enhancement strength');
ylabel('Suppression strength');
title('Micro Sacs');


%plot modulation timing
jit_amp = 0.001; %prevent points from occluding each other with some random jitter
f2 = figure(); 
subplot(2,1,1);hold on
curSUs = intersect(jbe_SUs,gsac_used_SUs);
plot(gsac_exctime(curSUs)+randn(length(curSUs),1)*jit_amp,gsac_inhtime(curSUs)+randn(length(curSUs),1)*jit_amp,'ro','markersize',mS);
curSUs = intersect(lem_SUs,gsac_used_SUs);
plot(gsac_exctime(curSUs)+randn(length(curSUs),1)*jit_amp,gsac_inhtime(curSUs)+randn(length(curSUs),1)*jit_amp,'o','markersize',mS);
% legend('JBE','LEM');
line(xr,xr,'color','k');
xlim(xr); ylim(xr);
xlabel('Enhancement peak time (s)');
ylabel('Suppression peak time (s)');
title('Guided Sacs');

subplot(2,1,2);hold on
curSUs = intersect(jbe_SUs,msac_used_SUs);
plot(msac_exctime(curSUs)+randn(length(curSUs),1)*jit_amp,msac_inhtime(curSUs)+randn(length(curSUs),1)*jit_amp,'ro','markersize',mS);
curSUs = intersect(lem_SUs,msac_used_SUs);
plot(msac_exctime(curSUs)+randn(length(curSUs),1)*jit_amp,msac_inhtime(curSUs)+randn(length(curSUs),1)*jit_amp,'o','markersize',mS);
xlim(xr); ylim(xr);
line(xr,xr,'color','k');
xlabel('Enhancement peak time (s)');
ylabel('Suppression peak time (s)');
title('Micro Sacs');


fig_width = 3.5; rel_height = 1.8;

figufy(f1);
fname = [fig_dir 'SUA_EnhSup_scatter.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'SUA_EnhSupTime_scatter.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);


%% PREF VS NONPREF ORI

gsac_used_SUs = find(N_gray_gsacs >= min_Nsacs & avg_rates >= min_rate & ...
    [all_SU_data(:).has_NP]);
gsac_used_SUs_NP = find([all_SU_NPdata(:).N_gsacs_gray] >= min_Nsacs & ...
    [all_SU_NPdata(:).avg_rates]/dt >= min_rate);

%find units that are usable at two different orientations
ind_set1 = [[all_SU_data(gsac_used_SUs).expt_num]' [all_SU_data(gsac_used_SUs).SU_numbers]'];
ind_set2 = [[all_SU_NPdata(gsac_used_SUs_NP).expt_num]' [all_SU_NPdata(gsac_used_SUs_NP).SU_numbers]'];
[C,IA,IB] = intersect(ind_set1,ind_set2,'rows');

all_gsac_gray = reshape([all_SU_data(gsac_used_SUs((IA))).gsac_gray_avg],[],length(IA))';
all_gsac_gray_NP = reshape([all_SU_NPdata(gsac_used_SUs_NP(IB)).gsac_gray_avg],[],length(IB))';

if sm_sigma > 0
    for ii = 1:size(all_gsac_gray,1)
        all_gsac_gray(ii,:) = jmm_smooth_1d_cor(all_gsac_gray(ii,:),sm_sigma);
        all_gsac_gray_NP(ii,:) = jmm_smooth_1d_cor(all_gsac_gray_NP(ii,:),sm_sigma);
    end
end

% close all
xl = [-0.15 0.4];
yl = [0.7 1.3];

f1 = figure(); hold on
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray),nanstd(all_gsac_gray)/sqrt(length(IA)),{'color','r'});
h2=shadedErrorBar(tlags,nanmean(all_gsac_gray_NP),nanstd(all_gsac_gray_NP)/sqrt(length(IA)),{'color','b'});
xlim(xl);
ylim(yl);
legend([h1.mainLine h2.mainLine],{'PREF','Non-PREF'});
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

fig_width = 3.5; rel_height = 0.8;

figufy(f1);
fname = [fig_dir 'SUA_Gsac_TA_PREFNONPREF.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% MUA GRAYBACK TA

mua_gsac_avg = [all_MU_data(:).gsac_gray_avg]';
mua_msac_avg = [all_MU_data(:).msac_gray_avg]';

if mua_sm_sigma > 0
    for ii = 1:size(mua_gsac_avg,1)
        mua_gsac_avg(ii,:) = jmm_smooth_1d_cor(mua_gsac_avg(ii,:),mua_sm_sigma);
        mua_msac_avg(ii,:) = jmm_smooth_1d_cor(mua_msac_avg(ii,:),mua_sm_sigma);
    end
end

xl = [-0.15 0.4];
yl = [0.6 1.5];

f1 = figure(); hold on
% h1=plot(tlags,nanmean(mua_gsac_avg(jbe_MUs,:)),'r','linewidth',2);
% h2=plot(tlags,nanmean(mua_gsac_avg(lem_MUs,:)),'b','linewidth',2);
h1=shadedErrorBar(tlags,nanmean(mua_gsac_avg(jbe_MUs,:)),nanstd(mua_gsac_avg(jbe_MUs,:)),{'color','r'});
h2=shadedErrorBar(tlags,nanmean(mua_gsac_avg(lem_MUs,:)),nanstd(mua_gsac_avg(lem_MUs,:)),{'color','b'});
% h2=shadedErrorBar(tlags,nanmean(mua_gsac_avg(lem_fov_MUs,:)),nanstd(mua_gsac_avg(lem_fov_MUs,:))/sqrt(length(lem_fov_MUs)),{'color','b'});
% h3=shadedErrorBar(tlags,nanmean(mua_gsac_avg(lem_parafov_MUs,:)),nanstd(mua_gsac_avg(lem_parafov_MUs,:))/sqrt(length(lem_parafov_MUs)),{'color','r'});
xlim(xl);ylim(yl);
% legend([h1 h2.mainLine h3.mainLine],{'JBE','LEM-fov','LEM-parafov'});
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
title('Gsac TA Grayback');

xl = [-0.15 0.4];
yl = [0.7 1.3];

f2 = figure(); hold on
% h1=plot(tlags,nanmean(mua_msac_avg(jbe_MUs,:)),'r','linewidth',2);
% h2=plot(tlags,nanmean(mua_msac_avg(lem_MUs,:)),'b','linewidth',2);
h1=shadedErrorBar(tlags,nanmean(mua_msac_avg(jbe_MUs,:)),nanstd(mua_msac_avg(jbe_MUs,:)),{'color','r'});
h2=shadedErrorBar(tlags,nanmean(mua_msac_avg(lem_MUs,:)),nanstd(mua_msac_avg(lem_MUs,:)),{'color','b'});
% h2=shadedErrorBar(tlags,nanmean(mua_msac_avg(lem_fov_MUs,:)),nanstd(mua_msac_avg(lem_fov_MUs,:))/sqrt(length(lem_fov_MUs)),{'color','b'});
% h3=shadedErrorBar(tlags,nanmean(mua_msac_avg(lem_parafov_MUs,:)),nanstd(mua_msac_avg(lem_parafov_MUs,:))/sqrt(length(lem_parafov_MUs)),{'color','r'});
xlim(xl);ylim(yl);
% legend([h1 h2.mainLine h3.mainLine],{'JBE','LEM-fov','LEM-parafov'});
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
title('Msac TA Grayback');


fig_width = 3.5; rel_height = 0.8;

figufy(f1);
fname = [fig_dir 'MUA_Gsac_TA_Gback.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'MUA_Msac_TA_Gback.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);



%% SU GSACS AND SIMSACS (IMAGE BACK) SU

gsac_used_SUs = find(N_im_gsacs >= min_Nsacs & avg_rates >= min_rate);
simsac_used_SUs = find(N_simsacs >= min_Nsacs & avg_rates >= min_rate);

all_gsac_im = reshape([all_SU_data(gsac_used_SUs).gsac_im_avg],[],length(gsac_used_SUs))';
all_simsac = reshape([all_SU_data(simsac_used_SUs).simsac_avg],[],length(simsac_used_SUs))';
if sm_sigma > 0
    for ii = 1:size(all_gsac_im,1)
        all_gsac_im(ii,:) = jmm_smooth_1d_cor(all_gsac_im(ii,:),sm_sigma);
        all_simsac(ii,:) = jmm_smooth_1d_cor(all_simsac(ii,:),sm_sigma);
    end
end

xl = [-0.15 0.4];
yl = [0.6 1.3];

f1 = figure(); hold on

curSUs = find(ismember(gsac_used_SUs,jbe_SUs));
h1=shadedErrorBar(tlags,nanmean(all_gsac_im(curSUs,:)),nanstd(all_gsac_im(curSUs,:))/sqrt(length(curSUs)),{'color','r'});
curSUs = find(ismember(gsac_used_SUs,lem_SUs));
h2=shadedErrorBar(tlags,nanmean(all_gsac_im(curSUs,:)),nanstd(all_gsac_im(curSUs,:))/sqrt(length(curSUs)),{'color','b'});
% curSUs = find(ismember(gsac_used_SUs,[jbe_SUs lem_SUs]));
% h1=shadedErrorBar(tlags,nanmean(all_gsac_im(curSUs,:)),nanstd(all_gsac_im(curSUs,:))/sqrt(length(curSUs)),{'color','r'});

curSUs = find(ismember(simsac_used_SUs,jbe_SUs));
h3=shadedErrorBar(tlags,nanmean(all_simsac(curSUs,:)),nanstd(all_simsac(curSUs,:))/sqrt(length(curSUs)),{'color','m'});
curSUs = find(ismember(simsac_used_SUs,lem_SUs));
h4=shadedErrorBar(tlags,nanmean(all_simsac(curSUs,:)),nanstd(all_simsac(curSUs,:))/sqrt(length(curSUs)),{'color','k'});
% curSUs = find(ismember(simsac_used_SUs,[lem_SUs jbe_SUs]));
% h4=shadedErrorBar(tlags,nanmean(all_simsac(curSUs,:)),nanstd(all_simsac(curSUs,:))/sqrt(length(curSUs)),{'color','k'});


xlim(xl);ylim(yl)
% legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'JBE-gsac','LEM-gsac','JBE-simsac','LEM-simsac'},'Location','Southeast');
line(xl,[1 1],'color','k');
line([0 0],yl,'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

fig_width = 3.5; rel_height = 0.8;

figufy(f1);
fname = [fig_dir 'SUA_Gsac_Simsac_TA.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);


%% MUA SIMSAC TA
mua_gsac_avg = [all_MU_data(:).gsac_im_avg]';
mua_simsac_avg = [all_MU_data(:).simsac_avg]';
mua_simmsac_avg = [all_MU_data(:).simmsac_avg]';

if mua_sm_sigma > 0
    for ii = 1:size(mua_gsac_avg,1)
        mua_gsac_avg(ii,:) = jmm_smooth_1d_cor(mua_gsac_avg(ii,:),mua_sm_sigma);
        mua_simsac_avg(ii,:) = jmm_smooth_1d_cor(mua_simsac_avg(ii,:),mua_sm_sigma);
        mua_simmsac_avg(ii,:) = jmm_smooth_1d_cor(mua_simmsac_avg(ii,:),mua_sm_sigma);
    end
end

xl = [-0.15 0.4];
f1 = figure(); hold on
h1=plot(tlags,nanmean(mua_gsac_avg(jbe_MUs,:)),'r','linewidth',2);
% h2=shadedErrorBar(tlags,nanmean(mua_gsac_avg(lem_MUs,:)),nanstd(mua_gsac_avg(lem_MUs,:))/sqrt(sum(~isnan(mua_gsac_avg(lem_MUs,1)))),{'color','b'});
h2=plot(tlags,nanmean(mua_gsac_avg(lem_MUs,:)),'b','linewidth',2);
h3=plot(tlags,nanmean(mua_simsac_avg(jbe_MUs,:)),'m','linewidth',2);
% h4=shadedErrorBar(tlags,nanmean(mua_simsac_avg(lem_MUs,:)),nanstd(mua_simsac_avg(lem_MUs,:))/sqrt(sum(~isnan(mua_simsac_avg(lem_MUs,1)))),{'color','r'});
h4=plot(tlags,nanmean(mua_simsac_avg(lem_MUs,:)),'k','linewidth',2);
% h4=shadedErrorBar(tlags,nanmean(mua_simsac_avg(lem_fov_MUs,:)),nanstd(mua_simsac_avg(lem_fov_MUs,:))/sqrt(sum(~isnan(mua_simsac_avg(lem_fov_MUs,1)))),{'color','r'});
xlim(xl);
legend([h1 h2 h3 h4],{'JBE-gsac','LEM-gsac','JBE-simsac','LEM-simsac'},'Location','Southeast');
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.65 1.4])


%for simulated microsacs
% f2 = figure(); hold on
% h1=plot(tlags,nanmean(mua_gsac_avg(jbe_MUs,:)),'k','linewidth',2);
% h2=shadedErrorBar(tlags,nanmean(mua_gsac_avg(lem_MUs,:)),nanstd(mua_gsac_avg(lem_MUs,:))/sqrt(sum(~isnan(mua_gsac_avg(lem_MUs,1)))),{'color','b'});
% h3=plot(tlags,nanmean(mua_simmsac_avg(jbe_MUs,:)),'m','linewidth',2);
% h4=shadedErrorBar(tlags,nanmean(mua_simmsac_avg(lem_MUs,:)),nanstd(mua_simmsac_avg(lem_MUs,:))/sqrt(sum(~isnan(mua_simmsac_avg(lem_MUs,1)))),{'color','r'});
% xlim(xl);
% legend([h1 h2.mainLine h3 h4.mainLine],{'JBE-gsac','LEM-gsac','JBE-simmsac','LEM-simmsac'},'Location','Southeast');
% line(xl,[1 1],'color','k');
% xlabel('Time (s)');
% ylabel('Relative rate');
% ylim([0.65 1.4])

fig_width = 3.5; rel_height = 0.8;

figufy(f1);
fname = [fig_dir 'MUA_Gsac_Simsac_TA.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);



%% COMPARE GSACS on gray and image backs

gsac_im_used_SUs = find(N_im_gsacs >= min_Nsacs & avg_rates >= min_rate);

all_gsac_gray = reshape([all_SU_data(gsac_im_used_SUs).gsac_gray_avg],[],length(gsac_im_used_SUs))';
all_gsac_im = reshape([all_SU_data(gsac_im_used_SUs).gsac_im_avg],[],length(gsac_im_used_SUs))';
if sm_sigma > 0
    for ii = 1:size(all_gsac_gray,1)
        all_gsac_gray(ii,:) = jmm_smooth_1d_cor(all_gsac_gray(ii,:),sm_sigma);
        all_gsac_im(ii,:) = jmm_smooth_1d_cor(all_gsac_im(ii,:),sm_sigma);
    end
end

xl = [-0.15 0.4];

f1 = figure(); hold on
curSUs = find(ismember(gsac_im_used_SUs,jbe_SUs));
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','r'});
curSUs = find(ismember(gsac_im_used_SUs,lem_SUs));
h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','b'});

curSUs = find(ismember(gsac_im_used_SUs,jbe_SUs));
h3=shadedErrorBar(tlags,nanmean(all_gsac_im(curSUs,:)),nanstd(all_gsac_im(curSUs,:))/sqrt(length(curSUs)),{'color','m'});
curSUs = find(ismember(gsac_im_used_SUs,lem_SUs));
h4=shadedErrorBar(tlags,nanmean(all_gsac_im(curSUs,:)),nanstd(all_gsac_im(curSUs,:))/sqrt(length(curSUs)),{'color','k'});

xlim(xl);
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'JBE-gray','LEM-gray','JBE-im','LEM-im'},'Location','Southeast');
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.65 1.3])


mua_gsac_gray_avg = [all_MU_data(:).gsac_gray_avg]';
mua_gsac_im_avg = [all_MU_data(:).gsac_im_avg]';
if mua_sm_sigma > 0
    for ii = 1:size(mua_gsac_gray_avg,1)
        mua_gsac_gray_avg(ii,:) = jmm_smooth_1d_cor(mua_gsac_gray_avg(ii,:),mua_sm_sigma);
        mua_gsac_im_avg(ii,:) = jmm_smooth_1d_cor(mua_gsac_im_avg(ii,:),mua_sm_sigma);
    end
end

f2 = figure(); hold on
h1=plot(tlags,nanmean(mua_gsac_gray_avg(jbe_MUs,:)),'r','linewidth',2);
% h2=shadedErrorBar(tlags,nanmean(mua_gsac_gray_avg(lem_MUs,:)),nanstd(mua_gsac_gray_avg(lem_MUs,:))/sqrt(sum(~isnan(mua_gsac_gray_avg(lem_MUs,1)))),{'color','b'});
h2=plot(tlags,nanmean(mua_gsac_gray_avg(lem_MUs,:)),'b','linewidth',2);
h3=plot(tlags,nanmean(mua_gsac_im_avg(jbe_MUs,:)),'m','linewidth',2);
% h4=shadedErrorBar(tlags,nanmean(mua_gsac_im_avg(lem_MUs,:)),nanstd(mua_gsac_im_avg(lem_MUs,:))/sqrt(sum(~isnan(mua_gsac_im_avg(lem_MUs,1)))),{'color','m'});
h4=plot(tlags,nanmean(mua_gsac_im_avg(lem_MUs,:)),'k','linewidth',2);
xlim(xl);
legend([h1 h2 h3 h4],{'JBE-gray','LEM-gray','JBE-im','LEM-im'},'Location','Southeast');
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.65 1.4])


fig_width = 3.5; rel_height = 0.8;

figufy(f1);
fname = [fig_dir 'SUA_Gsac_GrayIm_TA.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'MUA_Gsac_GrayIm_TA.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);
% 


%% COMPARE GSACS and Blank TAs
% gsac_blank_used_SUs = find(N_blanks >= min_Nsacs & avg_rates >= min_rate);
% 
% all_gsac_gray = reshape([all_SU_data(gsac_blank_used_SUs).gsac_gray_avg],[],length(gsac_blank_used_SUs))';
% all_blank = reshape([all_SU_data(gsac_blank_used_SUs).blank_avg],[],length(gsac_blank_used_SUs))';
% if sm_sigma > 0
%     for ii = 1:size(all_gsac_gray,1)
%         all_gsac_gray(ii,:) = jmm_smooth_1d_cor(all_gsac_gray(ii,:),sm_sigma);
%         all_blank(ii,:) = jmm_smooth_1d_cor(all_blank(ii,:),sm_sigma);
%     end
% end
% 
% xl = [-0.2 0.4];
% 
% f1 = figure(); hold on
% curSUs = find(ismember(gsac_blank_used_SUs,jbe_SUs));
% h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','r'});
% h2=shadedErrorBar(tlags,nanmean(all_blank(curSUs,:)),nanstd(all_blank(curSUs,:))/sqrt(length(curSUs)),{'color','b'});
% 
% curSUs = find(ismember(gsac_blank_used_SUs,lem_SUs));
% h3=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','k'});
% h4=shadedErrorBar(tlags,nanmean(all_blank(curSUs,:)),nanstd(all_blank(curSUs,:))/sqrt(length(curSUs)),{'color','m'});
% 
% xlim(xl);
% legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'JBE-gray','LEM-gray','JBE-im','LEM-im'},'Location','Southeast');
% line(xl,[1 1],'color','k');
% xlabel('Time (s)');
% ylabel('Relative rate');
% ylim([0.65 1.4])


mua_gsac_gray_avg = [all_MU_data(:).gsac_gray_avg]';
mua_blank_avg = [all_MU_data(:).blank_avg]';
if mua_sm_sigma > 0
    for ii = 1:size(mua_gsac_gray_avg,1)
        mua_gsac_gray_avg(ii,:) = jmm_smooth_1d_cor(mua_gsac_gray_avg(ii,:),mua_sm_sigma);
        mua_blank_avg(ii,:) = jmm_smooth_1d_cor(mua_blank_avg(ii,:),mua_sm_sigma);
    end
end

uset = find(~isnan(mua_blank_avg(:,1)));

xl = [-0.1 0.3];

f2 = figure(); hold on
h1=shadedErrorBar(tlags,nanmean(mua_gsac_gray_avg(uset,:)),nanstd(mua_gsac_gray_avg(uset,:))/sqrt(length(uset)),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(mua_blank_avg(uset,:)),nanstd(mua_blank_avg(uset,:))/sqrt(length(uset)),{'color','k'});
xlim(xl);
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
% ylim([0.65 1.4])


fig_width = 3.5; rel_height = 0.8;


figufy(f2);
fname = [fig_dir 'MUA_BLANK_TA.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);



%% LARGE VS SMALL MSACS
mua_small_avg = [all_MU_data(:).small_msac_avg]';
mua_large_avg = [all_MU_data(:).large_msac_avg]';
mua_gsac_avg = [all_MU_data(:).gsac_avg]';

if mua_sm_sigma > 0
    for ii = 1:size(mua_small_avg,1)
        mua_small_avg(ii,:) = jmm_smooth_1d_cor(mua_small_avg(ii,:),mua_sm_sigma);
        mua_large_avg(ii,:) = jmm_smooth_1d_cor(mua_large_avg(ii,:),mua_sm_sigma);
        mua_gsac_avg(ii,:) = jmm_smooth_1d_cor(mua_gsac_avg(ii,:),mua_sm_sigma);
    end
end

f1 = figure(); hold on
h1=plot(tlags,nanmean(mua_small_avg(jbe_MUs,:)),'k','linewidth',2);
h3=plot(tlags,nanmean(mua_large_avg(jbe_MUs,:)),'m','linewidth',2);
% h5=plot(tlags,nanmean(mua_gsac_avg(jbe_MUs,:)),'g','linewidth',2);
h2=shadedErrorBar(tlags,nanmean(mua_small_avg(lem_MUs,:)),nanstd(mua_small_avg(lem_MUs,:))/sqrt(sum(~isnan(mua_small_avg(lem_MUs,1)))),{'color','b'});
h4=shadedErrorBar(tlags,nanmean(mua_large_avg(lem_MUs,:)),nanstd(mua_large_avg(lem_MUs,:))/sqrt(sum(~isnan(mua_large_avg(lem_MUs,1)))),{'color','r'});
% h6=shadedErrorBar(tlags,nanmean(mua_gsac_avg(lem_MUs,:)),nanstd(mua_gsac_avg(lem_MUs,:))/sqrt(sum(~isnan(mua_gsac_avg(lem_MUs,1)))),{'color','c'});
xlim(xl);


fig_width = 3.5; rel_height = 0.8;

% figufy(f1);
% fname = [fig_dir 'SUA_Gsac_TA_Gback.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'SUA_Msac_TA_Gback.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

%% MSAC DIRECTION DEPENDENCE
mua_vert_avg = [all_MU_data(:).msac_vert_avg]';
mua_hori_avg = [all_MU_data(:).msac_hor_avg]';
if mua_sm_sigma > 0
    for ii = 1:size(mua_vert_avg,1)
        mua_vert_avg(ii,:) = jmm_smooth_1d_cor(mua_vert_avg(ii,:),mua_sm_sigma);
        mua_hori_avg(ii,:) = jmm_smooth_1d_cor(mua_hori_avg(ii,:),mua_sm_sigma);
    end
end
hori_MUs = find(MU_bar_oris == 0);
vert_MUs = find(MU_bar_oris == 90);

f1 = figure(); 
subplot(2,1,1)
hold on
h1=shadedErrorBar(tlags,nanmean(mua_vert_avg(hori_MUs,:)),nanstd(mua_vert_avg(hori_MUs,:))/sqrt(sum(~isnan(mua_vert_avg(hori_MUs,1)))),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(mua_hori_avg(hori_MUs,:)),nanstd(mua_hori_avg(hori_MUs,:))/sqrt(sum(~isnan(mua_hori_avg(hori_MUs,1)))),{'color','r'});
xlim(xl);
subplot(2,1,2)
hold on
h1=shadedErrorBar(tlags,nanmean(mua_vert_avg(vert_MUs,:)),nanstd(mua_vert_avg(vert_MUs,:))/sqrt(sum(~isnan(mua_vert_avg(vert_MUs,1)))),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(mua_hori_avg(vert_MUs,:)),nanstd(mua_hori_avg(vert_MUs,:))/sqrt(sum(~isnan(mua_hori_avg(vert_MUs,1)))),{'color','r'});
xlim(xl);


% fig_width = 3.5; rel_height = 0.8;

%% ANALYZE LAMINAR DEPENDENCIES
load('/home/james/Analysis/bruce/FINsac_mod/layer_boundaries/layer_classification.mat')
boundary_enums = [boundary_class(:).Expt_num];

gsac_gray_used_SUs = find(N_gray_gsacs >= min_Nsacs & avg_rates >= min_rate);

mua_gsac_gray_avg = [all_MU_data(:).gsac_gray_avg]';
sua_gsac_gray_avg = reshape([all_SU_data(:).gsac_gray_avg],[],length(all_SU_data))';
if mua_sm_sigma > 0
    for ii = 1:size(mua_gsac_gray_avg,1)
        mua_gsac_gray_avg(ii,:) = jmm_smooth_1d_cor(mua_gsac_gray_avg(ii,:),mua_sm_sigma);
    end
end
if sm_sigma > 0
    for ii = 1:size(sua_gsac_gray_avg,1)
        sua_gsac_gray_avg(ii,:) = jmm_smooth_1d_cor(sua_gsac_gray_avg(ii,:),sm_sigma);
    end
end

un_lem_expts = unique(MU_expt_nums(lem_MUs));
n_lem_expts = length(un_lem_expts);

gran_gsac = nan(n_lem_expts,length(tlags));
supra_gsac = nan(n_lem_expts,length(tlags));
infra_gsac = nan(n_lem_expts,length(tlags));
all_gran_supt = [];
all_supra_supt = [];
all_infra_supt = [];
all_gran_enht = [];
all_supra_enht = [];
all_infra_enht = [];
all_supt = [];
all_enht = [];
all_relloc = [];
all_gran_SUs = [];
all_infra_SUs = [];
all_supra_SUs = [];
for ee = 1:n_lem_expts
   cur_mua_set = find(MU_expt_nums == un_lem_expts(ee));
   cur_bound_info = find(boundary_enums == un_lem_expts(ee),1);
   cur_ub = boundary_class(cur_bound_info).ub;
   cur_lb = boundary_class(cur_bound_info).lb;
    
   gran_probes = (cur_ub+1):(cur_lb-1);
   supra_probes = 1:(cur_ub-1);
   infra_probes = (cur_lb+1):24;
   
   [~,cur_supt] = min(mua_gsac_gray_avg(cur_mua_set,:),[],2);
   [~,cur_enht] = max(mua_gsac_gray_avg(cur_mua_set,:),[],2);
   all_supt = cat(1,all_supt,cur_supt);
   all_enht = cat(1,all_enht,cur_enht);
   cur_relpos = (1:24) - cur_ub;
   all_relloc = cat(2,all_relloc,cur_relpos);
   if ~isempty(gran_probes)
       gran_gsac(ee,:) = mean(mua_gsac_gray_avg(cur_mua_set(gran_probes),:));
       all_gran_supt = cat(1,all_gran_supt,cur_supt(gran_probes));
       all_gran_enht = cat(1,all_gran_enht,cur_enht(gran_probes));
   end
   if ~isempty(supra_probes)
       supra_gsac(ee,:) = mean(mua_gsac_gray_avg(cur_mua_set(supra_probes),:));
       all_supra_supt = cat(1,all_supra_supt,cur_supt(supra_probes));
       all_supra_enht = cat(1,all_supra_enht,cur_enht(supra_probes));
   end
   if ~isempty(infra_probes)
       infra_gsac(ee,:) = mean(mua_gsac_gray_avg(cur_mua_set(infra_probes),:));
       all_infra_supt = cat(1,all_infra_supt,cur_supt(infra_probes));
       all_infra_enht = cat(1,all_infra_enht,cur_enht(infra_probes));
   end
   
   cur_SU_set = find([all_SU_data(:).expt_num] == un_lem_expts(ee));
   cur_SU_set = cur_SU_set(ismember(cur_SU_set,gsac_gray_used_SUs));
   cur_SU_probenums = [all_SU_data(cur_SU_set).probe_numbers];
   all_gran_SUs = cat(2,all_gran_SUs,cur_SU_set(ismember(cur_SU_probenums,gran_probes)));
   all_infra_SUs = cat(2,all_infra_SUs,cur_SU_set(ismember(cur_SU_probenums,infra_probes)));
   all_supra_SUs = cat(2,all_supra_SUs,cur_SU_set(ismember(cur_SU_probenums,supra_probes)));
end

xr = [-0.05 0.3];
yl = [0.7 1.3];

f1 = figure(); hold on
h2=shadedErrorBar(tlags,nanmean(supra_gsac),nanstd(supra_gsac)/sqrt(n_lem_expts),{'color','b'});
h1=shadedErrorBar(tlags,nanmean(gran_gsac),nanstd(gran_gsac)/sqrt(n_lem_expts),{'color','k'});
h3=shadedErrorBar(tlags,nanmean(infra_gsac),nanstd(infra_gsac)/sqrt(n_lem_expts),{'color','r'});
xlim(xr);
ylim(yl);
% line(mean(tlags(all_gran_supt)) + [0 0],yl,'color','k');
% line(mean(tlags(all_supra_supt)) + [0 0],yl,'color','b');
% line(mean(tlags(all_infra_supt)) + [0 0],yl,'color','r');
% line(mean(tlags(all_gran_enht)) + [0 0],yl,'color','k');
% line(mean(tlags(all_supra_enht)) + [0 0],yl,'color','b');
% line(mean(tlags(all_infra_enht)) + [0 0],yl,'color','r');

% f2 = figure(); hold on
% h1=shadedErrorBar(tlags,nanmean(sua_gsac_gray_avg(all_gran_SUs,:)),nanstd(sua_gsac_gray_avg(all_gran_SUs,:))/sqrt(length(all_gran_SUs)),{'color','k'});
% h2=shadedErrorBar(tlags,nanmean(sua_gsac_gray_avg(all_supra_SUs,:)),nanstd(sua_gsac_gray_avg(all_supra_SUs,:))/sqrt(length(all_supra_SUs)),{'color','b'});
% h3=shadedErrorBar(tlags,nanmean(sua_gsac_gray_avg(all_infra_SUs,:)),nanstd(sua_gsac_gray_avg(all_infra_SUs,:))/sqrt(length(all_infra_SUs)),{'color','r'});

[granSU_sup,granSU_supt] = min(sua_gsac_gray_avg(all_gran_SUs,:),[],2);
[infraSU_sup,infraSU_supt] = min(sua_gsac_gray_avg(all_infra_SUs,:),[],2);
[supraSU_sup,supraSU_supt] = min(sua_gsac_gray_avg(all_supra_SUs,:),[],2);

% % trange = tlags(tlags > 0.02 & tlags < 0.125);
trange = linspace(0.025,0.125,20);
gran_hist = histc(tlags(all_gran_supt),trange);
sup_hist = histc(tlags(all_supra_supt),trange);
infra_hist = histc(tlags(all_infra_supt),trange);

f2 = figure();
% subplot(2,1,1);
hold on
stairs(trange,gran_hist/sum(gran_hist),'k');
stairs(trange,sup_hist/sum(sup_hist),'b');
stairs(trange,infra_hist/sum(infra_hist),'r')
xlabel('Suppresion timing');
ylabel('Relative freq');
legend('Granular','Supra-gran','Infra-gran');
xlim([0.025 0.125]);

% % trange = tlags(tlags > 0.05 & tlags < 0.3);
% trange = linspace(0.,0.3,20);
% gran_hist = histc(tlags(all_gran_enht),trange);
% sup_hist = histc(tlags(all_supra_enht),trange);
% infra_hist = histc(tlags(all_infra_enht),trange);
% subplot(2,1,2);
% hold on
% stairs(trange,gran_hist/sum(gran_hist),'k');
% stairs(trange,sup_hist/sum(sup_hist),'b');
% stairs(trange,infra_hist/sum(infra_hist),'r')
% xlim([0 0.3]);



fig_width = 3.5; rel_height = 0.8;
figufy(f1);
fname = [fig_dir 'MUA_Lamdep_STA.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'MUA_Lamdep_suphist.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);


% 
% jit_amp = 0.001;
% f2 = figure(); hold on
% plot(tlags(all_gran_supt)+randn(size(all_gran_supt))*jit_amp,tlags(all_gran_enht)+randn(size(all_gran_supt))*jit_amp,'k.');
% plot(tlags(all_supra_supt)+randn(size(all_supra_supt))*jit_amp,tlags(all_supra_enht)+randn(size(all_supra_enht))*jit_amp,'.');
% plot(tlags(all_infra_supt)+randn(size(all_infra_supt))*jit_amp,tlags(all_infra_enht)+randn(size(all_infra_enht))*jit_amp,'r.');
% xlim([0 0.15]); ylim([0.1 0.3]);



%% LOOK AT ECCENTRICITY DEPENDENCE OF SAC MOD WITH MUA
close all
xl = [-0.2 0.4];

use_hor = true;

if use_hor
    use_mua_data = all_jbehor_mua_data;
else
    use_mua_data = all_jbever_mua_data;
end

lem_mua_gsac_avg = cat(2,all_lem_mua_data(:).gsac_gray_avg)';
lem_mua_msac_avg = cat(2,all_lem_mua_data(:).msac_gray_avg)';

% jbe_mua_gsac_avg = cat(2,use_mua_data(:).gsac_gray_avg)';
jbe_mua_gsac_avg = cat(2,all_jbehor_mua_data(:).gsac_gray_avg,all_jbever_mua_data(:).gsac_gray_avg)';
jbe_mua_gsac_avg = reshape(jbe_mua_gsac_avg,96,[],length(tlags));
jbe_mua_gsac_avg = squeeze(nanmean(jbe_mua_gsac_avg,2));

% jbe_mua_msac_avg = cat(2,use_mua_data(:).msac_gray_avg)';
jbe_mua_msac_avg = cat(2,all_jbehor_mua_data(:).msac_gray_avg,all_jbever_mua_data(:).msac_gray_avg)';
jbe_mua_msac_avg = reshape(jbe_mua_msac_avg,96,[],length(tlags));
jbe_mua_msac_avg = squeeze(nanmean(jbe_mua_msac_avg,2));

jbe_mua_gsac_avg_pos = cat(2,use_mua_data(:).gsac_pos_avg)';
% jbe_mua_gsac_avg_pos = cat(2,use_mua_data(:).gsac_outpos_avg)';
% jbe_mua_gsac_avg_pos = cat(2,use_mua_data(:).gsac_inpos_avg)';
jbe_mua_gsac_avg_pos = reshape(jbe_mua_gsac_avg_pos,96,[],length(tlags));
jbe_mua_gsac_avg_pos = squeeze(nanmean(jbe_mua_gsac_avg_pos,2));

jbe_mua_gsac_avg_neg = cat(2,use_mua_data(:).gsac_neg_avg)';
% jbe_mua_gsac_avg_neg = cat(2,use_mua_data(:).gsac_outneg_avg)';
% jbe_mua_gsac_avg_neg = cat(2,use_mua_data(:).gsac_inneg_avg)';
jbe_mua_gsac_avg_neg = reshape(jbe_mua_gsac_avg_neg,96,[],length(tlags));
jbe_mua_gsac_avg_neg = squeeze(nanmean(jbe_mua_gsac_avg_neg,2));

jbe_mua_msac_avg_tow = cat(2,use_mua_data(:).msac_towards_avg)';
jbe_mua_msac_avg_tow = reshape(jbe_mua_msac_avg_tow,96,[],length(tlags));
jbe_mua_msac_avg_tow = squeeze(nanmean(jbe_mua_msac_avg_tow,2));

jbe_mua_msac_avg_awa = cat(2,use_mua_data(:).msac_away_avg)';
jbe_mua_msac_avg_awa = reshape(jbe_mua_msac_avg_awa,96,[],length(tlags));
jbe_mua_msac_avg_awa = squeeze(nanmean(jbe_mua_msac_avg_awa,2));

load ~/Data/bruce/general_array_data/array_pos_data.mat
array_ecc = sqrt(interp_x.^2 + interp_y.^2);
[~,ecc_ord] = sort(array_ecc);

for ii = 1:length(tlags)
    ecc_corr(ii) = corr(array_ecc,jbe_mua_gsac_avg(:,ii),'type','spearman'); 
    ecc_corr_msac(ii) = corr(array_ecc,jbe_mua_msac_avg(:,ii),'type','spearman'); 
    ecc_corr_pos(ii) = corr(array_ecc,jbe_mua_gsac_avg_pos(:,ii),'type','spearman'); 
    ecc_corr_neg(ii) = corr(array_ecc,jbe_mua_gsac_avg_neg(:,ii),'type','spearman'); 
end

f0 = figure();hold on
shadedErrorBar(tlags,mean(jbe_mua_gsac_avg(ecc_ord(1:10),:)),std(jbe_mua_gsac_avg(ecc_ord(1:10),:))/sqrt(10));
shadedErrorBar(tlags,mean(jbe_mua_gsac_avg(ecc_ord(87:end),:)),std(jbe_mua_gsac_avg(ecc_ord(87:end),:))/sqrt(10),{'color','r'});
plot(tlags,mean(lem_mua_gsac_avg(lem_fov_MUs,:)),'b','linewidth',2);
plot(tlags,mean(lem_mua_gsac_avg(lem_parafov_MUs,:)),'g','linewidth',2);
xlim(xl);
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

f01 = figure();
plot(tlags,mean(jbe_mua_gsac_avg(ecc_ord(1:10),:))-mean(jbe_mua_gsac_avg(ecc_ord(87:end),:)),'k','linewidth',2);
line(xl,[0 0],'color','k');
xlim([-0.1 0.3]);
xlabel('Time (s)');
ylabel('Relative rate difference');

f02 = figure();hold on
shadedErrorBar(tlags,mean(jbe_mua_msac_avg(ecc_ord(1:10),:)),std(jbe_mua_msac_avg(ecc_ord(1:10),:))/sqrt(10));
shadedErrorBar(tlags,mean(jbe_mua_msac_avg(ecc_ord(87:end),:)),std(jbe_mua_msac_avg(ecc_ord(87:end),:))/sqrt(10),{'color','r'});
plot(tlags,mean(lem_mua_msac_avg(lem_fov_MUs,:)),'b','linewidth',2);
plot(tlags,mean(lem_mua_msac_avg(lem_parafov_MUs,:)),'g','linewidth',2);
xlim(xl);
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');

% f1 = figure();
% imagesc(tlags,1:96,jbe_mua_gsac_avg(ecc_ord,:));

f1 = figure();
shadedErrorBar(tlags,mean(jbe_mua_gsac_avg_pos),std(jbe_mua_gsac_avg_pos)/10);
hold on
shadedErrorBar(tlags,mean(jbe_mua_gsac_avg_neg),std(jbe_mua_gsac_avg_neg)/10,{'color','r'});
plot(tlags,mean(lem_mua_gsac_avg(lem_fov_MUs,:)),'b','linewidth',2);
plot(tlags,mean(lem_mua_gsac_avg(lem_parafov_MUs,:)),'g','linewidth',2);
xlim(xl);
% legend([h1 h2.mainLine h3.mainLine],{'JBE','LEM-fov','LEM-parafov'});
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.6 1.5])

f3= figure;
shadedErrorBar(tlags,mean(jbe_mua_msac_avg_tow),std(jbe_mua_msac_avg_tow)/10);
hold on
shadedErrorBar(tlags,mean(jbe_mua_msac_avg_awa),std(jbe_mua_msac_avg_awa)/10,{'color','r'});
% plot(tlags,mean(lem_mua_msac_avg(lem_fov_MUs,:)),'b','linewidth',2);
% plot(tlags,mean(lem_mua_msac_avg(lem_parafov_MUs,:)),'g','linewidth',2);
xlim(xl);
% legend([h1 h2.mainLine h3.mainLine],{'JBE','LEM-fov','LEM-parafov'});
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.7 1.4])

f2 = figure; hold on
% plot(tlags,ecc_corr_pos,tlags,ecc_corr_neg,'r');
plot(tlags,ecc_corr);
xlim(xl);
xlabel('Time (s)');
ylabel('Correlation');

% f3 = figure();
% subplot(2,1,1);
% imagesc(tlags,1:96,jbe_mua_gsac_avg_pos(ecc_ord,:));
% caxis([0.6 1.5])
% subplot(2,1,2);
% imagesc(tlags,1:96,jbe_mua_gsac_avg_neg(ecc_ord,:));
% caxis([0.6 1.5])
% 
% f4 = figure();
% subplot(2,1,1);
% imagesc(tlags,1:96,jbe_mua_msac_avg_awa(ecc_ord,:));
% caxis([0.6 1.5])
% subplot(2,1,2);
% imagesc(tlags,1:96,jbe_mua_msac_avg_tow(ecc_ord,:));
% caxis([0.6 1.5])


fig_width = 4; rel_height = 0.8;

% 
% figufy(f0);
% fname = [fig_dir 'MUA_FP_ecc.pdf'];
% exportfig(f0,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f0);

figufy(f01);
fname = [fig_dir 'MUA_FP_gsac_ecc_diff.pdf'];
exportfig(f01,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f01);

% figufy(f02);
% fname = [fig_dir 'MUA_FP_msac_ecc.pdf'];
% exportfig(f02,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f02);
% 
% figufy(f1);
% fname = [fig_dir 'MUA_FP_verpn.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

% figufy(f2);
% fname = [fig_dir 'MUA_FP_horcorr.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

% figufy(f3);
% fname = [fig_dir 'MUA_FP_vermsac.pdf'];
% exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f3);


%%
close all
% jbe_mua_gsac_avg = cat(2,all_jbehor_mua_data(:).gsac_gray_avg,all_jbever_mua_data(:).gsac_gray_avg)';
% jbe_mua_gsac_avg = cat(2,all_jbehor_mua_data(:).gsac_gray_avg)';
jbe_mua_gsac_avg = cat(2,all_jbever_mua_data(:).gsac_gray_avg)';
jbe_mua_gsac_avg = reshape(jbe_mua_gsac_avg,96,[],length(tlags));
jbe_mua_gsac_avg = squeeze(nanmean(jbe_mua_gsac_avg,2));

jbe_mua_msac_avg = cat(2,all_jbehor_mua_data(:).msac_gray_avg,all_jbever_mua_data(:).msac_gray_avg)';
jbe_mua_msac_avg = reshape(jbe_mua_msac_avg,96,[],length(tlags));
jbe_mua_msac_avg = squeeze(nanmean(jbe_mua_msac_avg,2));

% jbe_mua_msac_avg_tow = cat(2,all_jbehor_mua_data(:).msac_towards_avg,all_jbever_mua_data(:).msac_towards_avg)';
% jbe_mua_msac_avg_tow = cat(2,all_jbehor_mua_data(:).msac_towards_avg)';
% jbe_mua_msac_avg_tow = cat(2,all_jbever_mua_data(:).msac_towards_avg)';
jbe_mua_msac_avg_tow = cat(2,all_jbever_mua_data(:).gsac_pos_avg)';
% jbe_mua_msac_avg_tow = cat(2,all_jbehor_mua_data(:).gsac_pos_avg)';
jbe_mua_msac_avg_tow = reshape(jbe_mua_msac_avg_tow,96,[],length(tlags));
jbe_mua_msac_avg_tow = squeeze(nanmean(jbe_mua_msac_avg_tow,2));

% jbe_mua_msac_avg_awa = cat(2,all_jbehor_mua_data(:).msac_away_avg,all_jbever_mua_data(:).msac_away_avg)';
% jbe_mua_msac_avg_awa = cat(2,all_jbehor_mua_data(:).msac_away_avg)';
% jbe_mua_msac_avg_awa = cat(2,all_jbever_mua_data(:).msac_away_avg)';
jbe_mua_msac_avg_awa = cat(2,all_jbever_mua_data(:).gsac_neg_avg)';
% jbe_mua_msac_avg_awa = cat(2,all_jbehor_mua_data(:).gsac_neg_avg)';
jbe_mua_msac_avg_awa = reshape(jbe_mua_msac_avg_awa,96,[],length(tlags));
jbe_mua_msac_avg_awa = squeeze(nanmean(jbe_mua_msac_avg_awa,2));

load ~/Data/bruce/general_array_data/array_pos_data.mat
array_ecc = sqrt(interp_x.^2 + interp_y.^2);
[~,ecc_ord] = sort(array_ecc);

for ii = 1:length(tlags)
    ecc_corr(ii) = corr(array_ecc,jbe_mua_gsac_avg(:,ii),'type','spearman'); 
    ecc_corr_msac(ii) = corr(array_ecc,jbe_mua_msac_avg(:,ii),'type','spearman'); 
end

f1 = figure();
imagesc(tlags,1:96,jbe_mua_gsac_avg(ecc_ord,:));

f2 = figure();
shadedErrorBar(tlags,mean(jbe_mua_msac_avg_tow),std(jbe_mua_msac_avg_tow)/10);
hold on
shadedErrorBar(tlags,mean(jbe_mua_msac_avg_awa),std(jbe_mua_msac_avg_awa)/10,{'color','r'});

f3 = figure();
subplot(2,1,1);
imagesc(tlags,1:96,jbe_mua_msac_avg_awa(ecc_ord,:));
caxis([0.6 1.5])
subplot(2,1,2);
imagesc(tlags,1:96,jbe_mua_msac_avg_tow(ecc_ord,:));
caxis([0.6 1.5])
