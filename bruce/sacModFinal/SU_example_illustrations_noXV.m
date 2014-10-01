clear all
close all
clc

Expt_name = 'G093';
bar_ori = 0;
cc = 103;

fit_unCor = 0;
include_bursts = 0;

fig_dir = '/home/james/Analysis/bruce/FINsac_mod/figures/';
base_tname = 'sac_trig_avg_data';
base_sname = 'sacStimProcFin_noXV';
base_yname = 'sacTypeDep_noXV';

if include_bursts
    base_tname = strcat(base_tname,'_withbursts');
    base_sname = strcat(base_sname,'_withbursts');
    base_yname = strcat(base_yname,'_withbursts');
end

if Expt_name(1) == 'M'
    rec_type = 'LP';
    good_coils = [1 1]; %which coils are usable
    use_coils = [1 1]; %[L R] Use info from coils?
    n_probes = 24;
elseif Expt_name(1) == 'G'
    good_coils = [1 0]; %which coils are usable
    use_coils = [0 0]; %[L R] Use info from coils?
    n_probes = 96;
end

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';

%if using coil info
if any(use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end

et_mod_data_name = [et_dir et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_dir et_anal_name sprintf('_ori%d',bar_ori)];
load(et_anal_name,'et_params');

Expt_num = str2num(Expt_name(2:end));
sac_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];

%load trig avg data
tname = strcat(sac_dir,base_tname,sprintf('_ori%d',bar_ori));
load(tname);
%load stimulus proc data
sname = strcat(sac_dir,base_sname,sprintf('_ori%d',bar_ori));
if fit_unCor
    sname = strcat(sname,'_unCor');
end
load(sname);
%load type-dep data (STILL NEED TO ADD THIS)
yname = strcat(sac_dir,base_yname,sprintf('_ori%d',bar_ori));
if fit_unCor
    yname = strcat(yname,'_unCor');
end
load(yname);

su_range = (n_probes+1):length(sacStimProc);
su_ind = find(su_range == cc);
clear SU_data
sacStimProc = sacStimProc(su_range(su_ind));
trig_avg = sua_data(su_ind);
type_dep = sacTypeDep(su_range(su_ind));


cur_GQM = sacStimProc.ModData.rectGQM;
flen = cur_GQM.stim_params(1).stim_dims(1);
sp_dx = et_params.sp_dx;
use_nPix = et_params.use_nPix;
use_nPix_us = use_nPix*et_params.spatial_usfac;

tlags = trig_avg_params.lags*trig_avg_params.dt;
cid = sprintf('E%d_C%d_',Expt_num,cc);

%% PLOT MODEL FILTERS
close all
stim_dims = sacStimProc.ModData.rectGQM.stim_params(1).stim_dims;
pix_ax = (1:stim_dims(2))*sp_dx;
pix_ax = pix_ax - mean(pix_ax);
lag_ax = ((1:stim_dims(1))*dt - dt/2)*1e3;
pix_lim = [-0.3 0.3];
lag_lim = [0 140];
[fig_props] = plot_NMM_filters_1d(sacStimProc.ModData.rectGQM,pix_ax,lag_ax);

for ii = 1:fig_props.nmods
    subplot(fig_props.dims(1),fig_props.dims(2),ii)
    xlim(pix_lim);
    ylim(lag_lim);
end

% fig_width = 2*fig_props.dims(2);
% rel_height = 0.8*fig_props.dims(1)/fig_props.dims(2);
% fname = [fig_dir cid sprintf('ori%d_',bar_ori) 'stim_mod.pdf'];
% exportfig(fig_props.h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(fig_props.h);

%% STA figs
sac_xr = [-0.1 0.3];

close all
stim_dims = sacStimProc.ModData.rectGQM.stim_params(1).stim_dims;
pix_ax = (1:stim_dims(2))*sp_dx;
pix_ax = pix_ax - mean(pix_ax);
lag_ax = ((1:stim_dims(1))*dt - dt/2)*1e3;

ov_sta = sacStimProc.ov_phaseDep_sta;
[~,sta_peakloc] = max(sum(ov_sta.^2,2));
% sta_peakloc = 7;
cond_STA = sacStimProc.gsac_phaseDep_sta;
cond_STA = squeeze(cond_STA(:,sta_peakloc,:));

space_sm = 0.75;
time_sm = 0;

%smooth conditional STA in space
if space_sm > 0
    for iii = 1:size(cond_STA,2)
        cond_STA(:,iii) = jmm_smooth_1d_cor(cond_STA(:,iii),space_sm);
    end
end

%smooth conditional STA in time
if time_sm > 0
    for iii = 1:size(cond_STA,1)
        cond_STA(iii,:) = jmm_smooth_1d_cor(cond_STA(iii,:),time_sm);
    end
end

xl = [-0.3 0.3];

f1 = figure();
imagesc(slags*dt,pix_ax,cond_STA');
cam = max(abs(cond_STA(:)));
caxis([-cam cam]*0.8);
yl = ylim();
xlim(sac_xr)
xlabel('Time (s)');
ylabel('Rel Position (deg)');
ylim(xl);

f2 = figure();
imagesc(pix_ax,lag_ax,ov_sta);
cam = max(abs(ov_sta(:)));
caxis([-cam cam]);
cid = sprintf('E%d_C%d_',Expt_num,cc);
set(gca,'ydir','normal');
xlim(xl);
line(xl,lag_ax([sta_peakloc sta_peakloc]),'color','k');


%% MAIN TB FIGURES
close all
sac_xr = [-0.1 0.3];
gr = [-2 3];
TB_lambda = 2;

Gtick = sacStimProc.gsac_TBmod{TB_lambda}.Gtick;
Xtick = sacStimProc.gsac_TBmod{TB_lambda}.lagX*dt;

base_slags = find(Xtick <= 0);

[mmm,TBmmmloc] = minmax(sacStimProc.gsac_TBmod{TB_lambda}.an_sac_rate);

%PLOT TB RATE MAP
f1 = figure();
imagesc(Xtick,Gtick,sacStimProc.gsac_TBmod{TB_lambda}.TB_rate); 
set(gca,'ydir','normal'); 
caxis([0 max(sacStimProc.gsac_TBmod{TB_lambda}.TB_rate(:))]*0.4);
yl = ylim();
ylim(gr); yl = gr;
line(Xtick(TBmmmloc([1 1])),yl,'color','w');
line(Xtick(TBmmmloc([2 2])),yl,'color','w');
xlim(sac_xr)
xlabel('Time (s)');
ylabel('Generating signal');

%PLOT TB MARGINAL FIRING RATE
f2 = figure();
plot(Xtick,sacStimProc.gsac_TBmod{TB_lambda}.an_sac_rate/dt);
hold on; axis tight
yl = ylim();
line(Xtick(TBmmmloc([1 1])),yl,'color','k','linestyle','--');
line(Xtick(TBmmmloc([2 2])),yl,'color','k','linestyle','--');
xlim(sac_xr);
line(sac_xr,[0 0]+sacStimProc.ModData.unit_data.avg_rate,'color','k');
xlabel('Time (s)');
ylabel('Firing rate (Hz)');

%PLOT TB RATEFUN SLICES
base_out = sacStimProc.gsac_TBmod{TB_lambda}.marg_grate;
cur_gdistY = sacStimProc.gsac_TBmod{TB_lambda}.equi_gdist;
cur_gdistX = sacStimProc.gsac_TBmod{TB_lambda}.equi_gX;
f3 = figure(); hold on
plot(Gtick,sacStimProc.gsac_TBmod{TB_lambda}.TB_rate(:,TBmmmloc(1))/dt,'g');
plot(Gtick,sacStimProc.gsac_TBmod{TB_lambda}.TB_rate(:,TBmmmloc(2))/dt,'r');
plot(Gtick,base_out/dt,'k','linewidth',2);
xlim(gr);
yl = ylim();
plot(cur_gdistX,cur_gdistY/max(cur_gdistY)*yl(2)*0.75,'k--');
xlabel('Generating signal');
ylabel('Firing rate (Hz)');
axis tight;
xlim(Gtick([1 end]));
xlim(gr);
xlabel('Generating signal');
ylabel('Firing rate (Hz)');


cur_sac_offset = sacStimProc.gsac_TBmod{TB_lambda}.sac_offset/dt;
cur_sac_offset = bsxfun(@minus,cur_sac_offset,mean(cur_sac_offset(base_slags))); %normalize by pre-saccade values
%TB RESPONSE OFFSETs
f4 = figure();
plot(Xtick,sacStimProc.cur_sac_offset);
xlim(sac_xr)
yl = ylim();
line(Xtick(TBmmmloc([1 1])),yl,'color','g');
line(Xtick(TBmmmloc([2 2])),yl,'color','r');
xlabel('Time (s)');
ylabel('Offset (Hz)');

cur_sac_gain = sacStimProc.gsac_TBmod{TB_lambda}.sac_gain;
cur_sac_gain = bsxfun(@rdivide,cur_sac_gain,mean(cur_sac_gain(base_slags)));
%TB RESPONSE GAINS
f5 = figure();
plot(Xtick,sacStimProc.cur_sac_gain,'r');
xlim(sac_xr)
yl = ylim();
line(Xtick(TBmmmloc([1 1])),yl,'color','g');
line(Xtick(TBmmmloc([2 2])),yl,'color','r');
xlabel('Time (s)');
ylabel('Gain');

%% PLOT TB INFO AND INFO RATE
TB_lambda = 2;
base_slags = find(Xtick <= 0);

TB_SSI = sacStimProc.gsac_TBmod{TB_lambda}.sac_modinfo;
TB_info_rate = sacStimProc.gsac_TBmod{TB_lambda}.sac_modinfo.*sacStimProc.gsac_avg_rate/dt;
% avg_info = sacStimProc.gsac_TBmod{TB_lambda}.ovInfo;
base_info = mean(TB_SSI(base_slags));
base_info_rate = mean(TB_info_rate(base_slags));

f1 = figure();
plot(slags*dt,TB_SSI);
xlim(sac_xr)
line(sac_xr,[0 0] + base_info,'color','k');
xlabel('Time (s)');
ylabel('SSI (bits/spk)');

f2 = figure();
plot(slags*dt,TB_info_rate,'k');
xlim(sac_xr)
set(gca,'YaxisLocation','right');
line(sac_xr,[0 0] + base_info_rate,'color','k');
xlabel('Time (s)');
ylabel('SSI (bits/sec)');

%% GO MODEL
lambda_off = 2;
lambda_gain = 2;
TB_lambda = 2;

base_slags = find(slags <= 0);

GO_offset = sacStimProc.gsac_post_mod{lambda_off,lambda_gain}.mods(2).filtK;
GO_offset = bsxfun(@minus,GO_offset,mean(GO_offset(base_slags)));

% INTERNAL OFFSET
f1 = figure();
plot(slags*dt,GO_offset)
xlim(sac_xr)
yl = ylim(); ya = max(abs(yl)); ylim([-ya ya]);
xlabel('Time (s)');
ylabel('Offset');

GO_gain = 1+sacStimProc.gsac_post_mod{lambda_off,lambda_gain}.mods(3).filtK;
GO_gain = bsxfun(@rdivide,GO_gain,mean(GO_gain(base_slags)));

%GO MODEL INTERNAL GAIN
f2 = figure();
plot(slags*dt,GO_gain,'r');
xlim(sac_xr)
yl = ylim(); ya = max(abs(yl-1)); ylim([-ya ya]+1);
xlabel('Time (s)');
ylabel('Gain');
set(gca,'YaxisLocation','right');

GO_sac_offset = sacStimProc.gsac_post_mod{lambda_off,lambda_gain}.sac_offset/dt;
GO_sac_offset = bsxfun(@minus,GO_sac_offset,mean(GO_sac_offset(base_slags)));
TB_sac_offset = sacStimProc.gsac_TBmod{TB_lambda}.sac_offset/dt;
TB_sac_offset = bsxfun(@minus,TB_sac_offset,mean(TB_sac_offset(base_slags)));

%GO MODEL RESPONSE OFFSET
f3 = figure(); hold on
plot(slags*dt,GO_sac_offset)
plot(Xtick,TB_sac_offset,'r')
xlim(sac_xr)
yl = ylim(); ya = max(abs(yl)); ylim([-ya ya]);
xlabel('Time (s)');
ylabel('Offset');

GO_sac_gain = sacStimProc.gsac_post_mod{lambda_off,lambda_gain}.sac_offset/dt;
GO_sac_gain = bsxfun(@rdivide,GO_sac_gain,mean(GO_sac_gain(base_slags)));
TB_sac_gain = sacStimProc.gsac_TBmod{TB_lambda}.sac_offset/dt;
TB_sac_gain = bsxfun(@rdivide,TB_sac_gain,mean(TB_sac_gain(base_slags)));

%GO MODEL RESPONSE GAIN
f4 = figure(); hold on
plot(slags*dt,sacStimProc.gsac_post_mod{lambda_off,lambda_gain}.sac_gain)
plot(Xtick,sacStimProc.gsac_TBmod{lambda_off,lambda_gain}.sac_gain,'r')
xlim(sac_xr)
yl = ylim(); ya = max(abs(yl-1)); ylim([-ya ya]+1);
xlabel('Time (s)');
ylabel('Gain');
set(gca,'YaxisLocation','right');

%% PRE POST GAIN COMPARISON
lambda_off = 2;
lambda_gain = 1;
pre_lambda = 1;

base_slags = find(slags <= 0);

post_gain = 1+sacStimProc.gsac_post_mod{lambda_off,lambda_gain}.mods(3).filtK;
pre_gain = 1+sacStimProc.gsacPreGainMod{pre_lambda}.stim_kernel;

post_gain = bsxfun(@rdivide,post_gain,mean(post_gain(base_slags)));
pre_gain = bsxfun(@rdivide,pre_gain,mean(pre_gain(base_slags)));

cur_xr = [-0.05 0.2];
f1 = figure();
plot(slags*dt,post_gain);
hold on
plot(slags*dt,pre_gain,'k');
xlim(cur_xr)
xlabel('Time (s)');
ylabel('Gain');

%% COMPARE INDIVIDUAL GAINS AND TEMPORAL KERNELS

lambda_d2T = 2;
lambda_L2 = 1;
base_slags = find(slags <= 0);

stim_dims = sacStimProc.ModData.rectGQM.stim_params(1).stim_dims;
lag_ax = ((1:stim_dims(1))*dt - dt/2)*1e3;

stim_mod = sacStimProc.ModData.rectGQM;
rel_weights = stim_mod.rel_filt_weights;
[all_tkerns,avg_Ekern,avg_Ikern] = get_hilbert_tempkerns(stim_mod);
nmods = size(all_tkerns,2);
all_Ntkerns = bsxfun(@rdivide,all_tkerns,rel_weights);
fpost_gains = 1+reshape([sacStimProc.gsac_post_Fullmod{1,lambda_d2T,lambda_L2}.mods(3).filtK],length(slags),nmods);
fpost_gains(:,rel_weights == 0) = [];

%divide by pre-sac gains?
fpost_gains = bsxfun(@rdivide,fpost_gains,mean(fpost_gains(base_slags,:)));

f1 = figure();
subplot(2,1,1);
imagesc(lag_ax,1:nmods,all_Ntkerns');
% imagesc(lag_ax,1:nmods,all_tkerns');
ca = caxis();
subplot(2,1,2);
imagesc(slags*dt,1:nmods,fpost_gains');
ca = caxis(); cam = max(abs(ca-1)); caxis([-cam cam]+1);
% xlim([0 0.15]);
%%



