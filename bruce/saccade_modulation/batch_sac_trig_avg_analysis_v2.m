%%
close all
clear all
fig_dir = '/home/james/Analysis/bruce/saccade_modulation/';

%% LOAD JBE HORIZONTAL
sname = 'sacStimProc_v2';
tname = 'sac_trig_avg_data6';
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
rmfield_list = {};

all_hor_data = [];
all_hor_tdata = [];
all_jbehor_sua_exptnum = [];
all_jbehor_mua_exptnum = [];
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
    cd(save_dir)
    load(sname)
    
    ucells = arrayfun(@(x) length(x.ModData),sacStimProc) > 0;
    cur_data = sacStimProc(ucells);
    cur_data = rmfields(cur_data,rmfield_list);
    [cur_data.expt_num] = deal(Expt_num);
    [cur_data.bar_ori] = deal(0);
    
    all_hor_data = cat(2,all_hor_data,cur_data);
    cur_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,cur_data);

    tdata = load(tname);
    tdat_SU_numbers = arrayfun(@(x) x.unit_data.SU_numbers,tdata.sua_trig_avgs);
    uset = find(ismember(tdat_SU_numbers,cur_SU_numbers));
    cur_tdata = tdata.sua_trig_avgs(uset);
    [cur_tdata.expt_num] = deal(Expt_num);
    [cur_tdata.bar_ori] = deal(0);
    all_hor_tdata = cat(2,all_hor_tdata,cur_tdata);

    all_jbehor_sua_data(ee) = tdata.sua_data;
    all_jbehor_sua_exptnum = cat(1,all_jbehor_sua_exptnum,ones(size(tdata.sua_data.msac_avg,2),1)*Expt_num);
    all_jbehor_mua_data(ee) = tdata.mua_data;
    all_jbehor_mua_exptnum = cat(1,all_jbehor_mua_exptnum,ones(size(tdata.mua_data.msac_avg,2),1)*Expt_num);
end

%% LOAD JBE VERTICAL
sname = 'sacStimProc_v2_vbars';
tname = 'sac_trig_avg_data6_vbars';
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
rmfield_list = {};

all_ver_data = [];
all_ver_tdata = [];
all_jbever_sua_exptnum = [];
all_jbever_mua_exptnum = [];
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
    cd(save_dir)
    load(sname)
    
    ucells = arrayfun(@(x) length(x.ModData),sacStimProc) > 0;
    cur_data = sacStimProc(ucells);
    cur_data = rmfields(cur_data,rmfield_list);
    [cur_data.expt_num] = deal(Expt_num);
    [cur_data.bar_ori] = deal(90);
    
    all_ver_data = cat(2,all_ver_data,cur_data);
    cur_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,cur_data);

    tdata = load(tname);
    tdat_SU_numbers = arrayfun(@(x) x.unit_data.SU_numbers,tdata.sua_trig_avgs);
    uset = find(ismember(tdat_SU_numbers,cur_SU_numbers));
    cur_tdata = tdata.sua_trig_avgs(uset);
    [cur_tdata.expt_num] = deal(Expt_num);
    [cur_tdata.bar_ori] = deal(90);
    all_ver_tdata = cat(2,all_ver_tdata,cur_tdata);

    all_jbever_sua_data(ee) = tdata.sua_data;
    all_jbever_sua_exptnum = cat(1,all_jbever_sua_exptnum,ones(size(tdata.sua_data.msac_avg,2),1)*Expt_num);
    all_jbever_mua_data(ee) = tdata.mua_data;
    all_jbever_mua_exptnum = cat(1,all_jbever_mua_exptnum,ones(size(tdata.mua_data.msac_avg,2),1)*Expt_num);
end

%% COMBINE JBE DATA
hor_SU_number = arrayfun(@(x) x.ModData.unit_data.SU_number,all_hor_data);
hor_expt_num = [all_hor_data(:).expt_num];
hor_mod_info = [all_hor_data(:).gsac_spost_ov_modinfo];
hor_avg_rate = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_hor_data);
hor_mod_inforate = hor_mod_info.*hor_avg_rate;

ver_SU_number = arrayfun(@(x) x.ModData.unit_data.SU_number,all_ver_data);
ver_expt_num = [all_ver_data(:).expt_num];
ver_mod_info = [all_ver_data(:).gsac_spost_ov_modinfo];
ver_avg_rate = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_ver_data);
ver_mod_inforate = ver_mod_info.*ver_avg_rate;

all_id_vecs = [hor_expt_num' hor_SU_number'; ver_expt_num' ver_SU_number'];
all_info_rate = [hor_mod_inforate'; ver_mod_inforate'];
all_comb_data = [all_hor_data'; all_ver_data'];
all_comb_tdata = [all_hor_tdata'; all_ver_tdata'];
[C,IA,IC] = unique(all_id_vecs,'rows');

%for each array SU, use data from either the horizontal or vertical bar
%stimuli, depending on which had the better model (in terms of info rate)
all_jbe_data = [];
all_jbe_tdata = [];
all_jbe_nonpref_data = [];
all_jbe_nonpref_tdata = [];
all_jbe_pref_data = [];
all_jbe_pref_tdata = [];
for ii = 1:length(C)
    curset = find(all_id_vecs(:,1) == C(ii,1) & all_id_vecs(:,2) == C(ii,2));
    if length(curset) > 1
        [~,better] = max(all_info_rate(curset));
        worse = setdiff(curset,curset(better));
        curset = curset(better);
        
        all_jbe_nonpref_data = cat(1,all_jbe_nonpref_data,all_comb_data(worse));
        all_jbe_nonpref_tdata = cat(1,all_jbe_nonpref_tdata,all_comb_tdata(worse));
        all_jbe_pref_data = cat(1,all_jbe_pref_data,all_comb_data(curset));
        all_jbe_pref_tdata = cat(1,all_jbe_pref_tdata,all_comb_tdata(curset));
    end
    all_jbe_data = cat(1,all_jbe_data,all_comb_data(curset));
    all_jbe_tdata = cat(1,all_jbe_tdata,all_comb_tdata(curset));
end

[all_jbe_data.animal] = deal('jbe');

%% LOAD LEM DATA
sname = 'sacStimProc_v2';
tname = 'sac_trig_avg_data6';
% Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294'};
% bar_oris = [80 60 135 70 140 90 160 40];
Expt_list = {'M266','M270','M275','M277','M281','M287','M294'}; %NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
bar_oris = [80 60 135 70 140 90 40];
rmfield_list = {};

all_lem_data = [];
all_lem_tdata = [];
lem_sua_exptnum = [];
lem_mua_exptnum = [];
for ee = 1:length(Expt_list)
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    save_dir = ['~/Analysis/bruce/' Expt_name '/sac_mod'];
    cd(save_dir)
    load(sname)
    
    ucells = arrayfun(@(x) length(x.ModData),sacStimProc) > 0;
    cur_data = sacStimProc(ucells);
    cur_data = rmfields(cur_data,rmfield_list);
    [cur_data.expt_num] = deal(Expt_num);
    [cur_data.bar_ori] = deal(bar_oris(ee));
    
    all_lem_data = cat(1,all_lem_data,cur_data');
    cur_SU_numbers = arrayfun(@(x) x.ModData.unit_data.SU_number,cur_data);

    tdata = load(tname);
    tdat_SU_numbers = arrayfun(@(x) x.unit_data.SU_numbers,tdata.sua_trig_avgs);
    uset = find(ismember(tdat_SU_numbers,cur_SU_numbers));
    cur_tdata = tdata.sua_trig_avgs(uset)';
    [cur_tdata.expt_num] = deal(Expt_num);
    [cur_tdata.bar_ori] = deal(bar_oris(ee));
    all_lem_tdata = cat(1,all_lem_tdata,cur_tdata);

    all_lem_sua_data(ee) = tdata.sua_data;
    lem_sua_exptnum = cat(1,lem_sua_exptnum,ones(size(tdata.sua_data.msac_avg,2),1)*Expt_num);
    all_lem_mua_data(ee) = tdata.mua_data;
    lem_mua_exptnum = cat(1,lem_mua_exptnum,ones(size(tdata.mua_data.msac_avg,2),1)*Expt_num);
end

[all_lem_data.animal] = deal('lem');

%% COMBINE JBE AND LEM DATA
all_SU_data = [all_jbe_data; all_lem_data];
all_SU_tdata = [all_jbe_tdata; all_lem_tdata];
tlags = tdata.trig_avg_params.lags;
dt = tdata.trig_avg_params.dt;
tlags = tlags*dt;

%% SELECT USABLE CELLS
avg_rates = arrayfun(@(x) x.unit_data.avg_rates,all_SU_tdata)/dt;
tot_spikes = arrayfun(@(x) x.unit_data.tot_nspikes,all_SU_tdata);
rec_dur = arrayfun(@(x) x.unit_data.N_used_samps,all_SU_tdata)*dt/60; %in min
expt_nums = [all_SU_data(:).expt_num];

clust_iso_dist = arrayfun(@(x) x.ModData.unit_data.SU_isodist,all_SU_data);
clust_Lratio = arrayfun(@(x) x.ModData.unit_data.SU_Lratio,all_SU_data);
clust_refract = arrayfun(@(x) x.ModData.unit_data.SU_refract,all_SU_data);
clust_dprime = arrayfun(@(x) x.ModData.unit_data.SU_dprime,all_SU_data);

sm_sigma = 0.0075/dt;
% sm_sigma = 0;
min_rate = 5; %in Hz
min_Nsacs = 100;
jbe_SUs = find(strcmp('jbe',{all_SU_data(:).animal}));
lem_SUs = find(strcmp('lem',{all_SU_data(:).animal}));

lem_fov_expt_nums = [266 275];
lem_parafov_expt_nums = [270 277 281 287 289 294];
fov_SUs = find([all_SU_data(:).expt_num] < 200 | ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
parafov_SUs = find([all_SU_data(:).expt_num] > 200 & ~ismember([all_SU_data(:).expt_num],lem_fov_expt_nums));
lem_fov_SUs = intersect(fov_SUs,lem_SUs);
lem_parafov_SUs = intersect(parafov_SUs,lem_SUs);

lem_fov_MUs = find(ismember(lem_mua_exptnum,lem_fov_expt_nums));
lem_parafov_MUs = find(ismember(lem_mua_exptnum,lem_parafov_expt_nums));

%% GSACS AND MSACS MADE ON GRAY BACKGROUNDS
N_gray_gsacs = arrayfun(@(x) x.unit_data.N_gsacs_gray,all_SU_tdata);
N_gray_msacs = arrayfun(@(x) x.unit_data.N_msacs_gray,all_SU_tdata);

gsac_used_SUs = find(N_gray_gsacs >= min_Nsacs & avg_rates >= min_rate);
msac_used_SUs = find(N_gray_msacs >= min_Nsacs & avg_rates >= min_rate);

all_gsac_gray = reshape([all_SU_tdata(:).gsac_gray],[],length(all_SU_tdata))';
all_msac_gray = reshape([all_SU_tdata(:).msac_gray],[],length(all_SU_tdata))';
if sm_sigma > 0
    for ii = 1:size(all_gsac_gray,1)
        all_gsac_gray(ii,:) = jmm_smooth_1d_cor(all_gsac_gray(ii,:),sm_sigma);
        all_msac_gray(ii,:) = jmm_smooth_1d_cor(all_msac_gray(ii,:),sm_sigma);
    end
end

close all

xl = [-0.2 0.4];

f1 = figure(); hold on
curSUs = intersect(jbe_SUs,gsac_used_SUs);
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','r'});
curSUs = intersect(lem_SUs,gsac_used_SUs);
h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','b'});
xlim(xl);
legend([h1.mainLine h2.mainLine],{'JBE','LEM'});
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
title('Gsac TA Grayback');
ylim([0.7 1.3]);

f2 = figure(); hold on
curSUs = intersect(jbe_SUs,msac_used_SUs);
h1=shadedErrorBar(tlags,nanmean(all_msac_gray(curSUs,:)),nanstd(all_msac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','r'});
curSUs = intersect(lem_SUs,msac_used_SUs);
h2=shadedErrorBar(tlags,nanmean(all_msac_gray(curSUs,:)),nanstd(all_msac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','b'});
xlim(xl);
legend([h1.mainLine h2.mainLine],{'JBE','LEM'});
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
title('Msac TA Grayback');
ylim([0.7 1.3]);


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

%%
N_gray_gsacs = arrayfun(@(x) x.unit_data.N_gsacs_gray,all_SU_tdata);
N_gray_msacs = arrayfun(@(x) x.unit_data.N_msacs_gray,all_SU_tdata);

gsac_used_SUs = find(N_gray_gsacs >= min_Nsacs & avg_rates >= min_rate);
msac_used_SUs = find(N_gray_msacs >= min_Nsacs & avg_rates >= min_rate);

all_gsac_gray = reshape([all_SU_tdata(:).gsac_gray],[],length(all_SU_tdata))';
all_msac_gray = reshape([all_SU_tdata(:).msac_gray],[],length(all_SU_tdata))';
all_gsac_gray_SD = reshape([all_SU_tdata(:).gsac_gray_SD],[],length(all_SU_tdata))';
all_msac_gray_SD = reshape([all_SU_tdata(:).msac_gray_SD],[],length(all_SU_tdata))';

% for ss = 1:length(gsac_used_SUs)
%     subplot(2,1,1);
%     shadedErrorBar(tlags,all_gsac_gray(gsac_used_SUs(ss),:),all_gsac_gray_SD(gsac_used_SUs(ss),:));
%     subplot(2,1,2);
%     plot(tlags,(all_gsac_gray(gsac_used_SUs(ss),:)-1)./all_gsac_gray_SD(gsac_used_SUs(ss),:));
%     
%     pause
%     clf
%     
% end

poss_lagrange = find(tlags > 0 & tlags < 0.3);
[gsac_exc,gsac_excloc] = max(all_gsac_gray(:,poss_lagrange),[],2);
[gsac_inh,gsac_inhloc] = min(all_gsac_gray(:,poss_lagrange),[],2);
[msac_exc,msac_excloc] = max(all_msac_gray(:,poss_lagrange),[],2);
[msac_inh,msac_inhloc] = min(all_msac_gray(:,poss_lagrange),[],2);
gsac_Efact = gsac_exc - 1;
gsac_Sfact = 1-gsac_inh;
msac_Efact = msac_exc - 1;
msac_Sfact = 1 - msac_inh;
gsac_exctime = tlags(poss_lagrange(gsac_excloc));
msac_exctime = tlags(poss_lagrange(msac_excloc));
gsac_inhtime = tlags(poss_lagrange(gsac_inhloc));
msac_inhtime = tlags(poss_lagrange(msac_inhloc));

mS = 3;

f1 = figure(); 
subplot(2,1,1);hold on
curSUs = intersect(jbe_SUs,gsac_used_SUs);
plot(gsac_Efact(curSUs),gsac_Sfact(curSUs),'o','markersize',mS);
curSUs = intersect(lem_SUs,gsac_used_SUs);
plot(gsac_Efact(curSUs),gsac_Sfact(curSUs),'ro','markersize',mS);
legend('JBE','LEM');
line([0 1],[0 1],'color','k');
xlabel('Enhancement');
ylabel('Suppression');
title('Guided Sacs');

subplot(2,1,2);hold on
curSUs = intersect(jbe_SUs,msac_used_SUs);
plot(msac_Efact(curSUs),msac_Sfact(curSUs),'o','markersize',mS);
curSUs = intersect(lem_SUs,msac_used_SUs);
plot(msac_Efact(curSUs),msac_Sfact(curSUs),'ro','markersize',mS);
line([0 1],[0 1],'color','k');
xlabel('Enhancement');
ylabel('Suppression');
title('Micro Sacs');


f2 = figure(); 
subplot(2,1,1);hold on
curSUs = intersect(jbe_SUs,gsac_used_SUs);
plot(gsac_exctime(curSUs),gsac_inhtime(curSUs),'o','markersize',mS);
curSUs = intersect(lem_SUs,gsac_used_SUs);
plot(gsac_exctime(curSUs),gsac_inhtime(curSUs),'ro','markersize',mS);
legend('JBE','LEM');
line([0 0.3],[0 0.3],'color','k');
xlabel('Enhancement peak');
ylabel('Suppression peak');
title('Guided Sacs');

subplot(2,1,2);hold on
curSUs = intersect(jbe_SUs,msac_used_SUs);
plot(msac_exctime(curSUs),msac_inhtime(curSUs),'o','markersize',mS);
curSUs = intersect(lem_SUs,msac_used_SUs);
plot(msac_exctime(curSUs),msac_inhtime(curSUs),'ro','markersize',mS);
line([0 0.3],[0 0.3],'color','k');
xlabel('Enhancement peak');
ylabel('Suppression peak');
title('Micro Sacs');

% figure;
% subplot(2,1,1);
% hold on
% curSUs = intersect(jbe_SUs,msac_used_SUs);
% plot(gsac_Efact(curSUs),msac_Efact(curSUs),'o');
% curSUs = intersect(lem_SUs,msac_used_SUs);
% plot(gsac_Efact(curSUs),msac_Efact(curSUs),'ro');
% line([0 1],[0 1],'color','k');
% xlabel('Gsac');
% ylabel('Msac');
% title('Enhancement');
% subplot(2,1,2);
% hold on
% curSUs = intersect(jbe_SUs,msac_used_SUs);
% plot(msac_Efact(curSUs),msac_Sfact(curSUs),'o');
% curSUs = intersect(lem_SUs,msac_used_SUs);
% plot(msac_Efact(curSUs),msac_Sfact(curSUs),'ro');
% line([0 1],[0 1],'color','k');
% xlabel('Gsac');
% ylabel('Msac');
% title('Suppression');

% figure;
% subplot(2,1,1);
% hold on
% curSUs = intersect(jbe_SUs,msac_used_SUs);
% plot(gsac_exctime(curSUs),msac_exctime(curSUs),'o');
% curSUs = intersect(lem_SUs,msac_used_SUs);
% plot(gsac_exctime(curSUs),msac_exctime(curSUs),'ro');
% line([0 0.3],[0 0.3],'color','k');
% xlabel('Enhancement peak');
% ylabel('Suppression peak');
% title('Micro Sacs');
% subplot(2,1,2);
% hold on
% curSUs = intersect(jbe_SUs,msac_used_SUs);
% plot(gsac_inhtime(curSUs),msac_inhtime(curSUs),'o');
% curSUs = intersect(lem_SUs,msac_used_SUs);
% plot(gsac_inhtime(curSUs),msac_inhtime(curSUs),'ro');
% line([0 0.3],[0 0.3],'color','k');
% xlabel('Enhancement peak');
% ylabel('Suppression peak');
% title('Micro Sacs');


fig_width = 3.5; rel_height = 1.6;

figufy(f1);
fname = [fig_dir 'SUA_EnhSup_scatter.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'SUA_EnhSupLag_scatter.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);


%% MUA GRAYBACK TA

lem_mua_gsac_avg = cat(2,all_lem_mua_data(:).gsac_gray_avg)';
lem_mua_msac_avg = cat(2,all_lem_mua_data(:).msac_gray_avg)';
jbe_mua_gsac_avg = cat(2,all_jbehor_mua_data(:).gsac_gray_avg,all_jbever_mua_data(:).gsac_gray_avg)';
jbe_mua_msac_avg = cat(2,all_jbehor_mua_data(:).msac_gray_avg,all_jbever_mua_data(:).msac_gray_avg)';

f1 = figure(); hold on
h1=plot(tlags,nanmean(jbe_mua_gsac_avg),'k','linewidth',2);
h2=shadedErrorBar(tlags,nanmean(lem_mua_gsac_avg(lem_fov_MUs,:)),nanstd(lem_mua_gsac_avg(lem_fov_MUs,:))/sqrt(length(lem_fov_MUs)),{'color','b'});
h3=shadedErrorBar(tlags,nanmean(lem_mua_gsac_avg(lem_parafov_MUs,:)),nanstd(lem_mua_gsac_avg(lem_parafov_MUs,:))/sqrt(length(lem_parafov_MUs)),{'color','r'});
xlim(xl);
% legend([h1 h2.mainLine h3.mainLine],{'JBE','LEM-fov','LEM-parafov'});
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
title('Gsac TA Grayback');
ylim([0.7 1.4]);

f2 = figure(); hold on
h1=plot(tlags,nanmean(jbe_mua_msac_avg),'k','linewidth',2);
h2=shadedErrorBar(tlags,nanmean(lem_mua_msac_avg(lem_fov_MUs,:)),nanstd(lem_mua_msac_avg(lem_fov_MUs,:))/sqrt(length(lem_fov_MUs)),{'color','b'});
h3=shadedErrorBar(tlags,nanmean(lem_mua_msac_avg(lem_parafov_MUs,:)),nanstd(lem_mua_msac_avg(lem_parafov_MUs,:))/sqrt(length(lem_parafov_MUs)),{'color','r'});
xlim(xl);
% legend([h1 h2.mainLine h3.mainLine],{'JBE','LEM-fov','LEM-parafov'});
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
title('Msac TA Grayback');
ylim([0.7 1.4]);

% lem_mua_msac_avg = cat(2,all_lem_mua_data(:).msac_gray_avg)';
% jbe_mua_toward_msac_avg = cat(2,all_jbehor_mua_data(:).msac_towards_avg,all_jbever_mua_data(:).msac_towards_avg)';
% jbe_mua_away_msac_avg = cat(2,all_jbehor_mua_data(:).msac_away_avg,all_jbever_mua_data(:).msac_away_avg)';


fig_width = 3.5; rel_height = 0.8;
% 
% figufy(f1);
% fname = [fig_dir 'MUA_Gsac_TA_Gback.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 
% figufy(f2);
% fname = [fig_dir 'MUA_Msac_TA_Gback.pdf'];
% exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f2);

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

%% GSACS AND SIMSACS (IMAGE BACK)
N_im_gsacs = arrayfun(@(x) x.unit_data.N_gsacs_im,all_SU_tdata);
N_simsacs = arrayfun(@(x) x.unit_data.N_simsacs,all_SU_tdata);

gsac_used_SUs = find(N_im_gsacs >= min_Nsacs & avg_rates >= min_rate);
simsac_used_SUs = find(N_simsacs >= min_Nsacs & avg_rates >= min_rate);

all_gsac_im = reshape([all_SU_tdata(:).gsac_im],[],length(all_SU_tdata))';
all_simsac = reshape([all_SU_tdata(:).simsac],[],length(all_SU_tdata))';
if sm_sigma > 0
    for ii = 1:size(all_gsac_im,1)
        all_gsac_im(ii,:) = jmm_smooth_1d_cor(all_gsac_im(ii,:),sm_sigma);
        all_simsac(ii,:) = jmm_smooth_1d_cor(all_simsac(ii,:),sm_sigma);
    end
end

xl = [-0.25 0.5];

f1 = figure(); hold on

curSUs = intersect(jbe_SUs,gsac_used_SUs);
h1=shadedErrorBar(tlags,nanmean(all_gsac_im(curSUs,:)),nanstd(all_gsac_im(curSUs,:))/sqrt(length(curSUs)),{'color','r'});
curSUs = intersect(lem_SUs,gsac_used_SUs);
h2=shadedErrorBar(tlags,nanmean(all_gsac_im(curSUs,:)),nanstd(all_gsac_im(curSUs,:))/sqrt(length(curSUs)),{'color','b'});

curSUs = intersect(jbe_SUs,simsac_used_SUs);
h3=shadedErrorBar(tlags,nanmean(all_simsac(curSUs,:)),nanstd(all_simsac(curSUs,:))/sqrt(length(curSUs)),{'color','m'});
curSUs = intersect(lem_SUs,simsac_used_SUs);
h4=shadedErrorBar(tlags,nanmean(all_simsac(curSUs,:)),nanstd(all_simsac(curSUs,:))/sqrt(length(curSUs)),{'color','k'});
xlim(xl);
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'JBE-gsac','LEM-gsac','JBE-simsac','LEM-simsac'},'Location','Southeast');
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.65 1.4])

fig_width = 4.5; rel_height = 0.8;

figufy(f1);
fname = [fig_dir 'SUA_Gsac_Simsac_TA.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);


%% MUA SIMSAC TA

lem_mua_gsac_avg = cat(2,all_lem_mua_data(:).gsac_im_avg)';
lem_mua_simsac_avg = cat(2,all_lem_mua_data(:).simsac_avg)';
jbe_mua_gsac_avg = cat(2,all_jbehor_mua_data(:).gsac_im_avg,all_jbever_mua_data(:).gsac_im_avg)';
jbe_mua_simsac_avg = cat(2,all_jbehor_mua_data(:).simsac_avg,all_jbever_mua_data(:).simsac_avg)';

f1 = figure(); hold on
h1=plot(tlags,nanmean(jbe_mua_gsac_avg),'k','linewidth',2);
h2=shadedErrorBar(tlags,nanmean(lem_mua_gsac_avg),nanstd(lem_mua_gsac_avg)/sqrt(length(lem_mua_exptnum)),{'color','b'});
h3=plot(tlags,nanmean(jbe_mua_simsac_avg),'m','linewidth',2);
h4=shadedErrorBar(tlags,nanmean(lem_mua_simsac_avg),nanstd(lem_mua_simsac_avg)/sqrt(length(lem_mua_exptnum)),{'color','r'});
xlim(xl);
legend([h1 h2.mainLine h3 h4.mainLine],{'JBE-gsac','LEM-gsac','JBE-simsac','LEM-simsac'},'Location','Southeast');
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.65 1.4])

fig_width = 4.5; rel_height = 0.8;

figufy(f1);
fname = [fig_dir 'MUA_Gsac_Simsac_TA.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% LARGE VS SMALL MSACS
lem_mua_small_msac_avg = cat(2,all_lem_mua_data(:).small_msac_avg)';
lem_mua_large_msac_avg = cat(2,all_lem_mua_data(:).large_msac_avg)';
jbe_mua_small_msac_avg = cat(2,all_jbehor_mua_data(:).small_msac_avg,all_jbever_mua_data(:).small_msac_avg)';
jbe_mua_large_msac_avg = cat(2,all_jbehor_mua_data(:).large_msac_avg,all_jbever_mua_data(:).large_msac_avg)';

f1 = figure(); hold on
h1=plot(tlags,nanmean(jbe_mua_small_msac_avg),'k','linewidth',2);
h2=shadedErrorBar(tlags,nanmean(lem_mua_small_msac_avg),nanstd(lem_mua_small_msac_avg)/sqrt(length(lem_mua_exptnum)),{'color','b'});
h3=plot(tlags,nanmean(jbe_mua_large_msac_avg),'m','linewidth',2);
h4=shadedErrorBar(tlags,nanmean(lem_mua_large_msac_avg),nanstd(lem_mua_large_msac_avg)/sqrt(length(lem_mua_exptnum)),{'color','r'});
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


%% COMPARE GSACS on gray and image backs
N_gray_gsacs = arrayfun(@(x) x.unit_data.N_gsacs_gray,all_SU_tdata);
N_im_gsacs = arrayfun(@(x) x.unit_data.N_gsacs_im,all_SU_tdata);

gsac_gray_used_SUs = find(N_gray_gsacs >= min_Nsacs & avg_rates >= min_rate);
gsac_im_used_SUs = find(N_im_gsacs >= min_Nsacs & avg_rates >= min_rate);

all_gsac_gray = reshape([all_SU_tdata(:).gsac_gray],[],length(all_SU_tdata))';
all_gsac_im = reshape([all_SU_tdata(:).gsac_im],[],length(all_SU_tdata))';
if sm_sigma > 0
    for ii = 1:size(all_gsac_gray,1)
        all_gsac_gray(ii,:) = jmm_smooth_1d_cor(all_gsac_gray(ii,:),sm_sigma);
        all_gsac_im(ii,:) = jmm_smooth_1d_cor(all_gsac_im(ii,:),sm_sigma);
    end
end

xl = [-0.2 0.4];

f1 = figure(); hold on
curSUs = intersect(jbe_SUs,gsac_gray_used_SUs);
h1=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','r'});
curSUs = intersect(lem_SUs,gsac_gray_used_SUs);
h2=shadedErrorBar(tlags,nanmean(all_gsac_gray(curSUs,:)),nanstd(all_gsac_gray(curSUs,:))/sqrt(length(curSUs)),{'color','b'});

curSUs = intersect(jbe_SUs,gsac_im_used_SUs);
h3=shadedErrorBar(tlags,nanmean(all_gsac_im(curSUs,:)),nanstd(all_gsac_im(curSUs,:))/sqrt(length(curSUs)),{'color','k'});
curSUs = intersect(lem_SUs,gsac_im_used_SUs);
h4=shadedErrorBar(tlags,nanmean(all_gsac_im(curSUs,:)),nanstd(all_gsac_im(curSUs,:))/sqrt(length(curSUs)),{'color','m'});

xlim(xl);
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'JBE-gray','LEM-gray','JBE-im','LEM-im'},'Location','Southeast');
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.65 1.4])

lem_mua_gsac_gray_avg = cat(2,all_lem_mua_data(:).gsac_gray_avg)';
lem_mua_gsac_im_avg = cat(2,all_lem_mua_data(:).gsac_im_avg)';
jbe_mua_gsac_gray_avg = cat(2,all_jbehor_mua_data(:).gsac_gray_avg,all_jbever_mua_data(:).gsac_gray_avg)';
jbe_mua_gsac_im_avg = cat(2,all_jbehor_mua_data(:).gsac_im_avg,all_jbever_mua_data(:).gsac_im_avg)';

f2 = figure(); hold on
h1=plot(tlags,nanmean(jbe_mua_gsac_gray_avg),'r','linewidth',2);
h2=shadedErrorBar(tlags,nanmean(lem_mua_gsac_gray_avg),nanstd(lem_mua_gsac_gray_avg)/sqrt(length(lem_mua_exptnum)),{'color','b'});
h3=plot(tlags,nanmean(jbe_mua_gsac_im_avg),'k','linewidth',2);
h4=shadedErrorBar(tlags,nanmean(lem_mua_gsac_im_avg),nanstd(lem_mua_gsac_im_avg)/sqrt(length(lem_mua_exptnum)),{'color','m'});
xlim(xl);
legend([h1 h2.mainLine h3 h4.mainLine],{'JBE-gray','LEM-gray','JBE-im','LEM-im'},'Location','Southeast');
line(xl,[1 1],'color','k');
xlabel('Time (s)');
ylabel('Relative rate');
ylim([0.65 1.4])


fig_width = 4.5; rel_height = 0.8;

figufy(f1);
fname = [fig_dir 'SUA_Gsac_GrayIm_TA.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

figufy(f2);
fname = [fig_dir 'MUA_Gsac_GrayIm_TA.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

%% ANALYZE LAMINAR DEPENDENCIES
load('/home/james/Analysis/bruce/saccade_modulation/layer_boundaries/layer_classification.mat')
boundary_enums = [boundary_class(:).Expt_num];

lem_mua_gsac_gray_avg = cat(2,all_lem_mua_data(:).gsac_gray_avg)';
un_lem_expts = unique(lem_mua_exptnum);
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
all_relloc = [];
for ee = 1:n_lem_expts
   cur_mua_set = find(lem_mua_exptnum == un_lem_expts(ee));
   cur_bound_info = find(boundary_enums == un_lem_expts(ee));
   cur_ub = boundary_class(cur_bound_info).ub;
   cur_lb = boundary_class(cur_bound_info).lb;
    
   gran_probes = (cur_ub+1):(cur_lb-1);
   supra_probes = 1:(cur_ub-1);
   infra_probes = (cur_lb+1):24;
   
   [~,cur_supt] = min(lem_mua_gsac_gray_avg(cur_mua_set,:),[],2);
   all_supt = cat(1,all_supt,cur_supt);
   cur_relpos = (1:24) - cur_ub;
   all_relloc = cat(2,all_relloc,cur_relpos);
   if ~isempty(gran_probes)
       gran_gsac(ee,:) = mean(lem_mua_gsac_gray_avg(cur_mua_set(gran_probes),:));
       [~,supt] = min(lem_mua_gsac_gray_avg(cur_mua_set(gran_probes),:),[],2);
       all_gran_supt = cat(1,all_gran_supt,supt);
       [~,enht] = max(lem_mua_gsac_gray_avg(cur_mua_set(gran_probes),:),[],2);
       all_gran_enht = cat(1,all_gran_enht,enht);
   end
   if ~isempty(supra_probes)
       supra_gsac(ee,:) = mean(lem_mua_gsac_gray_avg(cur_mua_set(supra_probes),:));
       [~,supt] = min(lem_mua_gsac_gray_avg(cur_mua_set(supra_probes),:),[],2);
       all_supra_supt = cat(1,all_supra_supt,supt);
       [~,enht] = max(lem_mua_gsac_gray_avg(cur_mua_set(supra_probes),:),[],2);
       all_supra_enht = cat(1,all_supra_enht,enht);
   end
   if ~isempty(infra_probes)
       infra_gsac(ee,:) = mean(lem_mua_gsac_gray_avg(cur_mua_set(infra_probes),:));
       [~,supt] = min(lem_mua_gsac_gray_avg(cur_mua_set(infra_probes),:),[],2);
       all_infra_supt = cat(1,all_infra_supt,supt);
       [~,enht] = max(lem_mua_gsac_gray_avg(cur_mua_set(infra_probes),:),[],2);
       all_infra_enht = cat(1,all_infra_enht,enht);
   end
end

f1 = figure(); hold on
h1=shadedErrorBar(tlags,nanmean(gran_gsac),nanstd(gran_gsac)/sqrt(n_lem_expts),{'color','b'});
h2=shadedErrorBar(tlags,nanmean(supra_gsac),nanstd(supra_gsac)/sqrt(n_lem_expts),{'color','r'});
h3=shadedErrorBar(tlags,nanmean(infra_gsac),nanstd(infra_gsac)/sqrt(n_lem_expts),{'color','k'});


trange = tlags(tlags > 0.02 & tlags < 0.125);
gran_hist = hist(tlags(all_gran_supt),trange);
sup_hist = hist(tlags(all_supra_supt),trange);
infra_hist = hist(tlags(all_infra_supt),trange);

f2 = figure();
% subplot(2,1,1);
hold on
plot(trange,gran_hist/sum(gran_hist));
plot(trange,sup_hist/sum(sup_hist),'r');
plot(trange,infra_hist/sum(infra_hist),'k')
xlabel('Suppresion timing');
ylabel('Relative freq');
legend('Granular','Supra-gran','Infra-gran');

trange = tlags(tlags > 0.1 & tlags < 0.25);
gran_hist = hist(tlags(all_gran_enht),trange);
sup_hist = hist(tlags(all_supra_enht),trange);
infra_hist = hist(tlags(all_infra_enht),trange);
% subplot(2,1,2);
% hold on
% plot(trange,gran_hist);
% plot(trange,sup_hist,'r');
% plot(trange,infra_hist,'k')
