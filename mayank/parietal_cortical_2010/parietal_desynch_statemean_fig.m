clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\Code\smoothing\software')
addpath('G:\WC_Germany\hsmm_state_detection\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_individual

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];

raw_Fs = 2016;
dsf = 40;
Fsd = raw_Fs/dsf;
lf_lcf = 0.05;
lf_hcf = 20;
Fsd2 = raw_Fs/8;
% xl = [730 760];
xl = [1052 1070];

% d = 21;
% d = 15;
d=15
cdir = sess_data(d).directory;
cdir(1) = 'G';
cd(cdir)

load used_data lf4 wcv_minus_spike
% load hsmm_state_seq4_seg_lf_4_28_10_v1
load ./hsmm_state_seq4_seg_lf_4_5_2011.mat

[lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd,[lf_lcf lf_hcf]);
[state_means,state_t,obs_dist,obs_range] = ...
    locate_state_means(lf4_lf,hsmm4.windowSize,5,Fsd);

figure
imagesc(state_t,obs_range,log(obs_dist'));shading flat, caxis([-4 -0.4])
hold on
for i = 1:hsmm4.Nsegs
    tseg = (hsmm4.UDS_segs(i,1):hsmm4.UDS_segs(i,2))/Fsd;
    plot(tseg,hsmm4.state(1).meanfun{i},'w','linewidth',2)
    plot(tseg,hsmm4.state(2).meanfun{i},'k','linewidth',2)
end
ylim([-2.5 3])
yl = ylim;
line([xl(1) xl(1)],yl,'color','k')
line([xl(2) xl(2)],yl,'color','k')
for i = 1:size(desynch_times_lf4{d},1)
    line([desynch_times_lf4{d}(i,1) desynch_times_lf4{d}(i,1)],yl,'color','w')
    line([desynch_times_lf4{d}(i,2) desynch_times_lf4{d}(i,2)],yl,'color','w')
end

%%
[desynch_times,P_mp,P_lfp,f,t] = locate_desynch_times(wcv_minus_spike,lf4);
df = f(2)-f(1);
f_so = find(f > 0.2 & f < 1.5);
f_theta = find(f > 4);
max_so_lfp = max(10*log10(P_lfp(:,f_so)),[],2);
net_hf_lfp = zscore(trapz(10*log10(P_lfp(:,f_theta)),2)*df);

f_log = logspace(log10(f(2)),log10(f(end)),length(f));
P_lfp_log = interp1(f,P_lfp',f_log)';
f_range = find(f_log >= 0.1 & f_log <= 10);

Fig = figure;
clf
set(Fig,'PaperUnits','centimeters');
set(Fig, 'PaperSize', [30 20],'Paperposition',[0,0,(get(Fig,'PaperSize'))]);
% imagesc(t,f,10*log10(P_lfp'));shading flat;
imagesc(t,f_log(f_range),10*log10(P_lfp_log(:,f_range)'))
caxis([-40 0]);
% ylim([0 10])
set(gca,'yscale','log')
yl = ylim;
for i = 1:size(desynch_times_lf4{d},1)
    line([desynch_times_lf4{d}(i,1) desynch_times_lf4{d}(i,1)],yl,'color','w')
    line([desynch_times_lf4{d}(i,2) desynch_times_lf4{d}(i,2)],yl,'color','w')
end

figure
plot(t,max_so_lfp,'b'), hold on
xlim([t(1) t(end)])
line([t(1) t(end)],[-6 -6],'Color','k'), hold on
if size(desynch_times,1) > 0
    plot(desynch_times(:,1),ones(size(desynch_times(:,1)))*-5,'go')
    plot(desynch_times(:,2),ones(size(desynch_times(:,2)))*-5,'ro')
end
plot(t,net_hf_lfp,'r')
xlim([t(1) t(end)])
line([t(1) t(end)],[-2 -2],'Color','k')

%%
lf4_state_seq = hsmm_state_seq4;
figure
plot((1:length(lf4_lf))/Fsd,lf4_lf)
hold on
for i = 1:hsmm4.Nsegs
    tseg = (hsmm4.UDS_segs(i,1):hsmm4.UDS_segs(i,2))/Fsd;
    plot(tseg,hsmm4.state(1).meanfun{i},'r','linewidth',2)
    plot(tseg,hsmm4.state(2).meanfun{i},'k','linewidth',2)
    plot(tseg,lf4_state_seq{i},'g')
end
xlim(xl)
%%

