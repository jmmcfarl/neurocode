clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')
addpath('G:\Code\general\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];

raw_Fs = 2016;
dsf = 40;
Fsd = raw_Fs/dsf;
lf_lcf = 0.05;
lf_hcf = 2;
Fsd2 = raw_Fs/8;
hf_lcf = 20;
hf_hcf = 80;
hf_smooth = 0.15;

%%
d=find_struct_field_vals(sess_data,'name','2007-04-05_C');
cdir = sess_data(d).directory;
cdir(1) = 'G';
cd(cdir)

load used_data lf4 wcv_minus_spike
load hsmm_state_seq4_seg_lf_4_28_10_v3
load hsmm_state_seq4_seg_hf_4_28_10_v3
load hsmm_state_seq4_seg_lfhf_4_28_10_v3_3

lf_bb = 0.05;
hf_bb = 20;
[lf4_bb,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd2,[lf_bb hf_bb]);
[lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd2,[lf_lcf lf_hcf]);
[lf4_hf,t_axis,Fs] = get_hf_features(lf4,raw_Fs,Fsd2,[hf_lcf hf_hcf],hf_smooth);

figure
plot(t_axis,lf4_bb)
hold on
plot(t_axis,lf4_hf,'k')
plot(t_axis,lf4_lf,'r')
[new_seg_inds] = resample_uds_seg_inds(hmm4.UDS_segs,hmm4.Fs,252,length(lf4_lf));
for i = 1:hmm4.Nsegs
    tseg = (new_seg_inds(i,1):new_seg_inds(i,2))/252;
    plot(tseg,hsmm_bbstate_seq4{i},'c','linewidth',2)
    plot(tseg,hsmm_bbstate_seq4_hf{i}+2.,'r','linewidth',2)
    plot(tseg,hsmm_bbstate_seq4_lfhf{i}-1,'g','linewidth',2)
end
% xlim([920 960])
% xlim([626 645])
% ylim([-2 3])
xl1 = [1150 1170];
xl2 = [200 220];
%%
dmean = [mean(hmm4_lfhf.state(1).meanfun{1}) mean(hmm4_lfhf.state(1).meanfun2{1})];
umean = [mean(hmm4_lfhf.state(2).meanfun{1}) mean(hmm4_lfhf.state(2).meanfun2{1})];
dcovar = hmm4_lfhf.state(1).var;
ucovar = hmm4_lfhf.state(2).var;

% cur_ts = find(t_axis > 910 & t_axis < 970);

figure
xl = [];
yl = [];
% cl = [-0.8 1.4];
cl = [];
gsmooth = 4;
islog = 0;
hist3_plot(lf4_lf,lf4_hf,[200 200],xl,yl,cl,islog,gsmooth);
hold on
h = plot_gaussian_ellipsoid(umean, ucovar, 2)
set(h,'color','w'); 
h = plot_gaussian_ellipsoid(umean, ucovar, 1)
set(h,'color','w'); 
h2 = plot_gaussian_ellipsoid(dmean, dcovar, 2)
set(h2,'color','k'); 
h2 = plot_gaussian_ellipsoid(dmean, dcovar, 1)
set(h2,'color','k'); 