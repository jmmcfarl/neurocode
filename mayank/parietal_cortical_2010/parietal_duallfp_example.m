clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')
addpath('G:\Code\general\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8_3_03_10

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];

thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
sess_data = sess_data(thomas_el);


raw_Fs = 2016;
dsf = 40;
Fsd = raw_Fs/dsf;
lf_lcf = 0.05;
lf_hcf = 2;
Fsd2 = raw_Fs/16;
hf_lcf = 20;
hf_hcf = 80;
hf_smooth = 0.1;

%%
d=5;
cdir = sess_data(d).directory;
cdir(1) = 'G';
cd(cdir)

load used_data lf8 lf4
load hsmm_state_seq8_seg_lf_4_10_10
load hsmm_state_seq4_seg_lf_4_10_10
load hsmm_state_seq84_seg_lf_duallfp

lf_bb = 0.05;
hf_bb = 20;
[lf8_bb1,t_axis,Fs] = get_lf_features(lf8,raw_Fs,Fsd2,[lf_bb hf_bb]);
[lf8_bb2,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd2,[lf_bb hf_bb]);

figure
plot(t_axis,lf8_bb1)
hold on
plot(t_axis,lf8_bb2,'k')
for i = 1:hmm8.Nsegs
    trange = hmm8.UDS_segs(i,:)/Fsd;  
    tseg = linspace(trange(1),trange(2),length(hsmm_bbstate_seq8{i}));
    plot(tseg,hsmm_bbstate_seq8{i},'c','linewidth',2)
     trange = hmm4.UDS_segs(i,:)/Fsd;  
    tseg = linspace(trange(1),trange(2),length(hsmm_bbstate_seq4{i}));
   plot(tseg,hsmm_bbstate_seq4{i}+0.5,'r','linewidth',2)
     trange = hmm84.UDS_segs(i,:)/Fsd;  
    tseg = linspace(trange(1),trange(2),length(hsmm_bbstate_seq84{i}));
    plot(tseg,hsmm_bbstate_seq84{i}-0.5,'g','linewidth',2)
end
xlim([295 308])
ylim([-2 3])
%%
dmean = [mean(hmm84.state(1).meanfun{1}) mean(hmm84.state(1).meanfun2{1})];
umean = [mean(hmm84.state(2).meanfun{1}) mean(hmm84.state(2).meanfun2{1})];
dcovar = hmm84.state(1).var;
ucovar = hmm84.state(2).var;

lf_lcf = 0.05;
lf_hcf = 2;
[lf8_lf1,t_axis,Fs] = get_lf_features(lf8,raw_Fs,Fsd2,[lf_lcf lf_hcf]);
[lf8_lf2,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd2,[lf_lcf lf_hcf]);

xl = [];
yl = [-3 3];
cl = [-1 1.5];
% cl = [];
gsmooth = 4;
islog = 1;
hist3_plot(lf8_lf1,lf8_lf2,[200 200],xl,yl,cl,islog,gsmooth);
hold on
h = plot_gaussian_ellipsoid(umean, ucovar, 2)
set(h,'color','w'); 
h = plot_gaussian_ellipsoid(umean, ucovar, 1)
set(h,'color','w'); 
h2 = plot_gaussian_ellipsoid(dmean, dcovar, 2)
set(h2,'color','k'); 
h2 = plot_gaussian_ellipsoid(dmean, dcovar, 1)
set(h2,'color','k'); 