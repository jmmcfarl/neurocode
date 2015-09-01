clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_individual

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf4(interneurons) = [];

raw_Fs = 2016;
dsf = 40;
Fsd = raw_Fs/dsf;
lf_lcf = 0.05;
lf_hcf = 2;
Fsd2 = raw_Fs/16;

%%
% d = 8;
% d=9
d=27
% d=31 %skewed dist
% d=34;  
% d=38;

cdir = sess_data(d).directory;
cdir(1) = 'G';
cd(cdir)

load used_data lf4 wcv_minus_spike
% load hsmm_state_seq4_seg_lf_4_28_10_v3
load hsmm_state_seq4_seg_lf_4_5_2011

% load fixmean_state_seqs_4_10_10
% load fixmean_state_seqs4_4_10_10
% load fixmean_state_seqs4_zm_4_26_10_v2
load fm_state_seq4_lf_4_5_2011
fixmean_state_seq4 = fm_state_seqz4;

lf_lcf = 0.05;
lf_hcf = 20;
f_lcf = 2;
[lf4_bb,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd2,[lf_lcf lf_hcf]);
[lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd2,[lf_lcf f_lcf]);
% [wcv_bb,t_axis,Fs] = get_lf_features(wcv_minus_spike,raw_Fs,Fsd2,[lf_lcf lf_hcf]);
figure
plot(t_axis,lf4_bb)
hold on
% plot(t_axis,lf4_lf,'k')

for i = 1:hmm4.Nsegs
    trange = hmm4.UDS_segs(i,:)/Fsd;  
    tseg = linspace(trange(1),trange(2),length(hsmm_bbstate_seq4{i}));
    plot(tseg,hsmm_bbstate_seq4{i},'k','linewidth',2)
    tsegf = linspace(trange(1),trange(2),length(fixmean_state_seq4{i}));
    plot(tsegf,fixmean_state_seq4{i}-1.2,'r','linewidth',2)
        tseg2 = (hmm4.UDS_segs(i,1):hmm4.UDS_segs(i,2))/hmm4.Fs;
        plot(tseg2,hmm4.state(1).meanfun{i},'g')
        plot(tseg2,hmm4.state(2).meanfun{i},'c')

    line([tseg(1) tseg(end)],[fhmmz4.threshold fhmmz4.threshold],'color','k')
end
xl1 = [250 265]; 
xl2 = [312.5 338];
% xl2 = [760 830];
% xl2 = [759 815];
xl3 = [780 825];
xl4 = [1020 1061];
xlim(xl2)
% xlim([1420 1460])

%%
figure
x_axis = linspace(-2,3,1000);
% ksdensity(lf8_bb)
ksdensity(lf4_lf,x_axis)
hold on
y1 = fhmm4.priors(1)/sqrt(2*pi*fhmm4.state_vars(1))*exp(-(x_axis-fhmm4.state_means(1)).^2/(2*fhmm4.state_vars(1)));
y2 = fhmm4.priors(2)/sqrt(2*pi*fhmm4.state_vars(2))*exp(-(x_axis-fhmm4.state_means(2)).^2/(2*fhmm4.state_vars(2)));
plot(x_axis,y1,'r')
plot(x_axis,y2,'g')
ylim([0 0.6])
yl = ylim;
% load fixmean_state_seqs4_zm_4_26_10_v2
load ./fm_state_seq4_lf_4_5_2011.mat
line([fhmm4.threshold fhmm4.threshold],yl,'color','k')
% load fixmean_state_seqs4_4_26_10_v2
line([fhmmz4.threshold fhmmz4.threshold],yl,'color','k')
xlim([-2 2.5])
%%
[lf4_lf,t_axis,Fs] = get_lf_features(lf4,raw_Fs,Fsd,[lf_lcf f_lcf]);

[state_means,state_t,obs_dist,obs_range] = ...
    locate_state_means(lf4_lf,50,2.5,Fsd);
figure
imagesc(state_t,obs_range,log(obs_dist'));shading flat, caxis([-4 -0.4])
hold on
for i = 1:hmm4.Nsegs
    tseg = (hmm4.UDS_segs(i,1):hmm4.UDS_segs(i,2))/Fsd;
    plot(tseg,hsmm4.state(1).meanfun{i},'w','linewidth',2)
    plot(tseg,hsmm4.state(2).meanfun{i},'k','linewidth',2)
    line([tseg(1) tseg(end)],[fhmm4.threshold fhmm4.threshold],'color','k')
    line([tseg(1) tseg(end)],[fhmmz4.threshold fhmmz4.threshold],'color','k')
end
% xlim([800 1650])
% ylim([-1.9 3])
caxis([-3 -0.5])
ylim([-2 2.5])

yl = ylim;
line([xl1(1) xl1(1)],yl,'color','k')
line([xl1(2) xl1(2)],yl,'color','k')
line([xl2(1) xl2(1)],yl,'color','k')
line([xl2(2) xl2(2)],yl,'color','k')
line([xl3(1) xl3(1)],yl,'color','k')
line([xl3(2) xl3(2)],yl,'color','k')
line([xl4(1) xl4(1)],yl,'color','k')
line([xl4(2) xl4(2)],yl,'color','k')