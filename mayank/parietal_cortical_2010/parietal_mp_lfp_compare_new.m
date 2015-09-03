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
% desynch_times_mp(interneurons) = [];
% desynch_times_lf8(interneurons) = [];
% desynch_times_lf4(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);

% thom_pfc = [thom_pfc 43 44];
sess_data = sess_data(thom_pfc);
% desynch_times_mp = desynch_times_mp(thom_pfc);
% desynch_times_lf4 = desynch_times_lf4(thom_pfc);
% desynch_times_lf8 = desynch_times_lf8(thom_pfc);

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;


n = length(sess_data);
for d = 1:9
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
    
    load ./hsmm_state_seq_seg_lf_4_5_2011
    mp_hsmm_state_seq = hsmm_bbstate_seq;    
    mp_hmm_state_seq = hmm_bbstate_seq;
    load ./hsmm_state_seq_seg_lffm_4_5_2011
    mp_hsmm_state_seq_lffm = hsmm_bbstate_seq;
    load ./fm_state_seq_lf_4_11_2011
    mp_fm_state_seq = fm_state_seq;
    mp_fm_state_seqz = fm_state_seqz;
    mp_fm_state_seqx = fm_state_seqx;
    
    load ./hsmm_state_seq4_seg_lf_4_5_2011
    lf4_hsmm_state_seq = hsmm_bbstate_seq4;    
    lf4_hmm_state_seq = hmm_bbstate_seq4;
    load ./hsmm_state_seq4_seg_hf_4_5_2011
    lf4_hsmm_state_seq_hf = hsmm_bbstate_seq4_hf;
    load ./hsmm_state_seq4_seg_lfhf_4_5_2011
    lf4_hsmm_state_seq_lfhf = hsmm_bbstate_seq4_lfhf;
    load ./hsmm_state_seq4_seg_lffm_4_5_2011
    lf4_hsmm_state_seq_lffm = hsmm_bbstate_seq4;
    load ./fm_state_seq4_lf_4_11_2011
    lf4_fm_state_seq = fm_state_seq4;
    lf4_fm_state_seqz = fm_state_seqz4;
    lf4_fm_state_seqx = fm_state_seqx4;
    
    [hsmm_ham(d),hsmm_coinup(d),hsmm_coindown(d),hsmm_fp(d),hsmm_fn(d),hsmm_tep(d)] = compute_state_seq_seg_hamdist_varuds(mp_hsmm_state_seq,lf4_hsmm_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    [hsmmhf_ham(d),hsmmhf_coinup(d),hsmmhf_coindown(d),hsmmhf_fp(d),hsmmhf_fn(d),hsmmhf_tep(d)] = compute_state_seq_seg_hamdist_varuds(mp_hsmm_state_seq,lf4_hsmm_state_seq_hf,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    [hsmmlfhf_ham(d),hsmmlfhf_coinup(d),hsmmlfhf_coindown(d),hsmmlfhf_fp(d),hsmmlfhf_fn(d),hsmmlfhf_tep(d)] = compute_state_seq_seg_hamdist_varuds(mp_hsmm_state_seq,lf4_hsmm_state_seq_lfhf,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    [hmm_ham(d),hmm_coinup(d),hmm_coindown(d),hmm_fp(d),hmm_fn(d),hmm_tep(d)] = compute_state_seq_seg_hamdist_varuds(mp_hmm_state_seq,lf4_hmm_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    [hsmmfm_ham(d),hsmmfm_coinup(d),hsmmfm_coindown(d),hsmmfm_fp(d),hsmmfm_fn(d),hsmmfm_tep(d)] = compute_state_seq_seg_hamdist_varuds(mp_hsmm_state_seq_lffm,lf4_hsmm_state_seq_lffm,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    [fm_ham(d),fm_coinup(d),fm_coindown(d),fm_fp(d),fm_fn(d),fm_tep(d)] = compute_state_seq_seg_hamdist_varuds(mp_fm_state_seq,lf4_fm_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    [fmz_ham(d),fmz_coinup(d),fmz_coindown(d),fmz_fp(d),fmz_fn(d),fmz_tep(d)] = compute_state_seq_seg_hamdist_varuds(mp_fm_state_seqz,lf4_fm_state_seqz,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    [fmx_ham(d),fmx_coinup(d),fmx_coindown(d),fmx_fp(d),fmx_fn(d),fmx_tep(d)] = compute_state_seq_seg_hamdist_varuds(mp_fm_state_seqx,lf4_fm_state_seqx,hmm.UDS_segs,hmm4.UDS_segs,Fsd);

    load ./used_data wcv
    wcv_d = downsample(wcv,8);
    [new_seg_inds4] = resample_uds_seg_inds(hmm4.UDS_segs,hmm4.Fs,252,length(wcv_d));

    %%
      [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,mp_hsmm_state_seq);
    [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,lf4_hsmm_state_seq);
    [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
    [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
    hsmm_sfn(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
    hsmm_sfp(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);

%           [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,mp_hsmm_state_seq);
    [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,lf4_hsmm_state_seq_hf);
    [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
    [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
    hsmmhf_sfn(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
    hsmmhf_sfp(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);

    %           [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,mp_hsmm_state_seq);
    [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,lf4_hsmm_state_seq_lfhf);
    [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
    [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
    hsmmlfhf_sfn(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
    hsmmlfhf_sfp(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);

          [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,mp_hsmm_state_seq_lffm);
    [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,lf4_hsmm_state_seq_lffm);
    [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
    [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
    hsmmfm_sfn(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
    hsmmfm_sfp(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);

    [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,mp_hmm_state_seq);
    [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,lf4_hmm_state_seq);
    [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
    [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
    hmm_sfn(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
    hmm_sfp(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);
    
    [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seq);
    [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seq4);
    [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
    [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
    fm_sfn(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
    fm_sfp(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);

    [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seqz);
    [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seqz4);
    [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
    [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
    fmz_sfn(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
    fmz_sfp(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);

        [mp_uptrans,mp_downtrans,mp_upsegnums,mp_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seqx);
    [lf4_uptrans,lf4_downtrans,lf4_upsegnums,lf4_downsegnums] = compute_state_transitions_seg(new_seg_inds4,fm_state_seqx4);
    [corresp_lf4_upinds,corresp_lf4_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        mp_uptrans,mp_downtrans,lf4_uptrans,lf4_downtrans);
    [corresp_mp_upinds,corresp_mp_downinds] = greedy_find_corresponding_ncx_state_transitions_simp(...
        lf4_uptrans,lf4_downtrans,mp_uptrans,mp_downtrans);
    fmx_sfn(d) = sum(isnan(corresp_lf4_upinds))/length(mp_uptrans);
    fmx_sfp(d) = sum(isnan(corresp_mp_upinds))/length(lf4_uptrans);

end

cd G:\WC_Germany\parietal_cortical_2010
save parietal_mp_lfp_compare_data_4_11_2011 

hsmm_sep = (hsmm_sfn + hsmm_sfp)/2;
hsmmfm_sep = (hsmmfm_sfn + hsmmfm_sfp)/2;
hsmmhf_sep = (hsmmhf_sfn + hsmmhf_sfp)/2;
hsmmlfhf_sep = (hsmmlfhf_sfn + hsmmlfhf_sfp)/2;
hmm_sep = (hmm_sfn + hmm_sfp)/2;
fm_sep = (fm_sfn + fm_sfp)/2;
fmz_sep = (fmz_sfn + fmz_sfp)/2;
fmx_sep = (fmx_sfn + fmx_sfp)/2;

hsmm_coinIn = (hsmm_coinup+hsmm_coindown)/2;
hsmmfm_coinIn = (hsmmfm_coinup+hsmmfm_coindown)/2;
hsmmhf_coinIn = (hsmmhf_coinup+hsmmhf_coindown)/2;
hsmmlfhf_coinIn = (hsmmlfhf_coinup+hsmmlfhf_coindown)/2;
hmm_coinIn = (hmm_coinup+hmm_coindown)/2;
fm_coinIn = (fm_coinup+fm_coindown)/2;
fmz_coinIn = (fmz_coinup+fmz_coindown)/2;
fmx_coinIn = (fmx_coinup+fmx_coindown)/2;

hsmm_ep = (hsmm_fp + hsmm_fn)/2;
hsmmfm_ep = (hsmmfm_fp + hsmmfm_fn)/2;
hsmmhf_ep = (hsmmhf_fp + hsmmhf_fn)/2;
hsmmlfhf_ep = (hsmmlfhf_fp + hsmmlfhf_fn)/2;
hmm_ep = (hmm_fp + hmm_fn)/2;
fm_ep = (fm_fp + fm_fn)/2;
fmz_ep = (fmz_fp + fmz_fn)/2;
fmx_ep = (fmx_fp + fmx_fn)/2;

%%
cd G:\WC_Germany\parietal_cortical_2010\
load ./roc_data.mat
addpath('G:\Code\eplots\')
X = [nanmean(fp); nanstd(fp)/3];
Y = [nanmean(fn); nanstd(fn)/3];
figure
set(gca,'fontsize',14,'fontname','arial')
% eplot(X',Y','b',0.5);
% plot(nanmean(fp),nanmean(fn),'linewidth',2)
hold on
plot(hsmm_fp,hsmm_fn,'ro','markersize',10)
% plot(hmm_fp,hmm_fn,'co')
plot(fm_fp,fm_fn,'k*','markersize',10)
plot(fmz_fp,fmz_fn,'g*','markersize',10)
% plot(fmx_fp,fmx_fn,'c*')
ylim([0 0.25])
xlim([0 0.35])
X = [nanmean(hsmm_fp); nanstd(hsmm_fp)/3];
Y = [nanmean(hsmm_fn); nanstd(hsmm_fn)/3];
eplot(X',Y','r',2)
X = [nanmean(fm_fp); nanstd(fm_fp)/3];
Y = [nanmean(fm_fn); nanstd(fm_fn)/3];
eplot(X',Y','k',2)
X = [nanmean(fmz_fp); nanstd(fmz_fp)/3];
Y = [nanmean(fmz_fn); nanstd(fmz_fn)/3];
eplot(X',Y','g',2)
% X = [nanmean(hmm_fp); nanstd(hmm_fp)/3];
% Y = [nanmean(hmm_fn); nanstd(hmm_fn)/3];
% eplot(X',Y','c',2)
xlabel('False UP probability','fontsize',16,'fontname','arial')
ylabel('False DOWN probability','fontsize',16,'fontname','arial')

%%
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% plot((hsmm_fp-fm_fp)./fm_fp*100,(hsmm_fn-fm_fn)./fm_fn*100,'ro')
% plot((hsmm_fp-fmz_fp)./fmz_fp*100,(hsmm_fn-fmz_fn)./fmz_fn*100,'ko')
% plot((hsmm_fp-fmx_fp)./fmx_fp*100,(hsmm_fn-fmx_fn)./fmx_fn*100,'go')
% line([-100 500],[0 0],'color','k')
% line([0 0],[-100 300],'color','k')
% line([-300 300],[300 -300],'color','k')
% ylim([-100 300])
% xlim([-100 500])
% xlabel('False UP probability','fontsize',16,'fontname','arial')
% ylabel('False DOWN probability','fontsize',16,'fontname','arial')

%%
% figure
% set(gca,'fontsize',14,'fontname','arial')
% hold on
% plot((hsmm_sfp-fm_sfp)./fm_sfp*100,(hsmm_sfn-fm_sfn)./fm_sfn*100,'ro')
% plot((hsmm_sfp-fmz_sfp)./fmz_sfp*100,(hsmm_sfn-fmz_sfn)./fmz_sfn*100,'ko')
% plot((hsmm_sfp-fmx_sfp)./fmx_sfp*100,(hsmm_sfn-fmx_sfn)./fmx_sfn*100,'go')
% line([-100 500],[0 0],'color','k')
% line([0 0],[-100 300],'color','k')
% line([-300 300],[300 -300],'color','k')
% ylim([-100 300])
% xlim([-100 500])
% xlabel('False UP probability','fontsize',16,'fontname','arial')
% ylabel('False DOWN probability','fontsize',16,'fontname','arial')

%%
figure
set(gca,'fontsize',14,'fontname','arial')
hold on
plot(hsmm_sfp,hsmm_sfn,'ro','markersize',10)
% plot(hmm_sfp,hmm_sfn,'co')
% plot(hsmmfm_sfp,hsmmfm_sfn,'yo')
plot(fm_sfp,fm_sfn,'k*','markersize',10)
plot(fmz_sfp,fmz_sfn,'g*','markersize',10)
ylim([0 0.2])
xlim([0 0.08])
X = [nanmean(hsmm_sfp); nanstd(hsmm_sfp)/3];
Y = [nanmean(hsmm_sfn); nanstd(hsmm_sfn)/3];
eplot(X',Y','r',2)
% X = [nanmean(hmm_sfp); nanstd(hmm_sfp)/3];
% Y = [nanmean(hmm_sfn); nanstd(hmm_sfn)/3];
% eplot(X',Y','c',2)
% X = [nanmean(hsmmfm_sfp); nanstd(hsmmfm_sfp)/3];
% Y = [nanmean(hsmmfm_sfn); nanstd(hsmmfm_sfn)/3];
% eplot(X',Y','y',2)
X = [nanmean(fm_sfp); nanstd(fm_sfp)/3];
Y = [nanmean(fm_sfn); nanstd(fm_sfn)/3];
eplot(X',Y','k',2)
X = [nanmean(fmz_sfp); nanstd(fmz_sfp)/3];
Y = [nanmean(fmz_sfn); nanstd(fmz_sfn)/3];
eplot(X',Y','g',2)
xlabel('False state probability','fontsize',16,'fontname','arial')
ylabel('Missed state probability','fontsize',16,'fontname','arial')

%%
dhmm = -100*(hsmm_coinIn - hmm_coinIn)./hsmm_coinIn;
dhsmmfm = -100*(hsmm_coinIn - hsmmfm_coinIn)./hsmm_coinIn;
dhsmmhf = -100*(hsmm_coinIn - hsmmhf_coinIn)./hsmm_coinIn;
dhsmmlfhf = -100*(hsmm_coinIn - hsmmlfhf_coinIn)./hsmm_coinIn;
dfm = -100*(hsmm_coinIn - fm_coinIn)./hsmm_coinIn;
dfmz = -100*(hsmm_coinIn - fmz_coinIn)./hsmm_coinIn;

figure
boxplot([dhmm' dhsmmfm' dfm' dfmz'])
xl = xlim();
line(xl,[0 0],'color','k')

%%
dhmm = -100*(hsmm_ep - hmm_ep)./hsmm_ep;
dhsmmfm = -100*(hsmm_ep - hsmmfm_ep)./hsmm_ep;
dhsmmhf = -100*(hsmm_ep - hsmmhf_ep)./hsmm_ep;
dhsmmlfhf = -100*(hsmm_ep - hsmmlfhf_ep)./hsmm_ep;
dfm = -100*(hsmm_ep - fm_ep)./hsmm_ep;
dfmz = -100*(hsmm_ep - fmz_ep)./hsmm_ep;
dfmx = -100*(hsmm_ep - fmx_ep)./hsmm_ep;

figure
set(gca,'fontname','arial','fontsize',14)
boxplot([dhmm' dhsmmfm' dfm' dfmz'])
xl = xlim();
line(xl,[0 0],'color','k')
ylabel('Percent change','fontsize',16,'fontname','arial')

%%
dhmm = -100*(hsmm_sep - hmm_sep)./hsmm_sep;
dhsmmfm = -100*(hsmm_sep - hsmmfm_sep)./hsmm_sep;
dhsmmhf = -100*(hsmm_sep - hsmmhf_sep)./hsmm_sep;
dhsmmlfhf = -100*(hsmm_sep - hsmmlfhf_sep)./hsmm_sep;
dfm = -100*(hsmm_sep - fm_sep)./hsmm_sep;
dfmz = -100*(hsmm_sep - fmz_sep)./hsmm_sep;

figure
set(gca,'fontname','arial','fontsize',14)
boxplot([dhmm' dhsmmfm' dfm' dfmz'])
xl = xlim();
line(xl,[0 0],'color','k')
ylabel('Percent change','fontsize',16,'fontname','arial')
