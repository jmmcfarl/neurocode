
clear all
close all

%%
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir.mat

raw_Fs = 2016;
dsf = 2;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;

maxLag = round(2*Fsd);
backlag = 2*Fsd;
forwardlag = 2*Fsd;
lags = -maxLag:maxLag;
NLags = length(lags);

rate_sm = round(Fsd*0.05);

min_durs = 60; %minimum segment duration

params.Fs = Fsd;
params.fpass = [5 420];
params.tapers = [2 3];
lines = [50 100 150 200 250 300 350 400];
f_i = linspace(10,400,400);

%%
for d = 1:length(combined_dir)
    if ~isnan(hpc_mua(d))
        d
        cd(combined_dir{d})
        pwd
        load ./used_data
        if ctx_lfp(d) == 7
            lf8 = lf7;
        end
        [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
        wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);

%         lf2_rip = get_lf_features(lf2,raw_Fs,Fsd,[10 400]);
        lf3_rip = get_lf_features(lf3,raw_Fs,Fsd,[10 400]);
        lf7_rip = get_lf_features(lf7,raw_Fs,Fsd,[10 400]);
      
         load ./pa_hsmm_state_seq_combined_lf7
         mp_state_seq = hsmm_bbstate_seq;
         %         load ./pa_hsmm_state_seq7_combined
         load ./pa_hsmm_state_seq7_combined_lf7
         lfp_state_seq = hsmm_bbstate_seq7;
         
         [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,252,mp_state_seq);
         seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
        
         [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
         up_state_durs = (down_trans_inds-up_trans_inds)/Fsd;
        used_states = find(up_state_durs > 0.5);
        
         up_trans_inds = round(up_trans_inds*8/dsf);
         down_trans_inds = round(down_trans_inds*8/dsf);
         
%          ripS2 = mtspectrumc_weighted_avg(lf2_rip,up_trans_inds,...
%              down_trans_inds,params,f_i,lines,2);
         ripS3 = mtspectrumc_weighted_avg(lf3_rip,up_trans_inds(used_states),...
             down_trans_inds(used_states),params,f_i,lines,2);
         ripS7 = mtspectrumc_weighted_avg(lf7_rip,up_trans_inds(used_states),...
             down_trans_inds(used_states),params,f_i,lines,2);
         ripS3_d = mtspectrumc_weighted_avg(lf3_rip,down_trans_inds(1:end-1),...
             up_trans_inds(2:end),params,f_i,lines,2);
         ripS7_d = mtspectrumc_weighted_avg(lf7_rip,down_trans_inds(1:end-1),...
             up_trans_inds(2:end),params,f_i,lines,2);

        plot(f_i,ripS3,f_i,ripS7,'r')
        hold on
        plot(f_i,ripS3_d,'--',f_i,ripS7_d,'r--')
        pause
        close all
    end
end

