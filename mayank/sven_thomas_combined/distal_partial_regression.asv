
clear all
close all

%%
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
cd C:\WC_Germany\sven_thomas_combined\

load ./distal_dir
distal_dir = distal_dir(distal_usable);
mctx_lfp = ctx_lfp(distal_usable);
mhpc_lfp = hpc_lfp(distal_usable);
mhpc_mua = hpc_mua(distal_usable);

load ./distal_lec_dir.mat
distal_dir = [distal_dir distal_lec(usable_distal_lec)];
ctx_lfp = [mctx_lfp ctx_lfp(usable_distal_lec)];
hpc_lfp = [mhpc_lfp hpc_lfp(usable_distal_lec)];
hpc_mua = [mhpc_mua hpc_mua(usable_distal_lec)];
distal_mec = (1:length(distal_usable));
distal_lec = (length(distal_usable)+1):(length(distal_usable)+length(usable_distal_lec));
%%
raw_Fs = 2016;
dsf = 32;
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

int_width_range = [2 8]; %range of spike widths (in 32kHz samples) for interneuron spikes
pyr_width_range = [10 18]; %range of spike widths for pyramidal spikes
min_corr = 0.4; %minimum correlation coefficient with the corresponding spike template

%%
for d = 1:length(distal_dir)
    if ~isnan(hpc_mua(d))
        cd(distal_dir{d})
        pwd
        load ./used_data
        if ctx_lfp(d) == 7
            lf8 = lf7;
        end
        [lf8_lf,t_axis] = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
        wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
        if hpc_lfp(d) == 3
            load ./used_data lf3 lf5
            if ismember(d,old_data_inds)
                lf3 = lf3 + lf5; %redefine LF3 wrt gnd
            end
            hpc_lf = get_lf_features(lf3,raw_Fs,Fsd,[lcf hcf]);
            hpc_hf = get_hf_features(lf3,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
        else
            load ./used_data lf2
            hpc_lf = get_lf_features(lf2,raw_Fs,Fsd,[lcf hcf]);
            hpc_hf = get_hf_features(lf2,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
        end
        [desynch_times_lf8,desynch_inds,P_lf8,f,t] = locate_desynch_times_individual(lf8);
        
        if ~isnan(hpc_mua(d))
            %         load ./mua_data
            load ./sync_times.mat
            synct(length(wcv)+1:end) = [];
            synct_d = downsample(synct,dsf);
            %         synct_d = synct(1):s_fsd:synct(end);
            load ./mua_data3
            
            hpc_mua_times = mua_times{hpc_mua(d)};
            %         good_mua = find(mua_corr_pyr{hpc_mua(d)} > min_corr | mua_corr_int{hpc_mua(d)} > min_corr);
            %         hpc_mua_times_good = hpc_mua_times(good_mua);
            %         pyr_mua = find(mua_corr_pyr{hpc_mua(d)} > min_corr & mua_widths{hpc_mua(d)}' > pyr_width_range(1) & mua_widths{hpc_mua(d)}' < pyr_width_range(2));
            %         hpc_mua_times_pyr = hpc_mua_times(pyr_mua);
            %         int_mua = find(mua_corr_int{hpc_mua(d)} > min_corr & mua_widths{hpc_mua(d)}' > int_width_range(1) & mua_widths{hpc_mua(d)}' < int_width_range(2));
            %         hpc_mua_times_int = hpc_mua_times(int_mua);
            
            hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
            %         hpc_mua_times_good(hpc_mua_times_good > synct_d(end) | hpc_mua_times_good < synct_d(1)) = [];
            %         hpc_mua_times_pyr(hpc_mua_times_pyr > synct_d(end) | hpc_mua_times_pyr < synct_d(1)) = [];
            %         hpc_mua_times_int(hpc_mua_times_int > synct_d(end) | hpc_mua_times_int < synct_d(1)) = [];
            %
            hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
            mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
%             mua_rate = zscore(mua_rate');
            
            %         hpc_mua_rate_good =hist(hpc_mua_times_good,synct_d)*Fsd;
            %         hpc_mua_rate_good = jmm_smooth_1d_cor(hpc_mua_rate_good,rate_sm);
            %         mua_rate_good = zscore(hpc_mua_rate_good');
            %
            %         hpc_mua_rate_pyr =hist(hpc_mua_times_pyr,synct_d)*Fsd;
            %         hpc_mua_rate_pyr = jmm_smooth_1d_cor(hpc_mua_rate_pyr,rate_sm);
            %         mua_rate_pyr = zscore(hpc_mua_rate_pyr');
            %
            %         hpc_mua_rate_int =hist(hpc_mua_times_int,synct_d)*Fsd;
            %         hpc_mua_rate_int = jmm_smooth_1d_cor(hpc_mua_rate_int,rate_sm);
            %         mua_rate_int = zscore(hpc_mua_rate_int');
            
            %         [log_mua,offset(d)] = log_transform_sig(mua_rate);
            %         mua_rate = zscore(log_mua(:));
            if length(mua_rate) > length(t_axis)
                mua_rate = mua_rate(1:length(t_axis));
            end
        else
            mua_rate = nan(size(lf8_lf));
        end
        
        %     load ./pa_hsmm_state_seq_combined_fin.mat
        %     load ./pa_hsmm_state_seq7_combined_fin.mat
        load ./pa_hsmm_state_seq_combined_fin_nd.mat
        load ./pa_hsmm_state_seq7_combined_fin_nd.mat
        hsmm_bbstate_seq8 = hsmm_bbstate_seq7;
        
        bb_Fs = 252;
        temp_lf8 = get_lf_features(lf8,raw_Fs,252,[0.05 10]);
        [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,bb_Fs,length(temp_lf8));
        [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
        [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq8);
        
        
        up_trans_inds8 = round(up_trans_inds8/(dsf/8));
        up_trans_inds = round(up_trans_inds/(dsf/8));
        down_trans_inds8 = round(down_trans_inds8/(dsf/8));
        down_trans_inds = round(down_trans_inds/(dsf/8));
        
        up_trans_vec8 = zeros(size(lf8_lf));
        up_trans_vec8(up_trans_inds8) = 1;
        up_trans_vec = zeros(size(lf8_lf));
        up_trans_vec(up_trans_inds) = 1;
        down_trans_vec8 = zeros(size(lf8_lf));
        down_trans_vec8(down_trans_inds8) = 1;
        down_trans_vec = zeros(size(lf8_lf));
        down_trans_vec(down_trans_inds) = 1;
        
        Xwd = zeros(2*maxLag+1,length(lf8_lf));
        Xwu = zeros(2*maxLag+1,length(lf8_lf));
        for i = 1:maxLag
            Xwd(i,1:end-maxLag+i) = down_trans_vec(maxLag-i+1:end);
            Xwu(i,1:end-maxLag+i) = up_trans_vec(maxLag-i+1:end);
        end
        for i = (maxLag+1):(2*maxLag+1)
            Xwd(i,(i-maxLag+1):end) = down_trans_vec(1:end-i+maxLag);
            Xwu(i,(i-maxLag+1):end) = up_trans_vec(1:end-i+maxLag);
        end
        X8d = zeros(2*maxLag+1,length(lf8_lf));
        X8u = zeros(2*maxLag+1,length(lf8_lf));
        for i = 1:maxLag
            X8d(i,1:end-maxLag+i) = down_trans_vec8(maxLag-i+1:end);
            X8u(i,1:end-maxLag+i) = up_trans_vec8(maxLag-i+1:end);
        end
        for i = (maxLag+1):(2*maxLag+1)
            X8d(i,(i-maxLag+1):end) = down_trans_vec8(1:end-i+maxLag);
            X8u(i,(i-maxLag+1):end) = up_trans_vec8(1:end-i+maxLag);
        end
        
        %compute markers indicating segments of data to be used
        if ~isempty(desynch_times_lf8)
            desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,1)));
            desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,2)));
        else
            desynch_start = [];
            desynch_stop = [];
        end
        desynch_ind = zeros(size(lf8_lf));
        for i = 1:length(desynch_start)
            desynch_ind(desynch_start(i):desynch_stop(i)) = 1;
        end
        synch_starts = find(desynch_ind(1:end-1)==1 & desynch_ind(2:end)==0)+1;
        if desynch_ind(1) == 0
            synch_starts = [1; synch_starts];
        end
        synch_stops = find(desynch_ind(1:end-1)==0 & desynch_ind(2:end)==1)+1;
        if desynch_ind(end) == 0
            synch_stops = [synch_stops; length(lf8_lf)];
        end
        sMarkers = [synch_starts(:) synch_stops(:)];
        seg_durs = diff(sMarkers')'/Fsd;
        sMarkers(seg_durs < min_durs,:) = [];
        seg_durs(seg_durs < min_durs) = [];
        
        mua_rate = (mua_rate - mean(mua_rate(desynch_ind==0)))/std(mua_rate(desynch_ind==0));
        mua_rate = mua_rate(:);
        
        cnt = 0;
        clear B*
        
        if ~isnan(hpc_mua(d))
            %     if ~isnan(hpc_lfp(d))
            for i = 1:size(sMarkers,1)
                cnt = cnt+1;
                cl = length(sMarkers(i,1):sMarkers(i,2));
                [Bu(cnt,:),~,~,~,statsw(cnt,:)] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
                    [Xwu(:,sMarkers(i,1):sMarkers(i,2))' X8u(:,sMarkers(i,1):sMarkers(i,2))' ones(cl,1)]);
                [B8u(cnt,:),~,~,~,~] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
                    [X8u(:,sMarkers(i,1):sMarkers(i,2))' ones(cl,1)]);
                [Bwu(cnt,:),~,~,~,~] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
                    [Xwu(:,sMarkers(i,1):sMarkers(i,2))' ones(cl,1)]);
                [Bd(cnt,:),~,~,~,stats8(cnt,:)] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
                    [Xwd(:,sMarkers(i,1):sMarkers(i,2))' X8d(:,sMarkers(i,1):sMarkers(i,2))' ones(cl,1)]);
                [B8d(cnt,:),~,~,~,~] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
                    [X8d(:,sMarkers(i,1):sMarkers(i,2))' ones(cl,1)]);
                [Bwd(cnt,:),~,~,~,~] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
                    [Xwd(:,sMarkers(i,1):sMarkers(i,2))' ones(cl,1)]);
%                 [Bu(cnt,:),~,~,~,statsw(cnt,:)] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
%                     [Xwu(:,sMarkers(i,1):sMarkers(i,2))' X8u(:,sMarkers(i,1):sMarkers(i,2))']);
%                 [B8u(cnt,:),~,~,~,~] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
%                     [X8u(:,sMarkers(i,1):sMarkers(i,2))']);
%                 [Bwu(cnt,:),~,~,~,~] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
%                     [X8u(:,sMarkers(i,1):sMarkers(i,2))']);
%                 [Bd(cnt,:),~,~,~,stats8(cnt,:)] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
%                     [Xwd(:,sMarkers(i,1):sMarkers(i,2))' X8d(:,sMarkers(i,1):sMarkers(i,2))']);
%                 [B8d(cnt,:),~,~,~,~] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
%                     [X8d(:,sMarkers(i,1):sMarkers(i,2))']);
%                 [Bwd(cnt,:),~,R,~,~] = regress(mua_rate(sMarkers(i,1):sMarkers(i,2)),...
%                     [Xwd(:,sMarkers(i,1):sMarkers(i,2))']);
           end
            %
            
            %compute weighted averages (weighted by relative duration
            weighting = seg_durs/sum(seg_durs);
            
%             B_avg = sum(B.*repmat(weighting(:),1,size(B,2)),1);
            Bu_avg = sum(Bu.*repmat(weighting(:),1,size(Bu,2)),1);
            Bd_avg = sum(Bd.*repmat(weighting(:),1,size(Bd,2)),1);
            B8u_avg = sum(B8u.*repmat(weighting(:),1,size(B8u,2)),1);
            B8d_avg = sum(B8d.*repmat(weighting(:),1,size(B8d,2)),1);
            Bwu_avg = sum(Bwu.*repmat(weighting(:),1,size(Bwu,2)),1);
            Bwd_avg = sum(Bwd.*repmat(weighting(:),1,size(Bwd,2)),1);
            C3w_d(d,:) = Bd_avg(1:NLags);
            C3w_u(d,:) = Bu_avg(1:NLags);
            C3w_donly(d,:) = Bwd_avg(1:NLags);
            C3w_uonly(d,:) = Bwu_avg(1:NLags);
            C38_d(d,:) = Bd_avg(NLags+1:end-1);
            C38_u(d,:) = Bu_avg(NLags+1:end-1);
            C38_donly(d,:) = B8d_avg(1:NLags);
            C38_uonly(d,:) = B8u_avg(1:NLags);
            
            C3_d_back(d) = Bd_avg(end);
            C3_u_back(d) = Bu_avg(end);
            C3w_donly_back(d) = Bwd_avg(end);
            C3w_uonly_back(d) = Bwu_avg(end);
            C38_donly_back(d) = B8d_avg(end);
            C38_uonly_back(d) = B8u_avg(end);
%             C3w_dpart(d,:) = B_avg(1:NLags);
%             C3w_upart(d,:) = B_avg(NLags+1:2*NLags);
%             C38_upart(d,:) = B_avg(2*NLags+1:3*NLags);
%             C38_dpart(d,:) = B_avg(3*NLags+1:4*NLags);
                       
        else
%             C3w_dpart(d,:) = nan(1,NLags);
%             C3w_upart(d,:) = nan(1,NLags);
%             C38_upart(d,:) = nan(1,NLags);
%             C38_dpart(d,:) = nan(1,NLags);
            C3w_d(d,:) = nan(1,NLags);
            C3w_u(d,:) = nan(1,NLags);
            C38_d(d,:) = nan(1,NLags);
            C38_u(d,:) = nan(1,NLags);
            C3w_donly(d,:) = nan(1,NLags);
            C3w_uonly(d,:) = nan(1,NLags);
            C38_donly(d,:) = nan(1,NLags);
            C38_uonly(d,:) = nan(1,NLags);
        end
        %     clear mua_binned B B8 Bw
        clear mua_rate B*
    end
end

cd C:\WC_Germany\sven_thomas_combined\
save distal_partial_regress_mua_fin_nd C3* lags maxLag
% load ./mua_rate_data

%%
% uset = find(~isnan(hpc_mua));
% figure(1)
% shadedErrorBar(lags/Fsd,nanmean(C3w_upart(uset,:)),std(C3w_upart(uset,:))/sqrt(length(uset)),{'r'});
% hold on
% shadedErrorBar(lags/Fsd,nanmean(C38_upart(uset,:)),std(C38_upart(uset,:))/sqrt(length(uset)),{'k'});
% xlim([-1.5 1.5])
% xlabel('Time lag (s)','fontsize',14)
% ylabel('Regression Coefficient','fontsize',14)
% legend('MP','Ctx LFP')
% 
% figure(2)
% hold on
% shadedErrorBar(lags/Fsd,mean(C3w_dpart(uset,:)),std(C3w_dpart(uset,:))/sqrt(length(uset)),{'r'});
% shadedErrorBar(lags/Fsd,mean(C38_dpart(uset,:)),std(C38_dpart(uset,:))/sqrt(length(uset)),{'k'});
% xlim([-1.5 1.5])
% xlabel('Time lag (s)','fontsize',14)
% ylabel('Regression Coefficient','fontsize',14)
% legend('MP','Ctx LFP')

% load ./combined_partial_regress_mua_fin_nd
% load ./combined_dir_nd.mat
% uset = sort([l3mec l3lec]);
% all_cells = 1:61;
% l3mec = find(ismember(all_cells(uset),l3mec));
% l3lec = find(ismember(all_cells(uset),l3lec));
% combined_dir = combined_dir(uset);
% hpc_mua = hpc_mua(uset);
% l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));
% figure(1)
% shadedErrorBar(lags/Fsd,nanmean(C3w_upart(l3mec_m,:)),std(C3w_upart(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% hold on
% shadedErrorBar(lags/Fsd,nanmean(C38_upart(l3mec_m,:)),std(C38_upart(l3mec_m,:))/sqrt(length(l3mec_m)),{'b'});
% 
% figure(2)
% shadedErrorBar(lags/Fsd,nanmean(C3w_dpart(l3mec_m,:)),std(C3w_dpart(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% hold on
% shadedErrorBar(lags/Fsd,nanmean(C38_dpart(l3mec_m,:)),std(C38_dpart(l3mec_m,:))/sqrt(length(l3mec_m)),{'b'});


%%
cd C:\WC_Germany\sven_thomas_combined
load ./distal_dir
distal_dir = distal_dir(distal_usable);
mctx_lfp = ctx_lfp(distal_usable);
mhpc_lfp = hpc_lfp(distal_usable);
mhpc_mua = hpc_mua(distal_usable);

load ./distal_lec_dir.mat
distal_dir = [distal_dir distal_lec(usable_distal_lec)];
ctx_lfp = [mctx_lfp ctx_lfp(usable_distal_lec)];
hpc_lfp = [mhpc_lfp hpc_lfp(usable_distal_lec)];
hpc_mua = [mhpc_mua hpc_mua(usable_distal_lec)];
distal_mec = (1:length(distal_usable));
distal_lec = (length(distal_usable)+1):(length(distal_usable)+length(usable_distal_lec));

load ./distal_partial_regress_mua_fin_nd
uset = find(~isnan(hpc_mua));

figure(1)
subplot(2,1,1)
shadedErrorBar(lags/Fsd,mean(C3w_u(uset,:)),std(C3w_u(uset,:))/sqrt(length(uset)),{'r'});
title('Distal MUA','fontsize',16)
hold on
shadedErrorBar(lags/Fsd,mean(C38_u(uset,:)),std(C38_u(uset,:))/sqrt(length(uset)),{'k'});
xlim([-1.5 1.5])
% ylim([-0.6 0.6])
xlabel('Time lag (s)','fontsize',14)
ylabel('Regression Coefficient','fontsize',14)
% legend('MP','Ctx LFP')
ylim([-0.6 1])
set(gca,'fontsize',16,'fontname','arial')

figure(2)
subplot(2,1,1)
shadedErrorBar(lags/Fsd,mean(C3w_d(uset,:)),std(C3w_d(uset,:))/sqrt(length(uset)),{'r'});
title('Distal MUA','fontsize',16)
hold on
shadedErrorBar(lags/Fsd,mean(C38_d(uset,:)),std(C38_d(uset,:))/sqrt(length(uset)),{'k'});
xlim([-1.5 1.5])
% ylim([-0.6 0.7])
xlabel('Time lag (s)','fontsize',14)
ylabel('Regression Coefficient','fontsize',14)
% legend('MP','Ctx LFP')
ylim([-0.6 1])
set(gca,'fontsize',16,'fontname','arial')

%%
load ./combined_partial_regress_mua_fin_nd
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));
figure(1)
subplot(2,1,2)
title('Proximal MUA','fontsize',16)
hold on
shadedErrorBar(lags/Fsd,mean(C3w_u(l3mec_m,:)),std(C3w_u(l3mec_m,:))/sqrt(length(l3mec_m)),{'r'});
shadedErrorBar(lags/Fsd,mean(C38_u(l3mec_m,:)),std(C38_u(l3mec_m,:))/sqrt(length(l3mec_m)),{'k'});
ylim([-0.6 1])
xlim([-1.5 1.5])
xlabel('Time lag (s)','fontsize',14)
ylabel('Regression Coefficient','fontsize',14)
set(gca,'fontsize',16,'fontname','arial')

figure(2)
subplot(2,1,2)
title('Proximal MUA','fontsize',16)
hold on
shadedErrorBar(lags/Fsd,mean(C3w_d(l3mec_m,:)),std(C3w_d(l3mec_m,:))/sqrt(length(l3mec_m)),{'r'});
shadedErrorBar(lags/Fsd,mean(C38_d(l3mec_m,:)),std(C38_d(l3mec_m,:))/sqrt(length(l3mec_m)),{'k'});
ylim([-0.6 1])
xlim([-1.5 1.5])
xlabel('Time lag (s)','fontsize',14)
ylabel('Regression Coefficient','fontsize',14)
set(gca,'fontsize',16,'fontname','arial')

%%
cd C:\WC_Germany\sven_thomas_combined
load ./distal_dir
distal_dir = distal_dir(distal_usable);
mctx_lfp = ctx_lfp(distal_usable);
mhpc_lfp = hpc_lfp(distal_usable);
mhpc_mua = hpc_mua(distal_usable);

load ./distal_lec_dir.mat
distal_dir = [distal_dir distal_lec(usable_distal_lec)];
ctx_lfp = [mctx_lfp ctx_lfp(usable_distal_lec)];
hpc_lfp = [mhpc_lfp hpc_lfp(usable_distal_lec)];
hpc_mua = [mhpc_mua hpc_mua(usable_distal_lec)];
distal_mec = (1:length(distal_usable));
distal_lec = (length(distal_usable)+1):(length(distal_usable)+length(usable_distal_lec));

load ./distal_partial_regress_mua_fin_nd
uset = find(~isnan(hpc_mua));

figure(1)
subplot(2,1,1)
shadedErrorBar(lags/Fsd,mean(C3w_uonly(uset,:)),std(C3w_uonly(uset,:))/sqrt(length(uset)),{'r'});
title('Distal MUA','fontsize',16)
hold on
shadedErrorBar(lags/Fsd,mean(C38_uonly(uset,:)),std(C38_uonly(uset,:))/sqrt(length(uset)),{'k'});
xlim([-1.5 1.5])
% ylim([-0.6 0.6])
xlabel('Time lag (s)','fontsize',14)
ylabel('Regression Coefficient','fontsize',14)
% legend('MP','Ctx LFP')
ylim([-0.6 1.5])

figure(2)
subplot(2,1,1)
shadedErrorBar(lags/Fsd,mean(C3w_donly(uset,:)),std(C3w_donly(uset,:))/sqrt(length(uset)),{'r'});
title('Distal MUA','fontsize',16)
hold on
shadedErrorBar(lags/Fsd,mean(C38_donly(uset,:)),std(C38_donly(uset,:))/sqrt(length(uset)),{'k'});
xlim([-1.5 1.5])
% ylim([-0.6 0.7])
xlabel('Time lag (s)','fontsize',14)
ylabel('Regression Coefficient','fontsize',14)
% legend('MP','Ctx LFP')
ylim([-0.6 1.2])

load ./combined_partial_regress_mua_fin_nd
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec]);
all_cells = 1:61;
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));
figure(1)
subplot(2,1,2)
title('Proximal MUA','fontsize',16)
hold on
shadedErrorBar(lags/Fsd,mean(C3w_uonly(l3mec_m,:)),std(C3w_uonly(l3mec_m,:))/sqrt(length(l3mec_m)),{'r'});
shadedErrorBar(lags/Fsd,mean(C38_uonly(l3mec_m,:)),std(C38_uonly(l3mec_m,:))/sqrt(length(l3mec_m)),{'k'});
ylim([-0.6 1.5])
xlim([-1.5 1.5])
xlabel('Time lag (s)','fontsize',14)
ylabel('Regression Coefficient','fontsize',14)

figure(2)
subplot(2,1,2)
title('Proximal MUA','fontsize',16)
hold on
shadedErrorBar(lags/Fsd,mean(C3w_donly(l3mec_m,:)),std(C3w_donly(l3mec_m,:))/sqrt(length(l3mec_m)),{'r'});
shadedErrorBar(lags/Fsd,mean(C38_donly(l3mec_m,:)),std(C38_donly(l3mec_m,:))/sqrt(length(l3mec_m)),{'k'});
ylim([-0.6 1.2])
xlim([-1.5 1.5])
xlabel('Time lag (s)','fontsize',14)
ylabel('Regression Coefficient','fontsize',14)
