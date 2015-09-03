
clear all
close all
cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\persistent_2010\')
addpath('G:\WC_Germany\persistent_revised')
addpath('G:\Code\WC_anal\general\')

raw_Fs = 2016;
niqf = raw_Fs/2;
dsf = 8;
Fsd = raw_Fs/dsf;

lcf = 0.05;
hcf = 20;

lcf2 = 1.5;
hcf2 = 8;

lcf3 = 2;
hcf3 = 5;
% niqf2 = Fsd/2;
% [b,a] = butter(2,[lcf3 hcf3]/niqf2);

winSize = 4;

d = 23;

cdir = sess_data(d).directory;
cdir(1) = 'G';
disp(sprintf('session %d',d))
cd(cdir);

load used_data lf8 wcv_minus_spike

load hsmm_state_seq_seg_lf_indivdesynch_dp
mp_state_seq = hsmm_bbstate_seq;
load hsmm_state_seq8_seg_lf_indivdesynch_dp
lfp_state_seq = hsmm_bbstate_seq8;

used_state_seq = lfp_state_seq;

lf8_f = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
lf8_f2 = get_lf_features(lf8,raw_Fs,Fsd,[lcf2 hcf2]);
lf8_f3 = get_lf_features(lf8,raw_Fs,Fsd,[lcf3 hcf3]);
lf8_f3a = get_lf_features_acaus(lf8,raw_Fs,Fsd,[lcf3 hcf3]);
lf8_f3c = get_lf_features_caus(lf8,raw_Fs,Fsd,[lcf3 hcf3]);
t = (1:length(lf8_f))/Fsd;
lf8_f3ap = angle(hilbert(lf8_f3a));
lf8_f3p = angle(hilbert(lf8_f3));
lf8_f3cp = angle(hilbert(lf8_f3c));

% wcv_f2c = get_lf_features_caus(wcv_minus_spike,raw_Fs,Fsd,[lcf2 hcf2]);
% wcv_f2a = get_lf_features_acaus(wcv_minus_spike,raw_Fs,Fsd,[lcf2 hcf2]);

%extract the index values of LF8 up transitions
% down_inds = [];
up_inds = [];
UDS_segs = (hmm.UDS_segs-1)*5+1;
clfp_state_seq = ones(UDS_segs(1,1)-1,1);
for i = 1:length(lfp_state_seq)
    up_trans = UDS_segs(i,1)+find(used_state_seq{i}(1:end-1) == 1 & used_state_seq{i}(2:end) == 2);
    down_trans = UDS_segs(i,1)+find(used_state_seq{i}(1:end-1) == 2 & used_state_seq{i}(2:end) == 1);
    %     up_trans(up_trans < down_trans(1)) = [];
    %     down_trans(down_trans > up_trans(end)) = [];
    %     cdown_inds(1,:) = down_trans;
    %     cdown_inds(2,:) = up_trans;
    %      down_inds = [down_inds cdown_inds];
    down_trans(down_trans < up_trans(1)) = [];
    up_trans(up_trans > down_trans(end)) = [];
    cup_inds(1,:) = up_trans;
    cup_inds(2,:) = down_trans;
    up_inds = [up_inds cup_inds];
    clear cup*
    clfp_state_seq= [clfp_state_seq; used_state_seq{i}];
    if length(lfp_state_seq) > i
        clfp_state_seq= [clfp_state_seq; ones(UDS_segs(i+1,1)-UDS_segs(i,2),1)];
    end
    
end

[urlid,urltime,urlamp,urlshift,urltau,urlerror,ut_90,ut_10] = ...
    get_lfp_wcv_sigmoid_fit_ut_12_4(up_inds(1,:),up_inds(2,:),lf8_f,t,Fsd);
[drlid,drltime,drlamp,drlshift,drltau,drlerror,dt_90,dt_10] = ...
    get_lfp_wcv_sigmoid_fit_dt_12_4(up_inds(2,1:end-1),up_inds(1,1:end-1),-lf8_f,t,Fsd);

%%

n_ups = size(up_inds,2);
min_dur = 1.3;
min_pow = 0.5;
temp_mats = [];
temp_mat = [];
n_points = 3000;
n_point = 1000;
for i = 2:n_ups
    
    prev_downseg = dt_90(i-1):ut_10(i);
    if length(prev_downseg)/Fsd > min_dur
            prev_down = detrend(zscore(lf8_f2(prev_downseg)));
        L = length(prev_down);
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
        Y = fft(prev_down,NFFT)/L;
        Y = 2*abs(Y(1:NFFT/2+1));
        f = Fsd/2*linspace(0,1,NFFT/2+1);
        [pk,pkloc] = max(Y);
        pkloc = f(pkloc);
        if pk > min_pow & pkloc >= 2 & pkloc <= 4.5
            new_ax = 1:1/pkloc:length(prev_down);
            cur_filt = lf8_f3(prev_downseg);
            cur_phase = lf8_f3cp(prev_downseg);
%             rel_t = linspace(-t(ut_10(i)),0,length(prev_down));
            dt = length(prev_downseg)/Fsd;
            cur_filti = interp1(1:length(prev_down),cur_filt,new_ax);
            cur_phasei = interp1(1:length(prev_down),cur_phase,new_ax);
%             rel_ti = linspace(0-t(ut_10(i))*pkloc,0,length(cur_filti));
            rel_pow = pk/sum(Y);
            cur_seg = (up_inds(1,i) - round(winSize/2*Fsd)):(up_inds(2,i)+round(winSize/2*Fsd));
            cur_up = up_inds(1,i):up_inds(2,i);
            cur_seg(cur_seg < 1 | cur_seg > length(lf8_f)) = [];
            
            
%             plot(t(cur_seg),lf8_f(cur_seg)), hold on
%             plot(t(cur_seg),lf8_f2(cur_seg)/3-3,'r')
%             plot(t(cur_seg),lf8_f3(cur_seg)/3-3,'k')
%             plot(t(cur_seg),lf8_f3p(cur_seg)/3-5,'r')
%             
%             plot(t(ut_10(i)),lf8_f(ut_10(i)),'go')
%             plot(t(ut_90(i)),lf8_f(ut_90(i)),'ko')
%             plot(t(up_inds(1,i)),0,'ro')
%             plot(t(up_inds(2,i)),0,'ro')
%             plot(t(prev_downseg),prev_down/3-3,'b')
%             title(strcat('pow: ',num2str(pk),' freq: ',num2str(pkloc),' rel: ',num2str(rel_pow)))
%             pause
%             clf
            
            used_chunk = cur_phasei;
            extra = length(used_chunk)-n_points;
            if length(used_chunk) > n_points
               used_chunk(1:extra) = [];
            else
                enans = nan(1,-extra);
                used_chunk = [enans used_chunk];
            end
            temp_mats = [temp_mats; used_chunk];
            used_chunk = cur_phase';
            extra = length(used_chunk)-n_point;
            if length(used_chunk) > n_point
               used_chunk(1:extra) = [];
            else
                enans = nan(1,-extra);
                used_chunk = [enans used_chunk];
            end
            temp_mat = [temp_mat; used_chunk];
        end
    end
end


%%

n_ups = size(up_inds,2);
min_dur = 0.8;
min_pow = 0.5;
temp_mats = [];
temp_mat = [];
n_points = 3000;
n_point = 1000;
for i = 1:n_ups-1
    
    cur_upseg = ut_90(i):dt_10(i);
    if length(cur_upseg)/Fsd > min_dur
            cur_up = detrend(zscore(lf8_f2(cur_upseg)));
        L = length(cur_up);
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
        Y = fft(cur_up,NFFT)/L;
        Y = 2*abs(Y(1:NFFT/2+1));
        f = Fsd/2*linspace(0,1,NFFT/2+1);
        [pk,pkloc] = max(Y);
        pkloc = f(pkloc);
        if pk > min_pow & pkloc >= 2 & pkloc <= 4.5
            new_ax = 1:1/pkloc:length(cur_up);
            cur_filt = lf8_f3(cur_upseg);
            cur_phase = lf8_f3p(cur_upseg);
            dt = length(cur_upseg)/Fsd;
            cur_filti = interp1(1:length(cur_upseg),cur_filt,new_ax);
            cur_phasei = interp1(1:length(cur_upseg),cur_phase,new_ax);
%             rel_pow = pk/sum(Y);
%             cur_seg = (up_inds(1,i) - round(winSize/2*Fsd)):(up_inds(2,i)+round(winSize/2*Fsd));
%             cur_up = up_inds(1,i):up_inds(2,i);
%             cur_seg(cur_seg < 1 | cur_seg > length(lf8_f)) = [];
            
            
%             plot(t(cur_seg),lf8_f(cur_seg)), hold on
%             plot(t(cur_seg),lf8_f2(cur_seg)/3-3,'r')
%             plot(t(cur_seg),lf8_f3(cur_seg)/3-3,'k')
%             plot(t(cur_seg),lf8_f3p(cur_seg)/3-5,'r')
%             
%             plot(t(ut_10(i)),lf8_f(ut_10(i)),'go')
%             plot(t(ut_90(i)),lf8_f(ut_90(i)),'ko')
%             plot(t(up_inds(1,i)),0,'ro')
%             plot(t(up_inds(2,i)),0,'ro')
%             plot(t(prev_downseg),prev_down/3-3,'b')
%             title(strcat('pow: ',num2str(pk),' freq: ',num2str(pkloc),' rel: ',num2str(rel_pow)))
%             pause
%             clf
            
            used_chunk = cur_phasei;
            extra = length(used_chunk)-n_points;
            if length(used_chunk) > n_points
               used_chunk(extra+1:end) = [];
            else
                enans = nan(1,-extra);
                used_chunk = [used_chunk enans];
            end
            temp_mats = [temp_mats; used_chunk];
            used_chunk = cur_phase';
            extra = length(used_chunk)-n_point;
            if length(used_chunk) > n_point
               used_chunk(extra+1:end) = [];
            else
                enans = nan(1,-extra);
                used_chunk = [used_chunk enans];
            end
            temp_mat = [temp_mat; used_chunk];
        end
    end
end


%%
lookback = round(2*Fsd);
used_ups = ut_10;
used_ups(used_ups < lookback) = [];
lf8_cmat = zeros(length(used_ups),lookback);
for i = 1:length(used_ups)
    lf8_cmat(i,:) = lf8_f3cp(used_ups(i)-lookback+1:used_ups(i));
end

lookforward = round(2*Fsd);
used_ups = ut_90;
used_ups(used_ups > length(lf8_f3c) - lookforward) = [];
lf8_cmat = zeros(length(used_ups),lookforward);
for i = 1:length(used_ups)
    lf8_cmat(i,:) = lf8_f3ap(used_ups(i):used_ups(i)+lookback-1);
end

lookback = round(2*Fsd);
used_ups = up_inds(2,:);
used_ups(used_ups < lookback) = [];
lf8_cmat = zeros(length(used_ups),lookback);
for i = 1:length(used_ups)
    lf8_cmat(i,:) = lf8_f3c(used_ups(i)-lookback+1:used_ups(i));
end

lookforward = round(2*Fsd);
used_ups = up_inds(2,:);
used_ups(used_ups > length(lf8_f3c) - lookforward) = [];
lf8_cmat = zeros(length(used_ups),lookforward);
for i = 1:length(used_ups)
    lf8_cmat(i,:) = lf8_f3c(used_ups(i):used_ups(i)+lookback-1);
end

lookforward = round(2*Fsd);
used_ups = round(rand(size(up_inds,2),1)*length(lf8_f3c));
used_ups(used_ups > length(lf8_f3c) - lookforward) = [];
lf8_cmat = zeros(length(used_ups),lookforward);
for i = 1:length(used_ups)
    lf8_cmat(i,:) = lf8_f3c(used_ups(i):used_ups(i)+lookback-1);
end
