clear all
close all
%%
load G:\WC_Germany\overall_EC\overall_allcells_dir
addpath('G:\Code\WC_anal\general\')
addpath('G:\WC_Germany\Overall_EC\')
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\hsmm_state_detection\')

drive_letter = 'G';

%%
dsf = 16;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,40/niqf,'low');

params.Fs = Fsd;
params.err = [2 0.05];
params.tapers = [7 13];
W = 0.1;
params.fpass = [0 8];
winlength = 50;
winslide = 50;
movingwin = [winlength winslide];

f_i = linspace(W,10,250);

for d = 1:200
    
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf8 lf3 wcv_minus_spike
    
    %     if length(sess_data(d).gains) == 8
    %         lf8_gain = sess_data(d).gains(8);
    %         lf3_gain = sess_data(d).gains(3);
    %     else
    %         lf8_gain = sess_data(d).gains(4);
    %         lf3_gain = sess_data(d).gains(3);
    %     end
    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = downsample(wcv_d,dsf);
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_d,dsf);
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = downsample(lf3_d,dsf);
%     lf2_d = filtfilt(b,a,lf2);
%     lf2_d = downsample(lf2_d,dsf);
%     lf5_d = filtfilt(b,a,lf5);
%     lf5_d = downsample(lf5_d,dsf);
%     lf2_r = lf2_d - lf5_d;
 
    lf3_hf = get_hf_features(lf3,2016,Fsd,[20 80],0.05);

%     lf2_r = lf2-lf5;
%     lf2_r = downsample(filtfilt(b,a,lf2_r),dsf);
    
    t_axis = (1:length(lf3_d))/Fsd;
    load ec_hmm_state_seq8
    [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,Fsd,length(lf8_d));
%     seg_durs = diff(new_seg_inds,[],2)/Fsd;
%     too_long = find(seg_durs > 1500);
%     
%     if ~isempty(too_long)
%         disp('TOO LONG!')
%         pause
%     end
    
    
%     [Cw8,Phiw8,Sw8,Pww_88,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf8_d],W,params,new_seg_inds);
%     [Cw3,Phiw3,Sw3,Pww_33,~,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf3_d],W,params,new_seg_inds);
%     [C83,Phi83,S83,P88_33,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf3_d],W,params,new_seg_inds);
%     [Cw3h,Phiw3h,Sw3h,Pww_33h,~,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf3_hf],W,params,new_seg_inds);
%     [C83h,Phi83h,S83h,P88_33h,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf3_hf],W,params,new_seg_inds);

    [Cw8,Phiw8,Sw8,Pww_88,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_d lf8_d],movingwin,params,new_seg_inds);
    [Cw3,Phiw3,Sw3,Pww_33,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_d lf3_d],movingwin,params,new_seg_inds);
    [C83,Phi83,S83,P88_33,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_d lf3_d],movingwin,params,new_seg_inds);
    [Cw3h,Phiw3h,Sw3h,Pww_33h,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([wcv_d lf3_hf],movingwin,params,new_seg_inds);
    [C83h,Phi83h,S83h,P88_33h,f,ConfC,dof] = coherencyc_unequal_length_trials_jmm([lf8_d lf3_hf],movingwin,params,new_seg_inds);
  
    
    %     [Cw2r,Phiw2r,Sw2r,Pww_22r,~,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf2_r],W,params,new_seg_inds);
%     [C82r,Phi82r,S82r,P88_22r,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf2_r],W,params,new_seg_inds);
%     [Cw2,Phiw2,Sw2,Pww_22,~,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([wcv_d lf2_d],W,params,new_seg_inds);
%     [C82,Phi82,S82,P88_22,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf8_d lf2_d],W,params,new_seg_inds);
%     [C52r,Phi52r,S52r,P55_22r,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf5_d lf2_r],W,params,new_seg_inds);
%     [C53,Phi53,S53,P55_33,f,ConfC,dof] = coherencyc_unequal_length_trials_fixedW([lf5_d lf3_d],W,params,new_seg_inds);
    
    S38 = conj(S83);
%     S2r8 = conj(S82r);
%     S28 = conj(S82);
    S8w = conj(Sw8);
    S3w = conj(Sw3);
%     S2w = conj(Sw2);
%     S2rw = conj(Sw2r);
     S3h8 = conj(S83h);
    S3hw = conj(Sw3h);
   
    aCw8_cor_temp = atanh(Cw8) - 1/(dof-2);
    aCw8_cor(d,:) = interp1(f,aCw8_cor_temp,f_i);
    Pw8(d,:) = interp1(f,Phiw8,f_i);
    
    aCw3_cor_temp = atanh(Cw3) - 1/(dof-2);
    aCw3_cor(d,:) = interp1(f,aCw3_cor_temp,f_i);
    Pw3(d,:) = interp1(f,Phiw3,f_i);
 
    aCw3h_cor_temp = atanh(Cw3h) - 1/(dof-2);
    aCw3h_cor(d,:) = interp1(f,aCw3h_cor_temp,f_i);
    Pw3h(d,:) = interp1(f,Phiw3h,f_i);
    
%     aCw2_cor_temp = atanh(Cw2) - 1/(dof-2);
%     aCw2_cor(d,:) = interp1(f,aCw2_cor_temp,f_i);
%     Pw2(d,:) = interp1(f,Phiw2,f_i);
    
%     aCw2r_cor_temp = atanh(Cw2r) - 1/(dof-2);
%     aCw2r_cor(d,:) = interp1(f,aCw2r_cor_temp,f_i);
%     Pw2r(d,:) = interp1(f,Phiw2r,f_i);
    
    aC83_cor_temp = atanh(C83) - 1/(dof-2);
    aC83_cor(d,:) = interp1(f,aC83_cor_temp,f_i);
    P83(d,:) = interp1(f,Phi83,f_i);
    
    aC83h_cor_temp = atanh(C83h) - 1/(dof-2);
    aC83h_cor(d,:) = interp1(f,aC83h_cor_temp,f_i);
    P83h(d,:) = interp1(f,Phi83h,f_i);
    
%     aC82_cor_temp = atanh(C82) - 1/(dof-2);
%     aC82_cor(d,:) = interp1(f,aC82_cor_temp,f_i);
%     P82(d,:) = interp1(f,Phi82,f_i);
    
%     aC82r_cor_temp = atanh(C82r) - 1/(dof-2);
%     aC82r_cor(d,:) = interp1(f,aC82r_cor_temp,f_i);
%     P82r(d,:) = interp1(f,Phi82r,f_i);
    
%     aC52r_cor_temp = atanh(C52r) - 1/(dof-2);
%     aC52r_cor(d,:) = interp1(f,aC52r_cor_temp,f_i);
%     P52r(d,:) = interp1(f,Phi52r,f_i);
    
%     aC53_cor_temp = atanh(C53) - 1/(dof-2);
%     aC53_cor(d,:) = interp1(f,aC53_cor_temp,f_i);
%     P53(d,:) = interp1(f,Phi53,f_i);
    
    lPww(d,:) = interp1(f,log10(Pww_88(:,1)),f_i);
    lP88(d,:) = interp1(f,log10(Pww_88(:,2)),f_i);
%     lP55(d,:) = interp1(f,log10(P55_33(:,1)),f_i);
    lP33(d,:) = interp1(f,log10(Pww_33(:,2)),f_i);
    lP33h(d,:) = interp1(f,log10(Pww_33h(:,2)),f_i);
%     lP22(d,:) = interp1(f,log10(Pww_22(:,2)),f_i);
%     lP22r(d,:) = interp1(f,log10(Pww_22r(:,2)),f_i);
    
    num = abs(S3w.*Pww_88(:,2) - S38.*S8w).^2;
    denom = (Pww_33(:,2).*Pww_88(:,2) - abs(S38).^2).*(Pww_33(:,1).*Pww_88(:,2) - abs(Sw8).^2);
    partial_w3(d,:) = interp1(f,num./denom,f_i);
    
    num = abs(S3hw.*Pww_88(:,2) - S3h8.*S8w).^2;
    denom = (Pww_33h(:,2).*Pww_88(:,2) - abs(S3h8).^2).*(Pww_33h(:,1).*Pww_88(:,2) - abs(Sw8).^2);
    partial_w3h(d,:) = interp1(f,num./denom,f_i);

    %     num = abs(S2w.*Pww_88(:,2) - S28.*S8w).^2;
%     denom = (Pww_22(:,2).*Pww_88(:,2) - abs(S28).^2).*(Pww_22(:,1).*Pww_88(:,2) - abs(Sw8).^2);
%     partial_w2(d,:) = interp1(f,num./denom,f_i);
    
%     num = abs(S2rw.*Pww_88(:,2) - S2r8.*S8w).^2;
%     denom = (Pww_22r(:,2).*Pww_88(:,2) - abs(S2r8).^2).*(Pww_22r(:,1).*Pww_88(:,2) - abs(Sw8).^2);
%     partial_w2r(d,:) = interp1(f,num./denom,f_i);
    
    numer = abs(S8w.*Pww_33(:,2) - S83.*S3w).^2;
    denomer = (Pww_88(:,2).*Pww_33(:,2) - abs(S83).^2).*(Pww_88(:,1).*Pww_33(:,2) - abs(Sw3).^2);
    partial_w8(d,:) = interp1(f,numer./denomer,f_i);
    
    cd G:\WC_Germany\overall_EC\
    save overall_allcells_coherence_data_5_19_2011 f_i aC* P* lP* partial*
end

%%
cd G:\WC_Germany\overall_EC\
save overall_allcells_coherence_data_5_19_2011 f_i aC* P* lP* partial*

%%
used_freqs = find(f_i < 6);

frange = find(f_i > 0.25 & f_i < 0.6);
for ii = 1:size(P83,1)
    avg_P83(ii) = circ_mean(P83(ii,frange));
%     avg_P82r(ii) = circ_mean(P82r(ii,frange));
%     avg_P82(ii) = circ_mean(P82(ii,frange));
    avg_P53(ii) = circ_mean(P53(ii,frange));
%     avg_P52r(ii) = circ_mean(P52r(ii,frange));
end
correct_lf3phase = find(avg_P83 > -0.4);
max_C83 = max(aC83_cor,[],2);
% max_C82r = max(aC82r_cor,[],2);
% max_C82 = max(aC82_cor,[],2);
max_C53 = max(aC53_cor,[],2);
% max_C52r = max(aC52r_cor,[],2);

% max_P22r = nan(size(lP22r,1),1);
% max_P22r_loc = nan(size(lP22r,1),1);
% max_P22 = nan(size(lP22,1),1);
% max_P22_loc = nan(size(lP22,1),1);
max_P33 = nan(size(lP88,1),1);
max_P33_loc = nan(size(lP88,1),1);
max_P55 = nan(size(lP88,1),1);
max_P55_loc = nan(size(lP88,1),1);
max_P88 = nan(size(lP88,1),1);
max_P88_loc = nan(size(lP88,1),1);
for ii = 1:size(lP88,1)
    %     [max_P88(ii),max_P88_loc(ii)] = nanmax(lP88(ii,peakrange));
    %     [max_P33(ii),max_P33_loc(ii)] = nanmax(lP33(ii,peakrange));
    %     [max_P22r(ii),max_P22r_loc(ii)] = nanmax(lP22r(ii,peakrange));
%     if any(~isnan(lP22r(ii,:)))
%         smooth_spectrum = jmm_smooth_1d_cor(real(lP22r(ii,used_freqs)),3);
%         [temp,temploc] = findpeaks(smooth_spectrum,'npeaks',1);
%         if ~isempty(temp)
%             max_P22r(ii) = temp; max_P22r_loc(ii) = temploc;
%         end
%     end
%     if any(~isnan(lP22(ii,:)))
%         smooth_spectrum = jmm_smooth_1d_cor(real(lP22(ii,used_freqs)),3);
%         [temp,temploc] = findpeaks(smooth_spectrum,'npeaks',1);
%         if ~isempty(temp)
%             max_P22(ii) = temp; max_P22_loc(ii) = temploc;
%         end
%     end
    if any(~isnan(lP33(ii,:)))
        smooth_spectrum = jmm_smooth_1d_cor(real(lP33(ii,used_freqs)),3);
        [temp,temploc] = findpeaks(smooth_spectrum,'npeaks',1);
        if ~isempty(temp)
            max_P33(ii) = temp; max_P33_loc(ii) = temploc;
        end
    end
%     if any(~isnan(lP55(ii,:)))
%         smooth_spectrum = jmm_smooth_1d_cor(real(lP55(ii,used_freqs)),3);
%         [temp,temploc] = findpeaks(smooth_spectrum,'npeaks',1);
%         if ~isempty(temp)
%             max_P55(ii) = temp; max_P55_loc(ii) = temploc;
%         end
%     end
    if any(~isnan(lP88(ii,:)))
        smooth_spectrum = jmm_smooth_1d_cor(real(lP88(ii,used_freqs)),3);
        [temp,temploc] = findpeaks(smooth_spectrum,'npeaks',1);
        if ~isempty(temp)
            max_P88(ii)= temp; max_P88_loc(ii) = temploc;
        end
    end
end
good_P88 = find(~isnan(max_P88));
for i = 1:length(good_P88)
%     peak_P82r(good_P88(i)) = P82r(i,max_P88_loc(good_P88(i)));
    peak_P83(good_P88(i)) = P83(i,max_P88_loc(good_P88(i)));
end

good_P55 = find(~isnan(max_P88));
for i = 1:length(good_P55)
%     peak_P52r(good_P55(i)) = P52r(i,max_P88_loc(good_P88(i)));
    peak_P53(good_P55(i)) = P53(i,max_P88_loc(good_P88(i)));
end

good_P33 = find(~isnan(max_P33));
max_P33_loc(good_P33) = f_i(used_freqs(max_P33_loc(good_P33)));
max_P88_loc(good_P88) = f_i(used_freqs(max_P88_loc(good_P88)));
% max_P55_loc(good_P55) = f_i(used_freqs(max_P55_loc(good_P55)));

% good_P22r = find(~isnan(max_P22r));
% max_P22r_loc(good_P22r) = f_i(used_freqs(max_P22r_loc(good_P22r)));
% 
% good_P22 = find(~isnan(max_P22));
% max_P22_loc(good_P22) = f_i(used_freqs(max_P22_loc(good_P22)));

max_P33 = max_P33';
max_P88 = max_P88';
% max_P55 = max_P55';
% max_P22 = max_P22';
% max_P22r = max_P22r';

%% use agglomerative clustering with cylindrical distance metric to isolate
%% useable hippocampal LFPs

% good_lf2_cells = find(~isnan(peak_P82r) & ~isnan(max_P22r) & max_P22 > -5); %get rid of bad data and an outlier
good_lf3_cells = find(~isnan(peak_P83)); %get rid of bad data and an outlier

% Y = pdist(peak_P82r(good_lf2_cells)',@(Xi,Xj) abs(circ_dist2(Xi,Xj)));
% Y = pdist(peak_P82r(good_lf2_cells)',@(Xi,Xj) abs(circ_dist2(Xi,Xj)));
% Y = pdist(peak_P53(good_lf3_cells)',@(Xi,Xj) abs(circ_dist2(Xi,Xj)));
Y = pdist(peak_P83(good_lf3_cells)',@(Xi,Xj) abs(circ_dist2(Xi,Xj)));

% Y2 = pdist(max_P22r(good_lf2_cells)','minkowski',1);
% Y2 = pdist(max_C52r(good_cells),'minkowski',1);
% Y3 = Y/circ_std(peak_P82r(good_lf2_cells)) + Y2/std(max_P22r(good_lf2_cells)); %cylindrical '1-norm'
% Y3 = Y/circ_std(peak_P82r(good_cells)) + Y2/std(max_C52r(good_cells)); %cylindrical '1-norm'

% Z = linkage(Y,'average');
Z = linkage(Y,'complete');
% Z = linkage(Y3,'average');
% dendrogram(Z,'colorthreshold','default')
% title('LF2r dendrogram')

figure(1)
n_clusters = 2;
T = cluster(Z,n_clusters);
cmap = colormap(jet(n_clusters));
for i = 1:n_clusters
%     avg_cluster_pow(i) = mean(max_P22r(good_lf2_cells(T==i)));
    avg_cluster_pow(i) = mean(max_P33(good_lf3_cells(T==i)));
%     plot(max_P22(good_lf2_cells(T==i)),peak_P82r(good_lf2_cells(T==i)),'.','color',cmap(i,:)), hold on
    plot(max_P33(good_lf3_cells(T==i)),peak_P83(good_lf3_cells(T==i)),'.','color',cmap(i,:)), hold on
end
xlabel('Maximum LF2r-LF8 Coherence','fontsize',14)
ylabel('Relative LF2r - LF8 UDS phase','fontsize',14)

[~,good_clust] = max(avg_cluster_pow);
good_hipp_lf2r = good_lf2_cells(T==good_clust);

figure
plot(max_P22r,max_P33,'.')
line([-6 -0.5],[-6 -0.5],'color','k')
xlabel('LF2r Peak UDS power','fontsize',14)
ylabel('LF3 Peak UDS Power','fontsize',14)

figure
plot(max_P22r,max_P88,'.')
line([-6 -0.5],[-6 -0.5],'color','k')
xlabel('LF2r Peak UDS power','fontsize',14)
ylabel('LF8 Peak UDS Power','fontsize',14)

%%
good_lf3_cells = find(~isnan(peak_P83) & ~isnan(max_P33)); %get rid of bad data and an outlier

Y = pdist(peak_P83(good_lf3_cells)',@(Xi,Xj) abs(circ_dist2(Xi,Xj)));
Y2 = pdist(max_P33(good_lf3_cells)','minkowski',1);
Y3 = Y/circ_std(peak_P83(good_lf3_cells)) + Y2/std(max_P33(good_lf3_cells)); %cylindrical '1-norm'
% Y3 = Y + Y2; %cylindrical '1-norm'

% Z = linkage(Y,'complete');
Z = linkage(Y,'average');
% Z = linkage(Y3,'average');
% dendrogram(Z,0,'colorthreshold','default')

figure(2)
n_clusters = 2;
T = cluster(Z,n_clusters);
cmap = colormap(jet(n_clusters));
for i = 1:n_clusters
    avg_cluster_pow(i) = mean(max_P33(good_lf3_cells(T==i)));
    plot(max_P33(good_lf3_cells(T==i)),peak_P83(good_lf3_cells(T==i)),'.','color',cmap(i,:)), hold on
end
xlabel('Maximum LF3 UDS power','fontsize',14)
ylabel('Relative LF3 - LF8 UDS phase','fontsize',14)

[~,good_clust] = max(avg_cluster_pow);
good_hipp_lf3 = good_lf3_cells(T==good_clust);

%%
both_good = intersect(good_hipp_lf2r,good_hipp_lf3);
only_lf2r = setdiff(good_hipp_lf2r,good_hipp_lf3);
only_lf3 = setdiff(good_hipp_lf3,good_hipp_lf2r);
lf2r_better = both_good(max_P22r(both_good) > max_P33(both_good));
lf3_better = both_good(max_P33(both_good) > max_P22r(both_good));
use_lf2r = union(only_lf2r,lf2r_better);
use_lf3 = union(only_lf3,lf3_better);
used_hipp_lfp = union(use_lf2r,use_lf3);
best_hipp_lfp = zeros(size(peak_P83));
best_hipp_lfp(use_lf2r) = 2;
best_hipp_lfp(use_lf3) = 3;

figure(1)
hold on
plot(max_P22r(use_lf2r),peak_P82r(use_lf2r),'o')
figure(2)
hold on
plot(max_P33(use_lf3),peak_P83(use_lf3),'o')

%%
% save overall_allcells_coherence_clustered max_* avg_* good_hipp* best_hipp* use_* used_*
