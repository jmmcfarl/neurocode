clear all
close all

%%
cd C:\WC_Germany\sven_thomas_combined\
load .\distal_dir
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

raw_Fs = 2016;
params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 50;

lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;

rate_sm = round(raw_Fs*0.05);

freqs = linspace(0.05,80,5000);
%%
for d = 1:length(distal_dir)
    cd(distal_dir{d})
    pwd
%     if exist('./used_data2.mat','file')
%         load ./used_data2
%     else
        load ./used_data
%     end
    if ctx_lfp(d) == 7
        lf8 = lf7;
    elseif ctx_lfp(d) == 6
        lf8 = lf6;
    elseif ctx_lfp(d) == 5
        lf8 = lf5;
    end

%     [desynch_times,desynch_inds,P_lf8,f,t] = locate_desynch_times_individual(lf8);
    [desynch_times,desynch_inds,P_lf8,f,t] = locate_desynch_times_individual_v2(lf8,1);
    %compute markers indicating segments of data to be used
    t_axis = (1:length(lf8))/raw_Fs;
    if ~isempty(desynch_times)
        desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times(:,1)));
        desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times(:,2)));
    else
        desynch_start = [];
        desynch_stop = [];
    end
    desynch_ind = zeros(size(lf8));
    for i = 1:length(desynch_start)
        desynch_ind(desynch_start(i):desynch_stop(i)) = 1;
    end
    desynch_fract(d) = sum(desynch_ind)/length(desynch_ind);
    
    synch_starts = find(desynch_ind(1:end-1)==1 & desynch_ind(2:end)==0)+1;
    if desynch_ind(1) == 0
        synch_starts = [1; synch_starts];
    end
    synch_stops = find(desynch_ind(1:end-1)==0 & desynch_ind(2:end)==1)+1;
    if desynch_ind(end) == 0
        synch_stops = [synch_stops; length(lf8)];
    end
    sMarkers = [synch_starts(:) synch_stops(:)];
        
    if ~isnan(hpc_mua(d))
        load ./mua_data3
        load ./sync_times.mat
        synctt = synct;
        hpc_mua_times = mua_times{hpc_mua(d)};
        hpc_mua_times(hpc_mua_times < synctt(1) | hpc_mua_times > synctt(end)) = [];
        temp_mua_rate =hist(hpc_mua_times,synct)*raw_Fs;
        mua_rate = jmm_smooth_1d_cor(temp_mua_rate,rate_sm);
        mua_rate = zscore(mua_rate);
        if length(mua_rate) > length(t_axis)
            mua_rate = mua_rate(1:length(t_axis));
        end
    end
    
    %zscore signals
%     wcv_minus_spike = zscore(wcv_minus_spike);
%     lf8 = zscore(lf8);
%     lf7 = zscore(lf7);
%     lf6 = zscore(lf6);
%     lf5 = zscore(lf5);
%     lf4 = zscore(lf4);
%     lf3 = zscore(lf3);
%     lf2 = zscore(lf2);
    
    [S(d,1,:),f]=mtspectrumc_unequal_length_trials(wcv_minus_spike,[win win],params,sMarkers);
    [S(d,8,:),f]=mtspectrumc_unequal_length_trials(lf8,[win win],params,sMarkers);
    [S(d,7,:),f]=mtspectrumc_unequal_length_trials(lf7,[win win],params,sMarkers);
    [S(d,6,:),f]=mtspectrumc_unequal_length_trials(lf6,[win win],params,sMarkers);
    [S(d,5,:),f]=mtspectrumc_unequal_length_trials(lf5,[win win],params,sMarkers);
    [S(d,4,:),f]=mtspectrumc_unequal_length_trials(lf4,[win win],params,sMarkers);
    [S(d,3,:),f]=mtspectrumc_unequal_length_trials(lf3,[win win],params,sMarkers);
    [S(d,2,:),f]=mtspectrumc_unequal_length_trials(lf2,[win win],params,sMarkers);
    if ~isnan(hpc_mua(d))
    [Smua(d,:),f]=mtspectrumc_unequal_length_trials(mua_rate(:),[win win],params,sMarkers);
    else
        Smua(d,:) = nan(1,length(f));
    end
       
end

%%
n_recs = size(S,1);

uds_freqs = find(f > 0.2 & f < 1);
lS = real(10*log10(S));
[peak_S,peakf] = max(lS(:,:,uds_freqs),[],3);
peak_uds_freq = f(uds_freqs(peakf));

high_freqs = find(f >= 15);
hg_pow = trapz(lS(:,:,high_freqs),3);

for i = 1:n_recs
    if ~isnan(hpc_lfp(i))
        hpc_lS(i,:) = squeeze(lS(i,hpc_lfp(i),:));
    end
end

%%
cd C:\WC_Germany\sven_thomas_combined\
save distal_spectra_fin_nd S* f

%%
mp_S = squeeze(S(:,1,:));
for i = 1:size(S,1)
    ctx_S(i,:) = squeeze(S(i,ctx_lfp(i),:));
end
l3mec_m = find(~isnan(hpc_mua));
lmp_S = log10(mp_S);
lctx_S = log10(ctx_S);
lSmua = log10(Smua);
figure; set(gca,'fontname','arial','fontsize',14)
hold on
shadedErrorBar(f,mean(ctx_S(l3mec_m,:)),std(ctx_S(l3mec_m,:))/sqrt(length(l3mec_m)),{'color',[0.2 0.2 0.2]});
shadedErrorBar(f,mean(mp_S(l3mec_m,:)),std(mp_S(l3mec_m,:))/sqrt(length(l3mec_m)),{'r'});
% shadedErrorBar(f,mean(Smua(l3mec_m,:)),std(Smua(l3mec_m,:))/sqrt(length(l3mec_m)),{'color',[0.2 0.8 0.2]});
xlim([0 1.])
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Relative Power','fontsize',16)

%%
uds_freqs = find(f > 0.1 & f < 0.8);
[mp_p,mp_f] = max(mp_S(:,uds_freqs),[],2);
[lfp_p,lfp_f] = max(ctx_S(:,uds_freqs),[],2);
[mua_p,mua_f] = max(Smua(:,uds_freqs),[],2);
mp_m = mean(mp_S(:,uds_freqs),2);
lfp_m = mean(ctx_S(:,uds_freqs),2);
mua_m = mean(Smua(:,uds_freqs),2);
ctx_uds_freq = f(uds_freqs(lfp_f));
mp_uds_freq = f(uds_freqs(mp_f));
mua_uds_freq = f(uds_freqs(mua_f));

noise_freqs = find(f < 0.01);
ctx_offset = sum(ctx_S(:,noise_freqs),2);


%%
l3mec_m = find(~isnan(hpc_mua));
df = f(2)-f(1);
norm_mp_S = bsxfun(@rdivide,mp_S,trapz(f,mp_S,2));
norm_ctx_S = bsxfun(@rdivide,ctx_S,trapz(f,ctx_S,2));
norm_hpc_mua = bsxfun(@rdivide,Smua,trapz(f,Smua,2));
figure; set(gca,'fontname','arial','fontsize',14)
hold on
shadedErrorBar(f,mean(norm_ctx_S(l3mec_m,:)),std(norm_ctx_S(l3mec_m,:))/sqrt(length(l3mec_m)),{'color',[0.2 0.2 0.2]});
shadedErrorBar(f,mean(norm_mp_S(l3mec_m,:)),std(norm_mp_S(l3mec_m,:))/sqrt(length(l3mec_m)),{'r'});
shadedErrorBar(f,mean(norm_hpc_mua(l3mec_m,:)),std(norm_hpc_mua(l3mec_m,:))/sqrt(length(l3mec_m)),{'g'});
% shadedErrorBar(f,mean(norm_hpc_lfp(l3mec_mlfp,:)),std(norm_hpc_lfp(l3mec_mlfp,:))/sqrt(length(l3mec_mlfp)),{'c'});
xlim([0 1.])
xlabel('Frequency (Hz)','fontsize',16)
ylabel('Relative Power','fontsize',16)


