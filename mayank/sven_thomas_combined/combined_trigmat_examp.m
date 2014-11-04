clear all
close all

cur_dir_letter = 'C';
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\new_mec\')
cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir.mat
%%
uset = sort([l3mec l3lec]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);
lcf2 = 1/niqf;
hcf2 = 10/niqf;
[b2,a2] = butter(2,[lcf2 hcf2]);
lcf3 = 15/niqf;
hcf3 = 80/niqf;
[b3,a3] = butter(2,[lcf3 hcf3]);
dsf = 8;
Fsd = Fs/dsf;
pow_smooth = round(Fsd*0.05);

backlag = 4*Fsd;
forwardlag = 10*Fsd;
lags = (-backlag:forwardlag)/Fsd;

%%
d = 14;

cur_dir = combined_dir{d};
cd(cur_dir)
pwd

load ./used_data lf7 wcv_minus_spike
lf8 = lf7;
wcv_f = filtfilt(b,a,wcv_minus_spike);
lf8_f = filtfilt(b,a,lf8);
wcv_f = downsample(wcv_f,dsf)
lf8_f = downsample(lf8_f,dsf);

wcv_f = zscore(wcv_f);
lf8_f = zscore(lf8_f);

t_axis = (1:length(lf8_f))/Fsd;

%% extract up and down transition times for MP and LF8
load ./pa_hsmm_state_seq_combined_fin.mat
load ./pa_hsmm_state_seq7_combined_fin.mat
hsmm_bbstate_seq8 = hsmm_bbstate_seq7;

mp_state_seq_c =  hsmm_bbstate_seq;

[new_mp_seg_inds] = round(resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,hsmm_bbstate_seq));

mp_state_seq = nan(size(t_axis));

mp_utrans = [];
mp_dtrans = [];
for n = 1:hmm.Nsegs
    mp_state_seq(new_mp_seg_inds(n,1):new_mp_seg_inds(n,2)) = mp_state_seq_c{n};
    cur_mp_utrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 1 & mp_state_seq_c{n}(2:end) == 2);
    cur_mp_dtrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 2 & mp_state_seq_c{n}(2:end) == 1);
    cur_mp_dtrans(cur_mp_dtrans < cur_mp_utrans(1)) = [];
    cur_mp_utrans(cur_mp_utrans > cur_mp_dtrans(end)) = [];
    mp_utrans = [mp_utrans; cur_mp_utrans];
    mp_dtrans = [mp_dtrans; cur_mp_dtrans];
end
n_mp_ups = length(mp_utrans);

%% initialize
mp_utrig_mp_mat = nan(n_mp_ups,length(lags));
mp_utrig_lf8_mat = nan(n_mp_ups,length(lags));
mp_dtrig_mp_mat = nan(n_mp_ups,length(lags));
mp_dtrig_lf8_mat = nan(n_mp_ups,length(lags));

%     calculate mp utrigs
for i = 1:n_mp_ups
    if mp_utrans(i) > backlag && length(wcv_f) - mp_utrans(i) > forwardlag
        mp_utrig_mp_mat(i,:) = wcv_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
        mp_utrig_lf8_mat(i,:) = lf8_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
    end
end

%     calculate mp dtrigs
for i = 1:n_mp_ups
    if mp_dtrans(i) > backlag && length(wcv_f) - mp_dtrans(i) > forwardlag
        mp_dtrig_mp_mat(i,:) = wcv_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
        mp_dtrig_lf8_mat(i,:) = lf8_f(mp_dtrans(i)-backlag:mp_dtrans(i)+forwardlag);
    end
end
mp_updur = (mp_dtrans-mp_utrans)/Fsd;
mp_downdur = (mp_utrans(2:end)-mp_dtrans(1:end-1))/Fsd;

%%
cur_mp_mat = mp_utrig_mp_mat;
cur_lfp_mat = mp_utrig_lf8_mat;

%get rid of the very first row of the matrix 
bad_rows = 1;
bad_rows = [1 size(cur_mp_mat,1) size(cur_mp_mat,1)-1];
cur_mp_mat(bad_rows,:) = [];
cur_lfp_mat(bad_rows,:) = [];
mp_updur(bad_rows) = [];

[dummy,up_order] = sort(mp_updur);
%%
Fig = figure(1);
clf
set(gca,'fontname','arial','fontsize',14)
set(Fig,'PaperUnits','centimeters');
set(Fig, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
imagesc(lags,(1:length(mp_updur)),(cur_mp_mat(up_order,:)));shading flat;
hold on
plot(mp_updur(up_order),(1:length(mp_updur)),'r','linewidth',2)
line([0 0],[0 1],'Color','k')
caxis([-2.5 2.5]);colorbar
xlim([-2 10])
xlabel('Time (s)','FontSize',16,'fontname','arial')
ylabel('MP Up Duration Percentile','FontSize',16,'fontname','arial')

Fig2 = figure(2);
clf
set(gca,'fontname','arial','fontsize',14)
set(Fig2,'PaperUnits','centimeters');
set(Fig2, 'PaperSize', [30 20]);% paper size is in [width height] format
set(Fig2,'PaperPosition',[0,0,(get(Fig2,'PaperSize'))])
imagesc(lags,(1:length(mp_updur)),(cur_lfp_mat(up_order,:)));shading flat;
hold on
plot(mp_updur(up_order),(1:length(mp_updur)),'r','linewidth',2)
line([0 0],[0 1],'Color','k')
caxis([-1.5 2.5]);colorbar
xlim([-2 10])

trans_num_1 = 26;
begt_1 = 100.8;
endt_1 = 103.5;
trans_num_2 = 125;
begt_2 = 495.1;
endt_2 = 499.3;
trans_num_3 = 32;
begt_3 = 125.3;
endt_3 = 131.1;
trans_num_4 = 38;
begt_4 = 145.2;
endt_4 = 153.2;
