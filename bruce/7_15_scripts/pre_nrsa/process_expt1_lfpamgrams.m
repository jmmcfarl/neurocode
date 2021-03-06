clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12/G034/
load ./CellList.mat
load ./G034Expts.mat

dt = 118/1e4;
fst = 1/dt;
raw_Fs = 3e4;
dsf = 30;
Fs = raw_Fs/dsf;
n_scales = 30;
scales = logspace(log10(3),log10(40),n_scales);
wfreqs = scal2frq(scales,'cmor1-1',1/Fs);

load ./all_eyedata_expt1_34

load ./Expt1_newcompiled_data_fixeddelay_d1p25_34 full_t full_expt*
%%
% Expt_nu = [1 6 16 17 20 25 28];
Expt_nu = [1 2 17 18 19 23 24]; %expt 1 34
n_allunits = 96;

full_ampgrams = [];
for ee = 1:length(Expt_nu)
    
    clear ampgram
    for n = 1:n_allunits
        fprintf('Loading LFP %d of %d\n',n,n_allunits);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),n);
        [lfp_data,lfp_timestamps,Fs] = load_lfp_data(filename,dsf);
        temp = cwt(lfp_data,scales,'cmor1-1');
        ampgram(:,:,n) = abs(cwt(lfp_data,scales,'cmor1-1'))';       
    end
    cur_set = find(full_expt_vec == Expt_nu(ee));
    cur_interp_ampgrams = interp1(lfp_timestamps,ampgram,full_t(cur_set));
    full_ampgrams = [full_ampgrams; cur_interp_ampgrams];
end

%%
save expt1_ampgram_data full_ampgrams