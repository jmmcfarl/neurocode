clear all
load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update

dsf = 8;
params.Fs = 2016/dsf;
% params.err = 0;
params.fpass = [0 params.Fs/2];
Ktapers = 4;
NW = (Ktapers+1)/2;
params.tapers = [NW Ktapers];
parmams.pad = 4;

win = [50 10];
movingwin = [10 1];
tau = 10;
% f0 = [50 100 150 200 250 300 350 400];
% f0 = 49.9;
for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike


    down_w = downsample(wcv_minus_spike,dsf);
    down_8 = downsample(lf8,dsf);

% [S8,f]=mtspectrumsegc(down_8,win,params);
    
[datac,datafit,Amps,freqs]=rmlinesmovingwinc(down_8,movingwin,tau,params,0.1,'y');
%     data=rmlinesc(down_8,params,.05,'y',[]);
    
end
