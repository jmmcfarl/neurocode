clear all; close all;
%%%%%%% INITIAL PROCESSING OF Tim Data
dt   = 19.9920031987/1000;
% dt   = 20/1000;
sdim = 64;

rstimfiles = {'ara-64p-50h-2m-mn125-ct015.s0', 'wra-64p-50h-2m-mn125-ct015.s1', 'pra-64p-50h-2m-mn125-ct015.s2', 'nra-64p-50h-2m-mn125-ct015.s3', ...
    'nra-64p-50h-2m-mn125-ct025.s4', 'ara-64p-50h-2m-mn125-ct025.s5', 'wra-64p-50h-2m-mn125-ct025.s6', 'pra-64p-50h-2m-mn125-ct025.s7',...
    'pra-64p-50h-2m-mn125-ct035.s8', 'ara-64p-50h-2m-mn125-ct035.s9', 'nra-64p-50h-2m-mn125-ct035.s10', 'wra-64p-50h-2m-mn125-ct035.s11',...
    'nra-64p-50h-2m-mn125-ct045.s12', 'ara-64p-50h-2m-mn125-ct045.s13', 'pra-64p-50h-2m-mn125-ct045.s14', 'wra-64p-50h-2m-mn125-ct045.s15'};

stimfiles  = cellfun(@(x)x(1:26),rstimfiles,'UniformOutput',0);
nstims = length(stimfiles);

datdir = '~/Data/blanche/rec_74';
%stimfiles = dir([datdir,'/stim_data/*.stim']);
cd(datdir)
cd matlabdata/
save stimfiles74.mat stimfiles

nfiles=length(stimfiles); allstims =[]; rstimparts= cell(nfiles,1);
for istim =1:nfiles;
    % istim = 3;
    tstimfile = stimfiles{istim}; %current dataset
    %read in current raw data
    [fid,message]    = fopen([datdir,'/stim_data/',tstimfile],'r');
    rawdat = fread(fid,inf,'uchar');
    fclose(fid);
    temp(istim) = length(rawdat);
    
    rawdat(1:11) = []; %first 11 samples are no good
    
    NT = length(rawdat)/sdim/sdim;
    used_el = sdim^2*floor(NT);
    %     datmat = reshape(rawdat(1:used_el),sdim^2,floor(NT));
    
    rawdat = 2*((rawdat(1:used_el)/255)-0.5); %not sure why this transformation
    rstimparts{istim} = rawdat;
    allstims=[allstims;rawdat];
    
end;

% for i = 1:250
% imagesc(reshape(datmat(:,i),sdim,sdim));%drawnow();
% pause(.1)
% end
%% reshape, normalize & downsample the stimulus [-1,1]
dsr=2; %spatial downsampling factor
rep_length = 250; %number of samples per repeat
stimps74  = cellfun(@(x)reshape(x,sdim*sdim,rep_length)',rstimparts,'UniformOutput',0); %reshape the two minutes of data
dstimps74 = cellfun(@(x)dsamplebstim(x,dsr),stimps74,'UniformOutput',0);


%% Read timing data
fid  = fopen([datdir,'/stim_data/74_-_track_7c_reliability-a.din'],'r');
frt  = fread(fid,inf,'uint64'); fclose(fid);
ntot=length(frt);
rfrt = reshape(frt,2,ntot/2)';
rfrt(1:10,:)
ts0  = frt(1,1); %initial time stamp

tims      = (rfrt(:,1)-ts0)/1000; %zero initial time stamp
fonsets   = reshape(tims,4,size(tims,1)/4)';
% seqstarts = [0; fonsets( (1:23)*6000+1 , 1)]/1000 %start times of each 2min sequence (24 of them)

reps_per_seq = 24;
sequence_inds = ceil((1:size(fonsets,1))*nfiles/size(fonsets,1));
rep_inds = ceil((1:size(fonsets,1))*(nfiles*reps_per_seq)/size(fonsets,1));
% rep_inds = mod(rep_inds,reps_per_seq)+1;

seq_start_inds = find(diff(sequence_inds) == 1)+1;
seq_stop_inds = seq_start_inds-1;
seq_bound_inds = [([1 seq_start_inds])' ([seq_stop_inds length(rep_inds)])'];

rep_start_inds = find(diff(rep_inds) == 1)+1;
rep_stop_inds = rep_start_inds-1;
rep_bound_inds = [([1 rep_start_inds])' ([rep_stop_inds length(rep_inds)])'];

seq_bound_times = fonsets(seq_bound_inds)/1000;
rep_bound_times = fonsets(rep_bound_inds)/1000;


%% read the spikes
spikefiles = dir([datdir,'/spk_data/*.spk']);
ncells = length(spikefiles);
rspks74 = cell(ncells,1);
% Read spike times
for ispfile = 1:ncells
    fid = fopen([datdir,'/spk_data/',spikefiles(ispfile).name],'r');
    rawts = fread(fid,inf,'uint64'); fclose(fid);
    ts = (rawts - ts0)/1000; %subtract off initial time stamp
    rspks74{ispfile} = ts(ts > 0)/1000;
end

% bin_dt = 0.1;
% bin_ax = 0:bin_dt:2098;
% for ii = 1:ncells
%     cur_rate = hist(rspks74{ii},bin_ax)/bin_dt;
%     sm_rate(ii,:) = smooth(cur_rate,50);
% end
% z_sm_rate = zscore(sm_rate')';
% imagesc(bin_ax,1:ncells,z_sm_rate)
% xlim([0 2098])
% for i = 1:16
%     line([i*120 i*120],[0 ncells+1],'color','k','linewidth',2)
% end

%parse spike times into separate 2m recording sequences
seqlen = dt*6000; %number of time samples per sequence
spksegs74 = cell(ncells,nfiles);
spksegs74_trialrel = cell(ncells,nfiles);
for icell=1:ncells;
    tspks = rspks74{icell};
    for iep =1:nfiles;
        tstart = seq_bound_times(iep,1); tend = seq_bound_times(iep,2);
        % 		tstart = seqstarts(iep); tend=tstart+seqlen; %start and stop timestamps of current sequence
        spksegs74{icell,iep} = tspks( tspks>tstart & tspks < tend) - tstart; %spikes from current cell in current sequence
        spksegs74_trialrel{icell,iep} = spksegs74{icell,iep};
        for r = 1:reps_per_seq
            rstart = rep_bound_times(r,1); rend = rep_bound_times(r,2);
            cur_spks = find(spksegs74{icell,iep} > rstart & spksegs74{icell,iep} < rend);
            spksegs74_trialrel{icell,iep}(cur_spks) = spksegs74_trialrel{icell,iep}(cur_spks) - rstart;
        end
    end;
end;

%cell array containing all spike times (concatenated) for each cell
spks74 = cell(ncells,1);
for icell=1:ncells;
    spks =[];
    for iseq=1:nfiles; spks=[spks;spksegs74{icell,iseq}+(iseq-1)*seqlen]; end;
    spks74{icell}=spks;
end;

cell_ids = [0:6 8 10:20 22:27];

% figure;
% ecdf(spks74{1}); hold on; for icell =2:ncells; ecdf(spks74{icell}); end; hold off;
% vline((1:24)*6000*dt)
%
% ecdf(rspks74{1}); hold on; for icell =2:ncells; ecdf(rspks74{icell}); end; hold off;
% vline(seqstarts); axis tight
%
% plot(rspks74{1},1:length(rspks74{1})); hold on;
%    for icell =2:ncells; plot(rspks74{icell},1:length(rspks74{icell})); end; hold off;



% spksRec75b = [];
% for tcell = 1:ncells;
% 	tspks = [];
% 	for ifile =1:nfiles; tebo=ebounds{ifile};
% 		tspks=[tspks;spksegs{tcell,ifile}+tebo(1)];
% 	end;
% 	spksRec75b{tcell} = tspks;
% end;
%%

nspks = cellfun(@length,spks74)
cd(datdir)
cd matlabdata
save spksRec74.mat spikefiles spks74 dt stimfiles rep_bound_times seq_bound_times %save spike file names, and the concatenated arrays of spike times for each cell

save dstimpsRec74.mat dstimps74 %save cell array of spatially down-sampled stimulus sequences
save spksegsRec74.mat spksegs74; %save segmented spike times for each cell

%%
clear all
datdir = '~/Data/blanche/rec_74';
cd(datdir)
cd matlabdata
load ./spksRec74
load ./dstimpsRec74
load ./spksegsRec74
load ./stimfiles74
% dstimps74 = dstimps75;
% spksegs74 = spksegs75;

nconds   = length(dstimps74)
cnames   = [stimfiles]; cellfun(@disp,cnames)
keys = {'ar','wr','pr','nr'}
rids = cellfun(@(tkey)...
	find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
	keys,'UniformOutput',0); 
tns        = rids{4}(3:4); cellfun(@disp,cnames(tns))
tpn        = rids{3}(3:4); cellfun(@disp,cnames(tpn))
%%
usets = [1:16]; ctype='all';
clear STA
longstim = [];
for cond = 1:length(usets)
    cur_cond = usets(cond);
    longstim = [longstim; repmat(dstimps74{cur_cond},24,1)];
end
S0 = makeStimRows(longstim,flen,1);
ov_mean = mean(S0)';

for ccell =1:25;
    ccell
    flen = 6;
    stimlen = 6000;
    longspks = [];
    for cond = 1:length(usets)
        cur_cond = usets(cond);
        longspks = [longspks; (cond-1)*stimlen*dt+spksegs74{ccell,cur_cond}];
    end
    
    NT = size(longstim,1);
    spikebin_ids = 1+floor(longspks/dt);
    spikebins = hist(spikebin_ids,1:NT)';
    iisp = find(spikebins > 0);
    splen = length(iisp);
    nsp = sum(spikebins(iisp));
    SS = makeStimRows(longstim,flen,iisp);
    STA = (spikebins(iisp)'*SS)'/nsp;
    STA = STA-ov_mean;
    plotfilterbank(STA,32,1:1024)
    pause
    clf
end
%%



for cond = 1:16
%     cur_stim = repmat(dstimps74{cond},24,1);
        cur_stim = dstimps74{cond};
    cur_spike_times = spksegs74{ccell,cond};
    cur_spike_times = mod(cur_spike_times,dt*250);
    if length(cur_spike_times) > 10
        spikebin_ids = 1+floor(cur_spike_times/dt);
%         spikebin_ids = ceil(250*rand(size(cur_spike_times)));
%         spikebins = hist(spikebin_ids,1:6000)';
        spikebins = hist(spikebin_ids,1:250)';
        
        iisp = find(spikebins>0);
        splen = length(iisp);
        nsp = sum(spikebins(iisp));
        Msz = splen*size(cur_stim,2)*flen;  % size of "full" spike-triggered stimulus matrix
        %     SO = makeStimRows(cur_stim,flen,1);  % Convert stimulus to matrix where each row is one stim
        %     ov_mean = mean(SO);
        SS = makeStimRows(cur_stim, flen, iisp);
        rowlen = size(SS,2);
        STA(:,cond) = (spikebins(iisp)'*SS)'/nsp;
    end
end

usets = [tpn,tns]; ctype='all'; 
nspks = cellfun(@(x)length(x),spksegs74(ccell,:));
usets = find(nspks > 100);
% usets(usets > size(STA,2)) = [];
avg_sta = STA(:,usets)*nspks(usets)';
figure
plotfilterbank(avg_sta,32,1:1024);

%%
cd ~/Data/blanche/rec_74
stf74=load('stimfiles74.mat'); 
nconds   = length(dstimps74)
cnames   = [stf74.stimfiles]; cellfun(@disp,cnames)
keys = {'ar','wr','pr','nr'}
rids = cellfun(@(tkey)...
	find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
	keys,'UniformOutput',0); 

% tns        = rids{4}(3:4); cellfun(@disp,cnames(tns))
% tpn        = rids{3}(3:4); cellfun(@disp,cnames(tpn))
tconds = 14;
% tconds = [tpn,tns]; ctype='all'; 
% tconds = [rids{3}]; ctype='pink'; 

n_reps = 24;
selstim   = []; 
for icond=1:length(tconds); 
	tcond = tconds(icond)
    cur_stim = dstimps74{tcond};
%     cur_stim = (cur_stim - mean(cur_stim(:)))/std(cur_stim(:)); %z-score normalization
	selstim=[selstim;repmat(cur_stim,n_reps,1)];
end 
allspks  = spksegs74';  

ncells = size(allspks,2);
cellfun(@(x)size(x,1),dstimps74,'UniformOutput',0) 
stimlen = dt*6000;
nspks = cellfun(@length,allspks); 
cellfun(@disp,cnames(tconds))

aselspks = cell(ncells,1); 
for icell = 1:ncells; 
	selspks = []; 
    spk_ids = [];
	for icond=1:length(tconds); 
		tcond = tconds(icond)
		selspks=[selspks;((icond-1)*stimlen)+allspks{tcond,icell}];
        spk_ids = [spk_ids; ones(size(allspks{tcond,icell}))*icond];
	end;
	aselspks{icell} = selspks; 
end

%%
close all
ccell = 24;
for tconds = 1:16
    disp(stimfiles{tconds})

tempr = mod(allspks{tconds,ccell},250*dt);
% figure
raster(tempr)
title(stimfiles{tconds});

flen = 6;
selstim = repmat(dstimps74{tconds},n_reps,1);
NT   = size(selstim,1); SDIM  = size(selstim,2); NeK   = flen*SDIM;
X    = zeros(NT-flen+1,NeK);
for i = flen:NT; X(i-flen+1,1:NeK) = reshape(selstim((i-flen+1):i,:),1,NeK); end

tsbs = 1+floor(allspks{tconds,ccell}/dt);
spikebins = tsbs(tsbs>flen & tsbs<(size(X,1)+flen-1))-flen+1 ;

sta(tconds,:) = mean(X(spikebins,:)) - mean(X);
figure
plotfilterbank(sta(tconds,:)',32,1:1024);
pause
close all
end
nspks = cellfun(@(x)length(x),allspks(:,ccell));

tns        = rids{4}(3:4); cellfun(@disp,cnames(tns))
tpn        = rids{3}(3:4); cellfun(@disp,cnames(tpn))
usets = [tpn,tns]; ctype='all'; 
wsta = nspks'*sta;

% %% make a little stimulus movie
% tstim = dstimps74{4}; tsdim = sqrt(size(tstim,2))
% aviobj = avifile('natstimmovie2.avi','fps',15);
% for tslice =1:250;
% 	tslice
% 	plot2drfmat(reshape(tstim(tslice,:),tsdim,tsdim),[-0.75,0.75]);
% 	frame = getframe(gca); aviobj = addframe(aviobj,frame);
% end;
% aviobj = close(aviobj);
% 
% %%
% tstim = dstimps74{4}; tsdim = sqrt(size(tstim,2))
% for tslice =1:250;
% % 	tslice
% 	plot2drfmat(reshape(tstim(tslice,:),tsdim,tsdim),[-0.75,0.75]);
%     pause(.5)
% end;
% 

