clear all; dt   = 19.9920031987/1000; sdim = 64; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STIMULI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimfnames    = {'af3-64p-50h-2m-mn125-ct015.s0', 'wn3-64p-50h-2m-mn125-ct015.s1', 'pt3-64p-50h-2m-mn125-ct015.s2', 'ps3-64p-50h-2m-mn125-ct015.s3', 'ns3-64p-50h-2m-mn125-ct015.s4', 'pn3-64p-50h-2m-mn125-ct015.s5', 'af3-64p-50h-2m-mn125-ct025.s6', 'wn3-64p-50h-2m-mn125-ct025.s7', 'ns3-64p-50h-2m-mn125-ct025.s8', 'pn3-64p-50h-2m-mn125-ct025.s9', 'ps3-64p-50h-2m-mn125-ct025.s10', 'pt3-64p-50h-2m-mn125-ct025.s11', 'pn3-64p-50h-2m-mn125-ct035.s12', 'ps3-64p-50h-2m-mn125-ct035.s13', 'ns3-64p-50h-2m-mn125-ct035.s14', 'pt3-64p-50h-2m-mn125-ct035.s15', 'af3-64p-50h-2m-mn125-ct035.s16', 'wn3-64p-50h-2m-mn125-ct035.s17', 'pt3-64p-50h-2m-mn125-ct045.s18', 'pn3-64p-50h-2m-mn125-ct045.s19', 'ns3-64p-50h-2m-mn125-ct045.s20', 'af3-64p-50h-2m-mn125-ct045.s21', 'ps3-64p-50h-2m-mn125-ct045.s22', 'wn3-64p-50h-2m-mn125-ct045.s23'}; 
stimfnames30  = cellfun(@(x)x(1:26),stimfnames,'UniformOutput',0); 
nstims        = length(stimfnames30)
cellfun(@disp,stimfnames30)


datdir = '~/Data/blanche/stimuli'
%stimfiles = dir([datdir,'/stim_data/*.stim']); 
save stimfiles30.mat stimfnames30;  

rstimparts= cell(nstims,1); 
for istim =1:nstims; 
	tstimfile = [datdir,'/',stimfnames30{istim}]
	fid    = fopen(tstimfile,'r'); rawdat = fread(fid,inf,'uchar'); fclose(fid);	

	[min(rawdat),max(rawdat),mean(rawdat),var(rawdat)]
	rawdat            = 2*((rawdat/255)-0.5); 
    temp(istim) = length(rawdat);
	rstimparts{istim} = rawdat; 
end; 

%% reshape, normalize & downsample the stimulus [-1,1]
dsr=2; 
stimps30  = cellfun(@(x)reshape(x,sdim*sdim,6000)',rstimparts,'UniformOutput',0); 
dstimps30 = cellfun(@(x)dsamplebstim(x,dsr),stimps30,'UniformOutput',0); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TIMING of frames and SEquence START/ENDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datdir = '~/Data/blanche/rec_30'
fid  = fopen([datdir,'/stim_data/30_-_track_6_rf-mapping3.din'],'r');
frt  = fread(fid,inf,'uint64'); fclose(fid); ntot=length(frt); 
rfrt = reshape(frt,2,ntot/2)'; rfrt(1:10,:)
ts0  = frt(1,1);

tims      = (rfrt(:,1)-ts0)/1000; 
fonsets   = reshape(tims,4,size(tims,1)/4)'; 
seqstarts = [0; fonsets( (1:23)*6000+1 , 1)]/1000



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPIKES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the spikes
spikefiles = dir([datdir,'/spk_data/*.spk']); 
ncells     = length(spikefiles); 
rspks30    = cell(ncells,1); 

% Read spike times
for ispfile = 1:ncells
  fid = fopen([datdir,'/spk_data/',spikefiles(ispfile).name],'r');
  rawts = fread(fid,inf,'uint64'); fclose(fid);
  ts = (rawts - ts0)/1000;
  rspks30{ispfile} = ts(ts > 0)/1000;
end

%% cut them into the parts corresponding to stimulus presentation
seqlen = dt*6000;  
spksegs30 = cell(ncells,24);
for icell=1:ncells;
	tspks = rspks30{icell}; 
	for iep =1:24; 
		tstart = seqstarts(iep); tend=tstart+seqlen; 
		spksegs30{icell,iep} = tspks( tspks>tstart & tspks < tend) - tstart; 
	end; 
end; 

%% concatenate them as if they were presented without any gaps
spks30 = cell(ncells,1);
for icell=1:ncells; 
	spks =[]; 
	for iseq=1:24; spks=[spks;spksegs30{icell,iseq}+(iseq-1)*seqlen]; end; 
    spks30{icell}=spks; 
end; 


%% check alignment of firing episodes with stimulus segments
figure; 
ecdf(spks30{1}); hold on; for icell =2:ncells; ecdf(spks30{icell}); end; hold off; 
vline((1:24)*6000*dt)

figure; 
ecdf(rspks30{1}); hold on; for icell =2:ncells; ecdf(rspks30{icell}); end; hold off; 
vline(seqstarts); axis tight

plot(rspks30{1},1:length(rspks30{1})); hold on; 
   for icell =2:ncells; plot(rspks30{icell},1:length(rspks30{icell})); end; hold off;

   
save spksRec30.mat spikefiles rspks30 dt
save dstimpsRec30.mat dstimps30
save spksegsRec30.mat spksegs30; 

%%
cnames   = [stimfnames30]; cellfun(@disp,cnames)
nconds   = length(cnames);
keys = {'af','wn','pn','pt','ps','ns'}
rids = cellfun(@(tkey)...
    find(strcmp(cellfun(@(x)x(1:2),cnames,'UniformOutput',0),tkey)),...
    keys,'UniformOutput',0);

taf        = rids{1}; %cellfun(@disp,cnames(tpn))
twn        = rids{2};
tpn        = rids{3}; %cellfun(@disp,cnames(tpn))
tpt        = rids{4};
tps        = rids{5}; %cellfun(@disp,cnames(tps))
tns        = rids{6}; %cellfun(@disp,cnames(tns))
usets = [rids{3}(3:4) rids{4}(3:4) rids{5}(3:4) rids{6}(3:4)]; ctype='all'; 
ccell =2;
stimlen = 6000;
longspks = [];
longstim = [];
for cond = 1:length(usets)
   cur_cond = usets(cond);
   longstim = [longstim; dstimps30{cur_cond}];
   longspks = [longspks; (cond-1)*stimlen*dt+spksegs30{ccell,cur_cond}];
    
end

NT = size(longstim,1);
spikebin_ids = 1+floor(longspks/dt);
spikebins = hist(spikebin_ids,1:NT)';
iisp = find(spikebins > 0);
splen = length(iisp);
nsp = sum(spikebins(iisp));
SS = makeStimRows(longstim,flen,iisp);
STA = (spikebins(iisp)'*SS)'/nsp;
S0 = makeStimRows(longstim,flen,1);
ov_mean = mean(S0)';
STA = STA-ov_mean;
plotfilterbank(STA,32,1:1024)

ccell =10;
clear STA
flen = 6;
for cond = 1:24
        cur_stim = dstimps30{cond};
    cur_spike_times = spksegs30{ccell,cond};
    if length(cur_spike_times) > 10
    stimfnames30{cond}
        spikebin_ids = 1+floor(cur_spike_times/dt);
%         spikebin_ids = ceil(6000*rand(size(cur_spike_times)));
        spikebins = hist(spikebin_ids,1:6000)';
%         spikebins = hist(spikebin_ids,1:250)';
        
        iisp = find(spikebins>0);
        splen = length(iisp);
        nsp = sum(spikebins(iisp));
%         Msz = splen*size(cur_stim,2)*flen;  % size of "full" spike-triggered stimulus matrix
            SO = makeStimRows(cur_stim,flen,1);  % Convert stimulus to matrix where each row is one stim
            ov_mean = mean(SO);
        SS = makeStimRows(cur_stim, flen, iisp);
        rowlen = size(SS,2);
        STA(:,cond) = (spikebins(iisp)'*SS)'/nsp-ov_mean';
        plotfilterbank(STA(:,cond),32,1:1024)
%         figure
%         plotfilterbank(ov_mean',32,1:1024)
        pause
        close all
    end
end

% usets = [tpn,tns]; ctype='all'; 
nspks = cellfun(@(x)length(x),spksegs30(ccell,:));
usets = find(nspks > 100);
% usets(usets > size(STA,2)) = [];
avg_sta = STA(:,usets)*nspks(usets)';
figure
plotfilterbank(avg_sta,32,1:1024);




%%
% 
% %% make a little stimulus movie
% tstim = dstimps30{15}; tsdim = sqrt(size(tstim,2))
% 
% aviobj = avifile('natstimmovie.avi','fps',15);
% for tslice =1+(1:200);
% 	tslice
% 	plot2drfmat(reshape(tstim(tslice,:),tsdim,tsdim),[-3,3]);
% 	frame = getframe(gca); aviobj = addframe(aviobj,frame);
% end;
% aviobj = close(aviobj);
% 
% 
% 
% 
% ************
% 
% rfrt = reshape(frt,2,ntot/2)'; rfrt(1:10,:)
% ts0  = frt(1,1); 
% tims      = (rfrt(:,1)-ts0)/1000; 
% fonsets   = reshape(tims,4,size(tims,1)/4)'; 
% seqstarts = [0; fonsets( (1:23)*6000+1 , 1)]/1000
% %plot(fonsets(:,1))
% 
% %% read the spikes
% spikefiles = dir([datdir,'/spk_data/*.spk']); ncells = length(spikefiles)
% rspks30 = cell(ncells,1); 
% xspks = cell(ncells,1); 
% 
% 
% % Read spike times
% for ispfile = 1:ncells
%   fid = fopen([datdir,'/spk_data/',spikefiles(ispfile).name],'r');
%   rawts = fread(fid,inf,'uint64'); fclose(fid);
%   ts = (rawts - ts0)/1000;
%   x = ts(ts > 0)/1000;
%   rspks30{ispfile} = x(x>13.9)-13.9;   
% end
% 
% 
% %% total protocol len is nstim*stimlen + (n-1)*breaklen (5sec)
% seqlen = dt*6000;  
% spksegs30 = cell(ncells,nstims);
% for icell=1:ncells;
% 	tspks = rspks30{icell}; 
% 	for iep =1:nstims; 
% 		tstart = seqstarts(iep); tend=tstart+seqlen; 
% 		spksegs30{icell,iep} = tspks( tspks>tstart & tspks < tend) - tstart; 
% 	end; 
% end; 
% 
% 
% figure; 
% ecdf(rspks30{1}); hold on; for icell =2:ncells; ecdf(rspks30{icell}); end; hold off;  
% 
% tcell=10; ecdf(rspks30{tcell});
% axis tight; vline(seqstarts)
% nspks = cellfun(@length,rspks30)
% indlengs = cellfun(@length,spksegs30); [nspks,sum(indlengs,2)]

