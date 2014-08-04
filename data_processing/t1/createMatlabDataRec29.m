clear all; 
%%%%%%% INITIAL PROCESSING OF Tim Data
dt   = 19.9920031987/1000; 
sdim = 64; 

stimfnames = {'ns1-64p-50h-2m-mn125-ctNAT.s0', 'ns1-64p-50h-2m-mnNAT-ct035.s1', 'ns1-64p-50h-2m-mnNAT-ctNAT.s2', 'ns2-64p-50h-2m-mn125-ctNAT.s3', 'ns2-64p-50h-2m-mnNAT-ct035.s4', 'ns2-64p-50h-2m-mnNAT-ctNAT.s5', 'ns3-64p-50h-2m-mnNAT-ctNAT.s6'};
stimfnames29  = cellfun(@(x)x(1:26),stimfnames,'UniformOutput',0); 
nstims = length(stimfnames29); 
cellfun(@disp,stimfnames29)

datdir = '/Users/James/Data/blanche/rec_29/'
cd(datdir)
cd matlabdata
save stimfiles29.mat stimfnames29;  
cd(datdir)
cd ~/Data/blanche/stimuli/
rstimparts= cell(nstims,1); 
for istim =1:nstims; 
	tstimfile = [stimfnames29{istim}]
	fid    = fopen(tstimfile,'r'); rawdat = fread(fid,inf,'uchar'); fclose(fid);	
	[min(rawdat),max(rawdat),mean(rawdat),var(rawdat)]
	rawdat = 2*((rawdat/255)-0.5); 
	rstimparts{istim} = rawdat; 
end; 

%% reshape, normalize & downsample the stimulus [-1,1]
dsr=2; 
stimps29  = cellfun(@(x)reshape(x,sdim*sdim,6000)',rstimparts,'UniformOutput',0); 
dstimps29 = cellfun(@(x)dsamplebstim(x,dsr),stimps29,'UniformOutput',0); 



datdir = '/Users/timm/myrepo/workspace/DataRepository/blanche/rec_29'
%% Read timing data
fid  = fopen([datdir,'/stim_data/29_-_track_6_tracking.din'],'r');
frt  = fread(fid,inf,'uint64'); fclose(fid); ntot=length(frt); 

rfrt = reshape(frt,2,ntot/2)'; rfrt(1:10,:)
ts0  = frt(1,1)
tims      = (rfrt(:,1)-ts0)/1000; 
fonsets   = reshape(tims,4,size(tims,1)/4)'; 

seqstarts = [0; fonsets( (1:(nstims-1))*6000+1 , 1)]/1000


%% read the spikes
spikefiles = dir([datdir,'/spk_data/*.spk']); ncells = length(spikefiles)
rspks29 = cell(ncells,1); 
% Read spike times
for ispfile = 1:ncells
  fid = fopen([datdir,'/spk_data/',spikefiles(ispfile).name],'r');
  rawts = fread(fid,inf,'uint64'); fclose(fid);
  ts = (rawts - ts0)/1000;
  rspks29{ispfile} = ts(ts > 0)/1000;
end


%% total protocol len is nstim*stimlen + (n-1)*breaklen (5sec)
seqlen = dt*6000;  
spksegs29 = cell(ncells,nstims);
for icell=1:ncells;
	tspks = rspks29{icell}; 
	for iep =1:nstims; 
		tstart = seqstarts(iep); tend=tstart+seqlen; 
		spksegs29{icell,iep} = tspks( tspks>tstart & tspks < tend) - tstart; 
	end; 
end; 


figure; 
ecdf(rspks29{1}); hold on; for icell =2:ncells; ecdf(rspks29{icell}); end; hold off;  
axis tight; vline(seqstarts)
nspks = cellfun(@length,rspks29)

indlengs = cellfun(@length,spksegs29); [nspks,sum(indlengs,2)]


save spksRec29.mat spikefiles rspks29 dt
save dstimpsRec29.mat dstimps29
save spksegsRec29.mat spksegs29; 



%% make a little stimulus movie
tstim = dstimps29{3}; tsdim = sqrt(size(tstim,2))
aviobj = avifile('natstimmovie.avi','fps',15);
for tslice =1+(1:200);
	tslice
	plot2drfmat(reshape(tstim(tslice,:),tsdim,tsdim),[-3,3]);
	frame = getframe(gca); aviobj = addframe(aviobj,frame);
end;
aviobj = close(aviobj);


tstim = stimps29{7}; tsdim = sqrt(size(tstim,2))
aviobj = avifile('natstimmovie.avi','fps',15);
for tslice =1+(1:200);
	tslice
	plot2drfmat(reshape(tstim(tslice,:),tsdim,tsdim),[-3,3]);
	frame = getframe(gca); aviobj = addframe(aviobj,frame);
end;
aviobj = close(aviobj);

subplot(3,1,1); plot2drfmat(reshape(tstim(400,:),tsdim,tsdim)); axis square; 
subplot(3,1,2); plot2drfmat(reshape(tstim(1000,:),tsdim,tsdim));axis square; 
subplot(3,1,3); plot2drfmat(reshape(tstim(1300,:),tsdim,tsdim));axis square; 


