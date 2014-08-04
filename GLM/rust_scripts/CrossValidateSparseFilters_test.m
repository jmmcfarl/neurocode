nfold = 5; 
tuids = {'33-27','33-44','43-03','43-21','44-29','52-18'}
icell=4
tuid = tuids{icell}

cd ~/Data/rust/
dataname = sprintf('stimspikedat-%s.mat',tuid);     stcdat = load(dataname); 
modmatname = sprintf('rawSparseMods_jmm-%s.mat',tuid);  moddat = load(modmatname); 

rmodsLoc= moddat.rmodsLoc; %series of filter banks for localized filters
rmodsSTC= moddat.rmodsSTC; %series of filter banks for STC filters

stcstim = stcdat.stcstim; stcspks = stcdat.stcspks; 
flen = stcdat.pars.flen; foff = flen + stcdat.pars.hilen; lltolf =1e-4;


%%%% Training and Crossvalidation %%%%
stclen  = size(stcstim,1); 
partlen = ceil(stclen/nfold); 
partids = 1+floor(((1:stclen)-1)/partlen); %which CV part does each sample belong to
tabulate(partids)

%initialize a series of cell arrays, one for each CV iteration
[trllsLocs,xvllsLocs,trllsSTCs,xvllsSTCs,afmodsSTC,afmodsLoc] = deal(cell(nfold,1));

for ixvsamp = 1:nfold; disp(sprintf('*** xv samp %i',ixvsamp)); 
	%% crossval 
	xvids    = find(partids==ixvsamp); %find current CV chunk
    xvstart = xvids(1);
    xvend = xvids(end);
    xvstim   = stcstim(xvids,:);%CV stimulus matrix
    xvlen = size(xvstim,1); %number of samples in CV matrix
	xvspks   = stcspks((stcspks >= xvstart) & (stcspks<xvend))-xvstart+1; %recompute spike binning
	xvspksKH = xvspks(xvspks>foff & xvspks<xvlen)-foff+2; %shift spike binning

	%% put together rest for training
	trsamps = 1:nfold; trsamps(ixvsamp) =[]; trsamps
	atrspks =[]; atrspsksKH = []; trstimlen =0; trstim=[]; %initialize training data
	for itrsamp = 1:length(trsamps); %loop over training chunks
		sampid    = trsamps(itrsamp); %current chunk ID
		trsids    = find(partids==sampid); %sample inds for current chunk
        trsstart = trsids(1); trsend = trsids(end); 
		trstim  = [trstim;stcstim(trsstart:trsend,:)]; %concatenate stimulus chunk
		trspks    = stcspks((stcspks >= trsstart) & (stcspks<trsend))-trsstart +1;%shifted current spikes
		atrspks   = [atrspks; trspks + trstimlen]; %compile training spikes
		trstimlen = trstimlen + length(trsids); %increment length of training vector
    end
	trlen    = size(trstim,1);
	trspksKH = atrspks(atrspks>foff & atrspks< trstimlen)-foff+2;	
	
	%% fit models
	fmodsLoc  = cellfun(@(x)nfitNLHIc(x,trstim,trspksKH,lltolf),rmodsLoc(1:5),'UniformOutput',0);
	fmodsSTC  = cellfun(@(x)nfitNLHIc(x,trstim,trspksKH,lltolf),rmodsSTC(1:5),'UniformOutput',0);
	
	afmodsSTC{ixvsamp} = fmodsSTC; afmodsLoc{ixvsamp} = fmodsLoc; 
	w = cellfun(@(x)x.ll,fmodsLoc); trllsLocs{ixvsamp} = w; 
	x = cellfun(@(x)x.ll,fmodsSTC);  trllsSTCs{ixvsamp} = x;
	%% cross validated likelihoods
	y = cellfun(@(x)getLLGLM(x,xvstim,xvspksKH),fmodsLoc); xvllsLocs{ixvsamp} = y; 
	z = cellfun(@(x)getLLGLM(x,xvstim,xvspksKH),fmodsSTC); xvllsSTCs{ixvsamp} = z; 
	[w;y]
	[x;z]	
end; 


cellxvname = sprintf('Constr-XV--%s.mat',tuid)
save(cellxvname,'trllsLocs','xvllsLocs','trllsSTCs','xvllsSTCs','afmodsSTC','afmodsLoc','nfold'); 

