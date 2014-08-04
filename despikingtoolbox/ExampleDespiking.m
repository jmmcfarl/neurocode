%% An example of wideband signal despiking with preprocessing steps 
%  Load data
load('examplelfpdespiking.mat');

%%
%Here wb is the wideband signal sampled at 10kHz and cluster_class is the
%class assigned to the nth detected spike by the Quiroga sorting algorithm
%In this case there is only one identified waveform type

Fs = 10000;
spktimes = cluster_class(cluster_class(:,1)==1,2)*10;
plot((1:length(wb))/Fs,wb,spktimes/Fs,2000*ones(size(spktimes)),'.');
xlabel('time (seconds)');
ylabel('voltage');
title('Wideband signal before artifact removal and location of spikes');

%%
%First get rid of the large spikes which are artifacts of the reward system
%We do this by calling despikeLFP
%First for positive artifacts
sartifacts = find([0;diff(wb)] > 850 | ([0;diff(wb)] > 350 & wb > 2000));
sartifacts = sartifacts([0;diff(sartifacts)]~=1);
%assume the artifact runs from -5 before the peak artifact time to +20
sartifacts = sartifacts - 5; 

spktrainartifacts = zeros(length(wb),1);
spktrainartifacts(sartifacts) = 1;

% Same with negative artifacts
sartifacts = find([0;diff(wb)] < -850 | ([0;diff(wb)] < -350 & wb < -2000));
sartifacts = sartifacts([0;diff(sartifacts)]~=1);
%assume the artifact runs from -5 before the peak artifact time to +20
sartifacts = sartifacts - 5; 

spktrainartifacts2 = zeros(length(wb),1);
spktrainartifacts2(sartifacts) = 1;


opts.displaylevel = 2;
%Here the prior is set to ones (in the Fourier domain), corresponding to an
%assumption of whiteness in the time domain. This is fine for removing
%artifacts that shouldn't be correlated with the local field potential
prior = ones(size(wb));
resultsreward = despikeLFP(wb,[spktrainartifacts,spktrainartifacts2],eye(25),prior,opts);

%And plot
plot((1:length(wb))/Fs,wb,(1:length(wb))/Fs,resultsreward.z);
%Not perfect because of the clipping at +- 2048 (the 12-bit boundary) means 
%the artifacts aren't perfectly stereotyped. Still, better than nothing.
wbpost1 = resultsreward.z;


%%
%Now remove line noise peaks at 60, 180Hz
wbpost2 = delineSignal(wbpost1,10000,[60,180],2);

%%
%Find a good value for g according to the method shown in the appendix,
%fitting to a function with a modest number of free parameters
%Usually this takes some fiddling with the parameters
g = fitLFPpowerSpectrum(wbpost2,.01,250,10000);

%Export figure
%set(gcf,'Position',[520 475 458 423]);
%exportfig(gcf,'figs/suppfig3.eps','Color','cmyk');


%%
%Finally create the spike train
%Plot a random subset of spikes
snipidx = bsxfun(@plus,spktimes,-50:50);
snippets = wbpost2(snipidx);
plot(-50:50,snippets(rand(size(snippets,1),1)<.1,:))

%We'd like to assume that a spike lasts
%from -10 samples to +19 samples (3 ms) compared to its peak
%Since spktimes are the time of the peak of each spike, we subtract 15
%from spktimes to obtain the start times of the spikes
Bs = eye(30);

S = zeros(length(wbpost2),1);
S(spktimes - 10) = 1;

%%
%With the prep work done, finally despike the LFP
%This takes 20 seconds on my Windows computer, 7 seconds on the 64-bit
%Linux i7-920 computers we use for array data

%Full debugging info
opts.displaylevel = 2;
tic;
results = despikeLFP(wbpost2,S,Bs,g,opts);
toc;
%%
%Make sure to plot the obtained spike waveforms to see if the peak is in
%the right location. 
plot(results.phi);

%It is.

%%
%Look at a before and after
rg = (Fs*4:Fs*6);
subplot(2,2,1);
plot(rg/Fs,wbpost2(rg),rg/Fs,results.z(rg)+300,...
    spktimes(spktimes/Fs > 4 & spktimes/Fs < 6)/Fs,50*ones(nnz(spktimes/Fs > 4 & spktimes/Fs < 6),1),'.');
xlim([5,5.2])
ylim([-1750,250]);
box off;

xlabel('Time (s)')
ylabel('Voltage')

legend('WB before','WB after (offset)','spikes');
legend(gca,'Location','SouthWest');
legend boxoff;

title('Wideband signal, before and after');

%Looking good. Note that since the method relies on the spike times given
%to it, it doesn't remove the spike at 5.02 seconds, and that the method
%removes just the mean spike from every spike, so that some remnants are still
%apparent when spikes are not completely stereotyped, for example the second 
%spike at sample 5.13s
%%
%Show the same results in a different way
%Smooth with a butterworth filter with same specs as Plexon system
[z,p,k] = butter(4,170/5000);
[sos,gg] = zp2sos(z,p,k);    % Convert to SOS form
Hd = dfilt.df2tsos(sos,gg);    % Create a dfilt object
[B,A]=sos2tf(Hd.sosMatrix,Hd.Scalevalues);

[z,p,k] = butter(1,.7/5000,'high');
[sos,gg] = zp2sos(z,p,k);    % Convert to SOS form
Hd = dfilt.df2tsos(sos,gg);    % Create a dfilt object
[B2,A2]=sos2tf(Hd.sosMatrix,Hd.Scalevalues);

smoothb = filtfilt(B2,A2,filtfilt(B,A,wbpost2-mean(wbpost2)));
smootha = filtfilt(B2,A2,filtfilt(B,A,results.z));

subplot(2,2,3);


plot((1:length(wb))/Fs,smoothb,(1:length(wb))/Fs,smootha,find(S)/Fs,200*ones(size(spktimes)),'.');
xlim([5,5.2])
ylim([-500,250]);
box off;

xlabel('Time (s)')
ylabel('Voltage')

title('LFP signal, before and after');

%It seems like it's a very small difference between the LFP with and
%without the despiking. But notice how the difference is always of the same
%sign. One could infer spikes from the LFP from this downward dip with high
%accuracy, and thus conclude that the LFP signal influences the probability
%of spiking to larger degree than in reality. We'll confirm this later

%%
%Look at the STA of the wideband signal before and after despiking
subplot(2,2,2);
rg = -250:250;
plot(rg/Fs*1000,sta(wbpost2-mean(wbpost2),spktimes,min(rg),max(rg)),...
     rg/Fs*1000,sta(results.z,spktimes,min(rg),max(rg)))
ylim([-300,-100]);
xlim([-25,25])
box off;
%Pretty huge difference. The STA before despiking goes down to about -500 

xlabel('Time (ms)')
ylabel('Voltage')

title('STA of wideband signal before and after');


%%
%Finally the STA of the LFP before and after despiking

subplot(2,2,4);
rg = -300:300;
plot(rg/Fs*1000,sta(smoothb,spktimes,min(rg),max(rg)),...
     rg/Fs*1000,sta(smootha,spktimes,min(rg),max(rg)))
ylim([-95,-45]);
xlim([-25,25]);
box off;

xlabel('Time (ms)')
ylabel('Voltage')

title('STA of LFP before and after');

% export figure
%set(gcf,'Position',[135   359   738   546]);
%exportfig(gcf,'figs/suppfig1.eps','Color','cmyk');


%%
%Low-pass filtering helps reduce the spike artifact in the spike-triggered 
%LFP. But a linear filter can always undo the low-pass filtering, so that 
%the apparent predictability of spikes from LFPs with a linear model will 
%always be biased, regardless of the filter. This observation is complicated 
%by the downsampling usually used in LFP recording systems, which is a
%destructive operation and therefore could help remove spike artifacts.
%Let's see for ourselves

%Downsample signals to 500 Hz
R = 20;
lowLFPbefore = decimate(smoothb,20);
lowLFPafter  = decimate(smootha,20);
lowSpikes    = hist(spktimes/Fs,linspace(0,length(wb)/Fs,length(wb)/R+1))';

%Let's assume that spikes can be predicted by convolving the LFP
%spikes = conv(LFP,filter) + offset
%And solve for filter using least squares

%Functions to create time-lagged matrices
zeropad = @(x,n) [zeros((n-1)/2,size(x,2));x;zeros((n-1)/2,size(x,2))];
subsel  = @(x,y) x(y);
timelag = @(x,n) subsel(zeropad(x,n),bsxfun(@plus,(1:length(x))',(0:n-1)));

%We're assuming the filter goes from -50 ms to 50 ms
Xbefore = [timelag(lowLFPbefore,51),ones(size(lowLFPbefore))];
Xafter =  [timelag(lowLFPafter ,51),ones(size(lowLFPbefore))];

%Solve for the filter in both cases using least-squares on the first half
%of the data
kernel1 = Xbefore(1:end/2,:)\lowSpikes(1:end/2);
kernel2 =  Xafter(1:end/2,:)\lowSpikes(1:end/2);

%look at the filters
plot(-50:2:50,kernel1(1:end-1),-50:2:50,kernel2(1:end-1));
%Notice the large difference. Also notice that the filters are very
%high-pass because the spike trains are much higher frequency than the LFPs

%Look at the cross-validated coefficient coefficient in both cases
corr(Xbefore(end/2+1:end,:)*kernel1,lowSpikes(end/2+1:end))
corr( Xafter(end/2+1:end,:)*kernel2,lowSpikes(end/2+1:end))

%There's a three fold drop in cross-validated r, or equivalently a 
%ten-fold drop in explained variance. Clearly most of the predictability of
%spikes from the LFP in this case is purely artifactual.

%%
%We could also work in chunks if we wanted to. Here the data is too short
%to justify chunking, but nevertheless for reference here's how you can do that
%Assume a chunk size of 2^18
firstchunk = wbpost2(1000:1000+2^18-1);

%Find a good value for g for this small chunk
gchunk = fitLFPpowerSpectrum(firstchunk,.01,250,10000);

%%
%Now despike chunk-wise. The size of the g fed to the function determines 
%the size of chunks
opts.displaylevel = 1;
resultsc = despikeLFPbyChunks(wbpost2,S,Bs,gchunk,opts);

%%
%Now plot the recovered waveforms
subplot(1,2,1);plot(resultsc.phi');title('waveforms for each chunk');
subplot(1,2,2);plot(bsxfun(@minus,resultsc.phi',results.phi));title('difference b/w waveforms estimated in each chunk with the one waveform estimated without chunking');

%The differences are pretty marginal

%% As for z
plot((1:length(wb))/Fs,results.z-resultsc.z)

%Again the difference is marginal, considering the range of z is +-2000
%%
%Finally, show how you can express spikes in a non-Dirac basis
%This is useful if you have few spikes or high sampling rates so that
%the waveforms estimated by the model are noisy. 

%Clearly this is not the case here but just for the sake of showing an
%example, create a smooth basis undercomplete by a factor 2
Bs2 = conv2(Bs,fspecial('gauss',[1,7],1),'same');
Bs2 = Bs2(:,1:2:end);

plot(Bs2);

%This is of course a dumb basis here because the real waveform is highly
%peaked
%%
%And now despike in this new basis
opts.displaylevel = 2;
resultsb = despikeLFP(wbpost2,S,Bs2,g,opts);

%%
%And plot the resulting phi
plot((-10:19)*.1,Bs2*resultsb.phi,(-10:19)*.1,Bs*results.phi,(-10:19)*.1,Bs2*pinv(Bs2)*Bs*results.phi)

%The red line here is the original waveform projected onto the Bs2 basis in
%the optimal manner, and it is quite similar to waveform obtained directly
%in the Bs2 matrix, but not quite the same because now the model has less
%free parameters, and  assigns features to the spike versus the
%LFP in a slightly different manner

%%
%And now for a demonstration of the smoothness prior's role in the
%procedure

%Create a surrogate signal by randomizing the phase of the wideband signal,
%and adding that to a train of spike waveforms and a train of Gabors of a
%longer time scale than spike waveforms. This is a lazy way creating an
%LFP-like signal with a legitimate linear spike-LFP relationship

%Create the Gabor train
rg = (-300:300)';
sigma = 120;
wl = 240;
cg = 5;
gabor = -exp(-(rg-cg).^2/sigma^2).*cos((rg-cg)/wl*2*pi);
gabor = gabor - exp(-(rg-cg).^2/sigma^2)/mean(exp(-(rg-cg).^2/sigma^2))*mean(gabor);
st = zeros(length(wb),1);
st(spktimes+5) = 1;
gabortrain = conv(st,gabor,'same');

%Create the spike train
wftrain = conv(st,results.phi,'same');

%Randomize the phase of the wideband signal by convolving with white
%noise
a = fft(randn(length(wb),1));
a = a/mean(abs(a));
wbnoise = ifft(a.*fft(wb));

%Set the gain of the Gabor to 1/10th of the gain of the spike waveform
ggain = max(abs(results.phi))/10;
wbgab = wbnoise + wftrain + ggain*gabortrain;
wbgab = wbgab - mean(wbgab);

opts.displaylevel = 1;

%And now despike under an assumption of smoothness and no assumption of
%smoothness (g = 1, equivalent to a weight decay prior)

gabresultwith = despikeLFP(wbgab,S,Bs,g,opts);
gabresultwithout = despikeLFP(wbgab,S,Bs,ones(size(g)),opts);


%%
%Now show a with/without assumptions STA, with/without smoothing
smoothwith = filtfilt(B2,A2,filtfilt(B,A,gabresultwith.z));
smoothwithout = filtfilt(B2,A2,filtfilt(B,A,gabresultwithout.z));

st2 = zeros(size(st));
st2(find(st)-2-20-12) = 1;

rg = -400:400;

subplot(2,2,1);
plot(1000*rg/Fs,sta(gabresultwithout.z,spktimes,min(rg),max(rg)),...
     1000*rg/Fs,sta(gabresultwith.z,spktimes,min(rg),max(rg)));
xlabel('Time (ms)')
ylabel('Voltage')
xlim([-40,40]);

legend('without','with');
legend(gca,'Location','SouthWest');
legend boxoff;

title('STA of wideband signal');
box off;
 
subplot(2,2,3);
plot(1000*rg/Fs,sta(smoothwithout,spktimes,min(rg),max(rg)),...
     1000*rg/Fs,sta(smoothwith,spktimes,min(rg),max(rg))); 
xlabel('Time (ms)')
ylabel('Voltage')
xlim([-40,40]);

title('STA of LFP');
box off;

subplot(2,2,2);
plot((-10:19)/10,gabresultwithout.phi,(-10:19)/10,gabresultwith.phi); 
xlabel('Time (ms)')
ylabel('Voltage')
xlim([-1,1.9]);

title('Spike waveforms');
box off;

subplot(2,2,4);
plot((-10:19)/10,gabresultwith.phi-gabresultwithout.phi); 
xlabel('Time (ms)')
ylabel('Voltage')
xlim([-1,1.9]);

title('Difference of spike waveforms');
box off;

%The model does not properly interpolate at 0 without the assumption of 
%smoothness and ends up assigning the mean value around the spike time to
%the spike and not the LFP. You run into the same trouble with classic
%methods like removing the STA around the time of every spike. You can try
%and add some heuristics to deal with the problem (removing a mean and
%linear trend in the spike, for example), but it gets messy.

%Note: although the STA appears to have a non-zero mean, this is because
%we're only looking from -40 to +40 ms. In actuality, both gabresultwith.z
%and gabresultwithout.z have a mean of zero.

%%
%Export figure
%{
opts = struct('LineMode','fixed','LineWidth',1);
set(gcf,'Position',[135   359   738   546]);
exportfig(gcf,'figs/suppfig2.eps',opts,'Color','cmyk');
%}
