function [pred_rate, spks, G, gint, fgint] = NIMsimulate( nim, Xstim, XLin, Nreps)
%
% [pred_rate, spks, G, gint, fgint] = NIMsimulate( nim, Xstim, <XLin>, <Nreps> )
%
% Simulates spike trains of the specified model, also returns other useful outputs
%
% INPUTS:
%   nim: model structure
%   Robs: binned spikes
%   Xstim: time-embedded stimulus mat
%   <XLin>: Matrix specifying additional linear predictors
%
% OUTPUTS:
%   pred_rate: predicted firing rate without spike history term. (unitless, so divide by dt to get Hz).
%       Note that if there is a spike history term, it does not take into account, and insteads 'spks' 
%       output should be used to estimate rate (PSTH-style)
%   spks: spike times generated from the simulation. Multiple repeats are stored in one list, with
%       a '-1' separating each repeat.
%   G: generating function (output of the model before the spk NL)
%   gint: TxNmods matrix of the output of each subunit's stimulus filter (before applying upstream NL)
%   fgint: TxNmods matrix of the output of each subunit (after applying upstream NL)

if nargin < 4
    XLin = [];
end
if (nargin < 4) || isempty(Nreps)
	Nreps = 1;
end

%% Key parameters
NT = size(Xstim,1);

Nmods = length(nim.mods);
lin_dims = nim.stim_params.lin_dims;
spkhstlen = nim.spk_hist.spkhstlen;
dt = nim.stim_params.dt;

%% CREATE L2 REGULARIZATION MATRICES
L2_mats = create_L2_matrices(nim);

%% ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
theta = nim.spk_NL_params(1); %offset
G = theta + zeros(NT,1); %initialize overall generating function G

Kmat = [nim.mods(:).filtK];
gint = Xstim*Kmat; %subunit generating functions

fgint = nan(NT,Nmods);
for n = 1:Nmods
    
    %process subunit g's with upstream NLs
    if strcmp(nim.mods(n).NLtype,'nonpar')
        fgint(:,n) = piecelin_process(gint(:,n),nim.mods(n).NLy,nim.mods(n).NLx);
    elseif strcmp(nim.mods(n).NLtype,'quad')
        fgint(:,n) = gint(:,n).^2;
    elseif strcmp(nim.mods(n).NLtype,'lin')
        fgint(:,n) = gint(:,n);
    elseif strcmp(nim.mods(n).NLtype,'threshlin')
        fgint(:,n) = gint(:,n);
        fgint(fgint(:,n) < 0,n) = 0;
    else
        error('Invalid internal NL');
    end
    
    %multiply by weight and add to generating function
    G = G + fgint(:,n)*nim.mods(n).sign;
end

%add contribution from linear filter
if lin_dims > 0
    G = G + XLin*nim.kLin;
end

% Calculate predicted rate (without spike history)
if strcmp(nim.spk_NL_type,'logexp')
    max_gbeta = 50; %to prevent numerical overflow
    bgint = G*nim.spk_NL_params(2); %g*beta
    %expg = exp(bgint);
    too_large = (bgint > max_gbeta);
    pred_rate = nim.spk_NL_params(3)*log(1+exp(bgint)); %alpha*log(1+exp(gbeta))
    pred_rate(too_large) = nim.spk_NL_params(3)*bgint(too_large); %log(1+exp(x)) ~ x in limit of large x
		logexp = 1;
elseif strcmp(nim.spk_NL_type,'exp')
    %expg = exp(G);
    pred_rate = exp(G);
		logexp = 0;
else
    error('invalid spk nl');
		return
end


% Simulate Poisson spikes
spks = [];
if spkhstlen == 0
  % no spike history term, just need Poisson generator
	for rep = 1:Nreps
		spkstemp = find(rand(NT,1) < pred_rate)*dt - dt/2; 
	  spks = [spks' spkstemp' -1]';
	end	
else
	
	% then generating function affected by generated spikes
  Lh = nim.spk_hist.bin_edges(end); 
	h = zeros(1,Lh); % spike-history term
	for n = 1:nim.spk_hist.spkhstlen
		h(nim.spk_hist.bin_edges(n):(nim.spk_hist.bin_edges(n+1)-1)) = nim.spk_hist.coefs(n);
	end
	
	% Simulate over time (all reps at once)
	spkstemp = zeros(NT+Lh,Nreps);  % add buffer at beginning for spike history
	for t = 1:NT

		Gspkhist = ones(1,Nreps) * G(t) + h * spkstemp(Lh+t-(1:Lh),:);

		if logexp > 0
			r = nim.spk_NL_params(3)*log(1+exp(Gspkhist*nim.spk_NL_params(2)));
		else
			r = exp(Gspkhist);
		end
		spkstemp(t+Lh,:) = rand(1,Nreps) < r;

	end
	
	for n = 1:Nreps
		spks = [spks' ((find(spkstemp(:,n) > 0)'-Lh)*dt - dt/2) -1]';
	end
end

fprintf( 'Simulated %d repeats. Average firing rate = %0.2f Hz.\n', Nreps, (length(spks)-Nreps)/Nreps/(NT*dt) )
if Nreps == 1
	spks = spks(1:end-1); % take the -1 off the end if only one rep
end

 