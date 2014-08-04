function [LLs, rawLLs, penLLs] = NIMeval_LLreps( nim, XstimR, RobsR, Xlin )
%
% Usage: [LLs rawLLs penLLs] = NIMeval_LLreps( nim, XstimR, RobsR, <Xlin> )
%
% Computes the LL/spk in units of nats of the specified model for repeated trials. 
%
% INPUTS:
%   nim: model structure
%   XstimR: time-embedded stimulus mat
%   <XLin>: Matrix specifying additional linear predictors
%   RobsR: NT x Nreps matrix representing histogram of spike times
%
%   As another option, a list of spike times, with repeats separated by '-1' can 
%   be included in Robs place. Note extra code below to take care of that.
%
% OUTPUTS:
%   LL: log-likelihood (per spike) of the model, adjusted by the 'null model' (constant firing rate)
%   rawLLs: unadjusted LLs
%   penLLs: penalized log-likelihood (per spike)

if nargin < 4
	Xlin = [];
end

NT = size(XstimR,1);
[NT2 Nreps] = size(RobsR);
dt = nim.stim_params.dt;
T = NT*dt;

%% If passed spike list instead of RobsR (see above
if NT2 ~= NT
	% then assume this is a spike list instead, and convert to RobsR
	spksR = RobsR(:)';
	Nreps = length(find(spksR < 0));
	replocs = [0 find(spksR < 0)];
	RobsR = zeros(NT,Nreps);
	for rep = 1:Nreps
		rspks = spksR( (replocs(rep)+1):(replocs(rep+1)-1) );
		Robs = histc(rspks, 0:dt:T );
		% Remove last bin, because histc throws in an extra one
		RobsR(:,rep) = Robs(1:end-1);
	end	
end

%% Calculate log-likelihoods given RobsR
LLs = zeros(Nreps,1);
rawLLs = LLs; penLLs = LLs;

for rep = 1:Nreps
	% Extract spikes from that repeat
	
	[rawLLs(rep) penLLs(rep)] = NIMmodel_eval( nim, RobsR(:,rep), XstimR, Xlin );
	
	LLnull = log(sum(RobsR(:,rep))/T*dt) - 1;  % formula for LL/spk of null model
	LLs(rep) = rawLLs(rep)-LLnull;
end
