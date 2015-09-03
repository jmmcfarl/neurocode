
classdef NIM
    
    %Summary of this class goes here
    %   Detailed explanation goes here
    %
    % James McFarland, September 2015
    properties
        spkNL;          %struct defining the spiking NL function
        subunits;       %array of subunit objects
        stim_params;    %struct array of parameters characterizing the stimuli that the model acts on, must have a .dims field
        noise_dist;     %noise distribution class specifying the noise model
        spk_hist;       %class defining the spike-history filter properties
        fit_props;      %struct containing information about model fit evaluations
        init_props;     %struct containing details about model initialization
    end
    
    methods
        function nim = NIM(stim_params, NLtypes, mod_signs, varargin)
            %            nim = NIM(stim_params, NLtypes, mod_signs, <Xtargets>,<spkNL>,<noise_dist>)
            %constructor for class NIM
            %            INPUTS:
            %                 stim_params: struct array defining parameters for each stimulus the model acts on. Must specify the .dims field for each stim
            %                 NLtypes: string or cell array of strings specifying the
            %                 upstream NL type associated with each subunit. If it's a
            %                 single string, we use the same NL type throughout
            %                 mod_signs: vector specifying the weight associated with each subunit (typically +/- 1)
            %                 <Xtargets>: vector specifying the index of the stimulus each subunit acts on (defaults to ones)
            %                 <spkNL>: string specifying type of spkNL function
            %                 <noise_dist>: string specifying type of noise distribution
            
            nStims = length(stim_params); %number of stimuli
            assert(isfield(stim_params,'dims'),'need to input field dims to stim_params');
            nim.stim_params = stim_params;
            
            %COULD ADD CHECK STIM_PARAMS FUNCTION
            
            nSubs = length(mod_signs); %number of subunits
            
            %if NLtypes is specified as a single string, default to using
            %this NLtype for all subunits
            if ~iscell(NLtypes) && isstr(NLtypes)
                NLtypes = repmat({NLtypes},nSubs,1);
            else
                error('invalid format for NLtypes');
            end
            
            %set defaults
            Xtargets = ones(nSubs,1);
            spkNL = 'softplus';
            noise_dist = 'poisson';
            
            %parse input flags
            assert(mod(nargin - 3,2)==0,'input format for optional flags must be in pairs: ''argname'',''argval');
            j = 1; %initialize counter after required input args
            while j <= length(varargin)
                flag_name = varargin{j};
                flag_val = varargin{j+1};
                switch lower(flag_name)
                    case 'xtargets'
                        Xtargets = flag_val;
                        assert(all(ismember(Xtargets,1:nStims)),'invalid Xtargets specified');
                    case 'spknl'
                        spkNL = flag_val;
                        assert(isstr(spkNL),'spkNL must be a string');
                    case 'noise_dist'
                        noise_dist = lower(flag_val);
                        assert(isstr(noise_dist),'noise_dist must be a string');
                    otherwise
                        error('Invalid input flag');
                end
                j = j + 2;
            end
            
            %check and create spk NL function
            allowed_NLs = {'lin','rectlin','exp','softplus'}; %set of NL functions currently implemented
            assert(ismember(spkNL,allowed_NLs),'not an allowed spk NL type');
            nim.spkNL.type = spkNL;
            nim.spkNL.theta = 0; %initialize offset term
            %set default parameters for other spkNL parameters depending on
            %the type
            switch nim.spkNL.type
                case 'lin'
                    nim.spkNL.params = [1]; %defines beta in f(x; beta) = beta*x
                case 'rectlin'
                    nim.spkNL.params = [1]; %defines beta in f(x; beta) = beta*x iff x > 0
                case 'exp'
                    nim.spkNL.params = [1]; %defines beta in f(x; beta) = exp(beta*x)
                case 'softplus'
                    nim.spkNL.params = [1 1]; %defines beta, alpha in f(x; beta, alpha) = alpha*log(1+exp(beta*x))
            end
            
            %check and create noise distribution
            allowed_noise_dists = {'poisson','bernoulli','gaussian'}; %allowed noise distributions
            assert(ismember(noise_dist,allowed_noise_dists),'not an allowed noise distribution');
            nim.noise_dist = noise_dist;
            
            assert(length(Xtargets) == nSubs,'length of mod_signs and Xtargets must be equal');
            
            nim.init_props = rng(); %save state of RNG used for initializing filter weights
            for ii = 1:nSubs %loop initializing subunits
                stimD = prod(nim.stim_params(Xtargets(ii)).dims); %dimensionality of the current filter
                init_filt = randn(stimD,1)/sqrt(stimD); %initialize fitler coefs with gaussian noise
                nim.subunits(ii) = SUBUNIT(init_filt, mod_signs(ii), NLtypes{ii});
            end
        end
        
        
        function LL = internal_LL(nim,rPred,rObs)
            %internal evaluatation method for computing the total LL associated with the predicted rate rPred, given observed data rObs
            switch nim.noise_dist
                case 'poisson' %LL = Rlog(r) - r + C
                    LL = sum(rObs .* log2(rPred) -rPred);
                    
                case 'bernoulli' %LL = R*log(r) + (1-R)*log(1-r)
                    LL = sum(rObs.*log(rPred) + (1-rObs).*log(1-rPred));
                    
                case 'gaussian' %LL = (r-R)^2 + c
                    LL = sum((rPred - rObs).^2);
            end
        end
        
        function LL_deriv = internal_LL_deriv(nim,rPred,rObs)
            %internal method for computing the derivative of the LL wrt the
            %predicted rate at rPred, given rObs (as a vector over time)
            switch nim.noise_dist
                case 'poisson' %LL'[r] = R/r - 1
                    LL_deriv = rObs./rPred - 1;
                    
                case 'bernoulli' %LL'[r] = R/r - (1-R)/(1-r)
                    LL_deriv = rObs./rPred - (1-rObs)./(1-rPred);
                    
                case 'gaussian' %LL'[r] = 2*(r-R)
                    LL_deriv = 2*(rPred - rObs);
            end
            
        end
        
        function rate = apply_spkNL(nim,gen_signal)
            %apply the spkNL function to the input gen_signal. NOTE: the
            %offset term should already be added to gen_signal
            switch nim.spkNL.NLtype
                case 'lin' %F[x;beta] = beta*x
                    rate = gen_signal*nim.spkNL.params(1);
                    
                case 'rectlin' %F[x; beta] = beta*x iff x > 0; else 0
                    rate = gen_signal*nim.spkNL.params(1);
                    rate(rate < 0) = 0;
                    
                case 'exp' %F[x; beta] = exp(beta*x)
                    rate = exp(gen_signal*nim.spkNL.params(1));
                    
                case 'softplus' %F[x; beta, alpha] = alpha*log(1+exp(beta*x))
                    max_g = 50; %to prevent numerical overflow
                    gint = gen_signal*nim.spkNL.params(1);
                    rate = nim.spkNL.params(2)*log(1 + exp(gint));
                    rate(gint > max_g) = nim.spkNL.params(2)*gint(gint > max_g);
            end
        end
        
        function rate_deriv = apply_spkNL_deriv(nim,gen_signal)
            %apply the derivative of the spkNL to the input gen_signal.
            %Again, gen_signal should have the offset theta already added
            %in
            switch nim.spkNL.NLtype
                case 'lin' %F'[x; beta] = beta;
                    rate_deriv = nim.spkNL.params(1)*ones(size(gen_signal));
                    
                case 'rectlin' %F'[x; beta] = beta iff x > 0; else 0
                    rate_deriv = nim.spkNL.params(1)*(gen_signal > 0);
                    
                case 'exp' %F'[x; beta] = beta*exp(beta*x)
                    rate_deriv = nim.spkNL.params(1)*exp(nim.spkNL.params(1)*gen_signal);
                    
                case 'softplus' %F[x; beta, alpha] = alpha*beta*exp(beta*x)/(1+exp(beta*x))
                    max_g = 50; %to prevent numerical overflow
                    gint = gen_signal*nim.spkNL.params(1);
                    rate_deriv = nim.spkNL.params(1)*nim.spkNL.params(2)*exp(gint)./(1 + exp(gint));
                    rate_deriv(gint > max_g) = nim.spkNL.params(1)*nim.spkNL.params(2); %e^x/(1+e^x) => 1 for large x
            end
        end
        
    end
end