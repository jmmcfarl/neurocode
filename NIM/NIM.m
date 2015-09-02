
classdef NIM
    
    %Summary of this class goes here
    %   Detailed explanation goes here
    %
    % James McFarland, September 2015
    properties
        spkNL;          %spkNL object defining the spiking NL function
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

            nSubs = length(mod_signs);
            
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
                        noise_dist = flag_val;
                        assert(isstr(noise_dist),'noise_dist must be a string');
                    otherwise
                       error('Invalid input flag');
                end
                j = j + 2;
            end
            
            nim.spkNL = SPKNL(spkNL); %create spk NL function
            nim.noise_dist = NOISE_DIST(noise_dist); %create noise distribution model
            
            assert(length(Xtargets) == nSubs,'length of mod_signs and Xtargets must be equal');
            
            nim.init_props = rng(); %save state of RNG used for initializing filter weights
            for ii = 1:nSubs
                stimD = prod(nim.stim_params(Xtargets(ii)).dims); %dimensionality of the current filter
                init_filt = randn(stimD,1)/sqrt(stimD); %initialize fitler coefs with gaussian noise
                nim.subunits(ii) = SUBUNIT(init_filt, mod_signs(ii), NLtypes{ii});
            end
            
        end
    end
        
    
    
end