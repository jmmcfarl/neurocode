
classdef NIM
    
    % Class implementation of 'nonlinear-input model' (NIM). 
    %     Models consist of a set of linear-nonlinear processing
    %     'subunits', which act on different sets of predictors
    %     ('stimuli'). The summed output of these subunits is then
    %     transformed by a 'spiking nonlinearity' function to generate a
    %     predicted firing rate. Model parameters are optimized based on
    %     the penalized Log-likelihood. Several different noise models can
    %     be used.
    %   
    % Reference: 
    %         McFarland JM, Cui Y, Butts DA (2013) Inferring nonlinear neuronal computation based on physiologically plausible inputs. PLoS Computational Biology 9(7): e1003142
    % Created by James McFarland, September 2015
    
    %%
    properties
        spkNL;          %struct defining the spiking NL function
        subunits;       %array of subunit objects
        stim_params;    %struct array of parameters characterizing the stimuli that the model acts on, must have a .dims field
        noise_dist;     %noise distribution class specifying the noise model
        spk_hist;       %class defining the spike-history filter properties
        fit_props;      %struct containing information about model fit evaluations
        fit_hist;       %struct containing info about history of fitting
    end
    
    properties (Hidden)
        init_props;     %struct containing details about model initialization
        allowed_reg_types = {'nld2','d2xt','d2x','d2t','l2','l1'}; %set of allowed regularization types
        version = '0.0';
        min_pred_rate = 1e-50; %minimum predicted rate (for non-negative data) to avoid NAN LL values
        opt_check_FO = 1e-3; %threshold on first-order optimality for fit-checking
    end
    %%
    methods
        %% CONSTRUCTOR
        function nim = NIM(stim_params, NLtypes, mod_signs, varargin)
%         nim = NIM(stim_params, NLtypes, mod_signs, varargin) 
%         constructor for class NIM 
%            INPUTS:
%               stim_params: struct array defining parameters for each stimulus the model acts on.
%                   Must specify the .dims field for each stim 
%               NLtypes: string or cell array of strings specifying the upstream NL type associated
%                   with each subunit. If it's a single string, we use the same NL type throughout 
%               mod_signs: vector specifying the weight associated with each subunit (typically +/-1) 
%               optional_flags:
%                     ('Xtargets',Xtargs): vector specifying the index of the stimulus each subunit 
%                         acts on (defaults to ones) 
%                     ('spkNL',spkNL): string specifying type of spkNL function
%                     ('noise_dist',noise_dist): string specifying type of noise distribution 
%                     ('init_filts',init_filts): cell array of initial filter values 
%                     ('Ksign_cons',Ksign_cons): vector specifying any constraints on the filter
%                           coefs of each subunit. [-1 for neg; +1 for pos; nan for no cons]
%                     ('NLparams',NLparams): cell array of parameter values for the corresponding NL functions
%           OUTPUTS:
%                nim: initialized model object
                                     
            nStims = length(stim_params); %number of stimuli
            stim_params = NIM.check_stim_params(stim_params); %validate and format input stim_params
            nim.stim_params = stim_params;
            
            nSubs = length(mod_signs); %number of subunits
            %if NLtypes is specified as a single string, default to using
            %this NLtype for all subunits
            if ~iscell(NLtypes) && ischar(NLtypes); NLtypes = cellstr(NLtypes); end;
            if length(NLtypes) == 1 && nSubs > 1; NLtypes = repmat(NLtypes,nSubs,1); end
            
            %set defaults
            Xtargets = ones(nSubs,1);
            spkNL = 'softplus';
            noise_dist = 'poisson';
            init_filts = cell(nSubs,1);
            Ksign_cons = nan(nSubs,1);
            NLparams = cell(nSubs,1);
            
            %parse input flags
            assert(mod(nargin - 3,2)==0,'input format for optional flags must be in pairs: ''argname'',''argval');
            j = 1; %initialize counter after required input args
            while j <= length(varargin)
                switch lower(varargin{j})
                    case 'xtargets'
                        Xtargets = varargin{j+1};
                        assert(all(ismember(Xtargets,1:nStims)),'invalid Xtargets specified');
                        j = j + 2;
                    case 'spknl'
                        spkNL = lower(varargin{j+1});
                        assert(ischar(spkNL),'spkNL must be a string');
                        j = j + 2;
                    case 'noise_dist'
                        noise_dist = lower(varargin{j+1});
                        assert(ischar(noise_dist),'noise_dist must be a string');
                        j = j + 2;
                    case 'init_filts'
                        init_filts = varargin{j+1};
                        assert(iscell(init_filts),'init_filts must be a cell array');
                        j = j + 2;
                    case 'ksign_cons',
                        Ksign_cons = varargin{j+1};
                        assert(all(ismember(Ksign_cons,[-1 1 0])),'Ksign_cons must have values -1,0, or 1');
                        j = j + 2;
                    case 'nlparams'
                        NLparams = varargin{j+1};
                        assert(iscell(NLparams),'NLparams must be a cell array');
                        j = j + 2;
                    otherwise
                        error('Invalid input flag');
                end
            end
            
            %check and create spk NL function
            allowed_spkNLs = {'lin','rectquad','exp','softplus','logistic'}; %set of NL functions currently implemented
            assert(ismember(spkNL,allowed_spkNLs),'not an allowed spk NL type');
            nim.spkNL.type = spkNL;
            nim.spkNL.theta = 0; %initialize offset term
            %set default parameters for other spkNL parameters depending on
            %the type
            switch nim.spkNL.type
                case 'lin'
                    nim.spkNL.params = [1]; %defines beta in f(x; beta) = beta*x
                case 'rectquad'
                    nim.spkNL.params = [1]; 
                case 'exp'
                    nim.spkNL.params = [1]; %defines beta in f(x; beta) = exp(beta*x)
                case 'softplus'
                    nim.spkNL.params = [1 1]; %defines beta, alpha in f(x; beta, alpha) = alpha*log(1+exp(beta*x))
                case 'logistic'
                    nim.spkNL.params = [1]; %defines beta in f(x; beta) = 1/(1+exp(-beta*x))
            end
            
            %check and create noise distribution
            allowed_noise_dists = {'poisson','bernoulli','gaussian'}; %allowed noise distributions
            assert(ismember(noise_dist,allowed_noise_dists),'not an allowed noise distribution');
            nim.noise_dist = noise_dist;
            if strcmp(nim.noise_dist,'bernoulli') %make sure we have a logistic spk NL if using a bernoulli noise model
                assert(strcmp(nim.spkNL.type,'logistic'),'bernoulli noise model is only supported with a logistic spk NL function');
            end
            assert(length(Xtargets) == nSubs,'length of mod_signs and Xtargets must be equal');
            
            %initialize subunits
            nim.init_props = rng(); %save state of RNG used for initializing filter weights
            for ii = 1:nSubs %loop initializing subunits (start from last to initialize object array)
                stimD = prod(nim.stim_params(Xtargets(ii)).dims); %dimensionality of the current filter
                if isempty(init_filts{ii})
                    init_filt = randn(stimD,1)/stimD; %initialize fitler coefs with gaussian noise
                else
                    init_filt = init_filts{ii};
                end
                nim.subunits = cat(1,nim.subunits,SUBUNIT(init_filt, mod_signs(ii), NLtypes{ii},Xtargets(ii),NLparams{ii},Ksign_cons(ii)));
            end
            
            %initialize WITHOUT spike history term
            spk_hist.coefs = [];
            spk_hist.bin_edges = [];
            spk_hist.spkhstlen = 0;
            nim.spk_hist = spk_hist;
        end
        
        %% setting methods
        function nim = set_reg_params(nim, varargin)
%         nim = nim.set_reg_params(varargin)
%         set a desired set of regularization parameters to specified values, apply to specified set 
%         of subunits
%            INPUTS:
%               optional flags:
%               ('sub_inds',sub_inds): set of subunits to apply the new reg_params for
%               ('lambda_type',lambda_val): first input is a string specifying the type of
%                    regularization (e.g. 'd2t' for temporal smoothness). This must be followed by a
%                    scalar giving the associated lambda value
%            OUTPUTS:
%                nim: initialized model object

            sub_inds = 1:length(nim.subunits); %default is to apply the change to all subunits
            
            %INPUT PARSING
            j = 1;
            reg_types = {}; reg_vals = [];
            while j <= length(varargin)
                switch lower(varargin{j})
                    case 'sub_inds'
                        sub_inds =  varargin{j+1};
                        assert(all(ismember(sub_inds,1:length(nim.subunits))),'invalid target subunits specified');
                        j = j + 2;
                    case nim.allowed_reg_types
                        reg_types = cat(1,reg_types,lower(varargin{j}));
                        reg_vals = cat(1,reg_vals, varargin{j+1});
                        j = j + 2;
                    otherwise
                        error('Invalid input flag');
                end
            end
            
            if isempty(reg_vals)
                warning('No regularization values specified, no action taken');
            end
            for ii = 1:length(reg_vals)
                assert(reg_vals(ii) >= 0,'regularization hyperparameters must be non-negative');
                for jj = 1:length(sub_inds)
                    nim.subunits(sub_inds(jj)).reg_lambdas = setfield(nim.subunits(sub_inds(jj)).reg_lambdas,reg_types{ii},reg_vals(ii));
                end
            end
        end
        
        %%
        function nim = set_stim_params(nim, varargin)
%         nim = nim.set_stim_params(varargin)
%         set values of the stim_params struct for a desired 'Xtarget'
%           INPUTS:
%             optional flags:
%                ('xtarg',Xtarg): index of stimulus to apply the new stim_params for [default is 1]
%                ('dims', dims): dimensionality of stim: [Tdim, X1dim, X2dim] where Tdim is the
%                   number of temporal dimensions, etc.
%                ('boundary_conds',boundary_conds): boundary conditions on each stim dimension. Inf
%                    means free boundaryes, 0 means tied to 0, and -1 means periodic.
%                ('split_pts',split_pts): specify an 'internal boundary' over which you dont want to
%                   smooth. Must be a vector of the form: [direction split_ind boundary_cond], where
%                   direction is the index of the dimension of the stimulus you want to create a
%                   boundary (i.e. 1 for the temporal dimension), split_ind is the index value
%                   along that dimension after which you want to create the boundary, and
%                   boundary_cond specifies the boundary conditions you want to use.
%           OUPUTS: nim: new nim object

            Xtarg = 1; %default is to apply the change to stim 1
            allowed_flags = {'dims','boundary_conds','split_pts'}; %fields of the stim_params struct that we might want to set
            %INPUT PARSING
            j = 1; fields_to_set = {}; field_vals = {};
            while j <= length(varargin)
                switch lower(varargin{j})
                    case 'xtarg'
                        Xtarg = varargin{j+1};
                        j = j + 2;
                    case allowed_flags
                        fields_to_set = cat(1,fields_to_set,flag_name);
                        field_vals = cat(1,field_vals,varargin{j+1});
                        j = j + 2;
                    otherwise
                        error('Invalid input flag');
                end
            end
            
            if isempty(field_vals)
                warning('No stim_params values specified to change');
            end
            if Xtarg > length(nim.stim_params) %if we're creating a new stimulus
                dim_locs = find(ismember('dims',fields_to_set));
                assert(~isempty(dim_locs),'need to specify dims to initialize new stimulus');
                new_stim_params.dims = field_vals{dim_locs};
                new_stim_params = NIM.check_stim_params(new_stim_params); %initialize default params for new stimulus
                nim.stim_params(Xtarg) = new_stim_params;
            end
            for ii = 1:length(fields_to_set) %assign new field values
                nim.stim_params(Xtarg) = setfield(nim.stim_params(Xtarg),fields_to_set{ii},field_vals{ii});
            end
            while length(nim.stim_params(Xtarg).dims) < 3 %pad dims with 1s for book-keeping
                nim.stim_params(Xtarg).dims = cat(2,nim.stim_params(Xtarg).dims,1);
            end            
        end
        
        %% getting methods
        
        function lambdas = get_reg_lambdas(nim, varargin)
%         lambdas = nim.get_reg_lambdas(varargin)
%         get regularizatoin lambda values of specified type from a set of nim subunits
%           INPUTS:
%             optional flags:
%                ('sub_inds',sub_inds): vector specifying which subunits to extract lambda values from
%                lambda_type: string specifying the regularization type
%           OUTPUTS:
%              lambdas: [K,N] matrix of lambda values, K is the number of specified lambda_types and 
%                 N is the number of subunits
 
            sub_inds = 1:length(nim.subunits); %default is to grab reg values from all subunits
            %INPUT PARSING
            jj = 1;
            reg_types = {};
            while jj <= length(varargin)
                switch lower(varargin{jj})
                    case 'sub_inds'
                        sub_inds = varargin{jj+1};
                        assert(all(ismember(sub_inds,1:length(nim.subunits))),'invalid target subunits specified');
                        jj = jj + 2;
                    case nim.allowed_reg_types
                        reg_types = cat(1,reg_types,lower(varargin{jj}));
                        jj = jj + 1;
                    otherwise
                        error('Invalid input flag');
                end
            end
            
            lambdas = nan(length(reg_types),length(sub_inds));
            if isempty(reg_types)
                warning('No regularization type specified, returning nothing');
            end
            for ii = 1:length(reg_types)
                for jj = 1:length(sub_inds)
                    lambdas(ii,jj) = getfield(nim.subunits(sub_inds(jj)).reg_lambdas,reg_types{ii});
                end
            end
        end
        
        %%
        function filtKs = get_filtKs(nim,sub_inds)
%         filtKs = nim.get_filtKs(sub_inds)
%         get filters for specified set of subunits. 
%           INPUTS:
%             <sub_inds>: vector specifying which subunits to get filters from (default is all subs)
%           OUTPUTS:
%              filtKs: Cell array of filter coefs
            
            Nsubs = length(nim.subunits);
            if nargin < 2
                sub_inds = 1:Nsubs; %default is to grab filters for all subunits
            end
            filtKs = cell(Nsubs,1);
            for ii = sub_inds
                filtKs{ii} = nim.subunits(ii).get_filtK;
            end
        end
        
        %%
        function NLtypes = get_NLtypes(nim,sub_inds)
%          NLtypes = nim.get_NLtypes(<sub_inds>)
%          gets cell array of strings specifying NLtype of each subunit
%            INPUTS: 
%                 <sub_inds>: set of subunits to get NL types from (default is all)
%            OUTPUTS: NLtypes: cell array of NLtypes

            Nsubs = length(nim.subunits);
            if nargin < 2
                sub_inds = 1:Nsubs; %default is to grab filters for all subunits
            end
            NLtypes = cell(length(sub_inds),1);
            for ii = 1:length(sub_inds)
                NLtypes{ii} = nim.subunits(sub_inds(ii)).NLtype;
            end
        end
        
        %%
        function [filt_penalties,NL_penalties] = get_reg_pen(nim,Tmats)
%         [filt_penalties,NL_penalties] = nim.get_reg_pen(<Tmats>)
%          calculates the regularization penalties on each subunit, separately for filter and NL regularization
%            INPUTS: 
%                 <Tmats>: struct array of 'Tikhonov' regularization matrices
%            OUTPUTS: 
%                 filt_penalties: Kx1 vector of regularization penalties for each filter
%                 NL_penalties: Kx1 vector of regularization penalties for each upstream NL

            if nargin < 2 || isempty(Tmats) %if the Tmats are not precomputed and supplied, compute them here
                Tmats = make_Tikhonov_matrices(nim);
            end
            Nsubs = length(nim.subunits);
            Xtargs = [nim.subunits(:).Xtarg];
            filtKs = nim.get_filtKs();
            filt_penalties = zeros(1,Nsubs);
            for ii = 1:length(Tmats) %loop over the derivative regularization matrices
                cur_subs = find(Xtargs == Tmats(ii).Xtarg); %set of subunits acting on the stimulus given by this Tmat
                cur_penalties = sum((Tmats(ii).Tmat * cat(2,filtKs{cur_subs})).^2);
                cur_lambdas = nim.get_reg_lambdas(Tmats(ii).type,'sub_inds',cur_subs); %current lambdas
                filt_penalties(cur_subs) = filt_penalties(cur_subs) + cur_penalties.*cur_lambdas; %reg penalties for filters
            end
            l2_lambdas = nim.get_reg_lambdas('l2');
            if any(l2_lambdas > 0) %compute L2 penalties on the filter coefficients
                filt_penalties = filt_penalties + l2_lambdas.*cellfun(@(x) sum(x.^2),filtKs)';
            end
            nl_lambdas = nim.get_reg_lambdas('nld2');  %reg lambdas on the NL TB coefficients
            NL_penalties = zeros(1,Nsubs);
            if any(nl_lambdas > 0)
                Tmat = nim.make_NL_Tmat();
                nonpar_subs = find(strcmp(nim.get_NLtypes,'nonpar'))';
                for imod = nonpar_subs %compute the reg penalty on each subunit's NL
                    NL_penalties(imod) = nl_lambdas(imod)*sum((Tmat*nim.subunits(imod).TBy').^2);
                end
            end
        end
        
        %%
        function nim = add_subunits(nim,NLtypes,mod_signs,varargin)
%         nim = nim.add_subunits(NLtypes,mod_signs,varargin)
%         Add subunits to the model with specified properties. Can add multiple subunits in one
%         call. Default is to initialize regularization lambdas to be equal to an existing subunit
%         that acts on the same Xtarg (and has same NL type). Otherwise, they are initialized to 0.
%            INPUTS: 
%                 NLtypes: string or cell array of strings specifying upstream NL types
%                 mod_signs: vector of weights associated with each subunit (typically +/- 1)
%                 optional flags:
%                    ('xtargs',xtargs): specify vector of of Xtargets for each added subunit
%                    ('init_filts',init_filts): cell array of initial filter values for each subunit
%                    ('lambda_type',lambda_val): first input is a string specifying the type of
%                        regularization (e.g. 'd2t' for temporal smoothness). This must be followed by a
%                        scalar giving the associated lambda value. Applies these regularization
%                        lambdas to all added subunits
%            OUPUTS: nim: new nim object
            
            if ~iscell(NLtypes) && ischar(NLtypes);
                NLtypes = cellstr(NLtypes); %make NL types a cell array
            end
            nSubs = length(mod_signs); %number of subunits being added
            nStims = length(nim.stim_params);
            Xtargets = ones(nSubs,1); %default Xtargets to 1
            if length(NLtypes) == 1 && nSubs > 1
                NLtypes = repmat(NLtypes,nSubs,1); %if NLtypes is specified as a single string, assume we want this NL for all subunits
            end
            init_filts = cell(nSubs,1);
            
            %parse input flags
            j = 1; reg_types = {}; reg_vals = {};
            while j <= length(varargin)
                switch lower(varargin{j})
                    case 'xtargs'
                        Xtargets = varargin{j+1};
                        assert(all(ismember(Xtargets,1:nStims)),'invalid Xtargets specified');
                        j = j + 2;
                    case 'init_filts'
                        if ~iscell(varargin{j+1}) %if init_filts are specified as a matrix, make them a cell array
                            init_filts = cell(length(mod_signs),1);
                            for ii = 1:length(mod_signs)
                                init_filts{ii} = varargin{j+1}(:,ii);
                            end
                        else
                            init_filts = varargin{j+1};
                        end
                        j = j + 2;
                    case nim.allowed_reg_types
                        reg_types = cat(1,reg_types,lower(flag_name));
                        reg_vals = cat(1,reg_vals, varargin{j+1});
                        j = j + 2;
                    otherwise
                        error('Invalid input flag');
                end
            end
            
            assert(length(Xtargets) == nSubs,'length of mod_signs and Xtargets must be equal');
            %initialize subunits
            for ii = 1:nSubs %loop initializing subunits (start from last to initialize object array)
                stimD = prod(nim.stim_params(Xtargets(ii)).dims); %dimensionality of the current filter
                if isempty(init_filts{ii})
                    init_filt = randn(stimD,1)/stimD; %initialize fitler coefs with gaussian noise
                else
                    init_filt = init_filts{ii};
                end
                
               %use the regularization parameters from the most similar subunit if we have one,
               %otherwise use default init
               same_Xtarg = find(nim.subunits(:).Xtarg == Xtargets(ii),1); %find any existing subunits with this same Xtarget
               same_Xtarg_and_NL = same_Xtarg(strcmp(nim.get_NLtypes(same_Xtarg),NLtypes{ii})); %set that also have same NL type
               if ~isempty(same_Xtarg_and_NL) 
                    default_lambdas = nim.subunits(same_Xtarg_and_NL(1)).reg_lambdas;
               elseif ~isempty(same_Xtarg)
                    default_lambdas = nim.subunits(same_Xtarg(1)).reg_lambdas;
               else
                   default_lambdas = [];
               end
               
               nim.subunits = cat(1,nim.subunits,SUBUNIT(init_filt, mod_signs(ii), NLtypes{ii},Xtargets(ii))); %add new subunit
               if ~isempty(default_lambdas)
                   nim.subunits(end).reg_lambdas = default_lambdas;
               end
               for jj = 1:length(reg_vals) %add in user-specified regularization parameters
                   assert(reg_vals(jj) >= 0,'regularization hyperparameters must be non-negative');
                   nim.subunits(end).reg_lambdas = setfield(nim.subunits(end).reg_lambdas,reg_types{jj},reg_vals(jj));
               end
               
            end
        end
        
        %%
        function nim = init_spkhist(nim,n_bins,varargin)
%         nim = nim.init_spkhist(n_bins,varargin)
%         Adds a spike history term with specified parameters to an existing NIM.
%            INPUTS: 
%                 n_bins: number of coefficients in spike history filter
%                 optional flags:
%                    ('init_spacing',init_spacing): Initial spacing (in time bins) of piecewise constant 'basis functions'
%                    ('doubling_time',doubling_time): Make bin spacing logarithmically increasing, with given doubling time
%                    'negCon': If flag is included, constrain spike history filter coefs to be non-positive
%            OUPUTS: nim: new nim object
            
            % default inputs
            init_spacing = 1;
            doubling_time = n_bins;
            negCon = false;
            %parse input flags
            j = 1; %initialize counter after required input args
            while j <= length(varargin)
                switch lower( varargin{j})
                    case 'init_spacing'
                        init_spacing = varargin{j+1};
                        assert(init_spacing > 0,'invalid init_spacing');
                        j = j + 2;
                    case 'doubling_time'
                        doubling_time = varargin{j+1};
                        assert(doubling_time > 0,'invalid doubling time');
                        j = j + 2;
                    case 'negcon'
                        negCon = true;
                        j = j + 1;
                    otherwise
                        error('Invalid input flag');
                end
            end
            
            % COMPUTE RECTANGULAR BASIS FUNCTIONS 
            bin_edges = zeros(n_bins+1,1);
            inc = init_spacing; %bin increment size
            pos = 1; count = 0;
            for n = 1:n_bins+1
                bin_edges(n) = pos; %current bin edge loc
                pos = pos + inc; %move pointer by inc
                count = count + 1; %increase counter for doubling
                if count >= doubling_time %if reach doubling time, reset counter and double inc
                    count = 0; inc = inc * 2;
                end
            end
            
            % LOAD INTO A SPK_HIST STRUCT IN NIM
            nim.spk_hist.bin_edges = bin_edges;
            nim.spk_hist.coefs = zeros(n_bins,1); %init filter coefs to 0
            nim.spk_hist.negCon = negCon;
            nim.spk_hist.spkhstlen = n_bins;
        end
        
        %% 
        function nim = init_nonpar_NLs(nim, Xstims, varargin)
%         nim = nim.init_nonpar_NLs(Xstims,varargin)
%         Initializes the specified model subunits to have nonparametric (tent-basis) upstream NLs
%            INPUTS: 
%                 Xstims: cell array of stimuli
%                 optional flags:
%                    ('sub_inds',sub_inds): Index values of set of subunits to make nonpar (default is all)
%                    ('lambda_nld2',lambda_nld2): specify strength of smoothness regularization for the tent-basis coefs
%                    ('NLmon',NLmon): Set to +1 to constrain NL coefs to be monotonic increasing and
%                       -1 to make monotonic decreasing. 0 means no constraint. Default here is +1 (monotonic increasing)
%                    ('edge_p',edge_p): Scalar that determines the locations of the outermost tent-bases 
%                       relative to the underlying generating distribution
%                    ('n_bfs',n_bfs): Number of tent-basis functions to use 
%                    ('space_type',space_type): Use either 'equispace' for uniform bin spacing, or 'equipop' for 'equipopulated bins' 
%            OUPUTS: nim: new nim object
            
            Nsubs = length(nim.subunits);
            sub_inds = 1:Nsubs; %defualt to fitting all subunits 
            NLmon = 1; %default monotonic increasing TB-coefficients
            edge_p = 0.05; %relative to the generating distribution (pth percentile) where to put the outermost tent bases
            n_bfs = 25; %default number of tent basis functions
            space_type = 'equispace'; %default uninimrm tent basis spacing
            lambda_nld2 = 0; %default no smoothing on TB coefs
            
            j = 1;
            while j <= length(varargin)
                switch lower(varargin{j})
                    case 'sub_inds'
                        sub_inds = varargin{j+1};
                        assert(all(ismember(sub_inds,[1:Nsubs])),'invalid target subunits specified');
                        j = j + 2;
                    case 'nlmon'
                        NLmon = varargin{j+1};
                        assert(ismember(NLmon,[-1 0 1]),'NLmon must be a -1, 0 or 1');
                        j = j + 2;
                    case 'edge_p'
                        edge_p = varargin{j+1};
                        assert(edge_p > 0 & edge_p < 100,'edge_p must be between 0 and 100');
                        j = j + 2;
                    case 'n_bfs'
                        n_bfs = varargin{j+1};
                        assert(n_bfs > 0,'n_bfs must be greater than 0');
                        j = j + 2;
                    case 'space_type'
                        space_type = varargin{j+1};
                        assert(ismember(space_type,{'equispace','equipop'}),'unsupported bin spacing type');
                        j = j + 2;
                    case 'lambda_nld2'
                        lambda_nld2 = varargin{j+1};
                        assert(lambda_nld2 >= 0,'lambda must be >= 0');
                        j = j + 2;
                    otherwise
                        error('Invalid input flag');
                end
            end
            
            %store NL tent-basis parameters
            tb_params = struct('edge_p',edge_p,'n_bfs',n_bfs,'space_type',space_type);
            for ii = sub_inds %load the TB param struct into each subunit we're making nonpar
                nim.subunits(ii).TBparams = tb_params;
            end
            
            % Compute internal generating functions
            gint = nan(size(Xstims{1},1),Nsubs);
            for ii = 1:length(sub_inds)
                gint(:,sub_inds(ii)) = Xstims{nim.subunits(sub_inds(ii)).Xtarg} * nim.subunits(sub_inds(ii)).filtK;
            end
            
            prev_NL_types = nim.get_NLtypes(); %current NL types
            for imod = sub_inds
                nim.subunits(imod).NLtype = 'nonpar'; %set the subunit NL type to nonpar
                if strcmp(space_type,'equispace') %for equi-spaced bins
                    left_edge = my_prctile(gint(:,imod),edge_p);
                    right_edge = my_prctile(gint(:,imod),100-edge_p);
                    if left_edge == right_edge %if the data is constant over this range (e.g. with a 0 filter), just make the xrange unity
                        left_edge = right_edge - 0.5;
                        right_edge = right_edge + 0.5;
                    end
                    spacing = (right_edge - left_edge)/n_bfs;
                    %adjust the edge locations so one of the bins lands at 0
                    left_edge = ceil(left_edge/spacing)*spacing;
                    right_edge = floor(right_edge/spacing)*spacing;
                    TBx = linspace(left_edge,right_edge,n_bfs); %equispacing
                elseif strcmp(space_type,'equipop') %for equi-populated binning
                    if std(gint(:,imod)) == 0  % subunit with constant output
                        TBx = mean(gint(:,imod)) + linspace(-0.5,0.5,n_bfs); %do something sensible
                    else
                        TBx = my_prctile(gint(:,imod),linspace(edge_p,100-edge_p,n_bfs)); %equipopulated
                    end
                end
                %set nearest tent basis to 0 so we can keep it fixed during fitting
                [~,nearest] = min(abs(TBx));
                TBx(nearest) = 0;
                
                %initalize tent basis coefs
                switch prev_NL_types{imod}
                    case 'lin'
                        TBy = TBx;
                    case 'rectlin'
                        TBy = TBx; TBy(TBx < 0) = 0;
                    case 'rectquad'
                        TBy = TBx.^2; TBy(TBx < 0) = 0;
                    case 'quad'
                        TBy = TBx.^2;
                    case 'nonpar'
                        fprintf('upstream NL already set as nonparametric\n');
                    otherwise
                        error('Unsupported NL type');
                end
                nim.subunits(imod).TBy = TBy;
                nim.subunits(imod).TBx = TBx;
                nim.subunits(imod).reg_lambdas.nld2 = lambda_nld2;
                nim.subunits(imod).TBparams.NLmon = NLmon;
                nim.subunits(imod).TBy_deriv = nim.subunits(imod).get_TB_derivative(); %calculate derivative of Tent-basis coeffics
            end
        end
        
        %% MODEL EVAL
        
        function [LL, pred_rate, mod_internals, LL_data] = eval_model(nim, Robs, Xstims, varargin)
%           [LL, pred_rate, mod_internals, LL_data] = nim.eval_model(Robs, Xstims, <eval_inds>, varargin)
%            Evaluates the model on the supplied data
%                INPUTS:
%                   Robs: vector of observed data
%                   Xstims: cell array of stimuli
%                   <eval_inds>: optional vector of indices on which to evaluate the model
%                   optional flags:
%                       ('gain_funs',gain_funs): [TxK] matrix specifying gain at each timepoint for each subunit
%                OUTPUTS:
%                   LL: log-likelihood per spike
%                   pred_rate: predicted firing rates (in counts/bin)
%                   mod_internals: struct containing the internal components of the model prediction
%                      G: is the total generating signal (not including the constant offset theta). 
%                           This is the sum of subunit outputs (weighted by their subunit weights w)
%                      fgint: is the output of each subunit
%                      gint: is the output of each subunits linear filter
%                   LL_data: struct containing more detailed info about model performance:
%                      filt_pen: total regularization penalty on filter coefs 
%                      NL_pen total regularization penalty on filter upstream NLs
%                      nullLL: LL of constant-rate model
            
            Nsubs = length(nim.subunits); %number of subunits
            NT = length(Robs); %number of time points
            eval_inds = nan; %this default means evaluate on all data
            % PROCESS INPUTS
            gain_funs = []; %default has no gain_funs
            j = 1;
            while j <= length(varargin)
                flag_name = varargin{j};
                if ~ischar(flag_name)
                    eval_inds = flag_name;
                    j = j + 1;
                else
                    switch lower(flag_name)
                        case ~ischar(flag_name) %if the input is not a string, assume it's eval_inds
                            eval_inds = flag_name;
                            j = j + 1; %account for the fact that theres no associated input value
                        case 'gain_funs'
                            gain_funs = varargin{j+1};
                            j = j + 2;
                        otherwise
                            error('Invalid input flag');
                    end
                end
            end
            if size(Robs,2) > size(Robs,1); Robs = Robs'; end; %make Robs a column vector
            nim.check_inputs(Robs,Xstims,eval_inds,gain_funs); %make sure input format is correct
            if nim.spk_hist.spkhstlen > 0 % add in spike history term if needed
                Xspkhst = create_spkhist_Xmat( Robs, nim.spk_hist.bin_edges);
            end
            
            if ~isnan(eval_inds) %if specifying a subset of indices to train model params
                for nn = 1:length(Xstims)
                    Xstims{nn} = Xstims{nn}(eval_inds,:); %grab the subset of indices for each stimulus element
                end
                Robs = Robs(Uindx);
                if ~isempty(Xspkhst); Xspkhst = Xspkhst(eval_inds,:); end;
                if ~isempty(gain_funs); gain_funs = gain_funs(eval_inds,:); end;
            end
            [G, fgint, gint] = nim.process_stimulus(Xstims,1:Nsubs,gain_funs);
            if nim.spk_hist.spkhstlen > 0 % add in spike history term if needed
                G = G + Xspkhst*nim.spk_hist.coefs(:);
            end
            
            pred_rate = nim.apply_spkNL(G + nim.spkNL.theta); %apply spiking NL
            LL = nim.internal_LL(pred_rate,Robs); %compute LL
            Nspks = sum(Robs);
            LL = LL/Nspks; %normalize by spikes
            if nargout > 2 %if outputting model internals
                mod_internals.G = G;
                mod_internals.fgint = fgint;
                mod_internals.gint = gint;
            end
            if nargout > 3 %if we want more detailed model evaluation info, create an LL_data struct
                LL_data.LL = LL;
                [filt_penalties,NL_penalties] = nim.get_reg_pen(); %get regularization penalty for each subunit
                LL_data.filt_pen = sum(filt_penalties)/Nspks; %normalize by number of spikes
                LL_data.NL_pen = sum(NL_penalties)/Nspks;
                avg_rate = mean(Robs);
                null_prate = ones(NT,1)*avg_rate;
                nullLL = nim.internal_LL(null_prate,Robs)/Nspks;
                LL_data.nullLL = nullLL;
            end
        end
        %% display methods
        function [] = display_model(nim,Robs,Xstims,varargin)
%         [] = nim.display_model(<Robs>,<Xstims>,varargin)
%         Creates a display of the elements of a given NIM
%              INPUTS:
%                   <Robs>: observed spiking data. Needed if you want to utilize a spike-history
%                       filter. Otherwise set as empty
%                   <Xstims>: Stimulus cell array. Needed if you want to display the distributions of generating signals
%                   optional_flags:
%                         ('xtargs',xtargs): indices of stimuli for which we want to plot the filters
%                         'no_spknl': include this flag to suppress plotting of the spkNL 
%                         'no_spk_hist': include this flag to suppress plotting of spike history filter 
%                         ('gain_funs',gain_funs): if you want the computed generating signals to account for specified gain_funs
            
            if nargin < 2; Robs = []; end; 
            if nargin < 3; Xstims = []; end;
            
            Xtargs = [1:length(nim.stim_params)]; %default plot filters for all stimuli
            plot_spkNL = true;
            plot_spk_hist = true;
            gain_funs = [];
            j = 1; %initialize counter after required input args
            while j <= length(varargin)
                switch lower(varargin{j})
                    case 'xtargs'
                        Xtargs = varargin{j+1};
                        assert(all(ismember(Xtargs,1:length(nim.stim_params))),'invalid Xtargets specified');
                        j = j + 2;
                    case 'no_spknl'
                        plot_spkNL = false;
                        j = j + 1;
                    case 'gain_funs'
                        gain_funs = varargin{j+1};
                        j = j + 2;
                    case 'no_spk_hist'
                        plot_spk_hist = false;
                        j = j + 1;
                    otherwise
                        error('Invalid input flag');
                end
            end
            
            Nsubs = length(nim.subunits);
            spkhstlen = nim.spk_hist.spkhstlen;
            if spkhstlen > 0 && (plot_spk_hist || plot_spkNL)
                Xspkhst = create_spkhist_Xmat(Robs,nim.spk_hist.bin_edges);
            end
            n_hist_bins = 500; %internal parameter determining histogram resolution
            if ~isempty(Xstims)
                [G, ~, gint] = nim.process_stimulus(Xstims,1:Nsubs,gain_funs);
                G = G + nim.spkNL.theta; %add in constant term
                if spkhstlen > 0 %add in spike history filter output
                    G = G + Xspkhst*nim.spk_hist.coefs(:);
                end
            else
                G = []; gint = [];
            end
            
            % PLOT SPIKING NL FUNCTION
            if ~isempty(G) && plot_spkNL
                fig_handles.spk_nl = figure();
                n_bins = 1000; %bin resolution for G distribution
                [Gdist_y,Gdist_x] = hist(G,n_hist_bins); %histogram the generating signal
                
                %this is a hack to deal with cases where the threshold linear terms
                %create a min value of G
                if Gdist_y(1) > 2*Gdist_y(2)
                    Gdist_y(1) = 1.5*Gdist_y(2);
                end
                cur_xrange = Gdist_x([1 end]);
                
                cur_y = nim.apply_spkNL(Gdist_x);
                cur_y = cur_y/nim.stim_params(1).dt; %convert to correct firing rate units
                
                [ax,h1,h2] = plotyy(Gdist_x,cur_y,Gdist_x,Gdist_y);
                set(h1,'linewidth',1)
                yr = [min(cur_y) max(cur_y)];
                xlim(ax(1),cur_xrange)
                xlim(ax(2),cur_xrange);
                ylim(ax(1),yr);
                
                xlabel('Generating function')
                ylabel(ax(1),'Predicted firing rate','fontsize',14);
                ylabel(ax(2),'Probability','fontsize',14)
                set(ax(2),'ytick',[]);
                title('Spiking NL','fontsize',14)
            end
            
            if nim.spk_hist.spkhstlen > 0 && plot_spk_hist
                fig_handles.spk_hist = figure();
                subplot(2,1,1)
                stairs(nim.spk_hist.bin_edges(1:end-1)*nim.stim_params(1).dt,nim.spk_hist.coefs);
                xlim(nim.spk_hist.bin_edges([1 end])*nim.stim_params(1).dt)
                xl = xlim();
                line(xl,[0 0],'color','k','linestyle','--');
                xlabel('Time lag');
                ylabel('Spike history filter')
                title('Spike history term','fontsize',14)
                
                subplot(2,1,2)
                stairs(nim.spk_hist.bin_edges(1:end-1)*nim.stim_params(1).dt,nim.spk_hist.coefs);
                xlim(nim.spk_hist.bin_edges([1 end-1])*nim.stim_params(1).dt)
                set(gca,'xscale','log')
                xl = xlim();
                line(xl,[0 0],'color','k','linestyle','--');
                xlabel('Time lag');
                ylabel('Spike history filter')
                title('Spk Hist Log-axis','fontsize',14)
            end
            
            % CREATE FIGURE SHOWING INDIVIDUAL SUBUNITS
            for tt = Xtargs(Xtargs > 0) %loop over stimuli
                cur_subs = find([nim.subunits(:).Xtarg] == tt); %set of subunits acting on this stim
                
                fig_handles.stim_filts = figure();
                if nim.stim_params(tt).dims(3) > 1 %if 2-spatial-dimensional stim
                    n_columns = nim.stim_params(tt).dims(1) + 1;
                    n_rows = length(cur_subs);
                else
                    n_columns = max(round(sqrt(length(cur_subs)/2)),1);
                    n_rows = ceil(length(cur_subs)/n_columns);
                end
                nLags = nim.stim_params(tt).dims(1); %time lags
                dt = nim.stim_params(tt).dt; %time res
                nPix = squeeze(nim.stim_params(tt).dims(2:end)); %spatial dimensions
                %create filter time lag axis
                tax = (0:(nLags-1))*dt;
                tax = tax * 1000; % put in units of ms
                
                for imod = 1:length(cur_subs)
                    cur_sub = nim.subunits(cur_subs(imod));
                    
                    if nim.stim_params(tt).dims(3) == 1 %if < 2 spatial dimensions
                        %PLOT FILTER
                        subplot(n_rows,2*n_columns,(imod-1)*2+1);
                        if nPix == 1 %if temporal-only stim
                            %                             if isfield(thismod, 'keat_basis')
                            %                                 kblen = size(thismod.keat_basis,2);
                            %                                 tax = (0:kblen-1)*dt*1000;
                            %                                 plot(tax,thismod.filtK(:)'*thismod.keat_basis,'.-');
                            %                             else
                            plot(tax,cur_sub.filtK,'.-');
                            %                             end
                            xr = tax([1 end]);
                            line(xr,[0 0],'color','k','linestyle','--');
                            xlim(xr);
                            xlabel('Time lag')
                            ylabel('Filter coef');
                        elseif nPix(2) == 1
                            imagesc(1:nPix(1),tax,reshape(cur_sub.filtK,nLags,nPix(1)));
                            cl = max(abs(cur_sub.filtK));
                            caxis([-cl cl]);
                            %colormap(jet);
                            colormap(gray);
                            set(gca,'ydir','normal');
                            xlabel('Pixels')
                            ylabel('Time lags');
                        end
                        if strcmp(cur_sub.NLtype,'lin')
                            title('Linear stimulus filter','fontsize',14)
                        elseif cur_sub.weight > 0
                            title('Excitatory stimulus filter','fontsize',14);
                        elseif cur_sub.weight < 0
                            title('Suppressive stimulus filter','fontsize',14);
                        end
                    else %if 2-spatial dimensional stim
                        maxval = max(abs(cur_sub.filtK));
                        for jj = 1:nim.stim_params(tt).dims(1) %loop over time slices
                            subplot(n_rows,n_columns,(imod-1)*n_columns + jj);
                            cur_fdims = jj - 1 + (1:nim.stim_params(tt).dims(1):prod(nim.stim_params(tt).dims));
                            imagesc(1:nPix(1),1:nPix(2),reshape(cur_sub.filtK(cur_fdims),nim.stim_params(tt).dims(2:end)));
                            colormap(gray)
                            if strcmp(cur_sub.NLtype,'lin')
                                title(sprintf('Lin-input Lag %d',jj-1),'fontsize',10);
                            elseif cur_sub.weight > 0
                                title(sprintf('E-Input Lag %d',jj-1),'fontsize',10);
                            elseif cur_sub.weight < 0
                                title(sprintf('S-Input Lag %d',jj-1),'fontsize',10);
                            end
                            caxis([-maxval maxval]*0.85);
                        end
                    end
                    
                    %PLOT UPSTREAM NL
                    if nim.stim_params(tt).dims(3) == 1
                        subplot(n_rows,2*n_columns,(imod-1)*2+2);
                    else
                        subplot(n_rows,n_columns,(imod)*n_columns);
                    end
                    if ~isempty(gint) %if computing distribution of filtered stim
                        [gendist_y,gendist_x] = hist(gint(:,cur_subs(imod)),n_hist_bins);
                        
                        % Sometimes the gendistribution has a lot of zeros (dont want to screw up plot)
                        [a b] = sort(gendist_y);
                        if a(end) > a(end-1)*1.5
                            gendist_y(b(end)) = gendist_y(b(end-1))*1.5;
                        end
                    else
                        gendist_x = linspace(-3,3,n_hist_bins); %otherwise, just pick an arbitrary x-axis to plot the NL
                    end
                    if strcmp(cur_sub.NLtype,'nonpar')
                        cur_modx = cur_sub.TBx; cur_mody = cur_sub.TBy;
                    else
                        cur_modx = gendist_x; cur_mody = cur_sub.apply_NL(cur_modx);
                    end
                    cur_xrange = cur_modx([1 end]);
                    
                    if ~isempty(gint)
                        [ax,h1,h2] = plotyy(cur_modx,cur_mody,gendist_x,gendist_y);
                        if strcmp(cur_sub.NLtype,'nonpar')
                            set(h1,'Marker','o');
                        end
                        set(h1,'linewidth',1)
                        xlim(ax(1),cur_xrange)
                        xlim(ax(2),cur_xrange);
                        ylim(ax(1),[min(cur_mody) max(cur_mody)]);
                        set(ax(2),'ytick',[])
                        yl = ylim();
                        line([0 0],yl,'color','k','linestyle','--');
                        ylabel(ax(1),'Subunit output','fontsize',12);
                        ylabel(ax(2),'Probability','fontsize',12)
                    else
                        h = plot(cur_modx,cur_mody,'linewidth',1);
                        if strcmp(cur_sub.NLtype,'nonpar')
                            set(h,'Marker','o');
                        end
                        xlim(cur_xrange)
                        ylim([min(cur_mody) max(cur_mody)]);
                        ylabel('Subunit output','fontsize',12);
                    end
                    box off
                    xlabel('Internal generating function')
                    title('Upstream NL','fontsize',14)
                end
            end
        end
        
        %% fitting methods
        
        function nim = fit_filters(nim, Robs, Xstims, varargin)
%         nim = nim.fit_filters(Robs, Xstims, <train_inds>, varargin)
%         estimate filters of NIM model.
%         INPUTS:
%            Robs: vector of response observations (e.g. spike counts)
%            Xstims: cell array of stimuli
%            <train_inds>: index values of data on which to fit the model [default to all indices in provided data]
%            optional flags:
%                ('sub_inds',sub_inds): set of subunits whos filters we want to optimize [default is all]
%                ('gain_funs',gain_funs): matrix of multiplicative factors, one column for each subunit
%                ('optim_params',optim_params): struct of desired optimization parameters
%                'silent': include this flag to suppress the iterative optimization display
%                'hold_spkhist': include this flat to hold the spk NL filter constant
%         OUTPUTS:
%            new nim object with optimized subunit filters
% 
            Nsubs = length(nim.subunits); %number of subunits
            
            % PROCESS INPUTS
            fit_subs = 1:Nsubs; %defualt to fitting all subunits (plus -1 for spkHist filter)
            gain_funs = []; %default has no gain_funs
            train_inds = nan; %default nan means train on all data
            optim_params = []; %default has no user-specified optimization parameters
            fit_spk_hist = nim.spk_hist.spkhstlen > 0; %default is to fit the spkNL filter if it exists
            silent = false; %default is to display optimization output
            j = 1;
            while j <= length(varargin)
                flag_name = varargin{j};
                if ~ischar(flag_name)%if not a flag, it must be train_inds
                    train_inds = flag_name;
                    j = j + 1; %only one argument here
                else
                    switch lower(flag_name)
                        case 'sub_inds'
                            fit_subs = varargin{j+1};
                            assert(all(ismember(fit_subs,1:Nsubs)),'invalid target subunits specified');
                            j = j + 2;
                        case 'gain_funs'
                            gain_funs = varargin{j+1};
                            j = j + 2;
                        case 'optim_params'
                            optim_params = varargin{j+1};
                            assert(isstruct(optim_params),'optim_params must be a struct');
                            j = j + 2;
                        case 'silent'
                            silent = true;
                            j = j + 1;
                        case 'hold_spkhist'
                            fit_spk_hist = false;
                            j = j + 1;
                        otherwise
                            error('Invalid input flag');
                    end
                end
            end
            
            if size(Robs,2) > size(Robs,1); Robs = Robs'; end; %make Robs a column vector
            nim.check_inputs(Robs,Xstims,train_inds,gain_funs); %make sure input format is correct
            
            Nfit_subs = length(fit_subs); %number of targeted subunits
            non_fit_subs = setdiff([1:Nsubs],fit_subs); %elements of the model held constant
            spkhstlen = nim.spk_hist.spkhstlen; %length of spike history filter
            if fit_spk_hist; assert(spkhstlen > 0,'no spike history term initialized!'); end;
            if spkhstlen > 0 % create spike history Xmat IF NEEDED
                Xspkhst = create_spkhist_Xmat( Robs, nim.spk_hist.bin_edges);
            else
                Xspkhst = [];
            end
            if ~isnan(train_inds) %if specifying a subset of indices to train model params
                for nn = 1:length(Xstims)
                    Xstims{nn} = Xstims{nn}(train_inds,:); %grab the subset of indices for each stimulus element
                end
                Robs = Robs(Uindx);
                if ~isempty(Xspkhst); Xspkhst = Xspkhst(train_inds,:); end;
                if ~isempty(gain_funs); gain_funs = gain_funs(train_inds,:); end;
            end
            
            % PARSE INITIAL PARAMETERS
            init_params = [];
            lambda_L1 = zeros(size(init_params));
            sign_con = zeros(size(init_params));
            for imod = fit_subs
                cur_kern = nim.subunits(imod).filtK;
                if (nim.subunits(imod).Ksign_con ~= 0) %add sign constraints on the filters of this subunit if needed
                    sign_con(length(init_params)+(1:length(cur_kern))) = nim.subunits(imod).Ksign_con;
                end
                lambda_L1(length(init_params) + (1:length(cur_kern))) = nim.subunits(imod).reg_lambdas.l1;
                init_params = [init_params; cur_kern]; % add coefs to initial param vector
            end
            lambda_L1 = lambda_L1/sum(Robs); % since we are dealing with LL/spk
            Nfit_filt_params = length(init_params); %number of filter coefficients in param vector
            % Add in spike history coefs
            if fit_spk_hist
                init_params = [init_params; nim.spk_hist.coefs];
            end
            
            init_params(end+1) = nim.spkNL.theta; % add constant offset
            [nontarg_g] = nim.process_stimulus(Xstims,non_fit_subs,gain_funs);
            if ~fit_spk_hist && spkhstlen > 0 %add in spike history filter output, if we're not fitting it
                nontarg_g = nontarg_g + Xspkhst*nim.spk_hist.coefs(:);
            end
            
            % IDENTIFY ANY CONSTRAINTS
            use_con = 0;
            LB = -Inf*ones(size(init_params));
            UB = Inf*ones(size(init_params));
            % Constrain any of the filters to be positive or negative
            if any(sign_con ~= 0)
                LB(sign_con == 1) = 0;
                UB(sign_con == -1) = 0;
                use_con = 1;
            end
            if fit_spk_hist %if optimizing spk history term
                %negative constraint on spk history coefs
                if nim.spk_hist.negCon
                    spkhist_inds = Nfit_filt_params + (1:spkhstlen);
                    UB(spkhist_inds) = 0;
                    use_con = 1;
                end
            end
            
            % GENERATE REGULARIZATION MATRICES
            Tmats = nim.make_Tikhonov_matrices();
            
            fit_opts = struct('fit_spk_hist', fit_spk_hist, 'fit_subs',fit_subs); %put any additional fitting options into this struct
            %the function we want to optimize
            opt_fun = @(K) nim.internal_LL_filters(K,Robs,Xstims,Xspkhst,nontarg_g,gain_funs,Tmats,fit_opts);
            
            %determine which optimizer were going to use
            if max(lambda_L1) > 0
                assert(~use_con,'Can use L1 penalty with constraints');
                assert(exist('L1General2_PSSas','file') == 2,'Need Mark Schmidts optimization tools installed to use L1');
                optimizer = 'L1General_PSSas';
            else
                if ~use_con %if there are no constraints
                    if exist('minFunc','file') == 2
                        optimizer = 'minFunc';
                    else
                        optimizer = 'fminunc';
                    end
                else
                    if exist('minConf_TMP','file')==2 
                        optimizer = 'minConf_TMP';
                    else
                        optimizer = 'fmincon';
                    end
                end
            end
            optim_params = nim.set_optim_params(optimizer,optim_params,silent);
            if ~silent; fprintf('Running optimization using %s\n\n',optimizer); end;
            
            switch optimizer %run optimization
                case 'L1General2_PSSas'
                    [params] = L1General2_PSSas(opt_fun,init_params,lambda_L1,optim_params);
                case 'minFunc'
                    [params] = minFunc(opt_fun, init_params, optim_params);
                case 'fminunc'
                    [params] = fminunc(opt_fun, init_params, optim_params);
                case 'minConf_TMP'
                    [params] = minConf_TMP(opt_fun, init_params, LB, UB, optim_params);
                case 'fmincon'
                    [params] = fmincon(opt_fun, init_params, [], [], [], [], LB, UB, [], optim_params);
            end
            [~,penGrad] = opt_fun(params);
            first_order_optim = max(abs(penGrad));
            if first_order_optim > nim.opt_check_FO
                warning(sprintf('First-order optimality: %.3f, fit might not be converged!',first_order_optim));
            end
            
            % PARSE MODEL FIT
            nim.spkNL.theta = params(end); %set new offset parameter
            if fit_spk_hist
                nim.spk_hist.coefs = params(Nfit_filt_params + (1:spkhstlen));
            end
            kOffset = 0; %position counter for indexing param vector
            for ii = 1:Nfit_subs
                filtLen = length(nim.subunits(fit_subs(ii)).filtK);
                cur_kern = params((1:filtLen) + kOffset); %grab parameters corresponding to this subunit's filters
                nim.subunits(fit_subs(ii)).filtK = cur_kern(:); %assign new filter values
                kOffset = kOffset + filtLen;
            end
            
            [LL,~,mod_internals,LL_data] = nim.eval_model(Robs,Xstims,'gain_funs',gain_funs);
            nim = nim.set_subunit_scales(mod_internals.fgint); %update filter scales
            cur_fit_details = struct('fit_type','filter','LL',LL,'filt_pen',LL_data.filt_pen,...
                'NL_pen',LL_data.NL_pen,'FO_optim',first_order_optim);
            nim.fit_props = cur_fit_details; %store details of this fit
            nim.fit_hist = cat(1,nim.fit_hist,cur_fit_details);
        end
        
        %%
        function nim = fit_upstreamNLs(nim, Robs, Xstims, varargin)
%         nim = nim.fit_upstreamNLs(Robs, Xstims, <train_inds>, varargin)
%         Optimizes the upstream NLs (in terms of tent-basis functions) 
%         INPUTS:
%            Robs: vector of response observations (e.g. spike counts)
%            Xstims: cell array of stimuli
%            <train_inds>: index values of data on which to fit the model [default to all indices in provided data]
%            optional flags:
%                ('sub_inds',sub_inds): set of subunits whos filters we want to optimize [default is all]
%                ('gain_funs',gain_funs): matrix of multiplicative factors, one column for each subunit
%                ('optim_params',optim_params): struct of desired optimization parameters
%                'silent': include this flag to suppress the iterative optimization display
%                'hold_spkhist': include this flat to hold the spk NL filter constant
%                'no_rescaling': use this flag if you dont want to rescale the NLs after fitting
%        OUTPUTS:
%            nim: output model struct
            
            Nsubs = length(nim.subunits); %number of subunits
            NT = length(Robs); %number of time points
            
            % PROCESS INPUTS
            poss_targets = find(strcmp(nim.get_NLtypes,'nonpar'))'; %set of subunits with nonpar NLs
            fit_subs = poss_targets; %defualt to fitting all subunits that we can
            gain_funs = []; %default has no gain_funs
            train_inds = nan; %default nan means train on all data
            optim_params = []; %default has no user-specified optimization parameters
            silent = false; %default is show the optimization output
            fit_spk_hist = nim.spk_hist.spkhstlen > 0; %default is fit the spkNL filter if it exists
            rescale_NLs = true; %default is to rescale the y-axis of NLs after estimation
            
            j = 1;
            while j <= length(varargin)
                flag_name = varargin{j}; %if not a flag, it must be train_inds
                if ~ischar(flag_name)
                    train_inds = flag_name;
                    j = j + 1; %there's just one arg here
                else
                    switch lower(flag_name)
                        case 'sub_inds'
                            fit_subs = varargin{j+1};
                            assert(all(ismember(fit_subs,[poss_targets])),'specified target doesnt have non-parametric NL, or doesnt exist');
                            j = j + 2;
                        case 'gain_funs'
                            gain_funs = varargin{j+1};
                            j = j + 2;
                        case 'optim_params'
                            optim_params = varargin{j+1};
                            assert(isstruct(optim_params),'optim_params must be a struct');
                            j = j + 2;
                        case 'silent'
                            silent = true;
                            j = j + 1;
                        case 'no_rescaling'
                            rescale_NLs = false;
                            j = j + 1;
                        case 'hold_spkhist'
                            fit_spk_hist = false;
                            j = j + 1;
                        otherwise
                            error('Invalid input flag');
                    end
                end
            end
            
            if size(Robs,2) > size(Robs,1); Robs = Robs'; end; %make Robs a column vector
            nim.check_inputs(Robs,Xstims,train_inds,gain_funs); %make sure input format is correct
            
            Nfit_subs = length(fit_subs); %number of targeted subunits
            non_fit_subs = setdiff(1:Nsubs,fit_subs); %elements of the model held constant
            spkhstlen = nim.spk_hist.spkhstlen; %length of spike history filter
            if fit_spk_hist; assert(spkhstlen > 0,'no spike history term initialized!'); end;
            if fit_spk_hist
                Xspkhst = create_spkhist_Xmat( Robs, nim.spk_hist.bin_edges);
            else
                Xspkhst = [];
            end
            if ~isnan(train_inds) %if specifying a subset of indices to train model params
                for nn = 1:length(Xstims)
                    Xstims{nn} = Xstims{nn}(train_inds,:); %grab the subset of indices for each stimulus element
                end
                Robs = Robs(Uindx);
                if ~isempty(Xspkhst); Xspkhst = Xspkhst(train_inds,:); end;
                if ~isempty(gain_funs); gain_funs = gain_funs(train_inds,:); end;
            end
            
            n_TBs = arrayfun(@(x) length(x.TBx),nim.subunits(fit_subs));  %get the number of TBs for each subunit
            assert(length(unique(n_TBs)) == 1,'Have to have same number of tent-bases for each subunit');
            n_TBs = unique(n_TBs);
            
            nontarg_g = nim.process_stimulus(Xstims,non_fit_subs,gain_funs); %get output of nontarget subunits
            if ~fit_spk_hist && spkhstlen > 0 %add in spike history filter output, if we're not fitting it
                nontarg_g = nontarg_g + Xspkhst*nim.spk_hist.coefs(:);
            end
            
            % COMPUTE NEW X-MATRIX OUT OF TENT-BASIS OUTPUTS
            XNL = zeros(NT,Nfit_subs*n_TBs); %initialize X matrix which is for the NL BFs of each module
            for ii = 1:Nfit_subs %for each module
                tar = fit_subs(ii);
                gint = Xstims{nim.subunits(tar).Xtarg}*nim.subunits(tar).filtK;
                % The output of the current model's internal filter projected onto the tent basis representation
                if isempty(gain_funs)
                    tbf_out = nim.subunits(tar).weight * nim.subunits(tar).tb_rep(gint);
                else
                    tbf_out = nim.subunits(tar).weight * bsxfun(@times,nim.subunits(tar).tb_rep(gint),gain_funs(:,tar));
                end
                XNL(:,((ii-1)*n_TBs + 1):(ii*n_TBs)) = tbf_out; % assemble filtered NLBF outputs into X matrix
            end
            
            % CREATE INITIAL PARAMETER VECTOR
            % Compute initial fit parameters
            init_params = [];
            for imod = fit_subs
                init_params = [init_params; nim.subunits(imod).TBy']; %not incorporating the multiplier here because doing so messes with regularization
            end
            
            % Add in spike history coefs
            if fit_spk_hist
                init_params = [init_params; nim.spk_hist.coefs];
            end
            
            init_params(end+1) = nim.spkNL.theta; %add constant offset
            
            lambda_nl = nim.get_reg_lambdas('sub_inds',fit_subs,'nld2');
            if any(lambda_nl > 0)
                Tmat = nim.make_NL_Tmat;
            else
                Tmat = [];
            end
            
            % PROCESS CONSTRAINTS
            use_con = 0;
            LB = []; UB = []; A = []; Aeq = []; % initialize constraint parameters
            % Check for spike history coef constraints
            if fit_spk_hist
                % negative constraint on spk history coefs
                if nim.spk_hist.negCon
                    spkhist_inds = Nfit_subs*n_tbfs + (1:spkhstlen);
                    LB = -Inf*ones(size(initial_params));
                    UB = Inf*ones(size(initial_params));
                    UB(spkhist_inds) = 0;
                    use_con = 1;
                end
            end
            
            % Process NL monotonicity constraints, and constraints that the tent basis
            % centered at 0 should have coefficient of 0 (eliminate y-shift degeneracy)
            if any(arrayfun(@(x) x.TBparams.NLmon,nim.subunits(fit_subs)) ~= 0)
                zvec = zeros(1,length(init_params)); % indices of tent-bases centered at 0
                for ii = 1:Nfit_subs
                    cur_range = (ii-1)*n_TBs + (1:n_TBs);
                    % For monotonicity constraint
                    if nim.subunits(fit_subs(ii)).TBparams.NLmon ~= 0
                        for jj = 1:length(cur_range)-1 %create constraint matrix
                            cur_vec = zvec;
                            cur_vec(cur_range([jj jj + 1])) = nim.subunits(fit_subs(ii)).TBparams.NLmon*[1 -1];
                            A = cat(1,A,cur_vec);
                        end
                    end
                    
                    % Constrain the 0-coefficient to be 0
                    [~,zp] = find(nim.subunits(fit_subs(ii)).TBx == 0);
                    if isempty(zp)
                        error('Need one TB to be centered at 0')
                    end
                    cur_vec = zvec;
                    cur_vec(cur_range(zp)) = 1;
                    Aeq = cat(1,Aeq,cur_vec);
                end
                b = zeros(size(A,1),1);
                beq = zeros(size(Aeq,1),1);
                use_con = 1;
            end
            
            if ~use_con %if there are no constraints
                if exist('minFunc','file') == 2
                    optimizer = 'minFunc';
                else
                    optimizer = 'fminunc';
                end
            else
                optimizer = 'fmincon';
            end
            optim_params = nim.set_optim_params(optimizer,optim_params,silent);
            if ~silent; fprintf('Running optimization using %s\n\n',optimizer); end;
            
            fit_opts = struct('fit_spk_hist', fit_spk_hist, 'fit_subs',fit_subs); %put any additional fitting options into this struct
            opt_fun = @(K) nim.internal_LL_NLs(K, Robs, XNL, Xspkhst,nontarg_g, Tmat,fit_opts);
            
            switch optimizer %run optimization
                case 'L1General2_PSSas'
                    [params] = L1General2_PSSas(opt_fun,init_params,lambda_L1,optim_params);
                case 'minFunc'
                    [params] = minFunc(opt_fun, init_params, optim_params);
                case 'fminunc'
                    [params] = fminunc(opt_fun, init_params, optim_params);
                case 'minConf_TMP'
                    [params] = minConf_TMP(opt_fun, init_params, LB, UB, optim_params);
                case 'fmincon'
                    [params] = fmincon(opt_fun, init_params, A, b, Aeq, beq, LB, UB, [], optim_params);
            end
            [~,penGrad] = opt_fun(params);
            first_order_optim = max(abs(penGrad));
            if first_order_optim > nim.opt_check_FO
                warning(sprintf('First-order optimality %.3f, fit might not be converged!',first_order_optim));
            end
            
            nlmat = reshape(params(1:Nfit_subs*n_TBs),n_TBs,Nfit_subs); %take output K vector and restructure into a matrix of NLBF coefs, one for each module
            nlmat_resc = nlmat;
            for ii = 1:Nfit_subs;
                cur_pset = ((ii-1)*n_TBs+1) : (ii*n_TBs);
                thisnl = nlmat(:,ii); %NL coefs for current subunit
                cur_std = std(XNL(:,cur_pset)*thisnl);
                if rescale_NLs %rescale so that the std dev of the subunit output is conserved
                    thisnl = thisnl*nim.subunits(fit_subs(ii)).scale/cur_std;
                else
                    nim.subunits(fit_subs(ii)).scale = cur_std; %otherwise adjust the model output std dev
                end
                nim.subunits(fit_subs(ii)).TBy = thisnl';
                nlmat_resc(:,ii) = thisnl';
            end
            
            if fit_spk_hist
                nim.spk_hist.coefs = params((Nfit_subs*n_TBs+1):(Nfit_subs*n_TBs+spkhstlen));
            end
            
            % If rescaling the Nls, we need to resestimate the offset theta after scaling
            if rescale_NLs
                resc_nlvec = nlmat_resc(:);
                new_g_out = XNL*resc_nlvec;
                G = nontarg_g + new_g_out;
                if spkhstlen > 0
                    G = G + Xspkhst*nim.spk_hist.coefs;
                end
                init_theta = params(end);
                opts.Display = 'off';opts.GradObj = 'on'; opts.LargeScale = 'off';
                new_theta = fminunc( @(K) nim.internal_theta_opt(K,G,Robs), init_theta, opts);
                nim.spkNL.theta = new_theta;
            else
                nim.spkNL.theta = params(end);
            end
            
            [LL,~,mod_internals,LL_data] = nim.eval_model(Robs,Xstims,'gain_funs',gain_funs);
            nim = nim.set_subunit_scales(mod_internals.fgint); %update filter scales
            cur_fit_details = struct('fit_type','upstream_NLs','LL',LL,'filt_pen',LL_data.filt_pen,...
                'NL_pen',LL_data.NL_pen,'FO_optim',first_order_optim);
            nim.fit_props = cur_fit_details;
            nim.fit_hist = cat(1,nim.fit_hist,cur_fit_details);
        end
        
        %%
        function nim = fit_spkNL(nim, Robs, Xstims, varargin)
%         nim = nim.fit_spkNL(Robs, Xstims, <train_inds>, varargin)
%         Optimizes the parameters of the spkNL
%         INPUTS:
%            Robs: vector of response observations (e.g. spike counts)
%            Xstims: cell array of stimuli
%            <train_inds>: index values of data on which to fit the model [default to all indices in provided data]
%            optional flags:
%                ('gain_funs',gain_funs): matrix of multiplicative factors, one column for each subunit
%                ('optim_params',optim_params): struct of desired optimization parameters
%                'silent': include this flag to suppress the iterative optimization display
%                'hold_const': vector of parameter indices to hold constant
%        OUTPUTS:
%            nim: output model struct            
            
            % PROCESS INPUTS
            gain_funs = []; %default has no gain_funs
            train_inds = nan; %default nan means train on all data
            optim_params = []; %default has no user-specified optimization parameters
            silent = false; %default is show the optimization output
            hold_const = []; %default is fit all spk NL params
            j = 1;
            while j <= length(varargin)
                flag_name = varargin{j}; %if not a flag, it must be train_inds
                if ~ischar(flag_name)
                    train_inds = flag_name;
                    j = j + 1;
                else
                    switch lower(flag_name)
                        case 'gain_funs'
                            gain_funs = varargin{j+1};
                            j = j + 2;
                        case 'optim_params'
                            optim_params = varargin{j+1};
                            assert(isstruct(optim_params),'optim_params must be a struct');
                            j = j + 2;
                        case 'silent'
                            silent = true;
                        case 'hold_const'
                            hold_const = varargin{j+1};
                            j = j + 2;
                        otherwise
                            error('Invalid input flag');
                    end
                end
            end
            
            if size(Robs,2) > size(Robs,1); Robs = Robs'; end; %make Robs a column vector
            nim.check_inputs(Robs,Xstims,train_inds,gain_funs); %make sure input format is correct
            
            [~, ~, mod_internals] = nim.eval_model(Robs, Xstims, train_inds,'gain_funs',gain_funs);
            G = mod_internals.G;
            if ~isnan(train_inds) %if specifying a subset of indices to train model params
                Robs = Robs(train_inds);
            end
            
            init_params = [nim.spkNL.params nim.spkNL.theta]; %initialize parameters to fit (including the offset term theta)
            
            %BOUND CONSTRAINTS
            LB = -Inf*ones(size(init_params));
            UB = Inf*ones(size(init_params));
            if ismember(nim.spkNL.type,{'lin','rectlin','exp','softplus','logistic'})
                LB(1) = 0; %beta is non-negative
            end
            if ismember(nim.spkNL.type,{'softplus'})
                LB(2) = 0; %alpha is non-negative
            end
            %equality constraints
            Aeq = []; Beq = [];
            for i = 1:length(hold_const)
                Aeq = [Aeq; zeros(1,length(init_params))];
                Aeq(end,hold_const(i)) = 1;
                Beq = [Beq; init_params(hold_const(i))];
            end
            
            optimizer = 'fmincon';
            optim_params = nim.set_optim_params(optimizer,optim_params,silent);
            optim_params.GradObj = 'on';
            opt_fun = @(K) nim.internal_LL_spkNL(K, Robs, G);
            params = fmincon(opt_fun, init_params, [], [], Aeq, Beq, LB, UB, [], optim_params);
            [~,penGrad] = opt_fun(params);
            first_order_optim = max(abs(penGrad));
            if first_order_optim > nim.opt_check_FO
                warning(sprintf('First-order optimality: %.3f, fit might not be converged!',first_order_optim));
            end
            
            nim.spkNL.params = params(1:end-1);
            nim.spkNL.theta = params(end);
            
            [LL,~,mod_internals,LL_data] = nim.eval_model(Robs,Xstims,'gain_funs',gain_funs);
            nim = nim.set_subunit_scales(mod_internals.fgint); %update filter scales
            cur_fit_details = struct('fit_type','spkNL','LL',LL,'filt_pen',LL_data.filt_pen,...
                'NL_pen',LL_data.NL_pen,'FO_optim',first_order_optim);
            nim.fit_props = cur_fit_details;
            nim.fit_hist = cat(1,nim.fit_hist,cur_fit_details);
        end
    end
    
    methods (Hidden)
        %% internal methods
        
        function [penLL, penLLgrad] = internal_LL_filters(nim,params,Robs,Xstims,Xspkhst,nontarg_g,gain_funs,Tmats,fit_opts)
            %computes the penalized LL and its gradient wrt the filters for the given nim
            %with parameter vector params
            
            fit_subs = fit_opts.fit_subs;
            Nfit_subs = length(fit_subs); %number of targeted subs
            
            % USEFUL VALUES
            theta = params(end); % offset
            gint = nan(length(Robs),Nfit_subs); %initialize matrix for storing filter outputs
            filtLen = zeros(Nfit_subs,1); %store the length of each (target) sub's filter
            filtKs = cell(Nfit_subs,1); %store the filter coefs for all (target) subs)
            param_inds = cell(Nfit_subs,1); %this will store the index values of each subunit's filter coefs within the parameter vector
            Xtarg_set = [nim.subunits(fit_subs).Xtarg]; %vector of Xfit_subs for set of subunits being optimized
            un_Xtargs = unique(Xtarg_set); %set of unique Xfit_subs
            mod_NL_types = {nim.subunits(fit_subs).NLtype}; %NL types for each targeted subunit
            unique_NL_types = unique(mod_NL_types); %unique set of NL types being used
            mod_weights = [nim.subunits(fit_subs).weight]'; %signs of targeted subunits
            
            G = theta + nontarg_g; % initialize overall generating function G with the offset term and the contribution from nontarget subs
            
            NKtot = 0;  %init filter coef counter
            for ii = 1:Nfit_subs %loop over subunits, get filter coefs and their indices within the parameter vector
                filtLen(ii) = length(nim.subunits(fit_subs(ii)).filtK); % length of filter
                param_inds{ii} = NKtot + (1:filtLen(ii)); %set of param indices associated with this subunit's filters
                filtKs{ii} = params(param_inds{ii}); %store filter coefs
                NKtot = NKtot + filtLen(ii); %inc counter
            end
            for ii = 1:length(un_Xtargs) %loop over the unique Xtargs and compute the generating signals for all relevant filters
                cur_subs = find(Xtarg_set == un_Xtargs(ii)); %set of targeted subunits that act on this Xtarg
                gint(:,cur_subs) = Xstims{un_Xtargs(ii)} * cat(2,filtKs{cur_subs}); %apply filters to stimulus
            end
            
            fgint = gint; %init subunit outputs by filter outputs
            for ii = 1:length(unique_NL_types) %loop over unique subunit NL types and apply NLs to gint in batch
                cur_subs = find(strcmp(mod_NL_types,unique_NL_types{ii})); %set of subs with this NL type
                if strcmp(unique_NL_types{ii},'nonpar')
                    for jj = 1:length(cur_subs) %for TB NLs need to apply each subunit's NL individually
                        fgint(:,cur_subs(jj)) = nim.subunits(fit_subs(cur_subs(jj))).apply_NL(gint(:,cur_subs(jj)));
                    end
                elseif ~strcmp(unique_NL_types{ii},'lin') %if it's not just a linear sub, we have to do something
                    fgint(:,cur_subs) = nim.subunits(fit_subs(cur_subs(1))).apply_NL(gint(:,cur_subs)); %apply upstream NL to all subunits of this type
                end
            end
            
            % Multiply by weight (and multiplier, if appl) and add to generating function
            if isempty(gain_funs)
                G = G + fgint*mod_weights;
            else
                G = G + (fgint.*gain_funs(:,fit_subs))*mod_weights;
            end
            
            % Add contribution from spike history filter
            if fit_opts.fit_spk_hist
                G = G + Xspkhst*params(NKtot + (1:nim.spk_hist.spkhstlen));
            end
            
            pred_rate = nim.apply_spkNL(G);
            penLL = nim.internal_LL(pred_rate,Robs); %compute LL
            
            %residual = LL'[r].*F'[g]
            residual = nim.internal_LL_deriv(pred_rate,Robs) .* nim.apply_spkNL_deriv(G,pred_rate <= nim.min_pred_rate);
            
            penLLgrad = zeros(length(params),1); %initialize LL gradient
            penLLgrad(end) = sum(residual);      %Calculate derivatives with respect to constant term (theta)
            
            % Calculate derivative with respect to spk history filter
            if fit_opts.fit_spk_hist
                penLLgrad(NKtot+(1:nim.spk_hist.spkhstlen)) = residual'*Xspkhst;
            end
            
            for ii = 1:un_Xtargs %loop over unique Xfit_subs and compute LL grad wrt stim filters
                cur_sub_inds = find(Xtarg_set == un_Xtargs(ii)); %set of subunits with this Xtarget
                cur_NL_types = mod_NL_types(cur_sub_inds); %NL types of current subs
                cur_unique_NL_types = unique(cur_NL_types); %set of unique NL types
                
                if length(cur_sub_inds) == 1 && strcmp(cur_unique_NL_types,'lin') %if there's only a single linear subunit, this is a faster calc
                    if isempty(gain_funs)
                        penLLgrad(param_inds{cur_sub_inds}) = residual'*Xstims{un_Xtargs(ii)} * nim.subunits(cur_sub_inds).weight;
                    else
                        penLLgrad(param_inds{cur_sub_inds}) = (gain_funs.*residual)'*Xstims{un_Xtargs(ii)} * nim.subunits(cur_sub_inds).weight;
                    end
                else %otherwise, compute a matrix of upstream NL derivatives fpg
                    fpg = ones(length(residual),length(cur_sub_inds)); %initialize to linear NL derivative (all ones)
                    for jj = 1:length(cur_unique_NL_types) %loop over unique NL types
                        cur_sub_subinds = find(strcmp(cur_NL_types,cur_unique_NL_types{jj})); %indices of current subset of subunits
                        if strcmp(cur_unique_NL_types{jj},'nonpar')
                            for kk = 1:length(cur_sub_subinds) %if nonpar, need to apply each NL derivative individually
                                fpg(:,cur_sub_subinds(kk)) = nim.subunits(fit_subs(cur_sub_inds(cur_sub_subinds(kk)))).apply_NL_deriv(gint(:,cur_sub_inds(cur_sub_subinds(kk))));
                            end
                        else %otherwise we can apply the NL to all subunits at once
                            fpg(:,cur_sub_subinds) = nim.subunits(fit_subs(cur_sub_inds(cur_sub_subinds(1)))).apply_NL_deriv(gint(:,cur_sub_inds(cur_sub_subinds)));
                        end
                    end
                    target_params = cat(2,param_inds{cur_sub_inds}); %indices of filter coefs for current set of targeted subunits
                    %LL grad is residual * f'(.) *X *w, computed in parallel for all subunits targeting this Xtarg
                    if isempty(gain_funs)
                        penLLgrad(target_params) = bsxfun(@times,(bsxfun(@times,fpg,residual)'*Xstims{un_Xtargs(ii)}),mod_weights(cur_sub_inds))';
                    else
                        penLLgrad(target_params) = bsxfun(@times,(bsxfun(@times,fpg.*gain_funs(:,sub_inds),residual)'*Xstims{un_Xtargs(ii)}),mod_weights(cur_sub_inds))';
                    end
                end
            end
            
            net_penalties = zeros(size(fit_subs));
            net_pen_grads = zeros(length(params),1);
            for ii = 1:length(Tmats) %loop over the derivative regularization matrices
                cur_subs = find([nim.subunits(fit_subs).Xtarg] == Tmats(ii).Xtarg); %set of subunits acting on the stimulus given by this Tmat
                penalties = sum((Tmats(ii).Tmat * cat(2,filtKs{cur_subs})).^2);
                pen_grads = 2*(Tmats(ii).Tmat' * Tmats(ii).Tmat * cat(2,filtKs{cur_subs}));
                cur_lambdas = nim.get_reg_lambdas(Tmats(ii).type,'sub_inds',fit_subs(cur_subs)); %current lambdas
                net_penalties(cur_subs) = net_penalties(cur_subs) + penalties.*cur_lambdas;
                net_pen_grads(cat(2,param_inds{cur_subs})) = net_pen_grads(cat(2,param_inds{cur_subs})) + reshape(bsxfun(@times,pen_grads,cur_lambdas),[],1);
            end
            l2_lambdas = nim.get_reg_lambdas('l2');
            if any(l2_lambdas > 0)
                net_penalties = net_penalties + l2_lambdas.*cellfun(@(x) sum(x.^2),filtKs)';
                net_pen_grads(cat(2,param_inds{:})) = net_pen_grads(cat(2,param_inds{:})) + reshape(2*bsxfun(@times,l2_lambdas,cat(2,filtKs{:})),[],1);
            end
            
            penLL = penLL - sum(net_penalties);
            penLLgrad = penLLgrad - net_pen_grads;
            
            % CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
            Nspks = sum(Robs);
            penLL = -penLL/Nspks;
            penLLgrad = -penLLgrad/Nspks;
            
        end
        
        
        function [LL, LLgrad] = internal_LL_spkNL(nim,params, Robs, G)
            %computes the LL and its gradient for given set of spkNL parameters
                        
            % ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
            nim.spkNL.params = params(1:end-1);
            pred_rate = nim.apply_spkNL(G + params(end));
            
            LL = nim.internal_LL(pred_rate,Robs); %compute LL
            LL_deriv = nim.internal_LL_deriv(pred_rate,Robs);
            spkNL_grad = nim.spkNL_param_grad(params,G);
            LLgrad = sum(bsxfun(@times,spkNL_grad,LL_deriv));
            
            % CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
            Nspks = sum(Robs);
            LL = -LL/Nspks;
            LLgrad = -LLgrad/Nspks;
        end
        
        function [penLL, penLLgrad] = internal_LL_NLs(nim,params, Robs, XNL, Xspkhst, nontarg_g, Tmat,fit_opts)
            %computes the LL and its gradient for a given set of upstream NL parameters
            
            fit_subs = fit_opts.fit_subs;
            % Useful params
            Nfit_subs = length(fit_subs);
            n_TBs = length(nim.subunits(fit_subs(1)).TBx);
            spkhstlen = nim.spk_hist.spkhstlen;
            
            % ESTIMATE GENERATING FUNCTIONS (OVERALL AND INTERNAL)
            theta = params(end); %offset
            G = theta + nontarg_g;
            all_TBy = params(1:Nfit_subs*n_TBs);
            G = G + XNL*all_TBy;
            
            %add contribution from spike history filter
            if fit_opts.fit_spk_hist
                G = G + Xspkhst*params(Nfit_subs*n_TBs + (1:spkhstlen));
            end
            
            pred_rate = nim.apply_spkNL(G);            
            penLL = nim.internal_LL(pred_rate,Robs); %compute LL
            %residual = LL'[r].*F'[g]
            residual = nim.internal_LL_deriv(pred_rate,Robs) .* nim.apply_spkNL_deriv(G, pred_rate < nim.min_pred_rate);
            
            penLLgrad = zeros(length(params),1); %initialize LL gradient
            penLLgrad(1:Nfit_subs*n_TBs) = residual'*XNL;
            penLLgrad(end) = sum(residual);% Calculate derivatives with respect to constant term (theta)
            
            % Calculate derivative with respect to spk history filter
            if fit_opts.fit_spk_hist
                penLLgrad(Nfit_subs*n_TBs+(1:spkhstlen)) = residual'*Xspkhst;
            end
            
            % COMPUTE L2 PENALTIES AND GRADIENTS
            lambdas = nim.get_reg_lambdas('sub_inds',fit_subs,'nld2');
            if any(lambdas > 0)
                TBymat = reshape(all_TBy,n_TBs,[]);
                reg_penalties = lambdas.* sum((Tmat * TBymat).^2);
                pen_grads = 2*(Tmat' * Tmat * TBymat);
                pen_grads = reshape(bsxfun(@times,pen_grads,lambdas),[],1);
                penLL = penLL - sum(reg_penalties);
                penLLgrad(1:Nfit_subs*n_TBs) = penLLgrad(1:Nfit_subs*n_TBs) - pen_grads;
            end
            % CONVERT TO NEGATIVE LLS AND NORMALIZE BY NSPKS
            Nspks = sum(Robs);
            penLL = -penLL/Nspks;
            penLLgrad = -penLLgrad/Nspks;
        end
        
        
        function [LL,grad] = internal_theta_opt(nim,theta,G,Robs)
            %computes LL and its gradient for given additive offset term theta
            G = G + theta;
            pred_rate = nim.apply_spkNL(G);
            LL = nim.internal_LL(pred_rate,Robs);
            %residual = LL'[r].*F'[g]
            residual = nim.internal_LL_deriv(pred_rate,Robs) .* nim.apply_spkNL_deriv(G,pred_rate < nim.min_pred_rate);
            grad = sum(residual);
            Nspks = sum(Robs);
            LL=-LL/Nspks;
            grad=-grad'/Nspks;
        end
        
        
        function [G, fgint, gint] = process_stimulus(nim,Xstims,sub_inds,gain_funs)
%         [G, fgint, gint] = process_stimulus(nim,Xstims,sub_inds,gain_funs)
%         process the stimulus with the subunits specified in sub_inds
%             INPUTS:
%                 Xstims: stimulus as cell array
%                 sub_inds: set of subunits to process
%                 gain_funs: temporally modulated gain of each subunit
%             OUTPUTS:
%                 G: summed generating signal
%                 fgint: output of each subunit
%                 gint: output of each subunit filter
            
            NT = size(Xstims{1},1);
            if isempty(sub_inds);
                [G,fgint,gint] = deal(zeros(NT,1));
                return
            end
            Nsubs = length(sub_inds);
            Xtarg_set = [nim.subunits(sub_inds).Xtarg];
            un_Xtargs = unique(Xtarg_set); %set of Xtargets
            mod_NL_types = {nim.subunits(sub_inds).NLtype}; %NL types for each subunit
            unique_NL_types = unique(mod_NL_types); %unique set of NL types being used
            filtKs = cell(Nsubs,1);
            for ii = 1:Nsubs %loop over subunits, get filter coefs
                filtKs{ii} = nim.subunits(sub_inds(ii)).get_filtK();
            end
            gint = zeros(size(Xstims{1},1),Nsubs);
            for ii = 1:length(un_Xtargs) %loop over the unique Xtargs and compute the generating signals for all relevant filters
                cur_subs = find(Xtarg_set == un_Xtargs(ii)); %set of targeted subunits that act on this Xtarg
                gint(:,cur_subs) = Xstims{un_Xtargs(ii)} * cat(2,filtKs{cur_subs}); %apply filters to stimulus
            end
            fgint = gint; %init subunit outputs by filter outputs
            for ii = 1:length(unique_NL_types) %loop over unique subunit NL types and apply NLs to gint in batch
                if ~strcmp(unique_NL_types{ii},'lin') %if it's not just a linear sub, we have to do something
                    cur_subs = find(strcmp(mod_NL_types,unique_NL_types{ii})); %set of subs with this NL type
                    fgint(:,cur_subs) = nim.subunits(sub_inds(cur_subs(1))).apply_NL(gint(:,cur_subs)); %apply upstream NL to all subunits of this type
                end
            end
            if ~isempty(gain_funs)
                fgint = fgint.*gain_funs(:,sub_inds); %apply gain modulation if needed
            end
            G = fgint*[nim.subunits(sub_inds).weight]';
        end
        
        function LL = internal_LL(nim,rPred,Robs)
            %internal evaluatation method for computing the total LL associated with the predicted rate rPred, given observed data Robs
            switch nim.noise_dist
                case 'poisson' %LL = Rlog(r) - r + C
                    LL = sum(Robs .* log(rPred) -rPred);
                case 'bernoulli' %LL = R*log(r) + (1-R)*log(1-r)
                    LL = sum(Robs.*log(rPred) + (1-Robs).*log(1-rPred));
                case 'gaussian' %LL = (r-R)^2 + c
                    LL = sum((rPred - Robs).^2);
            end
        end
        
        function LL_deriv = internal_LL_deriv(nim,rPred,Robs)
            %computes the derivative of the LL wrt the predicted rate at rPred, given Robs (as a vector over time)
            switch nim.noise_dist
                case 'poisson' %LL'[r] = R/r - 1
                    LL_deriv = Robs./rPred - 1;                    
                case 'bernoulli' %LL'[r] = R/r - (1-R)/(1-r)
                    LL_deriv = Robs./rPred - (1-Robs)./(1-rPred);
                case 'gaussian' %LL'[r] = 2*(r-R)
                    LL_deriv = 2*(rPred - Robs);
            end
            
        end
        
        function rate = apply_spkNL(nim,gen_signal)
            %apply the spkNL function to the input gen_signal. NOTE: the offset term should already be added to gen_signal
            switch nim.spkNL.type
                case 'lin' %F[x;beta] = beta*x
                    rate = gen_signal*nim.spkNL.params(1);
                case 'rectquad' %F[x; beta] = (beta*x)^2 iff x > 0; else 0
                    rate = (gen_signal*nim.spkNL.params(1)).^2;
                    rate(gen_signal < 0) = 0;
                case 'exp' %F[x; beta] = exp(beta*x)
                    rate = exp(gen_signal*nim.spkNL.params(1));
                case 'softplus' %F[x; beta, alpha] = alpha*log(1+exp(beta*x))
                    max_g = 50; %to prevent numerical overflow
                    gint = gen_signal*nim.spkNL.params(1);
                    rate = nim.spkNL.params(2)*log(1 + exp(gint));
                    rate(gint > max_g) = nim.spkNL.params(2)*gint(gint > max_g);
                case 'logistic'
                    rate = 1./(1 + exp(-gen_signal*nim.spkNL.params(1)));
            end
            if ismember(nim.noise_dist,{'poisson','bernoulli'}) %cant allow rates == 0 because LL is undefined
               rate(rate < nim.min_pred_rate) = nim.min_pred_rate; 
            end
            if strcmp(nim.noise_dist,'bernoulli') %cant allow rates == 1 because LL is undefined
               rate(rate > (1 - nim.min_pred_rate)) = 1 - nim.min_pred_rate; 
            end
        end
        
        function rate_deriv = apply_spkNL_deriv(nim,gen_signal,thresholded_inds)
            %apply the derivative of the spkNL to the input gen_signal.
            %Again, gen_signal should have the offset theta already added
            %in
            if nargin < 3
                thresholded_inds = []; %this just specifies the index values where we've had to apply thresholding on the predicted rate to avoid Nan LLs
            end
            switch nim.spkNL.type
                case 'lin' %F'[x; beta] = beta;
                    rate_deriv = nim.spkNL.params(1)*ones(size(gen_signal));
                case 'rectquad' %F'[x; beta] = 2*beta*x iff x > 0; else 0
                    rate_deriv = 2*nim.spkNL.params(1)*gen_signal.*(gen_signal > 0);
                case 'exp' %F'[x; beta] = beta*exp(beta*x)
                    rate_deriv = nim.spkNL.params(1)*exp(nim.spkNL.params(1)*gen_signal);
                case 'softplus' %F[x; beta, alpha] = alpha*beta*exp(beta*x)/(1+exp(beta*x))
                    max_g = 50; %to prevent numerical overflow
                    gint = gen_signal*nim.spkNL.params(1);
                    rate_deriv = nim.spkNL.params(1)*nim.spkNL.params(2)*exp(gint)./(1 + exp(gint));
                    rate_deriv(gint > max_g) = nim.spkNL.params(1)*nim.spkNL.params(2); %e^x/(1+e^x) => 1 for large x
                case 'logistic'
                    rate_deriv = nim.spkNL.params(1)*exp(-gen_signal*nim.spkNL.params(1))./...
                        (1 + exp(-gen_signal*nim.spkNL.params(1))).^2; %e^(-x)/(1+e^(-x))^2
            end
            rate_deriv(thresholded_inds) = 0; %if thresholding the rate to avoid undefined LLs, set deriv to 0 at those points
        end
        
        function rate_grad = spkNL_param_grad(nim,params,x)
            %computes the gradient of the spkNL function with respect to its parameters (subroutine
            %for optimizing the spkNL params)
            
            rate_grad = zeros(length(x),length(params));
            switch nim.spkNL.type
                case 'lin'
                    rate_grad(:,1) = x; %dr/dbeta = x
                    rate_grad(:,2) = ones(size(x)); %dr/dtheta = 1
                case 'rectlin'
                    rate_grad(:,1) = gen_signal; %dr/dbeta = x iff x > 0
                    rate_grad(:,2) = ones(size(x)); %dr/dtheta = 1 iff x > 0
                    rate_grad(x < 0,:) = 0;
                case 'exp'
                    temp = exp(params(1)*(x + params(end)));
                    rate_grad(:,1) = (x + params(end)).*temp; %dr/dbeta = (x+theta)*exp(beta*(x+theta))
                    rate_grad(:,2) = params(1).*temp; %dr/dtheta = beta*exp(beta*(x+theta))
                case 'softplus'
                    temp = params(2)*exp(params(1)*(x + params(3)))./(1 + exp(params(1)*(x + params(3)))); %alpha*exp(beta*(x+theta))/(1 + exp(beta*(x+theta)))
                    rate_grad(:,1) = temp.*(x + params(3)); %dr/dbeta = temp*(x + theta)
                    rate_grad(:,2) = log(1 + exp(params(1)*(x + params(3)))); %dr/dalpha = log[]
                    rate_grad(:,3) = temp.*params(1); %dr/dtheta = temp*beta
                case 'logistic'
                    %                     nim.spkNL.params = [1]; %defines beta in f(x; beta) = 1/(1+exp(-beta*x))
                    error('undefined yet');
                otherwise
                    error('unsupported spkNL type');
            end            
        end
        
        function Tmats = make_Tikhonov_matrices(nim)
%         creates a struct containing the Tikhonov regularization matrices, given the stimulus
%         and regularization parameters specified in the nim
            Nstims = length(nim.stim_params); %number of unique stimuli
            Xtargs = [nim.subunits(:).Xtarg];
            
            deriv_reg_types = nim.allowed_reg_types(strncmp(nim.allowed_reg_types,'d',1)); %set of regularization types where we need a Tikhonov matrix
            cnt = 1;
            Tmats = [];
            for ii = 1:Nstims %for each stimulus
                cur_subs = find(Xtargs == ii); %get set of subunits acting on this stimuls
                for jj = 1:length(deriv_reg_types) %check each possible derivative regularization type
                    cur_lambdas = nim.get_reg_lambdas(deriv_reg_types{jj},'sub_inds',cur_subs);
                    if any(cur_lambdas > 0)
                        cur_Tmat = create_Tikhonov_matrix(nim.stim_params(ii),deriv_reg_types{jj});
                        Tmats(cnt).Tmat = cur_Tmat;
                        Tmats(cnt).Xtarg = ii;
                        Tmats(cnt).type = deriv_reg_types{jj};
                        cnt = cnt + 1;
                    end
                end
            end
        end
        
        function Tmat = make_NL_Tmat(nim)
            %make Tikhonov matrix for smoothness regularization of the TB NLs
            nonpar_set = find(strcmp(nim.get_NLtypes(),'nonpar'));
            assert(~isempty(nonpar_set),'no nonparametric NLs found');
            n_tbs = length(nim.subunits(nonpar_set(1)).TBx); %number of TBx (assume this is the same for all subunits)!
            et = ones(n_tbs,1);
            et([1 end]) = 0; %free boundaries
            Tmat = spdiags([et -2*et et], [-1 0 1], n_tbs, n_tbs)';
        end
        
        function [] = check_inputs(nim,Robs,Xstims,sub_inds,gain_funs)
            % checks the format of common inputs params
            if nargin < 4
                gain_funs = [];
            end
            if nargin < 5
                sub_inds = nan;
            end
            Nsubs = length(nim.subunits);
            for n = 1:Nsubs %check that stimulus dimensions match
                [NT,filtLen] = size(Xstims{nim.subunits(n).Xtarg}); %stimulus dimensions
                assert(filtLen == prod(nim.stim_params(nim.subunits(n).Xtarg).dims),'Xstim dims dont match stim_params');
            end
            assert(length(unique(cellfun(@(x) size(x,1),Xstims))) == 1,'Xstim elements need to have same size along first dimension');
            assert(size(Robs,2) == 1,'Robs must be a vector');
            assert(iscell(Xstims),'Xstims must for input as a cell array');
            if ~isempty(gain_funs)
                assert(size(gain_funs,1) == NT & size(gain_funs,2) == Nsubs,'format of gain_funs is incorrect');
            end
            if ~isnan(sub_inds)
                assert(min(sub_inds) > 0 & max(sub_inds) <= NT,'invalid data indices specified');
            end
        end
        
        function nim = set_subunit_scales(nim,fgint)
            %sets the 'scale' of each subunit based on the SD of its output fgint
            fgint_SDs = std(fgint);
            for ii = 1:length(nim.subunits)
                nim.subunits(ii).scale = fgint_SDs(ii);
            end
        end
        
    end
    
    methods (Static, Hidden)
        function stim_params = check_stim_params(stim_params)
%           internal function that checks stim_params struct formatting, and initializes default values if needed
            
            default_params.boundary_conds = [Inf 0 0];%[free boundary on first dim, and tied to 0 in other dims]
            default_params.split_pts = [];%no discontinuities in smoothing penalty
            default_params.dt = 1; %unitless measure of time scale
            default_params.dx = 1; %unitless measure of spatial scale
            default_params.name = 'unnamed'; %default name for the stimulus
            
            for ii = 1:length(stim_params)
                cur_stim_params = default_params; %start with default struct
                assert(isfield(stim_params(ii),'dims'),'need to specify the stimulus dimensions');
                spec_fields = fieldnames(stim_params(ii));
                for jj = 1:length(spec_fields) %override default values with user-specified ones
                    value = getfield(stim_params(ii),spec_fields{jj});
                    cur_stim_params = setfield(cur_stim_params,spec_fields{jj},value);
                end
                while length(cur_stim_params.dims) < 3 %pad dims with 1s for book-keeping
                    cur_stim_params.dims = cat(2,cur_stim_params.dims,1);
                end
                new_stim_params(ii) = cur_stim_params;
            end
            stim_params = new_stim_params;
        end
        
        function optim_params = set_optim_params(optimizer,input_params,silent)
            %internal function that checks stim_params struct formatting, and initializes default values for the given optimizer
            
            optim_params.maxIter = 500; %maximum number of iterations
            switch optimizer
                case 'fminunc'
                    optim_params.TolX = 1e-7; % termination tolerance on X
                    optim_params.TolFun = 1e-7; % termination tolerance on the function value
                    optim_params.LargeScale = 'off'; %dont use large-scale method
                    optim_params.HessUpdate = 'bfgs'; %update Hessian using BFGS
                    optim_params.GradObj = 'on'; %use gradient
                    if silent
                        optim_params.Display = 'off';
                    else
                        optim_params.Display = 'iter';
                    end
                case 'minFunc'
                    optim_params.optTol = 1e-5; %[minFunc] termination tolerance on first order optimality (max(abs(grad))
                    optim_params.progTol = 1e-8; %[minFunc] termination tolerance on function/parameter values
                    optim_params.Method = 'lbfgs'; %[minFunc] method
                    optim_params.verbose = 2; %display full iterative output
                    if silent
                        optim_params.Display = 'off';
                    else
                        optim_params.Display = 'iter';
                    end
                    
                case 'fmincon'
                    optim_params.Algorithm = 'active-set';
                    optim_params.GradObj = 'on';
                    optim_params.TolX = 1e-7;
                    optim_params.TolFun = 1e-7;
                    if silent
                        optim_params.Display = 'off';
                    else
                        optim_params.Display = 'iter';
                    end
                case 'L1GeneralPSSas'
                    optim_params.optTol = 1e-5;
                    optim_params.progTol = 1e-8;
                    optim_params.verbose = 2;
                    if silent
                        optim_params.verbose = 0;
                    else
                        optim_params.verbose = 2;
                    end
                case 'minConf_TMP'
                    optim_params.optTol = 1e-5;
                    optim_params.progTol = 1e-8;
                    optim_params.verbose = 2;
                    if silent
                        optim_params.verbose = 0;
                    else
                        optim_params.verbose = 2;
                    end
                    
                otherwise
                    error('unsupported optimizer');
                    
            end
            
            %load in specified parameters
            if ~isempty(input_params)
                spec_fields = fieldnames(input_params);
                for ii = 1:length(spec_fields)
                    value = getfield(input_params,spec_fields{ii});
                    optim_params = setfield(optim_params,spec_fields{ii},value);
                end
            end
        end
    end
end