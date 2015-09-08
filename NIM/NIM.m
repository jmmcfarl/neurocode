
classdef NIM
    
    %Summary of this class goes here
    %   Detailed explanation goes here
    %
    % James McFarland, September 2015
    
    %%
    properties
        spkNL;          %struct defining the spiking NL function
        subunits;       %array of subunit objects
        stim_params;    %struct array of parameters characterizing the stimuli that the model acts on, must have a .dims field
        noise_dist;     %noise distribution class specifying the noise model
        spk_hist;       %class defining the spike-history filter properties
        fit_props;      %struct containing information about model fit evaluations
        init_props;     %struct containing details about model initialization
    end
    
    properties (Hidden)
        allowed_reg_types = {'nld2','d2xt','d2x','d2t','l2','l1'}; %set of allowed regularization types
        version = '0.0';
    end
    %%
    methods
        %% CONSTRUCTOR
        function nim = NIM(stim_params, NLtypes, mod_signs, varargin)
            %             nim = NIM(stim_params, NLtypes, mod_signs, <Xtargets>,<spkNL>,<noise_dist>)
            %             constructor for class NIM
            %             INPUTS:
            %                   stim_params: struct array defining parameters for each stimulus the model acts on. Must specify the .dims field for each stim
            %                   NLtypes: string or cell array of strings specifying the
            %                     upstream NL type associated with each subunit. If it's a
            %                     single string, we use the same NL type throughout
            %                     mod_signs: vector specifying the weight associated with each subunit (typically +/- 1)
            %                     optional_flags:
            %                         <Xtargets>: vector specifying the index of the stimulus each subunit acts on (defaults to ones)
            %                         <spkNL>: string specifying type of spkNL function
            %                         <noise_dist>: string specifying type of noise distribution
            %                         <init_filts>: cell array of initial filter values
            
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
                    case 'init_filts'
                        init_filts = flag_val;
                    otherwise
                        error('Invalid input flag');
                end
                j = j + 2;
            end
            
            %check and create spk NL function
            allowed_spkNLs = {'lin','rectlin','exp','softplus','logistic'}; %set of NL functions currently implemented
            assert(ismember(spkNL,allowed_spkNLs),'not an allowed spk NL type');
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
                nim.subunits = cat(1,nim.subunits,SUBUNIT(init_filt, mod_signs(ii), NLtypes{ii},Xtargets(ii)));
            end
            
            %initialize WITHOUT spike history term
            spk_hist.coefs = [];
            spk_hist.bin_edges = [];
            spk_hist.spkhstlen = 0;
            nim.spk_hist = spk_hist;
            
        end
        
        %% setting methods
        function nim = set_reg_params(nim, varargin)
            %             set a desired set of regularization parameters to specified values, apply to specified set of subunits
            %             optional flags:
            %               'sub_inds': set of subunits to apply the new reg_params for
            %               'lambda_type': type of regularization
            
            sub_inds = 1:length(nim.subunits); %default is to apply the change to all subunits
            
            %INPUT PARSING
            j = 1;
            reg_types = {}; reg_vals = [];
            while j <= length(varargin)
                flag_name = varargin{j};
                flag_val = varargin{j+1};
                switch lower(flag_name)
                    case 'sub_inds'
                        sub_inds = flag_val;
                        assert(all(ismember(sub_inds,1:length(nim.subunits))),'invalid target subunits specified');
                    case nim.allowed_reg_types
                        reg_types = cat(1,reg_types,lower(flag_name));
                        reg_vals = cat(1,reg_vals,flag_val);
                    otherwise
                        error('Invalid input flag');
                end
                j = j + 2;
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
        
        function nim = set_stim_params(nim, varargin)
            %             adjust stimulus parameters for the desired Xtarg
            %             optional flags:
            %                 'xtarg': index of stimulus to apply the new stim_params for [default is 1]
            %                 'dims': dimensionality of stim: [Tdim, X1dim, X2dim] where Tdim is the number of temporal dimensions, etc.
            %                 'boundary_conds': boundary conditions on each stim dimension [Inf means free, 0 means tied, -1 means periodic]
            %                 'split_pts': set of stimulus indices to specify smoothing boundaries (NOT YET IMPLEMENTED)
            
            Xtarg = 1; %default is to apply the change to stim 1
            allowed_flags = {'dims','boundary_conds','split_pts'}; %fields of the stim_params struct that we might want to set
            %INPUT PARSING
            j = 1; fields_to_set = {}; field_vals = {};
            while j <= length(varargin)
                flag_name = varargin{j};
                flag_val = varargin{j+1};
                switch lower(flag_name)
                    case 'xtarg'
                        Xtarg = flag_val;
                    case allowed_flags
                        fields_to_set = cat(1,fields_to_set,flag_name);
                        field_vals = cat(1,field_vals,flag_val);
                    otherwise
                        error('Invalid input flag');
                end
                j = j + 2;
            end
            
            if isempty(field_vals)
                warning('No stim_params values specifid to change');
            end
            if Xtarg > length(nim.stim_params) %if we're creating a new stimulus
                assert(ismember('dims',fields_to_set),'need to specify dims to initialize new stimulus');
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
            %            get regularizatoin labmda values from nim subunits
            %            INPUT FLAGS:
            %                 sub_inds, followed by vector of indices, specifies which subunits to extract lambdas from
            %                 lambda_type: string of allowed regularization types
            %            OUTPUTS:
            %                 lambdas: [K,N] matrix of lambda values, K is the number of specified lambda_types and N is the number of subunits
            
            sub_inds = 1:length(nim.subunits); %default is to grab reg values from all subunits
            %INPUT PARSING
            jj = 1;
            reg_types = {};
            while jj <= length(varargin)
                flag_name = varargin{jj};
                switch lower(flag_name)
                    case 'sub_inds'
                        sub_inds = varargin{jj+1};
                        assert(all(ismember(sub_inds,1:length(nim.subunits))),'invalid target subunits specified');
                        jj = jj + 2;
                    case nim.allowed_reg_types
                        reg_types = cat(1,reg_types,lower(flag_name));
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
        
        
        function filtKs = get_filtKs(nim,sub_inds)
            %             filtKs = get_filtKs(sub_inds)
            %             get filters for specified set of subunits. Returns a cell array.
            %             Optional input argument sub_inds specifies subunits to get filters from
            
            Nsubs = length(nim.subunits);
            if nargin < 2
                sub_inds = 1:Nsubs; %default is to grab filters for all subunits
            end
            filtKs = cell(Nsubs,1);
            for ii = sub_inds
                filtKs{ii} = nim.subunits(ii).get_filtK;
            end
        end
        
        
        function NLtypes = get_NLtypes(nim)
            %NLtypes = get_NLtypes()
            %gets cell array of strings specifying NLtype of each subunit
            
            Nsubs = length(nim.subunits);
            NLtypes = cell(Nsubs,1);
            for ii = 1:Nsubs
                NLtypes{ii} = nim.subunits(ii).NLtype;
            end
        end
        
        
        function [filt_penalties,NL_penalties] = get_reg_pen(nim,Tmats)
            %  [filt_penalties,NL_penalties = get_reg_pen(<Tmats>)
            %calculates the regularization penalties on each subunit,
            %separately for filter and NL regularization
            
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
            %add subunits to the model with specified properties. Can add
            %multiple subunits at the same time
            %             INPUTS:
            %                 NLtypes: string or cell array of strings specifying the type of upstream NLs for each subunit
            %                 mod_signs: the weight associated with each subunit (typically +/- 1)
            %                 optional flags:
            %                     xtargs: vector of Xtargets for each added subunit
            %                     init_filts: cell array of initial filter values for each subunit
            %             OUPUTS:
            %                 nim: new model object
            
            if ~iscell(NLtypes) && ischar(NLtypes); NLtypes = cellstr(NLtypes); end;
            nSubs = length(mod_signs); %number of subunits being added
            nStims = length(nim.stim_params);
            Xtargets = ones(nSubs,1); %default Xtargets to 1
            if length(NLtypes) == 1 && nSubs > 1
                NLtypes = repmat(NLtypes,nSubs,1); %if NLtypes is specified as a single string, assume we want this NL for all subunits
            end
            init_filts = cell(nSubs,1);
            
            %parse input flags
            j = 1; %initialize counter after required input args
            while j <= length(varargin)
                flag_name = varargin{j};
                flag_val = varargin{j+1};
                switch lower(flag_name)
                    case 'xtargs'
                        Xtargets = flag_val;
                        assert(all(ismember(Xtargets,1:nStims)),'invalid Xtargets specified');
                    case 'init_filts'
                        if ~iscell(flag_val)
                            init_filts = cell(length(mod_signs),1);
                            for ii = 1:length(mod_signs)
                                init_filts{ii} = flag_val(:,ii);
                            end
                        else
                            init_filts = flag_val;
                        end
                    otherwise
                        error('Invalid input flag');
                end
                j = j + 2;
            end
            
            assert(length(Xtargets) == nSubs,'length of mod_signs and Xtargets must be equal');
            %initialize subunits
            for ii = 1:nSubs %loop initializing subunits (start from last to initialize object array)
                stimD = prod(nim.stim_params(Xtargets(ii)).dims); %dimensionality of the current filter
                if isempty(init_filts{ii})
                    init_filt = randn(stimD,1)/stimD; %initialize fitler coefs with gaussian noise
                else
                    if iscell(init_filts)
                    init_filt = init_filts{ii};
                    else
                        init_filt = init_filts(ii,:);
                    end
                end
                existing_sub = find(nim.subunits(:).Xtarg == Xtargets(ii),1); %find any existing subunits with this same Xtarget
                nim.subunits = cat(1,nim.subunits,SUBUNIT(init_filt, mod_signs(ii), NLtypes{ii},Xtargets(ii)));
                if ~isempty(existing_sub)
                    default_lambdas = nim.subunits(existing_sub).reg_lambdas;
                    nim.subunits(end).reg_lambdas = default_lambdas;
                end
            end
        end
        
        % initialize spike history filter
        function nim = init_spkhist(nim,n_bins,varargin)
            %
            % Usage: nim = NIMinit_spkhist(nim,n_bins,init_spacing,doubling_time,negCon)
            %
            % Adds a spike history term with specified parameters to an existing NIM.
            %
            % INPUTS:
            %     nim: input model structure
            %     n_bins: number of coefficients in spike history filter
            %     optional flags
            %       [init_spacing, XX]: Initial spacing (in time bins) of piecewise constant 'basis functions'
            %       [doubling_time, XX]: Make bin spacing logarithmically increasing, with given doubling time
            %       [negCon]: If flag is included, constrain spike history filter coefs to be non-positive
            % OUTPUTS:
            %     nim: output model
            
            % default inputs
            init_spacing = 1;
            doubling_time = n_bins;
            negCon = false;
            %parse input flags
            j = 1; %initialize counter after required input args
            while j <= length(varargin)
                flag_name = varargin{j};
                switch lower(flag_name)
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
            inc = init_spacing;
            pos = 1; count = 0;
            for n = 1:n_bins+1
                bin_edges(n) = pos;
                pos = pos + inc;  count = count + 1;
                if count >= doubling_time
                    count = 0; inc = inc * 2;
                end
            end
            
            % LOAD INTO A SPK_HIST STRUCT IN NIM
            nim.spk_hist.bin_edges = bin_edges;
            nim.spk_hist.coefs = zeros(n_bins,1); %init filter coefs to 0
            nim.spk_hist.negCon = negCon;
            nim.spk_hist.spkhstlen = n_bins;
        end
        
        %% initialize nonparametric subunit
        function nim = init_nonpar_NLs(nim, Xstims, varargin)
            %
            % Usage: nim = NMMinitialize_upstreamNLs( nim, Xstims, varargin)
            %
            % Initializes the specified model subunits to have nonparametric
            % (tent-basis) upstream NLs
            %
            % INPUTS:
            %     Xstim: time-embedded stim matrix
            %     <sub_inds>: vector specifying which subunits will be made to have nonparametric upstream NLs
            %     <lambda_nld2>: specifies strength of smoothness regularization for the tent-basis coefs
            %     <NLmon>: Set to +1 to constrain NL coefs to be monotonic increasing and -1 to make mono.
            %        decreasing. 0 means no constraint. Default here is +1 (monotonic increasing)
            %     <edge_p>: Determines the locations of the outermost tent-bases relative to the underlying generating distribution
            %     <n_bfs>: Number of tent-basis functions to use
            %     <space_type>: Use either 'equispace' for uniform bin spacing, or 'equipop' for 'equipopulated bins'
            %
            % OUTPUTS:
            %     nim: new model object
            
            Nsubs = length(nim.subunits);
            sub_inds = 1:Nsubs; %defualt to fitting all subunits (plus -1 for spkHist filter)
            NLmon = 1; %default monotonic increasing TB-coefficients
            edge_p = 0.05; %relative to the generating distribution (pth percentile) where to put the outermost tent bases
            n_bfs = 25; %default number of tent basis functions
            space_type = 'equispace'; %default uninimrm tent basis spacing
            lambda_nld2 = 0; %default no smoothing on TB coefs
            
            j = 1;
            while j <= length(varargin)
                flag_name = varargin{j};
                flag_val = varargin{j+1};
                switch lower(flag_name)
                    case 'sub_inds'
                        sub_inds = flag_val;
                        assert(all(ismember(sub_inds,[1:Nsubs])),'invalid target subunits specified');
                    case 'nlmon'
                        NLmon = flag_val;
                        assert(ismember(NLmon,[-1 0 1]),'NLmon must be a -1, 0 or 1');
                    case 'edge_p'
                        edge_p = flag_val;
                        assert(edge_p > 0 & edge_p < 100,'edge_p must be between 0 and 100');
                    case 'n_bfs'
                        n_bfs = flag_val;
                        assert(n_bfs > 0,'n_bfs must be greater than 0');
                    case 'space_type'
                        space_type = flag_val;
                        assert(ismember(space_type,{'equispace','equipop'}),'unsupported bin spacing type');
                    case 'lambda_nld2'
                        lambda_nld2 = flag_val;
                        assert(lambda_nld2 >= 0,'lambda must be >= 0');
                    otherwise
                        error('Invalid input flag');
                end
                j = j + 2;
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
                %set nearest tent basis to 0
                [~,nearest] = min(abs(TBx));
                TBx(nearest) = 0;
                
                %initalize tent basis coefs
                switch prev_NL_types{imod}
                    case 'lin'
                        TBy = TBx;
                    case 'rectlin'
                        TBy = TBx; TBy(TBy < 0) = 0;
                    case 'quad'
                        TBy = TBx.^2;
                    case 'nonpar'
                        fprintf('upstream NL already set as nonparametric\n');
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
            %           [LL, pred_rate, mod_internals, LL_data] = eval_model(nim, Robs, Xstims, <eval_inds>, varargin)
            %             Evaluates the mdoel on the supplied data
            %                 INPUTS:
            %                     Robs: vector of observed data
            %                     Xstims: cell array of stimuli
            %                     <eval_inds>: optional vector of indices on which to evaluate the model
            %                     optional flags:
            %                         gain_funs: matrix specifying temporal gain functions for each subunit
            %                 OUTPUTS:
            %                     LL: log-likelihood per spike
            %                     pred_rate: predicted firing rates (in counts/bin)
            %                     mod_internals: struct containing the internal components of the model prediction
            %                         G: is the total generating signal (not including
            %                           the constant offset theta). This is the sum of
            %                           subunit outputs (weighted by their subunit weights w)
            %                         fgint: is the output of each subunit
            %                         gint: is the output of each subunits linear filter
            
            Nsubs = length(nim.subunits); %number of subunits
            NT = length(Robs); %number of time points
            eval_inds = nan; %this default means evaluate on all data
            % PROCESS INPUTS
            gain_funs = []; %default has no gain_funs
            j = 1;
            while j <= length(varargin)
                flag_name = varargin{j};
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
        
        function [] = display_model(nim,varargin)
            % NIMdisplay_model(nim,<Xstim>,varargin)
            %
            % Creates a display of the elements of a given NIM
            % INPUTS:
            %     nim: model structure
            %     <Xstim>: provide the stimulus matrix if you want to display the distributions of generating signals
            % ADD IN GAIN FUNCTIONS?
            
            Xstims = [];
            Xtargs = [1:length(nim.stim_params)]; %default plot filters for all stimuli
            plot_spkNL = true;
            plot_spk_hist = true;
            gain_funs = [];
            j = 1; %initialize counter after required input args
            while j <= length(varargin)
                flag_name = varargin{j};
                if ~ischar(flag_name)
                    Xstims = flag_name;
                    assert(iscell(Xstims),'second argument must be Xstims if its not an option flag');
                    j = j+1; %account for the fact that theres no input value here
                else
                    flag_val = varargin{j+1};
                    switch lower(flag_name)
                        case 'xtargs'
                            Xtargs = flag_val;
                            assert(all(ismember(Xtargs,1:length(nim.stim_params))),'invalid Xtargets specified');
                        case 'plot_spknl'
                            plot_spkNL = flag_val;
                            assert(ismember(plot_spkNL,[0 1]),'plot_spkNL must be 0 or 1');
                        case 'gain_funs'
                            gain_funs = lower(flag_val);
                            assert(isstr(noise_dist),'noise_dist must be a string');
                        case 'plot_spk_hist'
                            plot_spk_hist = flag_val;
                            assert(ismember(plot_spk_hist,[0 1]),'plot_spk_hist must be a 0 or 1');
                        otherwise
                            error('Invalid input flag');
                    end
                    j = j + 2;
                end
            end
            
            Nsubs = length(nim.subunits);
            spkhstlen = nim.spk_hist.spkhstlen;
            if spkhstlen > 0
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
            %
            %             Usage: nim_out = NMMfit_filters( nim, Robs, Xstims, <train_inds>, varargin)
            %
            %             INPUTS:
            %                Robs: vector of observations (e.g. spike counts)
            %                Xstims: matrix or cell array of stimuli
            %                <train_inds>: index values on which to fit the model [default to all indices in provided data]
            %                   optional flags:
            %                      fit_subs: set of subunits whos filters we want to optimize [default is all]
            %                      gain_funs: matrix of multiplicative factors, one column for each subunit
            %                      optim_params: struct of desired optimization parameters
            %                      silent: include this flag to suppress the iterative optimization display
            %                      no_spkhist: include this flat to hold the spk NL filter constant
            %             OUTPUTS:
            %                new nim object with optimized subunit filters
            
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
                        case 'fit_subs'
                            fit_subs = varargin{j+1};
                            assert(all(ismember(fit_subs,[-1 1:Nsubs])),'invalid target subunits specified');
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
                        case 'no_spkhist'
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
            %             sign_con = [];
            for imod = fit_subs
                cur_kern = nim.subunits(imod).filtK;
                %                 if isfield(nim.mods(imod),'Kcon')
                %                     if nim.mods(imod).Kcon ~= 0
                %                         sign_con(length(initial_params)+(1:length(cur_kern))) = nim.mods(imod).Kcon;
                %                     end
                %                 end
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
            A = []; Aeq = []; %initialize constraint matrices
            LB = -Inf*ones(size(init_params));
            UB = Inf*ones(size(init_params));
            %             % Constrain any of the filters to be positive or negative
            %             if ~isempty(sign_con)
            %                 LB(sign_con == 1) = 0;
            %                 UB(sign_con == -1) = 0;
            %                 use_con = 1;
            %             end
            if spkhstlen > 0 && fit_spk_hist %if optimizing spk history term
                %negative constraint on spk history coefs
                if nim.spk_hist.negCon
                    spkhist_inds = Nfit_filt_params + (1:spkhstlen);
                    UB(spkhist_inds) = 0;
                    use_con = 1;
                end
            end
            
            beq = zeros(size(Aeq,1),1);
            b = zeros(size(A,1),1);
            
            % GENERATE REGULARIZATION MATRICES
            Tmats = nim.make_Tikhonov_matrices();
           
            fit_opts = struct('fit_spk_hist', fit_spk_hist, 'fit_subs',fit_subs); %put any additional fitting options into this struct
            %the function we want to optimize
            opt_fun = @(K) nim.internal_LL_filters(K,Robs,Xstims,Xspkhst,nontarg_g,Tmats,fit_opts);
            
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
                    [params] = fmincon(opt_fun, init_params, A, b, Aeq, beq, LB, UB, [], optim_params);
            end
            [penLL,penGrad] = opt_fun(params);
            first_order_optim = max(abs(penGrad));
            
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
            
            [LL,pred_rate,mod_internals,LL_data] = nim.eval_model(Robs,Xstims);
            nim = nim.set_subunit_scales(mod_internals.fgint); %update filter scales
            
            %         nim_out.LL_seq = cat(1,nim_out.LL_seq,LL);
            %         nim_out.penLL_seq = cat(1,nim_out.penLL_seq,penLL);
            %         nim_out.opt_history = cat(1,nim_out.opt_history,{'filt'});
        end
        
        %%
        function nim = fit_upstreamNLs(nim, Robs, Xstims, varargin)
            %
            %             Usage: nim = NMMfit_upstreamNLs( nim, Robs, Xstims, <train_inds>, varargin)
            %
            %             Optimizes the upstream NLs (in terms of tent-basis functions) (plus extra linear terms if desired) for
            %             given stimulus filters
            %
            %             INPUTS:
            %                   nim: model structure
            %                   Robs: binned spikes
            %                   Xstim: time-embedded stimulus mat
            %                   <train_inds>: indices of data to optimize on
            %                       Optional flags:
            %                       <fit_subs>: Vector of indices specifying which subunits to optimize.
            %                           (-1 = spk history filter) Default is to optimize all elements
            %                       <gain_funs>: matrix of gain values for each subunit
            %                       <no_rescaling>: use this flag if you dont want to rescale the NLs after fitting
            %                       <silent>: include this flag if you want to turn off the optimization display
            %                       <optim_params>: Struct of optimization parameters
            %                       <no_spkhist>: include this flag if you dont want to fit the spkNL filter
            %            OUTPUTS:
            %                   nim: output model struct
            
            Nsubs = length(nim.subunits); %number of subunits
            NT = length(Robs); %number of time points
            
            % PROCESS INPUTS
            poss_targets = find(strcmp(nim.get_NLtypes,'nonpar'))';
            fit_subs = poss_targets; %defualt to fitting all subunits
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
                        case 'fit_subs'
                            fit_subs = varargin{j+1};
                            assert(all(ismember(fit_subs,[poss_targets])),'specified target doesnt have non-parametric NL, or doesnt exist');
                        case 'gain_funs'
                            gain_funs = varargin{j+1};
                        case 'optim_params'
                            optim_params = varargin{j+1};
                            assert(isstruct(optim_params),'optim_params must be a struct');
                        case 'silent'
                            silent = true;
                        case 'no_rescaling'
                            rescale_NLs = false;
                        case 'no_spkhist'
                            fit_spk_hist = false;
                        otherwise
                            error('Invalid input flag');
                    end
                    j = j + 2;
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
                %                 if isempty(Gmults{tar})
                tbf_out = nim.subunits(tar).weight * nim.subunits(tar).tb_rep(gint);
                %                 else
                %                     %tbf_out = nim.mods(tar).sign * (Gmults{tar} .* tb_rep(gint(:,tar),nim.mods(tar).NLx));
                %                     tbf_out = nim.mods(tar).sign * bsxfun(@times, tb_rep(gint(:,tar),nim.mods(tar).NLx), Gmults{tar} );
                %                 end
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
            
            % IF RESCALING NLS, NEED TO RE-ESTIMATE OPTIMAL THETA
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
        end
        
    end
    
    methods (Hidden)
        %% internal methods
        
        function [penLL, penLLgrad] = internal_LL_filters(nim,params,Robs,Xstims,Xspkhst,nontarg_g,Tmats,fit_opts)
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
            %             if isempty(gain_funs)
            G = G + fgint*mod_weights;
            %             else
            %                 error('not implemented yet')
            %             end
            
            % Add contribution from spike history filter
            if fit_opts.fit_spk_hist
                G = G + Xspkhst*params(NKtot + (1:nim.spk_hist.spkhstlen));
            end
            
            pred_rate = nim.apply_spkNL(G);
            % Enforce minimum predicted firing rate to avoid nan LLs
            if ismember(nim.noise_dist,{'poisson','bernoulli'})
                min_pred_rate = 1e-50;
                if min(pred_rate) < min_pred_rate
                    pred_rate(pred_rate < min_pred_rate) = min_pred_rate; %minimum predicted rate
                end
            end
            
            penLL = nim.internal_LL(pred_rate,Robs); %compute LL
            
            %residual = LL'[r].*F'[g]
            residual = nim.internal_LL_deriv(pred_rate,Robs) .* nim.apply_spkNL_deriv(G);
            
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
                    penLLgrad(param_inds{cur_sub_inds}) = residual'*Xstims{un_Xtargs(ii)} * nim.subunits(cur_sub_inds).weight;
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
                    penLLgrad(target_params) = bsxfun(@times,(bsxfun(@times,fpg,residual)'*Xstims{un_Xtargs(ii)}),mod_weights(cur_sub_inds))';
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
        
        
        function [penLL, penLLgrad] = internal_LL_NLs(nim,params, Robs, XNL, Xspkhst, nontarg_g, Tmat,fit_opts)
            % DESCRIPTION HERE
            
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
            % Enforce minimum predicted firing rate to avoid nan LLs
            min_pred_rate = 1e-50;
            if min(pred_rate) < min_pred_rate
                pred_rate(pred_rate < min_pred_rate) = min_pred_rate; %minimum predicted rate
            end
            
            penLL = nim.internal_LL(pred_rate,Robs); %compute LL
            %residual = LL'[r].*F'[g]
            residual = nim.internal_LL_deriv(pred_rate,Robs) .* nim.apply_spkNL_deriv(G);
            
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
            %DESCRIPTION HERE
            G = G + theta;
            pred_rate = nim.apply_spkNL(G);
            %enforce minimum predicted firing rate to avoid nan LLs
            min_pred_rate = 1e-50;
            if min(pred_rate) < min_pred_rate
                pred_rate(pred_rate < min_pred_rate) = min_pred_rate; %minimum predicted rate
            end
            LL = nim.internal_LL(pred_rate,Robs);
            %residual = LL'[r].*F'[g]
            residual = nim.internal_LL_deriv(pred_rate,Robs) .* nim.apply_spkNL_deriv(G);
            grad = sum(residual);
            Nspks = sum(Robs);
            LL=-LL/Nspks;
            grad=-grad'/Nspks;
        end
        
        
        function [G, fgint, gint] = process_stimulus(nim,Xstims,sub_inds,gain_funs)
            %             [G, fgint, gint] = process_stimulus(nim,Xstims,sub_inds,gain_funs)
            %             process the stimulus with the subunits specified in sub_inds
            %                 INPUTS:
            %                     Xstims: stimulus as cell array
            %                     sub_inds: set of subunits to process
            %                     gain_funs: temporally modulated gain of each subunit
            %                 OUTPUTS:
            %                     G: summed generating signal
            %                     fgint: output of each subunit
            %                     gint: output of each subunit filter
            
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
                fgint = fgint.*gain_funs; %apply gain modulation if needed
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
            %internal method for computing the derivative of the LL wrt the
            %predicted rate at rPred, given Robs (as a vector over time)
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
            %apply the spkNL function to the input gen_signal. NOTE: the
            %offset term should already be added to gen_signal
            switch nim.spkNL.type
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
                case 'logistic'
                    rate = 1./(1 + exp(-gen_signal*nim.spkNL.params(1)));
            end
        end
        
        function rate_deriv = apply_spkNL_deriv(nim,gen_signal)
            %apply the derivative of the spkNL to the input gen_signal.
            %Again, gen_signal should have the offset theta already added
            %in
            switch nim.spkNL.type
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
                case 'logistic'
                    rate_deriv = nim.spkNL.params(1)*exp(-gen_signal*nim.spkNL.params(1))./...
                        (1 + exp(-gen_signal*nim.spkNL.params(1))).^2; %e^(-x)/(1+e^(-x))^2
            end
        end
        
        function Tmats = make_Tikhonov_matrices(nim)
            %             creates a struct containing the Tikhonov regularization matrices, given the stimulus
            %             and regularization parameters specified in the nim
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
            %             check_inputs(nim,Robs,Xstims,<gain_funs>,<sub_inds>)
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
            %sets the 'scale' of each subunits output (SD)
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
            %internal function that checks stim_params struct formatting,
            %and initializes default values for the given optimizer
            
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