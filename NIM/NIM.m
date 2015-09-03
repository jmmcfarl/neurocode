
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
%                   mod_signs: vector specifying the weight associated with each subunit (typically +/- 1)
%                   <Xtargets>: vector specifying the index of the stimulus each subunit acts on (defaults to ones)
%                   <spkNL>: string specifying type of spkNL function
%                   <noise_dist>: string specifying type of noise distribution
            
            nStims = length(stim_params); %number of stimuli
            stim_params = NIM.check_stim_params(stim_params); %validate and format input stim_params
            nim.stim_params = stim_params;
            
            nSubs = length(mod_signs); %number of subunits
            %if NLtypes is specified as a single string, default to using
            %this NLtype for all subunits
            if ~iscell(NLtypes) && isstr(NLtypes)
                NLtypes = repmat({NLtypes},nSubs,1);
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
            
            %initialize subunits
            nim.init_props = rng(); %save state of RNG used for initializing filter weights
            for ii = 1:nSubs %loop initializing subunits (start from last to initialize object array)
                stimD = prod(nim.stim_params(Xtargets(ii)).dims); %dimensionality of the current filter
                init_filt = randn(stimD,1)/stimD; %initialize fitler coefs with gaussian noise
                nim.subunits = cat(1,nim.subunits,SUBUNIT(init_filt, mod_signs(ii), NLtypes{ii},Xtargets(ii)));
            end
        end
        
        %% setting methods
        function nim = set_reg_params(nim, varargin)
            %             set a desired set of regularization parameters to specified values, apply to specified set of subunits
            %             optional flags:
            %                 'target_subs': set of subunits to apply the new reg_params for
            %                 'lambda_type': type of regularization
            
            target_subs = 1:length(nim.subunits); %default is to apply the change to all subunits
            
            %INPUT PARSING
            j = 1;
            reg_types = {}; reg_vals = [];
            while j <= length(varargin)
                flag_name = varargin{j};
                flag_val = varargin{j+1};
                switch lower(flag_name)
                    case 'target_subs'
                        target_subs = flag_val;
                        assert(all(ismember(target_subs,1:length(nim.subunits))),'invalid target subunits specified');
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
                for jj = 1:length(target_subs)
                    nim.subunits(target_subs(jj)).reg_lambdas = setfield(nim.subunits(target_subs(jj)).reg_lambdas,reg_types{ii},reg_vals(ii));
                end
            end
        end
        
        function nim = set_stim_params(nim, varargin)
            %             adjust stimulus parameters for the desired Xtarg
            %             optional flags:
            %                 'xtarg': index of stimulus to apply the new stim_params
            %                   for [default is 1]
            %                 'dims': dimensionality of stim: [Tdim, X1dim, X2dim] where Tdim is the number of temporal dimensions, etc.
            %                 'boundary_conds': boundary conditions on each stim dimension [Inf means free, 0 means tied, -1 means periodic]
            %                 'split_pts': set of stimulus indices to specify smoothing boundaries (NOT YET IMPLEMENTED)
            
            Xtarg = 1; %default is to apply the change to stim 1
            allowed_flags = {'dims','boundary_conds','split_pts'};
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
        end
        
        %% getting methods
        
        function lambdas = get_reg_lambdas(nim, varargin)
           %get regularizatoin labmda values from nim subunits
%            INPUT FLAGS:
%                 target_subs, followed by vector of indices, specifies which subunits to extract lambdas from
%                 lambda_type: string of allowed regularization types
%            OUTPUTS: 
%                 lambdas: [K,N] matrix of lambda values, K is the number of specified lambda_types and N is the number of subunits
           
           target_subs = 1:length(nim.subunits); %default is to grab reg values from all subunits
           %INPUT PARSING
           jj = 1;
           reg_types = {}; 
           while jj <= length(varargin)
               flag_name = varargin{jj};
               switch lower(flag_name)
                   case 'target_subs'
                       target_subs = varargin{jj+1};
                       assert(all(ismember(target_subs,1:length(nim.subunits))),'invalid target subunits specified');
                       jj = jj + 2;
                   case nim.allowed_reg_types
                       reg_types = cat(1,reg_types,lower(flag_name));
                       jj = jj + 1;
                   otherwise
                       error('Invalid input flag');
               end
           end
           
           lambdas = nan(length(reg_types),length(target_subs));
           if isempty(reg_types)
               warning('No regularization type specified, returning nothing');
           end
           for ii = 1:length(reg_types)
               for jj = 1:length(target_subs)
                   lambdas(ii,jj) = getfield(nim.subunits(target_subs(jj)).reg_lambdas,reg_types{ii});
               end
           end
        end
        
        %% display methods
        
        function h = display_filters(nim,Xtarg)
            %plot filters: [description goes here]
            h = figure();
            Xtargs = [nim.subunits(:).Xtarg];
            stim_mods = find(Xtargs == Xtarg);
            n_columns = max(round(sqrt(length(stim_mods)/2)),1);
            n_rows = ceil(length(stim_mods)/n_columns);
            
            nLags = nim.stim_params(Xtarg).dims(1);
            nPix = nim.stim_params(Xtarg).dims(2:3);
            
            pix_ax = 1:nPix;
            lag_ax = 1:nLags;
            for imod = 1:length(stim_mods)
                thismod = nim.subunits(stim_mods(imod));
                
                %PLOT FILTER
                subplot(n_rows,n_columns,(imod-1)+1);
                imagesc(pix_ax,lag_ax,reshape(thismod.filtK,nLags,nPix(1)));
                cl = max(abs(thismod.filtK));
                caxis([-cl cl]);
                %colormap(jet);
                colormap(gray);
                set(gca,'ydir','normal');
                xlabel('Pixels')
                ylabel('Time lags');
%                 if ~isempty(xrange)
%                     xlim(xrange);
%                 end
%                 if ~isempty(trange)
%                     ylim(trange);
%                 end
                NLtype = 'NP';
                if strcmp(thismod.NLtype,'lin')
                    NLtype = 'Lin';
                elseif strcmp(thismod.NLtype,'quad')
                    NLtype = 'Quad';
                end
                if thismod.weight == 1
                    NLsign = 'E';
                else
                    NLsign = 'I';
                end
                title(sprintf('%s %s-filt',NLtype,NLsign));
            end
            fig_props.dims = [n_rows n_columns];
            fig_props.nmods = length(stim_mods);
            
        end
        
        %% fitting methods
        
        function nim = fit_filters(nim, rObs, Xstims)
            %
            % Usage: nim_out = NMMfit_filters( nim, Robs, Xstims, <Gmults>, <Uindx>, <silent>, <desired_optim_params>, <regmat_custom>,<targets> )
            %
            % INPUTS:
            % OUTPUTS:
            
            % PROCESS INPUTS
            Nsubs = length(nim.subunits);
            targets = 1:Nsubs;
            % Index X-matrices and Robs
            %             RobsFULL = Robs;
            %             if ~isempty(Uindx)
            %                 for nn = 1:length(Xstims)
            %                     Xstims{nn} = Xstims{nn}(Uindx,:);
            %                 end
            %                 Robs = RobsFULL(Uindx);
            %             end
            %             if isempty(targets) %default is to optimize all model components
            %                 targets = 1:Nmods;
            %                 if spkhstlen > 0
            %                     targets = [targets -1]; %optimize spike hist filter
            %                 end
            %             elseif targets == -2
            %                 targets = [];
            %             end
            
            %             % Make sure Robs is a column vector
            %             Robs = Robs(:);
            %             spkhstlen = nim.spk_hist.spkhstlen;
            
            %INPUT CHECKING
            for n = 1:Nsubs %check that stimulus dimensions match
                [NT,filtLen] = size(Xstims{nim.subunits(n).Xtarg}); %stimulus dimensions
                assert(filtLen == prod(nim.stim_params(nim.subunits(n).Xtarg).dims),'Xstim dims dont match stim_params');
            end
            assert(all(ismember(targets,[1:Nsubs -1])),'Invalid target specified');
            if ismember(-1,targets); assert(spkhstlen > 0,'no spike history term initialized!'); end;
            
            Ntargets = sum(targets > 0); %number of targeted subunits
            non_targets = setdiff([1:Nsubs -1],targets); %elements of the model held constant
            
            % PARSE INITIAL PARAMETERS
            init_params = [];
            %             sign_con = [];
            for imod = targets(targets > 0)
                cur_kern = nim.subunits(imod).filtK;
                %                 if isfield(nim.mods(imod),'Kcon')
                %                     if nim.mods(imod).Kcon ~= 0
                %                         sign_con(length(initial_params)+(1:length(cur_kern))) = nim.mods(imod).Kcon;
                %                     end
                %                 end
                init_params = [init_params; cur_kern]; % add coefs to initial param vector
            end
            
            %             % Add in spike history coefs
            %             if ismember(-1,targets)
            %                 init_params = [init_params; nim.spk_hist.coefs];
            %             end
            
            init_params(end+1) = nim.spkNL.theta; % add constant offset
            
            %             % COMPUTE L1 PENALTY IF APPLICABLE
            %             lambda_L1 = zeros(size(initial_params));
            %             cnt = 0;
            %             for ii = 1:Ntargets
            %                 filtLen = length(nim.mods(targets(ii)).filtK);
            %                 cur_inds = (1:filtLen) + cnt;
            %                 lambda_L1(cur_inds) = nim.mods(targets(ii)).reg_params.lambda_L1;
            %                 cnt = cnt + filtLen;
            %             end
            %             lambda_L1 = lambda_L1/sum(Robs); % since we are dealing with LL/spk
            
            %             % PRECOMPUTE 'TENT-BASIS' DERIVATIVES OF UPSTREAM NLS IF NEEDED
            %             if any(strcmp('nonpar',{nim.mods(targets(targets > 0)).NLtype}))
            %                 for ii = 1:Ntargets
            %                     if strcmp(nim.mods(targets(ii)).NLtype,'nonpar')
            %                         NLx = nim.mods(targets(ii)).NLx;
            %                         NL = nim.mods(targets(ii)).NLy;
            %
            %                         % Compute derivative of non-linearity
            %                         fpr = zeros(1,length(NLx)-1);
            %                         for n = 1:length(fpr)
            %                             fpr(n) = (NL(n+1)-NL(n))/(NLx(n+1)-NLx(n));
            %                         end
            %                         fprimes{ii} = fpr;
            %                     else
            %                         fprimes{ii} = [];
            %                     end
            %                 end
            %             else
            %                 fprimes = [];
            %             end
            
            %             % CREATE SPIKE HISTORY Xmat IF NEEDED
            %             if nim.spk_hist.spkhstlen > 0
            %                 Xspkhst = create_spkhist_Xmat( RobsFULL, nim.spk_hist.bin_edges );
            %                 if ~isempty(Uindx)
            %                     Xspkhst = Xspkhst(Uindx,:);
            %                 end
            %             else
            %                 Xspkhst = [];
            %             end
            
            % COMPUTE NET OUPTUT OF ALL NON-TARGET PREDICTORS
            nontarg_g = zeros(NT,1);
            for imod = non_targets(non_targets > 0) %for all subunits that aren't targeted
                fgint = Xstims{nim.subunits(imod).Xtarg} * nim.subunits(imod).get_filtK; %apply filter to stim
                fgint = nim.subunits(imod).apply_NL(fgint); %apply upstream NL
                %                 if isempty(Gmults{imod})
                nontarg_g = nontarg_g + fgint*nim.subunits(imod).weight; %add to net non-target output
                %                 else
                %                     nontarg_g = nontarg_g + (fgint.*Gmults{imod}) * nim.mods(imod).sign;
                %                 end
            end
            
            %             if ismember(-1,non_targets) && spkhstlen > 0
            %                 nontarg_g = nontarg_g + Xspkhst*nim.spk_hist.coefs(:);
            %             end
            
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
            %             if spkhstlen > 0 && ismember(-1,targets) %if optimizing spk history term
            %                 %negative constraint on spk history coefs
            %                 if nim.spk_hist.negCon == 1
            %                     %spkhist_inds = (Ntargets*filtLen + 1):(Ntargets*filtLen + spkhstlen);
            %                     spkhist_inds = cnt + (1:spkhstlen);
            %                     UB(spkhist_inds) = 0;
            %                     use_con = 1;
            %                 end
            %             end
            
            beq = zeros(size(Aeq,1),1);
            b = zeros(size(A,1),1);
            
            %             % GENERATE REGULARIZATION MATRICES
            Tmats = nim.make_Tikhonov_matrices();
            
            %             if max(lambda_L1) > 0 && use_con == 1
            %                 disp('Cant use L1 with constrained optimization, aborting constraints');
            %                 use_con = 0;
            %             end
            
            opt_fun = @(K) nim.internal_LL_grad_filters(K,rObs,Xstims,targets,nontarg_g,Tmats);
            
            %
            optim_params.MaxFunEvals = 100*length(init_params);
            optim_params.MaxIter = 1e3;
            optim_params.Display = 'off';
            %             if silent == 0
            optim_params.Display = 'iter';
            %             end
            %             if use_con == 0 %if no constraints
            %
            %                 %if using L1 reg
            %                 if max(lambda_L1) > 0
            %                     if exist('L1General2_PSSas','file') == 2
            %                         optim_params.optTol = 1e-4;
            %                         optim_params.progTol = 1e-8;
            %                         if silent == 0
            %                             optim_params.verbose = 2;
            %                         else
            %                             optim_params.verbose = 0;
            %                         end
            %                         % Load in specified optimization parameters
            %                         if ~isempty(desired_optim_params)
            %                             spec_fields = fieldnames(desired_optim_params);
            %                             for i = 1:length(spec_fields)
            %                                 optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
            %                             end
            %                         end
            %
            %                         [params] = L1General2_PSSas(@(K) LLfit_filters_internal(nim, K, Robs, Xstims,Xspkhst,Gmults,L2_mats,targets,nt_gout,fprimes),...
            %                             initial_params,lambda_L1,optim_params);
            %                         %             [params] = L1General2_PSSas(@(K) NIM_fit_filters_internal(nim, K, Robs, Xstim,Xspkhst,XLin,L2_mats,targets,nt_gout,fprimes),...
            %                         %                 initial_params,lambda_L1);
            %                     else
            %                         error('Need to install Mark Schmidts L1 optimization toolbox for using L1');
            %                     end
            %                 else % if not using L1 reg
            
            if exist('minFunc','file') == 2 %try to use Mark Schmidt's optimizer
                % if using Mark Schmidt's optimization, some differences in option parameters
                optim_params.optTol = 1e-4;
                optim_params.progTol = 1e-8;
                optim_params.Method = 'lbfgs';
%                 if silent == 0
                    optim_params.verbose = 2;
%                 else
%                     optim_params.verbose = 0;
%                 end
%                 %load in specified optimization parameters
%                 if ~isempty(desired_optim_params)
%                     spec_fields = fieldnames(desired_optim_params);
%                     for i = 1:length(spec_fields)
%                         optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
%                     end
%                 end
                
                [params] = minFunc(opt_fun, init_params, optim_params);
                
            else %if using Matlab Optim toolbox:
                
                % Default optimization parameters
                optim_params.LargeScale = 'off';
                optim_params.TolFun = 1e-6;
                optim_params.TolX = 1e-6;
                optim_params.HessUpdate = 'bfgs';
                optim_params.GradObj = 'on';
                %
                %                         %load in specified optimization parameters
                %                         if ~isempty(desired_optim_params)
                %                             spec_fields = fieldnames(desired_optim_params);
                %                             for i = 1:length(spec_fields)
                %                                 optim_params = setfield(optim_params,spec_fields{i},getfield(desired_optim_params,spec_fields{i}));
                %                             end
                %                         end
                
                [params] = fminunc(opt_fun, init_params, optim_params);
                
            end
            %                 end
            %             else %if there are constraints
            
            %                 % Try to use Mark Schmidt constrained optimizer
            %                 if exist('minConf_TMP','file') == 2 && isempty(A) && isempty(Aeq)
            %                     % if using Mark Schmidt's optimization, some differences in option parameters
            %                     optim_params.optTol = 1e-4;
            %                     optim_params.progTol = 1e-6;
            %                     if silent == 0
            %                         optim_params.verbose = 2;
            %                     else
            %                         optim_params.verbose = 0;
            %                     end
            %                     [params] = minConf_TMP( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes ),...
            %                         initial_params, LB,UB,optim_params);
            %                 else
            %                     % otherwise resort to matlab's
            %                     optim_params.GradObj = 'on';
            %                     optim_params.LargeScale = 'off';
            %                     optim_params.Algorithm = 'active-set';
            %                     optim_params.optTol = 1e-4;
            %                     optim_params.progTol = 1e-6;
            %                     [params] = fmincon( @(K) LLfit_filters_internal( nim, K, Robs, Xstims, Xspkhst, Gmults, L2_mats, targets, nt_gout, fprimes ),...
            %                         initial_params, A,b,Aeq,beq,LB,UB,[],optim_params);
            %                 end
            
            %         end
            
            % PARSE MODEL FIT
            nim.spkNL.theta = params(end); %set new offset parameter
            %         if ismember(-1,targets)
            %             nim_out.spk_hist.coefs = params(cnt + (1:spkhstlen));
            %         end
            kOffset = 0; %position counter for indexing param vector
            for ii = 1:Ntargets
                filtLen = length(nim.subunits(targets(ii)).filtK);
                cur_kern = params((1:filtLen) + kOffset); %grab parameters corresponding to this subunit's filters
                nim.subunits(targets(ii)).filtK = cur_kern(:); %assign new filter values
                kOffset = kOffset + filtLen;
            end
            
            %         [LL, nullLL, ~, G, gint, fgint, penLL] = NMMeval_model( nim_out, Robs, Xstims, Gmults, [], regmat_custom );
            %         nim_out.LL_seq = cat(1,nim_out.LL_seq,LL);
            %         nim_out.penLL_seq = cat(1,nim_out.penLL_seq,penLL);
            %         nim_out.opt_history = cat(1,nim_out.opt_history,{'filt'});
        end
    end
    
    methods (Hidden)
        %% internal methods
        
        function [penLL, penLLgrad] = internal_LL_grad_filters(nim,params,rObs,Xstims,targets,nontarg_g, Tmats)
            %computes the penalized LL and its gradient wrt the filters for the given nim
            %with parameter vector params
            
            if ismember(-1,targets); opt_spk_NL = true; else opt_spk_NL = false; end; %are we optimizing the spk NL term?
            targets = targets(targets > 0);
            Ntargets = length(targets); %number of targeted subs
            
            % USEFUL VALUES
            theta = params(end); % offset
            gint = nan(length(rObs),Ntargets); %initialize matrix for storing filter outputs
            filtLen = zeros(Ntargets,1); %store the length of each (target) sub's filter
            filtKs = cell(Ntargets,1); %store the filter coefs for all (target) subs)
            param_inds = cell(Ntargets,1); %this will store the index values of each subunit's filter coefs within the parameter vector
            Xtarg_set = [nim.subunits(targets).Xtarg]; %vector of Xtargets for set of subunits being optimized
            un_Xtargs = unique(Xtarg_set); %set of unique Xtargets
            mod_NL_types = {nim.subunits(targets).NLtype}; %NL types for each targeted subunit
            unique_NL_types = unique(mod_NL_types); %unique set of NL types being used
            mod_weights = [nim.subunits(targets).weight]'; %signs of targeted subunits
            
            G = theta + nontarg_g; % initialize overall generating function G with the offset term and the contribution from nontarget subs
            
            NKtot = 0;  %init filter coef counter
            for ii = 1:Ntargets %loop over subunits, get filter coefs and their indices within the parameter vector
                filtLen(ii) = length(nim.subunits(targets(ii)).filtK); % length of filter
                param_inds{ii} = NKtot + (1:filtLen(ii)); %set of param indices associated with this subunit's filters
                filtKs{ii} = params(param_inds{ii}); %store filter coefs
                NKtot = NKtot + filtLen(ii); %inc counter
            end
            for ii = 1:length(un_Xtargs) %loop over the unique Xtargs and compute the generating signals for all relevant filters
                cur_mods = find(Xtarg_set == un_Xtargs(ii)); %set of targeted subunits that act on this Xtarg
                gint(:,cur_mods) = Xstims{un_Xtargs(ii)} * cat(2,filtKs{cur_mods}); %apply filters to stimulus
            end
            
            fgint = gint; %init subunit outputs by filter outputs
            for ii = 1:length(unique_NL_types) %loop over unique subunit NL types and apply NLs to gint in batch
                if ~strcmp(unique_NL_types{ii},'lin') %if it's not just a linear sub, we have to do something
                    cur_mods = find(strcmp(mod_NL_types,unique_NL_types{ii})); %set of subs with this NL type
                    fgint(:,cur_mods) = nim.subunits(targets(cur_mods(1))).apply_NL(gint(:,cur_mods)); %apply upstream NL to all subunits of this type
                end
            end
            
            % Multiply by weight (and multiplier, if appl) and add to generating function
            %             if isempty(gain_funs)
            G = G + fgint*mod_weights;
            %             else
            %                 error('not implemented yet')
            %             end
            
            %             % Add contribution from spike history filter
            %             if spkhstlen > 0 && ismember(-1,targets)
            %                 G = G + Xspkhst*params(NKtot + (1:spkhstlen));
            %             end
            
            pred_rate = nim.apply_spkNL(G);
            % Enforce minimum predicted firing rate to avoid nan LLs
            if ismember(nim.noise_dist,{'poisson','bernoulli'})
                min_pred_rate = 1e-50;
                if min(pred_rate) < min_pred_rate
                    pred_rate(pred_rate < min_pred_rate) = min_pred_rate; %minimum predicted rate
                end
            end
            
            penLL = nim.internal_LL(pred_rate,rObs); %compute LL
            
            %residual = LL'[r].*F'[g]
            residual = nim.internal_LL_deriv(pred_rate,rObs) .* nim.apply_spkNL_deriv(G);
            
            penLLgrad = zeros(length(params),1); %initialize LL gradient
            penLLgrad(end) = sum(residual);      %Calculate derivatives with respect to constant term (theta)
            
            %             % Calculate derivative with respect to spk history filter
            %             if spkhstlen > 0 && ismember(-1,targets)
            %                 LLgrad(NKtot+(1:spkhstlen)) = residual'*Xspkhst;
            %             end
            
            for ii = 1:un_Xtargs %loop over unique Xtargets and compute LL grad wrt stim filters
                cur_sub_inds = find(Xtarg_set == un_Xtargs(ii)); %set of subunits with this Xtarget
                cur_NL_types = mod_NL_types(cur_sub_inds); %NL types of current subs
                cur_unique_NL_types = unique(cur_NL_types); %set of unique NL types
                
                if length(cur_sub_inds) == 1 && strcmp(cur_unique_NL_types,'lin') %if there's only a single linear subunit, this is a faster calc
                    penLLgrad(param_inds{cur_sub_inds}) = residual'*Xstims{un_Xtargs(ii)} * nim.subunits(cur_sub_inds).weight;
                else %otherwise, compute a matrix of upstream NL derivatives fpg
                    fpg = ones(length(residual),length(cur_sub_inds)); %initialize to linear NL derivative (all ones)
                    for jj = 1:length(cur_unique_NL_types) %loop over unique NL types
                        cur_sub_subinds = find(strcmp(cur_NL_types,cur_unique_NL_types{jj})); %indices of current subset of subunits
                        fpg(:,cur_sub_subinds) = nim.subunits(targets(cur_sub_inds(cur_sub_subinds(1)))).apply_NL_deriv(gint(:,cur_sub_inds(cur_sub_subinds)));
                    end
                    target_params = cat(2,param_inds{cur_sub_inds}); %indices of filter coefs for current set of targeted subunits
                    %LL grad is residual * f'(.) *X *w, computed in parallel for all subunits targeting this Xtarg
                    penLLgrad(target_params) = bsxfun(@times,(bsxfun(@times,fpg,residual)'*Xstims{un_Xtargs(ii)}),mod_weights(cur_sub_inds))';
                end
            end
            
            net_penalties = zeros(size(targets));
            net_pen_grads = zeros(length(params),1);
            for ii = 1:length(Tmats) %loop over the derivative regularization matrices
                cur_subs = find([nim.subunits(targets).Xtarg] == Tmats(ii).Xtarg); %set of subunits acting on the stimulus given by this Tmat 
                penalties = sum((Tmats(ii).Tmat * cat(2,filtKs{cur_subs})).^2);
                pen_grads = 2*(Tmats(ii).Tmat' * Tmats(ii).Tmat * cat(2,filtKs{cur_subs}));
                cur_lambdas = nim.get_reg_lambdas(Tmats(ii).type,'target_subs',targets(cur_subs)); %current lambdas
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
            Nspks = sum(rObs);
            penLL = -penLL/Nspks;
            penLLgrad = -penLLgrad/Nspks;
        
        end
        
        
        function LL = internal_LL(nim,rPred,rObs)
            %internal evaluatation method for computing the total LL associated with the predicted rate rPred, given observed data rObs
            switch nim.noise_dist
                case 'poisson' %LL = Rlog(r) - r + C
                    LL = sum(rObs .* log(rPred) -rPred);
                    
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
                cur_mods = find(Xtargs == ii); %get set of subunits acting on this stimuls
                for jj = 1:length(deriv_reg_types) %check each possible derivative regularization type
                    cur_lambdas = nim.get_reg_lambdas(deriv_reg_types{jj},'target_subs',cur_mods);
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
        
    end
    
    methods (Static, Hidden)
        function stim_params = check_stim_params(stim_params)
            %           internal function that checks stim_params struct formatting, and initializes default values if needed
            
            default_split_pts = []; %no discontinuities in smoothing penalty
            default_boundary_conds = [Inf 0 0]; %[free boundary on first dim, and tied to 0 in other dims]
            for ii = 1:length(stim_params)
                assert(isfield(stim_params(ii),'dims'),'need to specify the stimulus dimensions');
                
                %assign defaults if they aren't specified
                if ~isfield(stim_params(ii),'boundary_conds') || isempty(stim_params(ii).boundary_conds)
                    stim_params(ii).boundary_conds = default_boundary_conds;
                end
                if ~isfield(stim_params(ii),'split_pts')
                    stim_params(ii).split_pts = default_split_pts;
                end
                
                %make sure .dims field a row vector
                if size(stim_params(ii).dims,1) > size(stim_params(ii).dims,2)
                    stim_params(ii).dims = stim_params(ii).dims';
                end
                %if the length of .dims is less than 3, pad with 1's
                while length(stim_params(ii).dims) < 3
                    stim_params(ii).dims = [stim_params(ii).dims 1];
                end
            end
        end
    end
end