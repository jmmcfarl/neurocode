
classdef SUBUNIT
    % Class implementing the subunits comprising an NIM model
    %
    % James McFarland, September 2015
    %%
    properties
        filtK;       % filter coefficients, [dx1] array where d is the dimensionality of the target stimulus
        NLtype;      % upstream nonlinearity type (string)
        NLparams;    % vector of parameters associated with the upstream NL function (for parametric functions)
        weight;      % subunit weight (typically +/- 1)
        Xtarg;       % index of stimulus the subunit filter acts on
        reg_lambdas; % struct of regularization hyperparameters
        Ksign_con;   %scalar defining any constraints on the filter coefs [-1 is negative con; +1 is positive con; 0 is no con]
        scale;       %SD of the subunit output derived from most-recent fit
        TBparams;    %struct of parameters associated with a 'nonparametric' NL
        TBy;         %tent-basis coefficients
        TBx;         %tent-basis center positions
    end
    properties (Hidden)
        TBy_deriv;   %internally stored derivative of tent-basis NL
    end
    
    %%
    methods
        function subunit = SUBUNIT(init_filt, weight, NLtype, Xtarg, NLparams,Ksign_con)
%         subunit = SUBUNIT(init_filt, weight, NLtype, <NLparams>, <Ksign_con>)
%         constructor for SUBUNIT class.
%             INPUTS:
%                 init_filt: vector of initial filter coefs
%                 weight: scalar weighting associated with subunit (typically +/-1)
%                 NLtype: string specifying the type of upstream NL
%                 Xtarg: scalar index specifying which stimulus element this subunit acts on
%                 <NLparams>: vector of parameters for the upstream NL
%                 <Ksign_con>: constraint on filter coefs [-1 = negative; +1 = positive; 0 is no con]
%             OUTPUTS: subunit: subunit object

            if (nargin < 4 || isempty(Xtarg)); Xtarg = 1; end %default Xtarget is 1
            if nargin < 5; NLparams = []; end;
            if (nargin < 6 || isempty(Ksign_con)); Ksign_con = 0; end; %default no constraints on filter coefs
            
            assert(length(weight) == 1,'weight must be scalar!');
            assert(ischar(NLtype),'NLtype must be a string');
            subunit.filtK = init_filt;
            subunit.weight = weight;
            if ~ismember(weight,[-1 1])
                warning('Best to initialize subunit weights to be +/- 1');
            end
            subunit.Xtarg = Xtarg;
            subunit.NLtype = lower(NLtype);
            allowed_NLs = {'lin','quad','rectlin','rectquad','softplus'}; %set of NL functions currently implemented
            assert(ismember(subunit.NLtype,allowed_NLs),'invalid NLtype!');
            
            %if using an NLtype that has parameters, check that input
            %parameter vector is the right size, or initialize to default
            %values
            switch subunit.NLtype
                case {'rectlin','rectquad'}
                    if isempty(NLparams) %if parameters are not specified
                        NLparams = [0]; %defines threshold in f(x) = (x-c) iff x >= c
                    else
                        assert(length(NLparams) == 1,'invalid NLparams vector');
                    end
                case 'softplus'
                    if isempty(NLparams)
                        NLparams = [1 0]; %defines beta and c in f(x) = log(1 + exp(beta*x + c))
                    else
                        assert(length(NLparams) == 2,'invalid NLparams vector');
                    end
            end
            subunit.NLparams = NLparams;            
            subunit.reg_lambdas = SUBUNIT.init_reg_lamdas();
            subunit.Ksign_con = Ksign_con;
        end
        
        function filtK = get_filtK(subunit)
            %get vector of filter coefs from the subunit
            filtK = subunit.filtK;
        end
        
        %to ADD: compute subunit output given the stimulus
        %         function sub_out = get_sub_out(subunit, Xstim)
        %
        %         end
        
        function sub_out = apply_NL(subunit,gen_signal)
            %apply subunit NL to the input generating signal
            switch subunit.NLtype                
                case 'lin' %f(x) = x
                    sub_out = gen_signal;
                    
                case 'quad' %f(x) = x^2
                    sub_out = gen_signal.^2;
                    
                case 'rectlin' %f(x;c) = (x-c) iff x >= c; else x = 0
                    sub_out = (gen_signal - subunit.NLparams(1));
                    sub_out(sub_out < subunit.NLparams(1)) = 0;
                    
                case 'rectquad' %f(x;c) = (x-c)^2 iff x >= c; else x = 0
                    sub_out = (gen_signal - subunit.NLparams(1)).^2;
                    sub_out(gen_signal < subunit.NLparams(1)) = 0;
                    
                case 'softplus' %f(x;beta, c) = log(1 + exp(beta*x + c))
                    max_g = 50; %to prevent numerical overflow
                    gint = gen_signal*subunit.NLparams(1) + subunit.NLparams(2); %gen_signal*beta + c
                    expg = exp(gint);
                    sub_out = log(1 + expg);
                    sub_out(gint > max_g) = gint(gint > max_g);
                    
                case 'nonpar' %f(x) piecewise constant with knot points TBx and coefficients TBy
                    sub_out = zeros(size(gen_signal));
                    %Data where X < TBx(1) are determined entirely by the first tent basis
                    left_edge = find(gen_signal < subunit.TBx(1));
                    sub_out(left_edge) = sub_out(left_edge) + subunit.TBy(1);
                    %similarly for the right edge
                    right_edge = find(gen_signal >= subunit.TBx(end));
                    sub_out(right_edge) = sub_out(right_edge) + subunit.TBy(end);
                    slopes = diff(subunit.TBy)./diff(subunit.TBx);
                    for j = 1:length(subunit.TBy)-1
                        cur_set = find(gen_signal >= subunit.TBx(j) & gen_signal < subunit.TBx(j+1));
                        sub_out(cur_set) = sub_out(cur_set) + subunit.TBy(j) + slopes(j)*(gen_signal(cur_set) - subunit.TBx(j));
                    end
            end
        end
        
        
        function sub_deriv = apply_NL_deriv(subunit,gen_signal)
            %apply derivative of subunit NL to input gen_signal
            switch subunit.NLtype
                
                case 'lin' %f'(x) = 1 (shouldnt ever need to use this within the optimization...)
                    sub_deriv = ones(size(gen_signal));
                    
                case 'quad' %f'(x) = 2x
                    sub_deriv = 2*gen_signal;
                    
                case 'rectlin' %f'(x) = 1 iff x >= c; else 0
                    sub_deriv = gen_signal >= subunit.NLparams(1);
                    
                case 'rectquad' %f'(x) = 2x iff x >= c; else 0
                    sub_deriv = 2*gen_signal;
                    sub_deriv(gen_signal < subunit.NLparams(1)) = 0;
                    
                case 'softplus' %f'(x) = beta*exp(beta*x + c)/(1 + exp(beta*x + c))
                    max_g = 50; %to prevent numerical overflow
                    gint = gen_signal*subunit.NLparams(1) + subunit.NLparams(2); %gen_signal*beta + c
                    sub_deriv = subunit.NLparams(1)*exp(gint)./(1 + exp(gint));
                    sub_deriv(gint > max_g) = 1; %for large gint, derivative goes to 1
                    
                case 'nonpar'
                    ypts = subunit.TBy_deriv; xedges = subunit.TBx;
                    sub_deriv = zeros(length(gen_signal),1);
                    for n = 1:length(xedges)-1
                        sub_deriv((gen_signal >= xedges(n)) & (gen_signal < xedges(n+1))) = ypts(n);
                    end
            end
        end
    end
    
    methods (Hidden)
        
        function fprime = get_TB_derivative(subunit)
            %calculate the derivative of the piecewise linear function wrt x
            fprime = zeros(1,length(subunit.TBx)-1);
            for n = 1:length(fprime)
                fprime(n) = (subunit.TBy(n+1)-subunit.TBy(n))/(subunit.TBx(n+1)-subunit.TBx(n));
            end
        end
        
        function gout = tb_rep(subunit,gin)
            %project the input signal gin onto the tent-basis functions associated with this subunit
            n_tbs =length(subunit.TBx); %number of tent-basis functions
            gout = zeros(length(gin),n_tbs);
            gout(:,1) = get_tentbasis_output(gin,subunit.TBx(1),[-Inf subunit.TBx(2)]);
            gout(:,end) = get_tentbasis_output(gin,subunit.TBx(end),[subunit.TBx(end-1) Inf]);
            for n = 2:n_tbs-1
                gout(:,n) = get_tentbasis_output(gin, subunit.TBx(n), [subunit.TBx(n-1) subunit.TBx(n+1)] );
            end
        end
    end
    
    methods (Static)
        function reg_lambdas = init_reg_lamdas()
            %creates reg_lambdas struct and sets default values to 0
            reg_lambdas.nld2 = 0; %second derivative of tent basis coefs
            reg_lambdas.d2xt = 0; %spatiotemporal laplacian
            reg_lambdas.d2x = 0; %2nd spatial deriv
            reg_lambdas.d2t = 0; %2nd temporal deriv
            reg_lambdas.l2 = 0; %L2 on filter coefs
            reg_lambdas.l1 = 0; %L1 on filter coefs
        end
    end
end

