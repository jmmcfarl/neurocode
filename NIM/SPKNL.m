
classdef SPKNL
    
    %Summary of this class goes here
    %   Detailed explanation goes here
    %   NOTE: the additive constant is not stored as a parameter of the
    %   spkNL function
    %
    % James McFarland, September 2015
    properties
        params;      % vector of parameters associated with the spkNL function
        NLtype;      % type of spkNL function (string)
    end
    
    methods
        function spkNL = SPKNL(NLtype, params)
            %constructor for SPKNL class.
            if nargin < 2
                params = [];
            end
            assert(isstr(NLtype),'NLtype must be a string!');
            spkNL.NLtype = lower(NLtype);
            allowed_NLs = {'lin','rectlin','exp','softplus'}; %set of NL functions currently implemented
            assert(ismember(NLtype,allowed_NLs),'invalid spkNL type!');
            
            %if using an NLtype that has parameters, check that input param
            %vector has right size, or initialize to default values
            switch NLtype
                case 'lin'
                    if isempty(params)
                        params = [1]; %defines beta in f(x; beta) = beta*x
                    else
                        assert(length(params) == 1,'invalid spkNL param vector');
                    end
                    
                case 'rectlin'
                    if isempty(params)
                        params = [1]; %defines beta in f(x; beta) = beta*x iff x > 0
                    else
                        assert(length(params) == 1,'invalid spkNL param vector');
                    end
                    
                case 'exp'
                    if isempty(params)
                        params = [1]; %defines beta in f(x; beta) = exp(beta*x)
                    else
                        assert(length(params) == 1,'invalid spkNL param vector');
                    end
                    
                case 'softplus'
                    if isempty(params)
                        params = [1 1]; %defines beta, alpha in f(x; beta, alpha) = alpha*log(1+exp(beta*x))
                    else
                        assert(length(params) == 2,'invalid spkNL param vector');
                    end
            end
            
            spkNL.params = params;
        end
        
        
        function rate = apply_spkNL(spkNL,gen_signal)
           %apply the spkNL function to the input gen_signal
           switch spkNL.NLtype
               case 'lin' %F[x;beta] = beta*x 
                   rate = gen_signal*spkNL.params(1);
                   
               case 'rectlin' %F[x; beta] = beta*x iff x > 0; else 0
                   rate = gen_signal*spkNL.params(1);
                   rate(rate < 0) = 0;
                   
               case 'exp' %F[x; beta] = exp(beta*x)
                   rate = exp(gen_signal*spkNL.params(1));
                   
               case 'softplus' %F[x; beta, alpha] = alpha*log(1+exp(beta*x))
                   max_g = 50; %to prevent numerical overflow
                   gint = gen_signal*spkNL.params(1);
                   rate = spkNL.params(2)*log(1 + exp(gint));
                   rate(gint > max_g) = spkNL.params(2)*gint(gint > max_g);
           end
        end
        
        function rate_deriv = apply_spkNL_deriv(spkNL,gen_signal)
           %apply the spkNL derivative to the input gen_signal 
           switch spkNL.NLtype
               case 'lin' %F'[x; beta] = beta;
                   rate_deriv = spkNL.params(1)*ones(size(gen_signal));
                   
               case 'rectlin' %F'[x; beta] = beta iff x > 0; else 0
                   rate_deriv = spkNL.params(1)*(gen_signal > 0);
                   
               case 'exp' %F'[x; beta] = beta*exp(beta*x)
                   rate_deriv = spkNL.params(1)*exp(spkNL.params(1)*gen_signal);
                   
               case 'softplus' %F[x; beta, alpha] = alpha*beta*exp(beta*x)/(1+exp(beta*x))
                   max_g = 50; %to prevent numerical overflow
                   gint = gen_signal*spkNL.params(1);
                   rate_deriv = spkNL.params(1)*spkNL.params(2)*exp(gint)./(1 + exp(gint));
                   rate_deriv(gint > max_g) = spkNL.params(1)*spkNL.params(2); %e^x/(1+e^x) => 1 for large x
           end
        end
        
    end
end