classdef NOISE_DIST
    
    %Summary of this class goes here
    %   Detailed explanation goes here
    %
    % James McFarland, September 2015
    properties
        noise_type;      % type of noise distribution (string)
    end
    
    
    methods
        function noise_dist = NOISE_DIST(noise_type)
            %constructor for NOISE_DIST class.
            assert(isstr(noise_type),'Noise type must be a string');
            noise_dist.noise_type = lower(noise_type);
            allowed_types = {'poisson','bernoulli','gaussian'}; %allowed noise distributions
            assert(ismember(noise_dist.noise_type,allowed_types),'not an allowed noise distribution');
        end
        
        function LL = get_LL(noise_dist,rPred,rObs)
            %evaluate the total LL associated with the predicted rate rPred, given observed data rObs, under the specified noise model
            switch noise_dist.noise_type
                case 'poisson' %LL = Rlog(r) - r + C
                    LL = sum(rObs .* log2(rPred) -rPred);
                    
                case 'bernoulli' %LL = R*log(r) + (1-R)*log(1-r)
                    LL = sum(rObs.*log(rPred) + (1-rObs).*log(1-rPred));
                    
                case 'gaussian' %LL = (r-R)^2 + c
                    LL = sum((rPred - rObs).^2);
            end
        end
        
        
        function LL_deriv = get_LL_deriv(noise_dist,rPred,rObs)
            %evaluate, elementwise, the derivative of the LL wrt the predicted rate at rPred, given rObs
            switch noise_dist.noise_type
                case 'poisson' %LL'[r] = R/r - 1
                    LL_deriv = rObs./rPred - 1;
                    
                case 'bernoulli' %LL'[r] = R/r - (1-R)/(1-r)
                    LL_deriv = rObs./rPred - (1-rObs)./(1-rPred);
                    
                case 'gaussian' %LL'[r] = 2*(r-R)
                    LL_deriv = 2*(rPred - rObs);
            end
            
        end
    end
end
