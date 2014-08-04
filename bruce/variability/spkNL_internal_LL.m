    function f = spkNL_internal_LL(K,G,Robs)
        beta = K(2)*(G+K(1));
        too_big = find(beta > 50);
        prate = log(1+exp(beta));
        prate(too_big) = beta(too_big);
        f = -nansum(Robs.*log(prate) - prate);
    end
