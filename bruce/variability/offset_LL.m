    function f = offset_LL(x,prate,Robs)
        cur_rate = prate;
        cur_rate(cur_rate < x) = x;
        f = -nansum(Robs.*log(cur_rate) - cur_rate);
    end
