function [log_sig,offset] = log_transform_sig(sig)
 
    offset = fmincon(@(X) get_lstat(X,sig),mean(sig)/10,[],[],[],[],0.01,1000);
    log_sig = log(sig+offset);

end

function lstat = get_lstat(offset,sig)

lsig = log(sig+offset);
[~,~,lstat] = lillietest(lsig);

end