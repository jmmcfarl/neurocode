function [spk_id,  sig_thresh, noise_sigma] = triggerSpikes(Vsig,thresh_sign,target_Nspks,sig_thresh)

if nargin < 4
    sig_thresh = [];
end

if thresh_sign == -1 %if triggering off negative peaks
    Vsig = -Vsig;
end
if size(Vsig,2) < size(Vsig,1)
    Vsig = Vsig';
end

sgn = diff(sign(diff(Vsig,1,2)),1,2);
id = find(sgn(1,:) < 0)+1;
noise_sigma = median(abs(Vsig))/0.6745; %Quiroga 2004 Neural Comp

if isempty(sig_thresh)
    n_extrema = length(id);
    if strcmp(target_Nspks,'median')
        sig_thresh = 4*noise_sigma; %Quiroga 2004 Neural Comp
    else
        prc = target_Nspks .* 100./n_extrema; % get 1000 spikes
        sig_thresh = prctile(Vsig(id),100-prc);
    end
end
id = id(Vsig(id) > sig_thresh);

%if the ISI is very  short, and the trigger channel does not go back to
%near zero (th/3) between two trigger points, then throw away the one
%witht the smaller triggger
isicheck = [20 3];
sid = find(diff(id) < isicheck(1));
if length(sid) < length(Vsig)*0.01 %must be a low trigger
    okid = [];
    for j = 1:length(sid)
        if min(Vsig(id(sid(j)):id(sid(j)+1))) < sig_thresh/isicheck(2)
            okid = [okid sid(j)];
        end
    end
    sid = setdiff(sid,okid);
    v = cat(1,Vsig(id(sid)),Vsig(id(sid+1)));
    [a,b] = min(v);
    xid = id(sid+b-1);
    %     fprintf('Removing %d double Triggers (from %d/%d maxima)\n',length(xid),length(id),n_extrema);
    spk_id = setdiff(id,xid);
else
    fprintf(sprintf('Too many double Triggers (%d/%d)\n',length(sid),length(Vsig)));
    spk_id = id; 
end

