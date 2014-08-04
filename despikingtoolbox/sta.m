function [ms,snips] = sta(y,spiketimes,ml,mr)
    ms = zeros(mr-ml+1,1);
    %pad y both ways
    pad = max([abs(ml),abs(mr)]);
    y = [zeros(pad,1);y;zeros(pad,1)];
    snips = zeros(mr-ml+1,length(spiketimes));
    for ii = 1:length(spiketimes)
        ms = ms + y(spiketimes(ii) + pad + (ml:mr)');
        snips(:,ii) = y(spiketimes(ii) + pad + (ml:mr)');
    end
    ms = ms/length(spiketimes);
end