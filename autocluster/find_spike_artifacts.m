function artifact_ids = find_spike_artifacts(Spikes,params)

%this code is directly copied from AllVPcs. Haven't looked at it in detail
%yet.

meanV = squeeze(mean(Spikes.V,3))';

allbid = [];
gid = 1:size(Spikes.V,1);
avar = sum(meanV(:,gid).^2);
bid = find(avar > prctile(avar,99) * 2);
while length(bid)
    allbid = [allbid bid];
    gid = setdiff(1:size(Spikes.V,1),allbid);
    avar = sum(meanV(:,gid).^2);
    bid = find(avar > prctile(avar,99) * 2);
    bid = gid(bid);
end

avar = sum(abs(diff(meanV(:,gid))));
bvar = smooth(avar,10);
bid = find(avar > prctile(avar,99) * 1.5);
while length(bid)
    allbid = [allbid bid];
    gid = setdiff(1:size(Spikes.V,1),allbid);
    avar = sum(abs(diff(meanV(:,gid))));
    bid = find(avar > prctile(avar,99) * 1.5);
    bid = gid(bid);
end
bid = allbid;

artifact_ids = bid;


