function Spikes = getSpikeSnippets(V,Vtime,spk_id,spk_pts,trig_ch)
% Spikes = getSpikeSnippets(V,spk_id,spk_pts)
if nargin < 5 || isempty(trig_ch)
    trig_ch = 1;
end
allid = bsxfun(@plus,spk_id,spk_pts');
% Spikes.V = reshape(V(:,allid),[size(allid) size(V,1) size(V,3)])';
Spikes.V = reshape(V(allid',:),[length(spk_id) length(spk_pts) size(V,2)]);
Spikes.times = Vtime(spk_id);
Spikes.trig_vals = V(spk_id,trig_ch);
