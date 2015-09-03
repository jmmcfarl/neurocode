function [state_seq] = get_fixedthresh_state_seq(hmm,emiss)

UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = hmm.UDS_segs(i,1):hmm.UDS_segs(i,2);
    UDS_samples = [UDS_samples hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)];
end

%set all instances where the observation is greater than threshold to the
%up state
for i = 1:hmm.Nsegs
    
    cur_data = emiss(UDS_seg_inds{i});
    state_seq{i} = ones(size(cur_data));
    state_seq{i}(cur_data > hmm.threshold) = 2;
    
end
