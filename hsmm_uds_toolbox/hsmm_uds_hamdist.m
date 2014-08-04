function [ham_dist] = hsmm_uds_hamdist(seq1,seq2,hmm1,hmm2,Fs,T)

new_seg_inds1 = resample_uds_seg_inds(hmm1.UDS_segs,hmm1.Fs,Fs,T);
new_seg_inds2 = resample_uds_seg_inds(hmm2.UDS_segs,hmm2.Fs,Fs,T);

used_inds1 = zeros(T,1);
used_inds2 = zeros(T,1);
ov_sig1 = zeros(T,1);
ov_sig2 = zeros(T,1);

n_segs = length(seq1);
for i = 1:n_segs
    cur_sec = new_seg_inds1(i,1):(new_seg_inds1(i,1)+length(seq1{i})-1);
    used_inds1(cur_sec) = 1;
    ov_sig1(cur_sec) = seq1{i};
end
n_segs = length(seq2);
for i = 1:n_segs
    cur_sec = new_seg_inds2(i,1):(new_seg_inds2(i,1)+length(seq2{i})-1);
    used_inds2(cur_sec) = 1;
    ov_sig2(cur_sec) = seq2{i};
end
used_inds = find(used_inds1 & used_inds2);

ham_dist = pdist([ov_sig1(used_inds)'; ov_sig2(used_inds)'],'hamming');