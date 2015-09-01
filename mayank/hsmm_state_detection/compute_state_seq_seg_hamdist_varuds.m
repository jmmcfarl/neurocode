function [ham_dist,coin_ups,coin_downs,fp,fn,tep] = compute_state_seq_seg_hamdist_varuds(seq1,seq2,uds1,uds2,Fs)

old_Fs = 50.4;
Fs_fac = round(Fs/old_Fs);
end_rec = round(max(max(uds1(:)),max(uds2(:)))*Fs_fac);
uds1 = floor((uds1(:,1)-1)*Fs_fac)+1;
uds2 = floor((uds2(:,1)-1)*Fs_fac)+1;
used_inds1 = zeros(end_rec,1);
used_inds2 = zeros(end_rec,1);
ov_sig1 = zeros(end_rec,1);
ov_sig2 = zeros(end_rec,1);

n_segs = length(seq1);
for i = 1:n_segs
    cur_sec = uds1(i,1):(uds1(i,1)+length(seq1{i})-1);
    used_inds1(cur_sec) = 1;
    ov_sig1(cur_sec) = seq1{i};
end
n_segs = length(seq2);
for i = 1:n_segs
    cur_sec = uds2(i,1):(uds2(i,1)+length(seq2{i})-1);
    used_inds2(cur_sec) = 1;
    ov_sig2(cur_sec) = seq2{i};
end
used_inds = find(used_inds1 & used_inds2);

ham_dist = pdist([ov_sig1(used_inds)'; ov_sig2(used_inds)'],'hamming');

upstates1 = find(ov_sig1(used_inds) == 2);
upstates2 = find(ov_sig2(used_inds) == 2);
downstates1 = find(ov_sig1(used_inds) == 1);
downstates2 = find(ov_sig2(used_inds) == 1);

coin_ups = length(intersect(upstates1,upstates2))/(0.5*(length(upstates1)+length(upstates2)));
coin_downs = length(intersect(downstates1,downstates2))/(0.5*(length(downstates1)+length(downstates2)));

fp = sum(ov_sig2(used_inds) == 2 & ov_sig1(used_inds) == 1)/sum(ov_sig1(used_inds)==1);
fn = sum(ov_sig2(used_inds) == 1 & ov_sig1(used_inds) == 2)/sum(ov_sig1(used_inds)==2);

p1 = sum(ov_sig1(used_inds)==1)/length(used_inds);
p2 = 1-p1;
tep = p1*fp + p2*fn;