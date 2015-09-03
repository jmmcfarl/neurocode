function [ham_dist] = compute_state_seq_seg_hamdist(seq1,seq2)

n_segs = length(seq1);
full_seq1 = [];
full_seq2 = [];
for i = 1:n_segs
    full_seq1 = [full_seq1; seq1{i}(:)];
    full_seq2 = [full_seq2; seq2{i}(:)];
end

ham_dist = pdist([full_seq1'; full_seq2'],'hamming');