function [comb_state_seq] = combined_lf_hf_state_seq(lf_state_seq,hf_state_seq,Fs)

max_diff = round(Fs*0.3);

up_trans_lf = find(lf_state_seq(1:end-1) == 1 & lf_state_seq(2:end)==2);
down_trans_lf = find(lf_state_seq(1:end-1)==2 & lf_state_seq(2:end)==1);
up_trans_hf = find(hf_state_seq(1:end-1) == 1 & hf_state_seq(2:end)==2);
down_trans_hf = find(hf_state_seq(1:end-1)==2 & hf_state_seq(2:end)==1);
up_trans_lf(up_trans_lf > down_trans_lf(end)) = [];
down_trans_lf(down_trans_lf < up_trans_lf(1)) = [];
up_trans_hf(up_trans_hf > down_trans_hf(end)) = [];
down_trans_hf(down_trans_hf < up_trans_hf(1)) = [];

comb_up_trans = up_trans_lf;
comb_down_trans = down_trans_lf;
for i = 2:length(up_trans_lf)-1
   
    [dummy,near_hf_up] = min(abs(up_trans_lf(i)-up_trans_hf));
    if up_trans_hf(near_hf_up) > down_trans_lf(i-1) & up_trans_hf(near_hf_up) < up_trans_lf(i+1)
    comb_up_trans(i) = up_trans_hf(near_hf_up);
    end
    
end

comb_up_trans(comb_up_trans > comb_down_trans(end)) = [];
comb_down_trans(comb_down_trans < comb_up_trans(1)) = [];

comb_state_seq = ones(size(lf_state_seq));
for i = 1:length(comb_up_trans)
    comb_state_seq(comb_up_trans(i):comb_down_trans(i)) = 2;
end
