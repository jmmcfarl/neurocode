function [new_state_seq] = fixedthresh_state_seq_smooth(hmm,emiss,state_seq)

min_dur = round(0.1*hmm.Fs);

for ns = 1:hmm.Nsegs
    new_state_seq{ns} = state_seq{ns};
    up_trans = 1+find(state_seq{ns}(1:end-1)==1 & state_seq{ns}(2:end) == 2);
    down_trans = 1+find(state_seq{ns}(1:end-1)==2 & state_seq{ns}(2:end) == 1);
    state_trans = sort([up_trans; down_trans]);
    for i = 2:length(state_trans)
       cur_dur = state_trans(i)-state_trans(i-1);
       if cur_dur < min_dur
          if ismember(state_trans(i),up_trans)
             new_state_seq{ns}(state_trans(i-1):state_trans(i)) = 2;
          else
             new_state_seq{ns}(state_trans(i-1):state_trans(i)) = 1;              
          end
       end
    end
end

