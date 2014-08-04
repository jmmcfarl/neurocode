function brfstim = tbrep(stim, tentx)

%stim is the input stimulus
%tentx is the tent basis boundary positions
%outputs brfstim which is (NTxNBFs)

NN=length(tentx);
brfstim = zeros(length(stim),NN);
for n = 1:NN
    if n == 1
      brfstim(:,n) = piece_proc_stim(stim, tentx(n), tentx(n+1) );
    else
      if n == NN
        brfstim(:,n) = piece_proc_stim(stim, tentx(n), tentx(n-1) );
      else
        brfstim(:,n) = piece_proc_stim(stim, tentx(n), [tentx(n-1) tentx(n+1)] );
      end
    end
end

end

