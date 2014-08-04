function proc_s = nlin_proc_stim( stim, NL, NLx )
%
% Usage: proc_s = nlin_proc_stim( stim, NL, NLx )

[NT NX] = size(stim);

proc_s = zeros(NT*NX,1); %initialize output of module NL

% not sure what this does.
if NT > 1
  stim = reshape(stim,NT*NX,1);
end

%the overall processed (filtered) stim is given by a sum over tent bases of
%the (filtered) stim passed through each tent basis function.
for j = 1:length(NL)
  if j == 1
    proc_s = proc_s + NL(j) * piece_proc_stim( stim, NLx(j), NLx(j+1) );
  else
    if j == length(NL)
      proc_s = proc_s + NL(j) * piece_proc_stim( stim, NLx(j), NLx(j-1) );
    else
      proc_s = proc_s + NL(j) * piece_proc_stim( stim, NLx(j), [NLx(j-1) NLx(j+1)] );
    end
  end
end

if NX > 1
  proc_s = reshape(proc_s,NT,NX);
end