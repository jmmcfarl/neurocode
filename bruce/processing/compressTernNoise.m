function stimComp = compressTernNoise(X)

[stimComp.NT,stimComp.K] = size(X);

stimComp.nzero = X~=0;
bitvals = X(stimComp.nzero);
bitvals(bitvals == -1) = 0;
stimComp.lvals = logical(bitvals);
