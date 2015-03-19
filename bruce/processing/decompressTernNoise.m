function stim_mat = decompressTernNoise(stimComp)

stim_mat = zeros(stimComp.NT,stimComp.K);

pix_data = ones(size(stimComp.lvals));
pix_data(~stimComp.lvals) = -1;
stim_mat(stimComp.nzero) = pix_data;


