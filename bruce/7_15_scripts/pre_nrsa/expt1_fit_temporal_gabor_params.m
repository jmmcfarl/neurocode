clear all
close all

%%
cd ~/Data/bruce/7_15_12/G029/
load ./expt1_eyecor_d1p25_nosac_v2 gabor*
load ./corrected_Xmat_flen8
flen = 8;
% load ./corrected_Xmat
% flen = 6;

fullX = fullX/std(fullX(:));
fullX_cor = fullX_cor/std(fullX_cor(:));
%%
NT = size(fullX,1);
sdim = length(xpatch_inds);
gabor_params = gabor_params_f{end};
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
for t = 1:96
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    
    gabor_filts(t,:) = gabor_emp1(:);
    gabor_filts2(t,:) = gabor_emp2(:);
end
klen = flen*sdim^2;
gabor_emp_filts = zeros(96,flen,klen);
gabor_emp_filts2 = zeros(96,flen,klen);
for f = 1:flen
    cur_kset = f:flen:klen;
    gabor_emp_filts(:,f,cur_kset) = gabor_filts;
    gabor_emp_filts2(:,f,cur_kset) = gabor_filts2;
end
gabor_emp_filts = reshape(gabor_emp_filts,96*flen,klen);
gabor_emp_filts2 = reshape(gabor_emp_filts2,96*flen,klen);

gabor_outs1 = fullX*gabor_emp_filts';
gabor_outs2 = fullX*gabor_emp_filts2';
gabor_outs1 = permute(reshape(gabor_outs1',96,flen,NT),[3 1 2]);
gabor_outs2 = permute(reshape(gabor_outs2',96,flen,NT),[3 1 2]);
gabor_outs = sqrt(gabor_outs1.^2 + gabor_outs2.^2);

gabor_outs1 = fullX_cor*gabor_emp_filts';
gabor_outs2 = fullX_cor*gabor_emp_filts2';
gabor_outs1 = permute(reshape(gabor_outs1',96,flen,NT),[3 1 2]);
gabor_outs2 = permute(reshape(gabor_outs2',96,flen,NT),[3 1 2]);
gabor_outs_cor = sqrt(gabor_outs1.^2 + gabor_outs2.^2);

%%
for t = 1:96
    spikebins = convert_to_spikebins(full_binned_spks(:,t));
    [fitp,grad] = GLMsolve_jmm(squeeze(gabor_outs(:,t,:)), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
    temp_gabor_weights(t,:) = fitp.k(1:end-1);
    temp_gabor_offset(t) = fitp.k(end);

    spikebins = convert_to_spikebins(full_binned_spks(:,t));
    [fitp,grad] = GLMsolve_jmm(squeeze(gabor_outs_cor(:,t,:)), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
    temp_gabor_weights_cor(t,:) = fitp.k(1:end-1);
    temp_gabor_offset_cor(t) = fitp.k(end);
end
    
%%
save temp_gabor_params_flen8 temp_gabor*