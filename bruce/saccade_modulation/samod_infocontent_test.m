[Xinds,Tinds] = meshgrid(1:use_nPix_us,1:flen);
ulag = 6;
uk =find(Tinds == ulag);

T = tabulate(cur_Robs);
un_resp = T(:,1);
marg_prob = T(:,3)/100;

poss_p = [-1 0 1];

cur_stims = all_Xmat_shift(:,uk);


for cur_pix = 1:use_nPix_us;

cur_pixvals = cur_stims(:,cur_pix);

pcond_rdist = nan(3,length(un_resp));
npix_counts = nan(3,1);
for vv = 1:length(poss_p)
   cset = find(cur_pixvals == poss_p(vv));
   npix_counts(vv) = length(cset);
   pcond_rdist(vv,:) = hist(cur_Robs(cset),un_resp);
   pcond_rdist(vv,:) = pcond_rdist(vv,:);
end
pcond_rdist = bsxfun(@rdivide,pcond_rdist,npix_counts);
pix_prob = npix_counts/sum(npix_counts);
pix_info(cur_pix) = pix_prob'*nansum(pcond_rdist.*log2(bsxfun(@rdivide,pcond_rdist,marg_prob')),2);
end

%%
sac_pix_info = nan(n_sac_bins,use_nPix_us);
sac_avg_rate = nan(n_sac_bins,1);
for nn = 1:n_sac_bins
    cur_sacinds = find(Xsac(cc_uinds,nn) == 1);
    
    cur_stims = all_Xmat_shift(cur_sacinds,uk);
    cur_resp = cur_Robs(cur_sacinds);
    marg_prob = hist(cur_resp,un_resp)'/length(cur_resp);
    for cur_pix = 1:use_nPix_us;
        
        cur_pixvals = cur_stims(:,cur_pix);
        
        pcond_rdist = nan(3,length(un_resp));
        npix_counts = nan(3,1);
        for vv = 1:length(poss_p)
            cset = find(cur_pixvals == poss_p(vv));
            npix_counts(vv) = length(cset);
            pcond_rdist(vv,:) = hist(cur_resp(cset),un_resp);
            pcond_rdist(vv,:) = pcond_rdist(vv,:);
        end
        pcond_rdist = bsxfun(@rdivide,pcond_rdist,npix_counts);
        pix_prob = npix_counts/sum(npix_counts);
        sac_pix_info(nn,cur_pix) = pix_prob'*nansum(pcond_rdist.*log2(bsxfun(@rdivide,pcond_rdist,marg_prob')),2);
    end
    sac_avg_rate(nn) = mean(cur_resp);
end
sac_pix_info_PS = bsxfun(@rdivide,sac_pix_info,sac_avg_rate);