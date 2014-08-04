function [tr_inds,xv_inds] = create_xv_set(NT,nfold,nparts)

partlen = floor(NT/nparts);
nxvparts = nparts/nfold;

%boundaries of parts
pbounds  = [(0:nparts-1)*partlen+1;(1:nparts)*partlen]';
cur_perm = randperm(nparts);

for j = 1:nfold
    cur_xv_parts = sort(cur_perm((j-1)*nxvparts+(1:nxvparts)));
    
    xv_inds{j} = [];
    for i = 1:nxvparts
        xv_inds{j} = [xv_inds{j} pbounds(cur_xv_parts(i),1):pbounds(cur_xv_parts(i),2)];
    end
    tr_inds{j} = setdiff(1:NT,xv_inds{j});
end
