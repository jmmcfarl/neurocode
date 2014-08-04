function [f] = plot2d_mod_dualeye(mod)


f = figure;

pids = mod.mods(1).pids;
sdim = sqrt(mod.mods(1).fsdim);
if round(sdim) ~= sdim
    error('non-square pixel array')
end
% k_mat = get_k_mat(mod);
% if strcmp(mod.basis,'white')
%     allks = (k_mat'*mod.pix_conv_mat)';
% elseif strcmp(mod.basis,'pix')
%     allks = k_mat;
% end
allks = get_pix_mat(mod);
n_ks = size(allks,1);
if nargin < 3
    targets = 1:size(allks,2);
end

%figure;
nslices = size(allks,1)/length(pids)/2;
nfilts  = length(targets);
n_cols = nslices;
zmin    = min(min(allks(:,targets))); zmax=max(max(allks(:,targets)));
disp('*** plotting normalized over all filters & times')
for ifilt = 1:nfilts;
    curmod = mod.mods(targets(ifilt));
    tcell_left = pad2dslices(allks(1:n_ks/2,targets(ifilt)),pids,sdim);
    tcell_right = pad2dslices(allks(n_ks/2+1:end,targets(ifilt)),pids,sdim);
    zmin    = min(allks(:,targets(ifilt)));
    zmax=max(allks(:,targets(ifilt)));
    zlarge = max(abs(zmin),abs(zmax));
%     zlarge = max((zmax));
%     zsmall = min(zmin);
%     if zlarge == 0
%         zlarge = 1;
%     end
    for islice = 1:nslices;
        subplot(2*nfilts,n_cols,(ifilt-1)*2*n_cols+islice);
        plot2drfmat(flipud(tcell_left{islice}),[-zlarge,zlarge]);
        subplot(2*nfilts,n_cols,(ifilt-1)*2*n_cols+n_cols+islice);
        plot2drfmat(flipud(tcell_right{islice}),[-zlarge,zlarge]);        
    end        
end

