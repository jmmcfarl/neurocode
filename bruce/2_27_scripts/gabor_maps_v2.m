clear all close all
addpath(genpath('~/Data/bruce/2_27_12'))
cd ~/Data/bruce/2_27_12/stimrecon/
addpath(genpath('~/James_scripts/'));

%%
Pix2Deg = 0.018837;
% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
Nyd = Nyp/dsfrac;
Nxd = Nxp/dsfrac;
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));


%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
% load fixation_stim_righteye
load fixation_data
% load fixation_image_patches
load fixation_image_patches
% X = X_right;
X = X_left;

NT = size(X,1);
SDIM = size(X,2);

%%

kern_len = SDIM-1;
[XX,YY] = meshgrid(-kern_len/2:kern_len/2,-kern_len/2:kern_len/2);
orientations = linspace(0,pi,9);
orientations(end) = [];
lambdas = 1./[1 1.5 2 2.75 3.75]*Fsd;
phases = [0 pi/2];
gabor_bank = zeros(length(orientations),length(lambdas),length(phases),size(XX,1),size(XX,2));
for i = 1:length(orientations)
    for j = 1:length(lambdas)
        for k = 1:length(phases)
            gabor_bank(i,j,k,:,:) = get_gabor_template(XX,YY,0,0,orientations(i),lambdas(j),phases(k));
        end
    end
end

%%

X_gabor_energy = zeros([length(orientations) length(lambdas) size(X)]);

% for i = 1:NT
%     fprintf('Fixation %d of %d\n',i,NT);
%     for j = 1:length(orientations)
%         for k = 1:length(lambdas)
%             phase1 = conv2(squeeze(X(i,:,:)),squeeze(gabor_bank(j,k,1,:,:)),'same');            
%             phase2 = conv2(squeeze(X(i,:,:)),squeeze(gabor_bank(j,k,2,:,:)),'same');
%             X_gabor_energy(j,k,i,:,:) = sqrt(phase1.^2+phase2.^2);
%         end
%     end
% end

resh_X = reshape(X,NT,SDIM^2);
delta_x = -kern_len/2:kern_len/2;
delta_y = -kern_len/2:kern_len/2;
for j = 1:length(orientations)
    fprintf('Orientation %d of %d\n',j,length(orientations));
    for k = 1:length(lambdas)
        cur_mask1 = squeeze(gabor_bank(j,k,1,:,:));
        cur_mask2 = squeeze(gabor_bank(j,k,2,:,:));
        for ii = 1:length(delta_x)
            for jj = 1:length(delta_y)
                smask1 = dist_shift2d(cur_mask1,delta_x(ii),2,0);
                smask1 = dist_shift2d(smask1,delta_y(jj),1,0);
                smask2 = dist_shift2d(cur_mask2,delta_x(ii),2,0);
                smask2 = dist_shift2d(smask2,delta_y(jj),1,0);
                smask1 = smask1(:);
                smask2 = smask2(:);
                gabor_out1 = resh_X*smask1;
                gabor_out2 = resh_X*smask2;
                X_gabor_energy(j,k,:,jj,ii) = sqrt(gabor_out1.^2+gabor_out2.^2);
            end
        end
    end
end

X_gabor_energy_norm = bsxfun(@minus,X_gabor_energy,mean(X_gabor_energy,3));
X_gabor_energy_norm = bsxfun(@rdivide,X_gabor_energy_norm,std(X_gabor_energy_norm,[],3));
clear X_gabor_energy

%%
spk_trg_avgs = zeros(10,length(orientations),length(lambdas),size(X,2),size(X,3));
for i = 1:10
    fprintf('Cell %d of %d\n',i,10);
    temp = bsxfun(@times,X_gabor_energy_norm,reshape(spk_cnts(used_fixs,i),[1 1 NT 1 1]));
    spk_trg_avgs(i,:,:,:,:) = squeeze(sum(temp,3))/sum(spk_cnts(used_fixs,i));
end
clear temp

%%
for cellid = 1:10
    fprintf('Cell %d of %d\n',cellid,10);
    cur_set = squeeze(spk_trg_avgs(cellid,:,:,:,:));
    [maxval,maxloc] = max(cur_set(:));
    [peak_or(cellid),peak_lam(cellid),peak_i(cellid),peak_j(cellid)] = ind2sub(size(cur_set),maxloc);
    or_tuning(cellid,:) = squeeze(cur_set(:,peak_lam(cellid),peak_i(cellid),peak_j(cellid)));
end

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
save gabor_map_v3_stas spk_trg_avgs orientations lambdas peak_* or_tuning* *patch_inds xax yax
%%
close all
cd /Users/James/James_scripts/bruce/gabor_models
for cellid = 1:10;
    
    cur_set = squeeze(spk_trg_avgs(cellid,:,:,:,:));
    minval = min(cur_set(:));
    maxval = max(cur_set(:));
    figure
    for i = 1:length(orientations)
        for j = 1:length(lambdas)
            subplot(length(orientations),length(lambdas),(i-1)*length(lambdas)+j)
            imagesc(xax(xpatch_inds),yax(ypatch_inds),squeeze(cur_set(i,j,:,:)));
%             axis square
            set(gca,'ydir','normal');
            caxis([minval maxval]);
            if i == peak_or(cellid) & j == peak_lam(cellid)
                hold on
%                 plot(xax(xpatch_inds(peak_j(cellid))),yax(ypatch_inds(peak_i(cellid))),'w*','markersize',8)
            end
        end
    end
    fname = sprintf('Cell%d_STAmap',cellid);
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [10 15]);
    print('-dpdf','-painters',fname);close all

end


%%
SDIM = size(X,2);
nphases = 10;
phases = linspace(0,(nphases-1)*pi/nphases,nphases);
gabor_bank = zeros(10,length(phases),size(XX,1),size(XX,2));
for i = 1:10
    for k = 1:length(phases)
        gabor_bank(i,k,:,:) = get_gabor_template(XX,YY,0,0,orientations(peak_or(i)),lambdas(peak_lam(i)),phases(k));
    end
end

NT = size(X,1);
nzpad = 30;
zmat = zeros(SDIM+2*nzpad,SDIM+2*nzpad);
X_gabor_amp = zeros([10 length(phases) NT]);
for j = 1:10
    for k = 1:length(phases)
        cur_mask = zmat;
        cur_mask(nzpad+(-kern_len/2:kern_len/2)+peak_i(j),nzpad+(-kern_len/2:kern_len/2)+peak_j(j)) = squeeze(gabor_bank(j,k,:,:));
        cur_mask = cur_mask(nzpad+(1:SDIM),nzpad+(1:SDIM));
        cur_mask = cur_mask(:);
        X_gabor_amp(j,k,:) = reshape(X,size(X,1),SDIM^2)*cur_mask;
    end
end

X_gabor_amp_norm = bsxfun(@minus,X_gabor_amp,mean(X_gabor_amp,3));
X_gabor_amp_norm = bsxfun(@rdivide,X_gabor_amp_norm,std(X_gabor_amp_norm,[],3));
clear X_gabor_amp

spk_trg_avgs_phase = zeros(10,length(phases));
for i = 1:10
    temp = bsxfun(@times,squeeze(X_gabor_amp_norm(i,:,:)),reshape(spk_cnts(used_fixs,i),[1 NT]));
    spk_trg_avgs_phase(i,:) = squeeze(sum(temp,2))/sum(spk_cnts(used_fixs,i));
end
clear temp

%%
norientations = 10;
cur_orientations = linspace(0,(norientations-1)*pi/norientations,norientations);
gabor_bank = zeros(10,length(norientations),length(phases),size(XX,1),size(XX,2));
for i = 1:10
    for j = 1:norientations
        for k = 1:nphases
            gabor_bank(i,j,k,:,:) = get_gabor_template(XX,YY,0,0,cur_orientations(j),lambdas(peak_lam(i)),phases(k));
        end
    end
end

NT = size(X,1);
nzpad = 30;
zmat = zeros(SDIM+2*nzpad,SDIM+2*nzpad);
X_gabor_amp = zeros([10 length(cur_orientations) length(phases) NT]);
for i = 1:10
    for j = 1:length(cur_orientations)
        for k = 1:length(phases)
            cur_mask = zmat;
            cur_mask(nzpad+(-kern_len/2:kern_len/2)+peak_i(i),nzpad+(-kern_len/2:kern_len/2)+peak_j(i)) = squeeze(gabor_bank(i,j,k,:,:));
            cur_mask = cur_mask(nzpad+(1:SDIM),nzpad+(1:SDIM));
            cur_mask = cur_mask(:);
            X_gabor_amp(i,j,k,:) = reshape(X,size(X,1),SDIM^2)*cur_mask;
        end
    end
end

X_gabor_amp_norm = bsxfun(@minus,X_gabor_amp,mean(X_gabor_amp,4));
X_gabor_amp_norm = bsxfun(@rdivide,X_gabor_amp_norm,std(X_gabor_amp_norm,[],4));

X_gabor_energy = squeeze(mean(abs(X_gabor_amp),3));
X_gabor_enery_norm = bsxfun(@minus,X_gabor_energy,mean(X_gabor_energy,3));
X_gabor_enery_norm = bsxfun(@rdivide,X_gabor_enery_norm,std(X_gabor_enery_norm,[],3));
clear X_gabor_amp X_gabor_energy

spk_trg_avgs_orientationphase = zeros(10,norientations,nphases);
spk_trg_avgs_orientation = zeros(10,norientations);
for i = 1:10
    temp = bsxfun(@times,squeeze(X_gabor_amp_norm(i,:,:,:)),reshape(spk_cnts(used_fixs,i),[1 1 NT]));
    spk_trg_avgs_orientationphase(i,:,:) = squeeze(sum(temp,3))/sum(spk_cnts(used_fixs,i));
    temp = bsxfun(@times,squeeze(X_gabor_enery_norm(i,:,:)),reshape(spk_cnts(used_fixs,i),[1 NT]));
    spk_trg_avgs_orientation(i,:) = squeeze(sum(temp,2))/sum(spk_cnts(used_fixs,i));
end
clear temp

%%

cd ~/James_scripts/bruce/
for cellid = 1:10;
    
    f1 = figure;
    
    subplot(2,2,1);
    imagesc(xax(xpatch_inds),yax(ypatch_inds),squeeze(spk_trg_avgs(cellid,peak_or(cellid),peak_lam(cellid),:,:)));
    set(gca,'ydir','normal');
    axis square
    size1 = get(gca,'Position');
    colorbar
    set(gca,'Position',size1);
    xlabel('Horizontal position (degrees)','fontsize',14)
    ylabel('Vertical position (degrees)','fontsize',14);
    
    subplot(2,2,2);
    imagesc(phases,cur_orientations,squeeze(spk_trg_avgs_orientationphase(cellid,:,:)));
    set(gca,'ydir','normal');
    size1 = get(gca,'Position');
    colorbar
     set(gca,'Position',size1);
   axis square
    xlabel('Gabor phase (degrees)','fontsize',14);
    ylabel('Gabor orientation (degrees)','fontsize',14);
    
    subplot(2,2,3);
    cmap = jet(nphases);
    hold on
    for i = 1:nphases
        plot(cur_orientations*180/pi,spk_trg_avgs_orientationphase(cellid,:,i),'color',cmap(i,:));
    end
    plot(cur_orientations*180/pi,spk_trg_avgs_orientation(cellid,:),'k','linewidth',2)
    axis tight
    xlabel('Gabor orientation (degrees)','fontsize',14);
    ylabel('Spike-triggered average stimulus (Z)','fontsize',14);
    
    
    subplot(2,2,4)
    hold on
    cmap = jet(norientations);
    for i = 1:norientations
        plot(phases*180/pi,squeeze(spk_trg_avgs_orientationphase(cellid,i,:)),'color',cmap(i,:));
    end
    plot(phases*180/pi,spk_trg_avgs_phase(cellid,:),'k','linewidth',2)
    axis tight
    xlabel('Gabor phase (degrees)','fontsize',14);
    ylabel('Spike-triggered average stimulus (Z)','fontsize',14);
    % axis square
    
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [11 8.5]);
    fname = sprintf('Gabor_properties_cell%d',cellid);
    print('-dpdf','-painters',fname);close all
    
end

