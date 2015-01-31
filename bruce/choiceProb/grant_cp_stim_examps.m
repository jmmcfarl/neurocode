close all
clear all
fig_dir = '/Users/james/Desktop/';

n_images = 5;
for nn = 1:n_images

n_pix = 20;
sparse = 0.25;

bp = randi(2,n_pix,2);
nz = rand(n_pix,2) < sparse;
rls = zeros(n_pix,2);
rls(nz) = bp(nz);
rls(rls==2) = -1;

[XX,YY] = meshgrid(1:n_pix);
ZZ = zeros(n_pix,n_pix);
ZZ = ZZ + repmat(rls(:,1),1,n_pix);
ZZ = ZZ + repmat(rls(:,2)',n_pix,1);

fname = [fig_dir sprintf('temp_imag%d',nn)];
f1 = figure();
imagesc(ZZ); colormap(gray);
caxis([-1.5 1.5]);
set(gca,'xtick',[],'ytick',[]);
axis square;

fig_width = 4; rel_height = 1;
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);


% pause
% close 

end