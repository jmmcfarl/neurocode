clear all
close all

outdirectory = '~/James_scripts/data_processing/Images/object_analysis';
obj_directory = '~/James_scripts/data_processing/Images/Object_images/';
white = 255;
black = 0;
resc_bounds = 2.5;
max_rot = 30;
rscl_factor = 1;
Pix2Deg = 0.018837;
 
cd '/Users/James/James_scripts/data_processing/Images'
load obj_params

cd(obj_directory)

dd = dir('./');
dd = dd(4:end);
n_objs = length(dd);

n_images = 50;
for ii = 1:n_objs
    
    %generate noise background
    Nx = 1280; Ny = 1024;
    dim = [Ny Nx];
    noise_back = spatialPattern(dim,-2);
    f = highpassfilter(dim,0.2*Pix2Deg,2); %cutoff freq at 2 cycle/degree, assuming 60 pix/degree
    noise_back = real(ifft2(fft2(noise_back).*f));
    noise_std = 60;
    noise_mean = (white+black)/2;
    noise_back = noise_back/std(noise_back(:))*noise_std + noise_mean;
    noise_back(noise_back > white) = white;
    noise_back(noise_back < black) = black;
    new_im = noise_back;
    
    cd(obj_directory)
    filename = dd(ii).name;
    disp(filename);
    imdata = imread(filename);
    imdata = rgb2gray(imdata);
    imdata = double(imdata); % convert to double format
    imdataorg = imdata;    
    
    %identify object pixels
    outmask = (imdata >= white - obj_params(ii,1));
    mask = ones(size(imdata));
    mask(outmask) = 0;
    [B,L,N,A] = bwboundaries(mask);
    
    bound_lengths = cellfun(@(x) length(x),B);
    if strcmp(filename,'lobster.jpg')
       bad_set = find(bound_lengths < 30);
    else
        bad_set = [];
    end
    [~,ord] = sort(bound_lengths,'descend');
    B = B(ord);
    new_L = L;
    for i = 1:length(unique(L(:)))-1
        new_L(L==ord(i)) = i;
    end
    L = new_L;
    A = A(ord,ord);
        
    if obj_params(ii,4) == 1
        enclosed_set = find(A(:,1));
    else
        enclosed_set = [];
    end
    enclosed_set(ismember(enclosed_set,bad_set)) = [];
    non_enclosed_set = setdiff(1:length(B),enclosed_set);
    L(ismember(L,non_enclosed_set)) = 1;
    obj_boundary = B{1};
    for j = 1:length(enclosed_set)
        obj_boundary = [obj_boundary; B{enclosed_set(j)}];
    end
    
    %get rid of pixels too close to object boundary
    in_obj = find(L==1);
    [sz_y,sz_x] = size(imdata);
    [X,Y] = meshgrid(1:sz_x,1:sz_y);
    obj_pos = [Y(in_obj) X(in_obj)];
    [D,I] = pdist2(obj_boundary,obj_pos,'euclidean','Smallest',1);
    boundary_pix = in_obj(D <= obj_params(ii,3));
    L(L > 1) = 0;
    L_old = L;
    L(boundary_pix) = 0;
    L_old(boundary_pix) = 2;
    
    %rescale object image
    imdata = imresize(imdata,rscl_factor);
    L = round(imresize(L,rscl_factor));
    in_obj = find(L==1);
    
%     %rescale object intensity distribution
%     bounds = prctile(imdata(in_obj),[resc_bounds (100-resc_bounds)]);
%     resc_factor = (white-black)/diff(bounds);
% %     resc_factor = min(resc_factor,max_resc_factor);
%     resc_imdata = imdata*resc_factor;
%     resc_imdata = (resc_imdata - bounds(1)*resc_factor + black);
%     resc_imdata(resc_imdata < black) = black;
%     resc_imdata(resc_imdata > white) = white;
%     imdata = resc_imdata;
    
    %
    imdata2 = imdataorg/white;
    J = imadjust(imdata2(in_obj),stretchlim(imdata2(in_obj),[.025 .975]),[0 1],obj_params(ii,2));
    imdata(in_obj) = J*white;
    imdata(L~=1) = white;
    
    %store stats of object intensity dist
    obj_mean(ii) = mean(imdata(in_obj));
    obj_std(ii) = std(imdata(in_obj));
    obj_area(ii) = length(find(L==1));
        
    
    
    x_shift = 0;
    y_shift = 0;
    [Nyo,Nxo] = size(L);
    bk_L = zeros(size(noise_back));
    y_inds = (round(Ny/2-Nyo/2):round(Ny/2+Nyo/2)-1)+y_shift;
    x_inds = (round(Nx/2-Nxo/2):round(Nx/2+Nxo/2)-1)+x_shift;
    used_y = find(y_inds > 0 & y_inds <= Ny);
    used_x = find(x_inds > 0 & x_inds <= Nx);
    bk_L(y_inds(used_y),x_inds(used_x)) = L(used_y,used_x);
    new_im(bk_L==1) = imdata(L==1);
        
    cd(outdirectory)
    
%     figure
%     new_im = flipud(new_im);
%     imagesc(new_im); colormap gray
%     set(gca,'ydir','normal')
%     
%     cur_fname = sprintf('objIM%.4d_2',ii);
%     print('-dpng',cur_fname);
%     close all
%     
%     figure
%     imagesc(imdata); colormap gray
%     hold on
%     plot(obj_boundary(:,2),obj_boundary(:,1),'r.','linewidth',2)
%     
%     cur_fname = sprintf('objmask%.4d',ii);
%     print('-dpng',cur_fname);
%     close all 
    
end