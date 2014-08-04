clear all
close all

outdirectory = '~/James_scripts/data_processing/Images/object_test_ims';
obj_directory = '~/James_scripts/data_processing/Images/Object_images/';
white = 255;
black = 0;
resc_bounds = 10;
max_rot = 30;
rscl_factor = 1;
Pix2Deg = 0.018837;
max_x_shift = round(3/Pix2Deg);
max_y_shift = round(3/Pix2Deg);
avg_obj = 2.5;

min_uncovered = 0.75;
texture_prob = 0.1; %probability of texture image

cd '/Users/James/James_scripts/data_processing/Images'
load obj_params

cd(obj_directory)

dd = dir('./');
dd = dd(4:end);
n_objs = length(dd);

n_images = 343;
cnt = 1;

while cnt < n_images
    
    %generate noise background
    Nx = 1280; Ny = 1024;
    dim = [Ny Nx];
    noise_back = spatialPattern(dim,-2);
    f = highpassfilter(dim,0.2*Pix2Deg,2); %cutoff freq at 2 cycle/degree, assuming 60 pix/degree
    noise_back = real(ifft2(fft2(noise_back).*f));
    noise_std = 61;
    noise_mean = (white+black)/2;
    noise_back = noise_back/std(noise_back(:))*noise_std + noise_mean;
    noise_back(noise_back > white) = white;
    noise_back(noise_back < black) = black;
    good_image = 0;
    
    %determine whether it's a texture patch or object image
    if rand > texture_prob
        n_used_objs = 2;
        n_textures = 0;
    else
        n_textures = 2;
        n_used_objs = 0;
    end
    
    %iterate until we have an accepted image
    while good_image == 0
        
        new_im = noise_back; %initialize to noise background
                
        used_set = [zeros(1,n_textures)];
        temp = randperm(n_objs);
        used_set = [used_set temp(1:n_used_objs)];
        
        n_items = n_used_objs + n_textures;
        ov_L = zeros(size(new_im));
        
        for n = 1:n_items
            
            if used_set(n) > 0 %it current item is an object
                cd(obj_directory)
                filename = dd(used_set(n)).name;
                disp(filename);
                imdata = imread(filename);
                imdata = rgb2gray(imdata);
                imdata = double(imdata); % convert to double format
                item_info(n).item_type = filename(1:end-4);
                
                %random rotation
                offset = 1;
                cur_rot = rand*2*max_rot-max_rot;
                imdata_r = imrotate(imdata+offset,cur_rot);
                imdata_r = imdata_r - offset;
                imdata_r(imdata_r == -1) = white;
                imdata = imdata_r;
                item_info(n).rot_angle = cur_rot;
                
                %identify object pixels
                outmask = (imdata >= white - obj_params(used_set(n),1));
                mask = ones(size(imdata));
                mask(outmask) = 0;
                [B,L,N,A] = bwboundaries(mask);
                
                bound_lengths = cellfun(@(x) length(x),B);
                if strcmp(filename,'lobster.jpg') %this is a hack just for the lobster image
                    bad_set = find(bound_lengths < 30);
                else
                    bad_set = [];
                end
                
                %sort boundaries by length
                [~,ord] = sort(bound_lengths,'descend');
                B = B(ord);
                new_L = L;
                for i = 1:length(unique(L(:)))-1
                    new_L(L==ord(i)) = i;
                end
                L = new_L;
                A = A(ord,ord);
                
                %find boundaries that are enclosed by the parent object
                if obj_params(used_set(n),4) == 1
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
                boundary_pix = in_obj(D <= obj_params(used_set(n),3));
                L(L > 1) = 0;
                L_old = L;
                L(boundary_pix) = 0;
                L_old(boundary_pix) = 2;
                
                %rescale object image
                imdata = imresize(imdata,rscl_factor);
                L = round(imresize(L,rscl_factor));
                in_obj = find(L==1);
                item_info(n).pix_area = length(in_obj);
                
                %     %rescale object intensity distribution
                %     bounds = prctile(imdata(in_obj),[resc_bounds (100-resc_bounds)]);
                %     resc_factor = (white-black)/diff(bounds);
                % %     resc_factor = min(resc_factor,max_resc_factor);
                %     resc_imdata = imdata*resc_factor;
                %     resc_imdata = (resc_imdata - bounds(1)*resc_factor + black);
                %     resc_imdata(resc_imdata < black) = black;
                %     resc_imdata(resc_imdata > white) = white;
                %     imdata = resc_imdata;
                
                % transform intensity distribution
                imdata2 = imdata/white;
                J = imadjust(imdata2(in_obj),stretchlim(imdata2(in_obj),[.025 .975]),[0 1],obj_params(used_set(n),2));
                imdata(in_obj) = J*white;
                imdata(L~=1) = white;
                
                item_info(n).mean_int = mean(imdata(L==1));
                item_info(n).std_int = std(imdata(L==1));
                
                
            else %for texture patches
                
                text_size = 400;
                [texture_patch,comp_gratings] = generate_texture_patch(text_size);
                texture_patch = texture_patch/std(texture_patch(:));
                texture_patch = texture_patch*noise_std + noise_mean;
                texture_patch(texture_patch < black) = black;
                texture_patch(texture_patch > white) = white;
                imdata = texture_patch;
                
                L = zeros(text_size);
                [X,Y] = meshgrid(0.5:text_size);
                X = X - text_size/2;
                Y = Y - text_size/2;
                
                rad_mean = 2.5/Pix2Deg;
                rad_std = 0.5/Pix2Deg;
                min_rad = 0.5/Pix2Deg;
                rand_rad = max(min_rad,randn*rad_std+rad_mean);
                L(X.^2+Y.^2 <= rand_rad.^2) = 1;
                item_info(n).pix_area = pi*rand_rad^2;
                item_info(n).mean_int = noise_mean;
                item_info(n).std_int = noise_std;
                item_info(n).item_type = 'texture_patch';
                item_info(n).texture_info = comp_gratings;
                
            end
            
            x_shift = round(2*rand*max_x_shift-max_x_shift);
            y_shift = round(2*rand*max_y_shift-max_y_shift);
            item_info(n).x_cent = x_shift;
            item_info(n).y_cent = y_shift;
            [Nyo,Nxo] = size(L);
            bk_L = zeros(size(noise_back));
            y_inds = (round(Ny/2-Nyo/2):round(Ny/2+Nyo/2)-1)+y_shift;
            x_inds = (round(Nx/2-Nxo/2):round(Nx/2+Nxo/2)-1)+x_shift;
            used_y = find(y_inds > 0 & y_inds <= Ny);
            used_x = find(x_inds > 0 & x_inds <= Nx);
            bk_L(y_inds(used_y),x_inds(used_x)) = L(used_y,used_x);
            new_im(bk_L==1) = imdata(L==1);
            
            ov_L(bk_L==1) = n;
        end
        
        obj_ids = [];
        overlap = ones(n_items,1);
        for i = 1:n_items
            if ~strcmp(item_info(i).item_type,'texture_patch')
                overlap(i) = length(find(ov_L==i))/item_info(i).pix_area;
                obj_ids = [obj_ids i];
            end
        end
        if min(overlap) < min_uncovered
            good_image = 0;
        else
            good_image = 1;
        end
        if good_image == 1
            [XX,YY] = meshgrid(1:Nx,1:Ny);
            XX = XX - Nx/2; YY = YY - Ny/2;
            cent_window = find(abs(XX) < 3/Pix2Deg & abs(YY) < 3/Pix2Deg);
            temp = ismember(ov_L(cent_window),obj_ids);
            obj_area(cnt) = sum(temp(:))/length(cent_window);
            temp = find(ov_L(cent_window) ~= 0);
            item_area(cnt) = length(temp)/length(cent_window);
            
            %         figure
            %         new_im = flipud(new_im);
            %         imagesc(new_im); colormap gray
            %         set(gca,'ydir','normal')
            %         hold on
            %         plot(Nx/2,Ny/2,'b+','linewidth',2,'markersize',8)
            %         x_cent = round(0.35/Pix2Deg);
            %         y_cent = round(-0.45/Pix2Deg);
            %         x_width = round(0.7/Pix2Deg);
            %         y_width = round(0.7/Pix2Deg);
            %         x_cent = x_cent - x_width/2;
            %         y_cent = y_cent - y_width/2;
            %         rectangle('Position',[Nx/2+x_cent Ny/2+y_cent x_width y_width],'linewidth',1,'edgecolor','r')
            %         x_cent = 0;
            %         y_cent = 0;
            %         x_width = round(10/Pix2Deg);
            %         y_width = round(10/Pix2Deg);
            %         x_cent = x_cent - x_width/2;
            %         y_cent = y_cent - y_width/2;
            %         rectangle('Position',[Nx/2+x_cent Ny/2+y_cent x_width y_width],'linewidth',1,'edgecolor','w')
            %         axis off
            %         cd(outdirectory)
            %         cur_fname = sprintf('objIM%.4d',cnt);
            %         print('-dpng',cur_fname);
            %         close all
            
            cd(outdirectory)
            imdatapng = new_im./white;
            cur_fname = sprintf('4%.4d',cnt);
            imwrite(imdatapng,strcat(outdirectory, '/', cur_fname, '.png'),'png');
            
            info_name = sprintf('info4%.4d',cnt);   
            save(info_name,'ov_L','item_info');
            
            cnt = cnt +1;
        else
            disp('BAD IMAGE')
            
        end
    end
end

%%

