up_sdim = 33;
f_ax = (0:-1:(-flen+1))*0.01;
up_bax = un_bar_pos(1):0.04:un_bar_pos(end);


close all
for cc = tr_set
    
    cc
    cur_n_mods = length(upquad_fit(cc).mods);
    
    k_mat_i = get_k_mat(upquad_fit(cc));
    k_mat_m = get_k_mat(it_quad_fit{1}(cc));
    k_mat_f = get_k_mat(it_quad_fit{end}(cc));
    
    pre_out = X_z*k_mat_i;
    pre_out(:,2:end) = pre_out(:,2:end).^2;
    
    post_out = X_sh*k_mat_f;
    post_out(:,2:end) = post_out(:,2:end).^2;
    
    pre_std(cc,1:cur_n_mods) = std(pre_out);
    post_std(cc,1:cur_n_mods) = std(post_out);
    
    disp('Pre std')
    pre_std(cc,:)
    disp('Post std')
    post_std(cc,:)
    
    k_mag_i(cc,1:cur_n_mods) = sqrt(sum(k_mat_i.^2));
    k_mag_f(cc,1:cur_n_mods) = sqrt(sum(k_mat_f.^2));
    
    
    ca = max(abs(k_mat_f));
    
    subplot(3,cur_n_mods,1)
    imagesc(up_bax,f_ax,flipud(reshape(k_mat_i(:,1),flen,up_sdim)));
    caxis([-ca(1) ca(1)]);
    title('Pre linear','fontsize',14)
    xlabel('Bar position (deg)','fontsize',14)
    ylabel('Time lag (s)','fontsize',14)
    
    subplot(3,cur_n_mods,cur_n_mods+1)
    imagesc(up_bax,f_ax,flipud(reshape(k_mat_m(:,1),flen,up_sdim)));
    caxis([-ca(1) ca(1)]);
    title('Post linear','fontsize',14)
    xlabel('Bar position (deg)','fontsize',14)
    ylabel('Time lag (s)','fontsize',14)
    
    
    subplot(3,cur_n_mods,2*cur_n_mods+1)
    imagesc(up_bax,f_ax,flipud(reshape(k_mat_f(:,1),flen,up_sdim)));
    caxis([-ca(1) ca(1)]);
    title('Post linear','fontsize',14)
    xlabel('Bar position (deg)','fontsize',14)
    ylabel('Time lag (s)','fontsize',14)

    subplot(3,cur_n_mods,2)
    imagesc(up_bax,f_ax,flipud(reshape(k_mat_i(:,2),flen,up_sdim)));
    caxis([-ca(2) ca(2)]);
    title('Pre quad 1','fontsize',14)
    xlabel('Bar position (deg)','fontsize',14)
    ylabel('Time lag (s)','fontsize',14)
   
    subplot(3,cur_n_mods,cur_n_mods+2)
    imagesc(up_bax,f_ax,flipud(reshape(k_mat_m(:,2),flen,up_sdim)));
    caxis([-ca(2) ca(2)]);
    title('Post quad 1','fontsize',14)
    xlabel('Bar position (deg)','fontsize',14)
    ylabel('Time lag (s)','fontsize',14)
    
    subplot(3,cur_n_mods,2*cur_n_mods+2)
    imagesc(up_bax,f_ax,flipud(reshape(k_mat_f(:,2),flen,up_sdim)));
    caxis([-ca(2) ca(2)]);
    title('Post quad 1','fontsize',14)
    xlabel('Bar position (deg)','fontsize',14)
    ylabel('Time lag (s)','fontsize',14)
    
    if cur_n_mods == 3
        subplot(3,cur_n_mods,3)
        imagesc(up_bax,f_ax,flipud(reshape(k_mat_i(:,3),flen,up_sdim)));
        caxis([-ca(3) ca(3)]);
        title('Pre quad 2','fontsize',14)
        xlabel('Bar position (deg)','fontsize',14)
        ylabel('Time lag (s)','fontsize',14)
        
        subplot(3,cur_n_mods,6)
        imagesc(up_bax,f_ax,flipud(reshape(k_mat_m(:,3),flen,up_sdim)));
        caxis([-ca(3) ca(3)]);
        title('Post quad 2','fontsize',14)
        xlabel('Bar position (deg)','fontsize',14)
        ylabel('Time lag (s)','fontsize',14)
        
        subplot(3,cur_n_mods,9)
        imagesc(up_bax,f_ax,flipud(reshape(k_mat_f(:,3),flen,up_sdim)));
        caxis([-ca(3) ca(3)]);
        title('Post quad 2','fontsize',14)
        xlabel('Bar position (deg)','fontsize',14)
        ylabel('Time lag (s)','fontsize',14)

    end
    pause
    clf
end
