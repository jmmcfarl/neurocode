% for i = 1:100
i = 17;
i
    %     cur_t = sac_times(i)+0.05;
    %     cur_t2 = sac_times(i+1)-0.025;
    cur_t = eye_tx(find(eye_speed < 2.5 & eye_tx' > sac_times(i),1,'first'));
    cur_t2 = eye_tx(find(eye_speed < 2.5 & eye_tx' < sac_times(i+1),1,'last'));
    
    % cur_data = find(eye_tx > 101.55 & eye_tx < 101.95);
    cur_data = find(eye_tx > cur_t & eye_tx < cur_t2);
    
    if ~any(eye_speed(cur_data) > 10) & length(cur_data) > 0.2*Eyedt
        
        % reye_xo = mean(reyepos(cur_data,1));
        % reye_yo = mean(reyepos(cur_data,2));
        % leye_xo = mean(leyepos(cur_data,1));
        % leye_yo = mean(leyepos(cur_data,2));
        
        reye_xo = reyepos(cur_data(1),1);
        reye_yo = reyepos(cur_data(1),2);
        leye_xo = leyepos(cur_data(1),1);
        leye_yo = leyepos(cur_data(1),2);
        
        subplot(3,1,1)
        % scatter(reyepos(cur_data,1)-reye_xo,reyepos(cur_data,2)-reye_yo,ones(size(cur_data))*20,'r')
        % hold on
        % scatter(leyepos(cur_data,1)-leye_xo,leyepos(cur_data,2)-leye_yo,ones(size(cur_data))*20,'b')
        set(gca,'fontsize',14)
        plot(reyepos(cur_data,1)-reye_xo,reyepos(cur_data,2)-reye_yo,'r-')
        hold on
        plot(leyepos(cur_data,1)-leye_xo,leyepos(cur_data,2)-leye_yo,'b-')
        legend('Right eye','Left eye')
        xlabel('Relative horizontal position (degrees)','fontsize',14)
        ylabel('Relative vertical position (degrees)','fontsize',14)
        subplot(3,1,2)
        % scatter(disparity(cur_data,1),disparity(cur_data,2),ones(size(cur_data))*20,'k')
        plot(disparity(cur_data,1),disparity(cur_data,2),'k-')
        xlabel('Horizontal disparity (degrees)','fontsize',14)
        ylabel('Vertical disparity (degrees)','fontsize',14)
        subplot(3,1,3)
        plot(eye_tx(cur_data),leyepos(cur_data,1)-leye_xo,'b')
        hold on
        plot(eye_tx(cur_data),leyepos(cur_data,2)-leye_yo,'b--')
        plot(eye_tx(cur_data),reyepos(cur_data,1)-reye_xo,'r')
        plot(eye_tx(cur_data),reyepos(cur_data,2)-reye_yo,'r--')
        axis tight
        legend('Left Hor','Left Vert','Right Hor','Right Vert')
        xlabel('Time (s)','fontsize',14)
        ylabel('Relative position (degrees)','fontsize',14)
        % subplot(3,1,1)
        % title('Right Eye','fontsize',14)
        % xlabel('Horizontal position (degrees)','fontsize',12)
        % ylabel('Vertical position (degrees)','fontsize',12)
        % subplot(3,1,2)
        % title('Left Eye','fontsize',14)
        % xlabel('Horizontal position (degrees)','fontsize',12)
        % ylabel('Vertical position (degrees)','fontsize',12)
        % subplot(3,1,3)
        % title('Disparity','fontsize',14)
        % xlabel('Horizontal position (degrees)','fontsize',12)
        % ylabel('Vertical position (degrees)','fontsize',12)
        
        
%         pause(1)
%         clf
%     end
end