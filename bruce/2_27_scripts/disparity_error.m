function tot_err = disparity_error(K,reye_pos,leye_pos,scale_pen)

cur_r(:,1) = K(1)*reye_pos(:,1) + K(2)*reye_pos(:,2);
cur_r(:,2) = K(3)*reye_pos(:,1) + K(4)*reye_pos(:,2);

cur_l(:,1) = K(5)*leye_pos(:,1) + K(6)*leye_pos(:,2);
cur_l(:,2) = K(7)*leye_pos(:,1) + K(8)*leye_pos(:,2);

cur_disp = cur_r - cur_l;

disp_err = mean(abs(cur_disp(:,1))+abs(cur_disp(:,2)));

Al = [K(1) K(2); K(3) K(4)];
Ar = [K(5) K(6); K(7) K(8)];
scale1 = (1-det(Al))^2;
scale2 = (1-det(Ar))^2;

tot_err = disp_err + scale_pen*scale1 + scale_pen*scale2;