win_dur = round(25*Fsd);
Nwins = ceil(length(t_axis)/win_dur);

pers_mp_ups = up_trans_inds(t2_ups{d});
rpers_mp_ups = up_trans_inds(rt2_ups{d});
pers_mp_downs = down_trans_inds(t2_downs{d})+1;
rpers_mp_downs = down_trans_inds(rt2_downs{d})+1;

close all
plot(t_axis,mp_state_vec+1.2)
hold on
plot(t_axis,lfp_state_vec,'r')
plot(t_axis(pers_mp_ups),mp_state_vec(pers_mp_ups)+1.2,'ko');
plot(t_axis(rpers_mp_ups),mp_state_vec(rpers_mp_ups)+1.2,'go');
plot(t_axis(pers_mp_downs),mp_state_vec(pers_mp_downs)+1.2,'ko');
plot(t_axis(rpers_mp_downs),mp_state_vec(rpers_mp_downs)+1.2,'go');

for ii = 1:length(up_trans_inds)
    if ~isnan(corresp_lf8_upinds{d}(ii))
   line(t_axis([up_trans_inds(ii) up_trans_inds8(corresp_lf8_upinds{d}(ii))]),[1.8 0.6],'color','k','linewidth',2) 
    end
    if ~isnan(corresp_lf8_downinds{d}(ii))
   line(t_axis([down_trans_inds(ii) down_trans_inds8(corresp_lf8_downinds{d}(ii))]),[1.8 0.6],'color','m','linewidth',2) 
    end
end

ylim([-0.2 2.4])
for ii = 1:Nwins
    cur_inds = (ii-1)*win_dur + (1:win_dur);
    xlim(t_axis(cur_inds([1 end])));
    pause
    
end


