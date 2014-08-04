function [] = moment_compare_all_cells(first_mom, second_mom)

cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\cortical_dir

[dummy,layer_23] = find_sess_data_fields(sess_data,'layer','23');
[dummy,layer_56] = find_sess_data_fields(sess_data,'layer','56');
[dummy,layer_5] = find_sess_data_fields(sess_data,'layer','5');
[dummy,parietal] = find_sess_data_fields(sess_data,'region','parietal');
[dummy,prefrontal] = find_sess_data_fields(sess_data,'region','prefrontal');
[dummy,frontal] = find_sess_data_fields(sess_data,'region','frontal');

l_23_pyr_par = find(layer_23 & parietal);
l_56_pyr_par = find(layer_56 & parietal);
l_23_pyr_fro = find(layer_23 & frontal);
l_5_pyr_fro = find(layer_5 & frontal);
l_23_pyr_pre = find(layer_23 & prefrontal);
l_5_pyr_pre = find(layer_5 & prefrontal);

plot(first_mom(l_23_pyr_par),second_mom(l_23_pyr_par),'o','linewidth',2)
hold on
plot(first_mom(l_56_pyr_par),second_mom(l_56_pyr_par),'ro','linewidth',2)
plot(first_mom(l_23_pyr_fro),second_mom(l_23_pyr_fro),'go','linewidth',2)
plot(first_mom(l_5_pyr_fro),second_mom(l_5_pyr_fro),'ko','linewidth',2)
plot(first_mom(l_23_pyr_pre),second_mom(l_23_pyr_pre),'co','linewidth',2)
plot(first_mom(l_5_pyr_pre),second_mom(l_5_pyr_pre),'o','Color',[0.8 0.3 0],'linewidth',2)
legend('L-23-Par','L-56-Par','L-23-Fro','L-5-Fro','L-23-Pre','L-5-Pre')


range = linspace(-0.5, 0.5, 100);
figure
stairs(range,hist(first_mom(l_23_pyr_par),range)/length(l_23_pyr_par))
hold on
stairs(range,hist(first_mom(l_56_pyr_par),range)/length(l_56_pyr_par),'Color','r')
stairs(range,hist(first_mom(l_23_pyr_fro),range)/length(l_23_pyr_fro),'Color','g')
stairs(range,hist(first_mom(l_5_pyr_fro),range)/length(l_5_pyr_fro),'Color','k')
stairs(range,hist(first_mom(l_23_pyr_pre),range)/length(l_23_pyr_pre),'Color','c')
stairs(range,hist(first_mom(l_5_pyr_pre),range)/length(l_5_pyr_pre),'Color',[0.8 0.3 0])

