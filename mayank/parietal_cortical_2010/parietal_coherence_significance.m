clear all
close all
cd F:\WC_Germany\parietal_cortical_2010\
load state_dep_lf5spectra_no_trans 
n = length(sess_data);

n_tests = 7;
p = 0.05;
bandwidth = 2*2/0.6;

n_up_greater = zeros(n,length(f));
n_down_greater = zeros(n,length(f));

for d = 1:n

    up_bias = 1/(dof_wup(d)-2);
    down_bias = 1/(dof_wdown(d)-2);
    
    y_up(d,:) = atanh(Cw8_wup(d,:)) - up_bias;
    y_down(d,:) = atanh(Cw8_wdown(d,:)) - down_bias;
    x_up(d,:) = squeeze(log(S_up(d,:,1)) + psi(dof_wup(d)/2) - log(dof_wup(d)/2));
    x_down(d,:) = squeeze(log(S_down(d,:,1)) + psi(dof_wdown(d)/2) - log(dof_wdown(d)/2));
    pow_diff(d,:) = x_up(d,:)-x_down(d,:);
    %     plot(f,y_up), hold on
%     plot(f,y_down,'r')
    
    deltay = (y_up(d,:) - y_down(d,:))/sqrt(up_bias + down_bias);
    
%     plot(f,deltay)
%     line([f(1) f(end)],[norminv(p/2/n_tests) norminv(p/2/n_tests)])
%     line([f(1) f(end)],[norminv(1-p/2/n_tests) norminv(1-p/2/n_tests)])
%     xlim([0 45])
%     
%     pause
%     clf
    
    n_up_greater(d,:) = deltay > norminv(1-p/2);
    n_down_greater(d,:) = deltay < norminv(p/2);
    
    up_start = find(n_up_greater(d,2:end) & ~n_up_greater(d,1:end-1));
    if n_up_greater(d,1)
        up_start = [1 up_start];
    end
    up_stop = find(~n_up_greater(d,2:end) & n_up_greater(d,1:end-1));
     if n_up_greater(d,end)
        up_stop = [up_stop length(f)];
     end
    n_ups = length(up_start);
    is_up = zeros(size(f));
    for i = 1:n_ups
        if up_stop(i)-up_start(i) > bandwidth
            is_up(up_start(i):up_stop(i)) = 1;
        end
    end
    n_up_greater(d,:) = is_up;
 
    down_start = find(n_down_greater(d,2:end) & ~n_down_greater(d,1:end-1));
    if n_down_greater(d,1)
        down_start = [1 down_start];
    end
    down_stop = find(~n_down_greater(d,2:end) & n_down_greater(d,1:end-1));
    if n_down_greater(d,end)
        down_stop = [down_stop length(f)];
    end
    n_downs = length(down_start);
    is_down = zeros(size(f));
    for i = 1:n_downs
        if down_stop(i)-down_start(i) > bandwidth
            is_down(down_start(i):down_stop(i)) = 1;
        end
    end
    n_down_greater(d,:) = is_down;

end

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
superficial = find_struct_field_vals(sess_data,'layer','23');

sup_par = 1:11;
deep_par = 12:21;
sup_fro = intersect(frontal,superficial);
sup_pre = intersect(prefrontal,superficial);
sup_pfc = union(sup_fro,sup_pre);
deep_pfc = setdiff(22:n,sup_pfc);

figure
plot(f,mean(n_up_greater(sup_par,:))), hold on
plot(f,mean(n_down_greater(sup_par,:)),'r')
xlabel('Frequency (hz)','fontsize',16)
ylabel('Fraction of data','fontsize',16)
legend('Up state greater','Down state greater')
title('superficial parietal')

figure
plot(f,mean(n_up_greater(deep_par,:))), hold on
plot(f,mean(n_down_greater(deep_par,:)),'r')
xlabel('Frequency (hz)','fontsize',16)
ylabel('Fraction of data','fontsize',16)
legend('Up state greater','Down state greater')
title('deep parietal')

figure
plot(f,mean(n_up_greater(sup_pfc,:))), hold on
plot(f,mean(n_down_greater(sup_pfc,:)),'r')
xlabel('Frequency (hz)','fontsize',16)
ylabel('Fraction of data','fontsize',16)
legend('Up state greater','Down state greater')
title('superficial pfc')

figure
plot(f,mean(n_up_greater(deep_pfc,:))), hold on
plot(f,mean(n_down_greater(deep_pfc,:)),'r')
xlabel('Frequency (hz)','fontsize',16)
ylabel('Fraction of data','fontsize',16)
legend('Up state greater','Down state greater')
title('deep pfc')

%%
figure
h=errorbar(f,mean(y_up(sup_par,:)),std(y_up(sup_par,:))/sqrt(length(sup_par)));
hold on
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(y_down(sup_par,:)),std(y_down(sup_par,:))/sqrt(length(sup_par)),'r');
errorbar_tick(h,.01,'units')
legend('Up state','Down state')
xlim([f(1) f(end)])
ylim([0 0.25])
xlabel('Frequency','fontsize',16)
ylabel('Coherence','fontsize',16)

figure
h=errorbar(f,mean(y_up(deep_par,:)),std(y_up(deep_par,:))/sqrt(length(deep_par)));
hold on
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(y_down(deep_par,:)),std(y_down(deep_par,:))/sqrt(length(deep_par)),'r');
errorbar_tick(h,.01,'units')
legend('Up state','Down state')
xlim([f(1) f(end)])
ylim([0 0.25])
xlabel('Frequency','fontsize',16)
ylabel('Coherence','fontsize',16)

figure
h=errorbar(f,mean(y_up(sup_pfc,:)),std(y_up(sup_pfc,:))/sqrt(length(sup_pfc)));
hold on
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(y_down(sup_pfc,:)),std(y_down(sup_pfc,:))/sqrt(length(sup_pfc)),'r');
errorbar_tick(h,.01,'units')
legend('Up state','Down state')
xlim([f(1) f(end)])
ylim([0 0.25])
xlabel('Frequency','fontsize',16)
ylabel('Coherence','fontsize',16)

figure
h=errorbar(f,mean(y_up(deep_pfc,:)),std(y_up(deep_pfc,:))/sqrt(length(deep_pfc)));
hold on
errorbar_tick(h,.01,'units')
h=errorbar(f,mean(y_down(deep_pfc,:)),std(y_down(deep_pfc,:))/sqrt(length(deep_pfc)),'r');
errorbar_tick(h,.01,'units')
legend('Up state','Down state')
xlim([f(1) f(end)])
ylim([0 0.25])
xlabel('Frequency','fontsize',16)
ylabel('Coherence','fontsize',16)

%%
%if you want to plot the up and down state power differences in the MP
figure
h=errorbar(f,mean(pow_diff(sup_par,:)),std(pow_diff(sup_par,:))/sqrt(length(sup_par)));
errorbar_tick(h,.01,'units')
xlim([f(1) f(end)])
% ylim([0 0.25])
xlabel('Frequency','fontsize',16)
ylabel('Coherence','fontsize',16)

figure
h=errorbar(f,mean(pow_diff(deep_par,:)),std(pow_diff(deep_par,:))/sqrt(length(deep_par)));
errorbar_tick(h,.01,'units')
xlim([f(1) f(end)])
% ylim([0 0.25])
xlabel('Frequency','fontsize',16)
ylabel('Coherence','fontsize',16)
% 
% figure
% h=errorbar(f,mean(pow_diff(sup_pfc,:)),std(pow_diff(sup_pfc,:))/sqrt(length(sup_pfc)));
% errorbar_tick(h,.01,'units')
% xlim([f(1) f(end)])
% % ylim([0 0.25])
% xlabel('Frequency','fontsize',16)
% ylabel('Coherence','fontsize',16)
% 
% figure
% h=errorbar(f,mean(pow_diff(deep_pfc,:)),std(pow_diff(deep_pfc,:))/sqrt(length(deep_pfc)));
% errorbar_tick(h,.01,'units')
% xlim([f(1) f(end)])
% % ylim([0 0.25])
% xlabel('Frequency','fontsize',16)
% ylabel('Coherence','fontsize',16)

%%
% later_pfc = [];
% later_par = [];
% earlier_pfc = [];
% earlier_par = [];
% for i = 1:n
%     if sess_data(i).name(4) == '7'
%         if strcmp(sess_data(i).region, 'parietal')
%             later_par = [later_par i];
%         else
%             later_pfc = [later_pfc i];
%         end
%     else
%         if strcmp(sess_data(i).region, 'parietal')
%             earlier_par = [earlier_par i];
%         else
%             earlier_pfc = [earlier_pfc i];
%         end
%     end
% end
% 
% figure
% h=errorbar(f,mean(y_up(later_pfc,:)),std(y_up(later_pfc,:))/sqrt(length(later_pfc)));
% hold on
% errorbar_tick(h,.01,'units')
% h=errorbar(f,mean(y_down(later_pfc,:)),std(y_down(later_pfc,:))/sqrt(length(later_pfc)),'r');
% errorbar_tick(h,.01,'units')
% legend('Up state','Down state')
% xlim([f(1) f(end)])
% ylim([0 0.25])
% xlabel('Frequency','fontsize',16)
% ylabel('Coherence','fontsize',16)
% 
% figure
% h=errorbar(f,mean(y_up(later_par,:)),std(y_up(later_par,:))/sqrt(length(later_par)));
% hold on
% errorbar_tick(h,.01,'units')
% h=errorbar(f,mean(y_down(later_par,:)),std(y_down(later_par,:))/sqrt(length(later_par)),'r');
% errorbar_tick(h,.01,'units')
% legend('Up state','Down state')
% xlim([f(1) f(end)])
% ylim([0 0.25])
% xlabel('Frequency','fontsize',16)
% ylabel('Coherence','fontsize',16)
% 
% figure
% h=errorbar(f,mean(y_up(earlier_pfc,:)),std(y_up(earlier_pfc,:))/sqrt(length(earlier_pfc)));
% hold on
% errorbar_tick(h,.01,'units')
% h=errorbar(f,mean(y_down(earlier_pfc,:)),std(y_down(earlier_pfc,:))/sqrt(length(earlier_pfc)),'r');
% errorbar_tick(h,.01,'units')
% legend('Up state','Down state')
% xlim([f(1) f(end)])
% ylim([0 0.25])
% xlabel('Frequency','fontsize',16)
% ylabel('Coherence','fontsize',16)
% 
% figure
% h=errorbar(f,mean(y_up(earlier_par,:)),std(y_up(earlier_par,:))/sqrt(length(earlier_par)));
% hold on
% errorbar_tick(h,.01,'units')
% h=errorbar(f,mean(y_down(earlier_par,:)),std(y_down(earlier_par,:))/sqrt(length(earlier_par)),'r');
% errorbar_tick(h,.01,'units')
% legend('Up state','Down state')
% xlim([f(1) f(end)])
% ylim([0 0.25])
% xlabel('Frequency','fontsize',16)
% ylabel('Coherence','fontsize',16)

% figure
% h=errorbar(f,mean(pow_diff(later,:)),std(pow_diff(later,:))/sqrt(length(later)));
% errorbar_tick(h,.01,'units')
% xlim([f(1) f(end)])
% % ylim([0 0.25])
% xlabel('Frequency','fontsize',16)
% ylabel('Coherence','fontsize',16)
% 
% figure
% h=errorbar(f,mean(pow_diff(earlier,:)),std(pow_diff(earlier,:))/sqrt(length(earlier)));
% errorbar_tick(h,.01,'units')
% xlim([f(1) f(end)])
% % ylim([0 0.25])
% xlabel('Frequency','fontsize',16)
% ylabel('Coherence','fontsize',16)
