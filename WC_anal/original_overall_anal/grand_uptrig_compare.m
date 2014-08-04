clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_data
load C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_midpts

%% lump LEC cells
cell_type(cell_type==5) = 4;
cell_type(cell_type==6)=4;
cell_type(cell_type==7) = 5;
cell_type(cell_type==8) = 6;

%%

diff_types = unique(cell_type);
time_vec = lf8_dtrig_lf2_mid;
time_range = -1200:200;
for i = 1:length(diff_types)
    
    y(i,:) = hist(1000*lags(time_vec(cell_type==diff_types(i))),time_range);
    y(i,:) = y(i,:)/sum(y(i,:));
    y(i,:) = cumsum(y(i,:));
    
end
total_f = size(y,1);
clear cell_type diff_types time_vec

load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cor_dir.mat')
load C:\WC_Germany\overall_calcs\HC_Cort\trig_avgs\trig_avg_data
load C:\WC_Germany\overall_calcs\HC_Cort\trig_avgs\trig_avg_midpts
diff_types = unique(cell_type);

time_vec = lf8_dtrig_lf2_mid;
for i = 1:length(diff_types)
    
    y(i+total_f,:) = hist(1000*lags(time_vec(cell_type==diff_types(i))),time_range);
    y(i+total_f,:) = y(i+total_f,:)/sum(y(i+total_f,:));
    y(i+total_f,:) = cumsum(y(i+total_f,:));
    
end

%% Get rid of CA pyrs
y(total_f+2,:) = [];
y(total_f+3,:) = [];

%%
num_cell_types = size(y,1);
cmap = colormap(jet(num_cell_types));
for i = 1:num_cell_types
plot(time_range,y(i,:),'color',cmap(i,:),'linewidth',2)
hold on
end
ylim([0 1])
%%
legend('MEC 3','MEC 2','MEC 5','LEC','MEC 3atr','MEC 2atr','Cort','CA1 int','DG','Location','NorthWest')
