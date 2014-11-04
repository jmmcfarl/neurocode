clear all
cd C:\WC_Germany\sven_thomas_combined
load ./border_mesh_data.mat
%%

dx=0.05;
dy=0.05;

x = data(:,1);
y = data(:,2);
x_edge=[floor(min(x)):dx:ceil(max(x))];
x_edge=[floor(min(x)):dx:-4];
y_edge=[floor(min(y)):dy:ceil(max(y))];

edge_alpha = 0.1;
face_alpha = 0.1;
% z=data(:,3);
% [X,Y]=meshgrid(x_edge,y_edge);
% Z=griddata(x,y,z,X,Y);
% h=mesh(X,Z,Y,'edgecolor','black','facecolor','black');
% set(h,'edgealpha',edge_alpha,'facealpha',face_alpha)
% hold on
% z=data(:,4);
% [X,Y]=meshgrid(x_edge,y_edge);
% Z=griddata(x,y,z,X,Y);
% h=mesh(X,Z,Y,'edgecolor','black','facecolor','black')
% set(h,'edgealpha',edge_alpha,'facealpha',face_alpha)
% % 

z=data(:,4);
missing = find(isnan(z));
z(missing) = data(missing,5);
[X,Y]=meshgrid(x_edge,y_edge);
Z=griddata(x,y,z,X,Y);
h=mesh(X,Z,Y,'edgecolor','black','facecolor','black')
hold on
set(h,'edgealpha',edge_alpha,'facealpha',face_alpha)
% 


% z=data(:,5);
% [X,Y]=meshgrid(x_edge,y_edge);
% use = X < -4;
% Z=griddata(x,y,z,X,Y);
% h = mesh(X,Z,Y,'edgecolor','red','facecolor','red')
% set(h,'edgealpha',edge_alpha,'facealpha',face_alpha)
% z=data(:,6);
% [X,Y]=meshgrid(x_edge,y_edge);
% Z=griddata(x,y,z,X,Y);
% h = mesh(X,Z,Y,'edgecolor','red','facecolor','red')
% set(h,'edgealpha',edge_alpha,'facealpha',face_alpha)

%%
hold on
ylim([2.5 4.5])
 xlim([-5.3 -3.4])
  zlim([-5.5 -2])
load ./anatomical_data.mat
load ./combined_dir_nd.mat
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
l3mec_np(l3mec_np==62) =[];
l3mec_np(l3mec_np == 59) = [];
l3mec_np(l3mec_np == 60) = [];
scatter3(anatomy(l3mec,3),anatomy(l3mec,5),anatomy(l3mec,4),'r','linewidth',2);
scatter3(anatomy(l3lec,3),anatomy(l3lec,5),anatomy(l3lec,4),'b','linewidth',2);
% scatter3(anatomy(l3mec_np,3),anatomy(l3mec_np,5),anatomy(l3mec_np,4),'k','linewidth',2);
% scatter3(anatomy(l3lec_np,3),anatomy(l3lec_np,5),anatomy(l3lec_np,4),'g','linewidth',2);
xlabel('Rostrocaudal','fontsize',16)
zlabel('Dorsoventral','fontsize',16)
ylabel('Mediolateral','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')

%%
load ./combined_core_analysis_fin_nd_np

int = 62;
offset = 4;
scale = 1200;
scatter3(anatomy(l3mec,3),anatomy(l3mec,5),anatomy(l3mec,4),scale*fract_rt2_ups(l3mec)+offset,'r','linewidth',2);
hold on
scatter3(anatomy(l3lec,3),anatomy(l3lec,5),anatomy(l3lec,4),scale*fract_rt2_ups(l3lec)+offset,'b','linewidth',2);
% scatter3(anatomy(l3mec_np,3),anatomy(l3mec_np,5),anatomy(l3mec_np,4),scale*fract_rt2_ups(l3mec_np)+offset,'k','linewidth',2);
% scatter3(anatomy(l3lec_np,3),anatomy(l3lec_np,5),anatomy(l3lec_np,4),scale*fract_rt2_ups(l3lec_np)+offset,'g','linewidth',2);
% scatter3(anatomy(int,3),anatomy(int,4),anatomy(int,5),scale*fract_rt2_ups(int)+offset,'c','linewidth',2);
xlabel('Rostrocaudal','fontsize',16)
zlabel('Dorsoventral','fontsize',16)
ylabel('Mediolateral','fontsize',16)
set(gca,'fontsize',12,'fontname','arial')
% ylim([2.5 4.5])
%  xlim([-5.3 -3.4])
%   zlim([-5.5 -2])

int = 62;
offset = 4;
scale = 1200;
figure
scatter3(anatomy(l3mec,3),anatomy(l3mec,5),anatomy(l3mec,4),scale*fract_rt2_downs(l3mec)+offset,'r','linewidth',2);
hold on
scatter3(anatomy(l3lec,3),anatomy(l3lec,5),anatomy(l3lec,4),scale*fract_rt2_downs(l3lec)+offset,'b','linewidth',2);
% scatter3(anatomy(l3mec_np,3),anatomy(l3mec_np,5),anatomy(l3mec_np,4),scale*fract_rt2_ups(l3mec_np)+offset,'k','linewidth',2);
% scatter3(anatomy(l3lec_np,3),anatomy(l3lec_np,5),anatomy(l3lec_np,4),scale*fract_rt2_ups(l3lec_np)+offset,'g','linewidth',2);
% scatter3(anatomy(int,3),anatomy(int,4),anatomy(int,5),scale*fract_rt2_ups(int)+offset,'c','linewidth',2);
xlabel('Rostrocaudal','fontsize',16)
zlabel('Dorsoventral','fontsize',16)
ylabel('Mediolateral','fontsize',16)
set(gca,'fontsize',12,'fontname','arial')

% scatter(anatomy(l3mec,3),anatomy(l3mec,4),scale*fract_rt2_ups(l3mec)+offset,'r','linewidth',2);
% hold on
% scatter(anatomy(l3lec,3),anatomy(l3lec,4),scale*fract_rt2_ups(l3lec)+offset,'b','linewidth',2);

%%
figure
plot(anatomy(l3mec,1),100*fract_rt2_ups(l3mec),'ro')
hold on
plot(anatomy(l3lec,1),100*fract_rt2_ups(l3lec),'o')
set(gca,'fontsize',14,'fontname','arial')
xlabel('Depth from pia (um)','fontsize',16)
ylabel('Percentage of persistent up states','fontsize',16)


