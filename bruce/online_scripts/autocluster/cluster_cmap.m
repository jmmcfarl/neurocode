function cmap = cluster_cmap(n_colors)

cmap(5,:) = [33 64 154]; %darkblue
cmap(1,:) = [8 135 67]; %darkgreen
cmap(2,:) = [241 140 34]; %orange
cmap(3,:) = [150 100 155]; %purple
cmap(4,:) = [71 195 211]; %lightblue
cmap(6,:) = [173 209 54]; %lightgreen
cmap(7,:) = [238 132 181]; %pink
cmap(8,:) = [255 222 23]; %yellow
cmap(9,:) = [127.5 0 0]; %darkred

cmap = cmap(1:n_colors,:)/256;
