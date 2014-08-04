function [] = visualize_movie(X,sdim,bounds,framerate)

dt = 1/framerate;

figure
for tslice = bounds(1):bounds(2)
	plot2drfmat(reshape(X(tslice,:),sdim(1),sdim(2)),[-1.,1.]); 
    pause(dt)
end

