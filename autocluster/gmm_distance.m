function distance  = gmm_distance(G)

try
    D = mahal(G,G.mu);
    distance = sqrt(2./((1./D(1,2))+(1./D(2,1))));
catch %in case fit failed
    distance = 0;
end
