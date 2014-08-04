function cur_dist = new_cosine_dist(x,y)

cur_dist = 1-abs(dot(x,y))/(norm(x)*norm(y));