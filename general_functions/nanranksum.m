function p = nanranksum(data1,data2)

p = ranksum(data1(~isnan(data1)),data2(~isnan(data2)));