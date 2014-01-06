function [d] = distanciaEuclidiana(r , y)

k = 0;
N = length(r);

for i = 1 : N
    k = k + r(i)*y(i);
end

d=k;
	
	

