function intVector = BinToInt(msg,m,n)
conv1 = zeros(n,1);
base2 = zeros(1,m);
    
for j = 1 : m
    base2(1,j) = 2^(m-j);
end
for i = 1 : n
    conv1(i,1) = base2*msg(i,:)';
end

intVector = conv1;
    