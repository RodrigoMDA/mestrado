% Função que converte um vetor de numeros inteiros/double para binário.
% A função do matlab, dec2bin, retorna uma string, esta 0's e 1's do tipo
% double.

function binVector = IntToBin(msg, m, n)
conv1 = dec2bin(msg,m);
conv2 = zeros(n,m);

for i = 1 : n
    for j = 1 : m
        conv2(i,j) = conv1(i,j) - '0';
    end
end
binVector = conv2;
