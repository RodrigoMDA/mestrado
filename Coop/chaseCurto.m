function [Z] = chaseCurto(posConf,r,Y,t,coder)

m = 5;
n = 2^m - 1;
k = 27;
M = n*m;

decoder = fec.rsdec(coder);

PTESTE = 2^t;

% Criando os 2^t padroes de teste 
pTeste = zeros(PTESTE,M);
posicoes = sort(posConf, 'descend');

comb = combn([0 1],t);

for i = 1 : PTESTE
    for j = 1:t
        h = posicoes(j);
        pTeste(i,h) = comb(i,t+1-j);
    end
end

% ----- %
% Gerar o conjunto de palavras código candidatas %
Yt = zeros(PTESTE,M);
lambda = zeros(PTESTE,M);

% Até a linha 74 gera o conjunto de palavras-código candidatas, as testa e
% gera o conjunto lambda.

for i = 1 : PTESTE
    Yt(i,:) = mod(Y+pTeste(i,:),2);
    % Soma módulo-2 da saida do demodulador Y, com o padrão de teste
    
    % Até a linha 40: Conversão de Yt, para simbolos m-arios
    Ypal = reshape(Yt(i,:),m,n); 
    Ycol = Ypal';
    
    Ydeci = BinToInt(Ycol,m,n);
    % Fim da conversão
    
    [decoded,cnumerr,ccode] = decode(decoder,Ydeci);
    % Decodificação RS
    
    % cnumerr informa se o decodificador conseguiu decodificar a palavra
    % Ydeci. Quando cnumerr = -1, não foi possível decodificar.
    if (cnumerr == -1) 
        for j = 1 : M
            lambda(i,j) = -1;
        end
    else
        % Caso o decodificador tenha decodificado: 
        ccodeBin = IntToBin(ccode,m,n);
        lambda(i,:) = reshape(ccodeBin',1,M); 
              
    end
    
end

% Seleciona a palavra código adequada dentre as palavaras do conjunto
% Lambda

minimo = 100000;
C = zeros(1,M);

for i = 1 : PTESTE
    if (lambda(i,1) ~= -1)
       C = lambda(i,:);
       dist = distanciaEuclidiana(r, C);
       
       if (dist < minimo)
           minimo = dist;
           Z = C;
       end
    end
end