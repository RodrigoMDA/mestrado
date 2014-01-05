% Transmissão digital - Modulação BPSK - Canal AWGN + Desvanecimento
% Rayleigh - Detecção coerente
% Código RS(31,27) - Soft-Decoded-> Algoritmo de Chase 

clear all;
close all;
clc;
% --------- %
L = 10^5;           % Palavras de k simbolos que serão geradas. 
% --------- %
% Parâmetros do código RS
m = 5;
p = 2^m;
n = p - 1;
k = 27;
tx = k/n;
% ---------- %
% Parâmetros gerais da simulação.
M = n*m;
bits_total = L*m*k;
EbdB = 10;
ber = zeros(1, length(EbdB));         
% --------- %

for i = 1 :length(EbdB)
    
    EcdB = EbdB(i) + 10*log10(tx);
    Es = 10^(EcdB/10);
    cont = 0;                       
    
    for j = 1 : L
        
        msgTx = randi([0 n],k,1);  % Gera vetor (1,k) de numeros inteiros
        TxBin = IntToBin(msgTx,m,k);
        msgTxBin = reshape(TxBin',1,k*m);
        
        coder = fec.rsenc(n,k);     % Cria um objeto codificador RS(n,k) 
        code = encode(coder, msgTx);% Codifica a mensagem
        
        codeBin = IntToBin(code,m,n);
        
        modulacao = zeros(n,m);
        msg_ruido = zeros(n,m);
        msg_dem = zeros(n,m);
        h = zeros(n,m);
        ruido = zeros(n,m);
        sgma = 1;
        
        % Broadcasting
        for ii = 1 : n
            for jj = 1 : m
        
                if codeBin(ii,jj) == 0          % Modulação BPSK
                    modulacao(ii,jj) = sqrt(Es);
                else
                    modulacao(ii,jj) = -sqrt(Es);
                end
                           
                h(ii,jj) = (1/sqrt(2))*(randn(1) + 1i*randn(1));        % Desvanecimento              
                ruido(ii,jj) = sgma*(1/sqrt(2))*(randn(1) + 1i*randn(1));   % Ruido AWGN
              
                msg_ruido(ii,jj) = h(ii,jj)*modulacao(ii,jj) + ruido(ii,jj);
        
                equ = msg_ruido(ii,jj)/h(ii,jj);
                
                if  (real(equ) > 0)        % Hard decision
                    msg_dem(ii,jj) = 0;
                else
                    msg_dem(ii,jj) = 1;
                end
            end
        end
        % ------------------------------------ %
        % Pré-Processamento
        h_conj = conj(h);
        yd = real(h_conj.*msg_ruido);
        % ------------------------------------ %
        
        r = reshape(yd',1,M);
        Y = reshape(msg_dem',1,M);
        t = 5;
        [valMenosConf, posConf] = getNElements(abs(r),t);
        Z = chaseCurto(posConf,r,Y,t,coder);
                
        msgRx = Z(1:k*m);
        
        for ii = 1 : k*m
            if msgRx(ii) ~= msgTxBin(ii)
                cont = cont + 1;
            end
            
        end
        
    end
    ber(i) = cont / bits_total;
    display(ber(i));
end
