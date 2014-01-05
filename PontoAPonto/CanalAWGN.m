% Transmissão digital utilizando modulação BPSK num canal com ruído AWGN.
% Note: o mapeamento ficou -sqrt(E) => 1. 
% Rodrigo Melo Donato, dezembro de 2013.

clear all;
close all;
clc;
% --------- %
L = 10^5;           % Palavras de k simbolos que serão geradas. 
% --------- %
EbdB = 2:1:8;
ber = zeros(1, length(EbdB));         
% --------- %

for i = 1 :length(EbdB)
    
    EcdB = EbdB(i);
    Es = 10^(EcdB/10);

    cont = 0;                       
    msg = rand(1,L)>0.5;
    s = zeros(1,L);
    
    for j = 1 : L
        if msg(j) == 1
            s(j) = -sqrt(Es);
        else
            s(j) = sqrt(Es);
        end
    end
        msg_ruido = zeros(1,L);
        msg_dem = zeros(1,L);
        ruido = zeros(1,L);
        sgma = 1;
        
        % Broadcasting
        for ii = 1 : L
             ruido(ii) = sgma*(1/sqrt(2))*(randn(1) + 1i*randn(1));   % Ruido AWGN
             msg_ruido(ii) = s(ii) + ruido(ii);
             y = msg_ruido(ii);
 
            if  (real(y) > 0)        % Hard decision
                msg_dem(ii) = 0;
            else
                msg_dem(ii) = 1;
            end
            
        end
        % Contagem de erros
        for ii = 1 : L
            if msg(ii) ~= msg_dem(ii)
                cont = cont + 1;
            end
            
        end
        
    ber(i) = cont / L;    
end

semilogy(EbdB, ber, 'cd-','LineWidth',2);
hold on

    
