clc
clear all
close all
% ---------- %
L = 10^5;             
% ---------- %
% Parâmetros do código RS
m = 5;
p = 2^m;
n = p - 1;
k = 27;
tx = k/n;  
% ---------- %
% Parâmetros gerais da simulação, arquivo de saída...

M = n*m;
bits_total = L*m*k;

%Probabilidade de erro de bits (BER)
ber=zeros(1,length(EtdB));
contadores = zeros(6,length(EtdB));
EbdB = 8: 2 : 12;

%Cria um arquivo txt de saída
arq=fopen('ber_Chase_paralelo-simetrico_8_12.txt','w'); 
    fprintf(arq,'EtdB \t ber \t qtdBits=%f \n', bits_total); 
fclose(arq);
% ---------- %

for i = 1 : length(EtdB)
    
    EcdB = EbdB(i) + log10(tx); % Energia por bit codificado em dB
    E = 10^(EcdB/10);        % Energia por símbolo, BPSK. 
    Es = E/2;           %Energia por bit da fonte 
    Er = E/2;           %Energia por bit do relay
    sgma=1;              %Variância - Link S-D / S-R / R-D
    
    % Estimativas de erros (contagem)
    num_err = 0;
    % Contagem das escolhas do decodificador suave
    rd_ch_win = 0;
    sd_ch_win = 0;
    rd_ch_aln = 0;
    sd_ch_aln = 0;
    rd_lrp_win = 0;
    sd_lrp_win = 0;
%--------------------------------------------------------------------------
%Gera mensagem no gf e codifica
    for j = 1 : L
        msgTx = randi([0 n],k,1);  % Gera vetor (1,k) de numeros inteiros
        TxBin = IntToBin(msgTx,m,k);
        msgTxBin = reshape(TxBin',1,k*m);
        
        coder = fec.rsenc(n,k);     % Cria um objeto codificador RS(n,k) 
        code = encode(coder, msgTx);% Codifica a mensagem
        
        codeBin = IntToBin(code,m,n);
             
%--------------------------------------------------------------------------
%Broadcasting
        modulacao_srd = zeros(n,m);
        
%S-R:
        msg_ruido_sr = zeros(n,m);
        h_sr = zeros(n,m);
        ruido_sr = zeros(n,m);
        
        msg_ruido_sd = zeros(n,m);
        h_sd = zeros(n,m);
        ruido_sd = zeros(n,m);

        for kk = 1 : n
            for ll = 1 : m
                if codeBin(kk,ll) == 1 
                    modulacao_srd(kk,ll) = -sqrt(Es);
                else
                    modulacao_srd(kk,ll) = sqrt(Es);
                end
                           
                h_sr(kk,ll) = (1/sqrt(2))*(randn(1) + 1i*randn(1));
                ruido_sr(kk,ll) = sgma*(1/sqrt(2))*(randn(1) + 1i*randn(1));
                msg_ruido_sr(kk,ll) = h_sr(kk,ll)*modulacao_srd(kk,ll) + ruido_sr(kk,ll);
           
%S-D:  
                h_sd(kk,ll) = (1/sqrt(2))*(randn(1) + 1i*randn(1));
                ruido_sd(kk,ll) = sgma*(1/sqrt(2))*(randn(1) + 1i*randn(1));
                msg_ruido_sd(kk,ll) = h_sd(kk,ll)*modulacao_srd(kk,ll) + ruido_sd(kk,ll);
                               
            end
        end
%--------------------------------------------------------------------------
%Segundo Time Slot:
%R-D:
        msg_ruido_rd = zeros(n,m);
        h_rd = zeros(n,m);
        ruido_rd = zeros(n,m);
        beta = zeros(n,m);
        rd = zeros(n,m);
        h_srd = zeros(n,m);
        
        for kk = 1 : n
            for ll = 1 : m
                
                h_rd(kk,ll) = (1/sqrt(2))*(randn(1) + 1i*randn(1));
                ruido_rd(kk,ll) = sgma*(1/sqrt(2))*(randn(1) + 1i*randn(1));
                
                beta(kk,ll) = sqrt( Er/( (abs(h_rd(kk,ll)^2))*Es + (sgma^2) ) );
                rd(kk,ll)=beta(kk,ll)*msg_ruido_sr(kk,ll);
                msg_ruido_rd(kk,ll) = h_rd(kk,ll)*rd(kk,ll) + ruido_rd(kk,ll);
                
                %%%%Calculo do h_srd %%%%
                h_srd(kk,ll) = sqrt( (Es*Er)/(Es*(abs(h_sr(kk,ll)^2)+ (sgma^2)) ))*h_sr(kk,ll)*h_rd(kk,ll);
                               
            end
        end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------        
%%%%Combinação no Destino%%%%

        hard_sd = zeros(n,m);
        hard_rd = zeros(n,m);
        
        for kk = 1 : n
            for ll = 1 : m
                equ_sd = msg_ruido_sd(kk,ll)/h_sd(kk,ll);
                if (real(equ_sd(kk,ll)) < 0)
                    hard_sd(kk,ll) = 0;
                else
                    hard_sd(kk,ll) = 1;
                end
                equ_rd = msg_ruido_rd(kk,ll);
                if (real(equ_rd(kk,ll)) < 0)
                    hard_rd(kk,ll) = 0;
                else
                    hard_rd(kk,ll) = 1;
                end
            end
        end
        % ------------------------------------ %
        % Pré-Processamento SD 
        hsd_conj = conj(h_sd);
        yd_sd = real(hsd_conj.*msg_ruido_sd);
        % ------------------------------------ %
        % Pré-Processamento RD
        hrd_conj = conj(h_srd);
        yd_rd = real(hrd_conj.*msg_ruido);
        % ------------------------------------ %
        r_sd = reshape(yd_sd',1,M);
        Y_sd = reshape(hard_sd',1,M);
        t = 4;
        [valMenosConf1, posConf1] = getNElements(abs(r_sd),t);
        Z_sd = chaseCurto(posConf1,r_sd,Y_sd,t,coder);
        
        r_rd = reshape(yd',1,M);
        Y_rd = reshape(msg_dem',1,M);
        t = 4;
        [valMenosConf2, posConf2] = getNElements(abs(r_rd),t);
        Z_rd = chaseCurto(posConf2,r_rd,Y_rd,t,coder);
        
        if (flagsd && flagrd)
            
            if Wsd < Wrd
                sd_ch_win = sd_ch_win +1;
                msgf_sd = Xest_sd(1:k*m);
                for kk = 1 : length(msgf_sd)
                    if (msgf_sd(kk) ~= msg_bin_line(kk))
                        num_err = num_err + 1;
                    end
                end
            else
                rd_ch_win = rd_ch_win +1;
                msgf_rd = Xest_rd(1:k*m);
                for kk = 1 : length(msgf_rd)
                    if (msgf_rd(kk) ~= msg_bin_line(kk))
                        num_err = num_err + 1;
                    end
                end
            end
        end
        if ((flagsd == 0) && (flagrd == 0))

            if Wsd > Wrd
                sd_lrp_win = sd_lrp_win +1;
                msgf_sd = Xest_sd(1:k*m);
                for kk = 1 : length(msgf_sd)
                    if (msgf_sd(kk) ~= msg_bin_line(kk))
                        num_err = num_err + 1;
                    end
                end
            else
                rd_lrp_win = rd_lrp_win +1;
                msgf_rd = Xest_rd(1:k*m);
                for kk = 1 : length(msgf_rd)
                    if (msgf_rd(kk) ~= msg_bin_line(kk))
                        num_err = num_err + 1;
                    end
                end
            end
        else
            if flagsd
                sd_ch_aln = sd_ch_aln +1;
                msgf_sd = Xest_sd(1:k*m);
                for kk = 1 : length(msgf_sd)
                    if (msgf_sd(kk) ~= msg_bin_line(kk))
                        num_err = num_err + 1;
                    end
                end
            else  
               rd_ch_aln = rd_ch_aln +1;
               msgf_rd = Xest_rd(1:k*m);
               for kk = 1 : length(msgf_rd)
                    if (msgf_rd(kk) ~= msg_bin_line(kk))
                        num_err = num_err + 1;
                    end
               end
            end
        end
        
        
    end
%------------Cálculo da BER (Probabilidades erro de bits)------------%
contadores(1,i) = sd_ch_win;contadores(2,i)=rd_ch_win;
contadores(3,i) = sd_ch_aln;contadores(4,i)=rd_ch_aln;
contadores(5,i) = sd_lrp_win;contadores(6,i)=rd_lrp_win;

ber(i)=num_err/bits_total;

%Imprimir as informações no arquivo txt     
arq=fopen('ber_Chase_paralelo-simetrico_8_12.txt','a'); 
    fprintf(arq,'%1.2f\t %1.4e\n', EtdB(i), ber(i)); 
fclose(arq);


end
