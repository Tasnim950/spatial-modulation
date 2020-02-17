clear all 
clc

Nt = 4;                       % # of transmit antennas
Nr = 4;                       % # of receive antennas
M = 8;                        % # of signal symbols
Nu = 3;                       % # of active transmit antenna
cmb = nchoosek(Nt, Nu);       % total # of antenna combinations

q = log2(M);                  % # of bits per signal symbol
N = floor(log2(cmb));         % # of bits per spatial symbol
Nc = 2.^N;                    % # of available combinations for transmission
TotalBits = N + q;

Txindx = nchoosek(1:Nt, Nu);  % antenna combination matrix

snr = 0:2:20;
snrL = 10.^(snr./10);         % linear snr
K = length(snr);

loops = 10e6;
ber = zeros(1,K);

Sym = qammod(0:M-1, M, 'gray', 'UnitAveragePower', true);

for ii = 1:K
    
    nv = 1./snrL(ii)./2;       % noise variance
    run_loops = loops;
    
    for jj = 1:loops
        %% Generating signal and spatial symbols
        
        BitsIn = randsrc(1, TotalBits, [0 1; 0.5 0.5]);
        
        CbNo = bi2de(BitsIn(1:N), 'left-msb')+1;
        SymNo = bi2de(BitsIn(N+1:end), 'left-msb')+1;
        
        s = qammod(SymNo, M, 'gray', 'UnitAveragePower', true);
        
        %% Channel Model
        
        H = sqrt(0.5).* (randn(Nr, Nt)+randn(Nr, Nt).*1i);
        nf = sum(sum(H.*conj(H))) /Nt /Nr;           % normalization factor
        H = H./sqrt(nf);
        
        noise = sqrt(nv*Nu).* (randn(Nr, 1)+randn(Nr, 1).*1i);
        
        h = sum(H(:,Txindx(CbNo,:)), 2);
        
        y = h*s + noise;                             % recieved signal
        
        %% maximum liklihood estimation for GSM
        
        minV = inf;
         for n = 1: Nc
             hc = sum(H(:,Txindx(n,:)),2);
           for k = 1:M
                temp = y - hc.*Sym(k);
                d = temp'*temp;
                if d < minV
                    minV = d;
                    EstSym = k;
                    EstCb = n;
                end
           end
           
        end
        
        
        %% GSM Demapping
        
        DemCb = de2bi(EstCb-1, N, 'left-msb');
        DemSym = qamdemod(Sym(EstSym), M, 'gray', 'OutputType', 'bit', 'UnitAveragePower', true);
        
        BitsOut = [DemCb DemSym.'];
        
        %% BER calculation
        
        ber(ii) = ber(ii) + sum(abs(BitsIn-BitsOut));
        if ber(ii) > 2000
            run_loops=jj;
 
            txt=['snr = ' num2str(snr(ii)), ' ber = ' num2str(ber(ii)./run_loops./TotalBits)];
            disp(txt);
 
            break;
         end
     end
  ber(ii)= ber(ii)./run_loops./TotalBits;
    
end

figure

semilogy(snr, ber, 'r+--', 'linewidth', 0.5);
grid on
 
xlabel('SNR (dB)');
ylabel('Average BER');
title('Generalised Spatial Modualtion');
