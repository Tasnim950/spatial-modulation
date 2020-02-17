clear all
clc

 
Nt = 4;                                                 % # of transmit antennas
Nr = 4;                                                 % # of receive antennas              
M = 4;                                                  % # of symbols in the modulation set
 
q = log2(M);                                            % # of bits per signal symbol
N = log2(Nt);                                           % # of bits per spatial symbol
TotalBits = log2(M*Nt);
 
 
snr = 0:2:20;                                           % dB SNR
snrL = 10.^(snr./10);
K = length(snr);
 
loops = 10e5;
ber = zeros(1,K);
 
Sym = qammod(0:M-1, M, 'gray', 'UnitAveragePower', true);
 

for ii = 1 : K
 
    nv = 1./snrL(ii)./2;                                  % noise variance
    run_loops = loops;
    
    for jj = 1 : loops
        %% Generating signal and spatial symbols
 
        BitsIn = randsrc(1,TotalBits, [0 1; 0.5 0.5]);    % generate 0 1 values with 0.5 probability
 
        AntNo = bi2de(BitsIn(1:N), 'left-msb') + 1;       % spatial symbol index
        SymNo = bi2de(BitsIn(N+1:end), 'left-msb');       % signal symbol index
 
        s = qammod(SymNo, M, 'gray', 'UnitAveragePower', true);
 
        %% Channel Model 

        H = sqrt(0.5).*(randn(Nr, Nt) + randn(Nr, Nt).*1i);  % channel
        nf = sum(sum(H.*conj(H))) / Nt / Nr;               % normalizing factor
        H = H ./ sqrt(nf);
 
        noise = sqrt(nv).* (randn(Nr, 1) + randn(Nr, 1).*1i);
   
        y = H(:, AntNo)* s + noise;                        % recieved signal
     
        %% Maximum likelihood estimation for SM
        minV = inf;
        for i1 = 1:Nt
            h_est = H(:,i1);                               % spatial estimation
            for mm = 1:M
                s_est = Sym(mm);                           % signal estimation
                
                temp = (y - h_est.*s_est);
                d = temp' * temp;
%                 d = norm(y - h_est.*s_est)^2;
                
                if d < minV
                    minV = d;
                    EstAnt = i1;
                    EstSym = mm;
                end
            end
        end
 
       
 
      %% SM Demapping 
        DemAnt = de2bi(EstAnt-1, N, 'left-msb');
        DemSym = qamdemod(Sym(EstSym), M, 'gray', 'OutputType', 'bit', 'UnitAveragePower', true);
 
        BitsOut = [DemAnt DemSym.'];

        %% BER calculation 
        ber(ii) = ber(ii) + sum(abs(BitsIn - BitsOut));
 
        if ber(ii) > 2000
            run_loops = jj;
 
            txt = ['snr = ' num2str(snr(ii)), ' ber = ' num2str(ber(ii)./ run_loops./ TotalBits)];
            disp(txt);
            break;
        end
    end
 
     ber(ii) = ber(ii) ./ run_loops ./ TotalBits;
  
end
    
    figure
    semilogy(snr, ber, 'md--', 'linewidth', 1);
    
    grid on
 
    xlabel('SNR (dB)');
 
    ylabel('Average BER');
    
    title('Spatial Modualtion');
