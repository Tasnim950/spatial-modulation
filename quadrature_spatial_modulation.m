clear all
clc

Nt = 4;                         % # of transmit antennas
Nr = 4;                         % # of receive antennas
M = 4;                          % # of signal symbols

q = log2(M);                    % # of bits per symbol
N = log2(Nt);                   % # of bits per antenna 
TotalBits = log2(M*Nt^2);

snr = 0:2.5:25;
snrL = 10.^(snr./10);           % linear snr
K = length(snr);

loops = 10e5;

ber = zeros(1, K);

Sym = qammod(0:M-1, M, 'gray', 'UnitAveragePower', true);


for ii = 1:K
    
    nv = 1./snrL(ii)./2;        % noise variance
    
    run_loops = loops;
    
    for jj = 1:loops
        %% Generating signal and spatial symbols
        
        BitsIn = randsrc(1, TotalBits, [0 1; 0.5 0.5]);        % generate 0 1 values with 0.5 probability

        i1 = bi2de(BitsIn(1:N), 'left-msb') + 1;               % in-phase spatial index
        i2 = bi2de(BitsIn(N+1:2*N), 'left-msb') + 1;           % quadrature spatial index
        SymNo = bi2de(BitsIn(2*N+1:end), 'left-msb') + 1;      % signal index
     
        s = qammod(SymNo, M, 'gray', 'UnitAveragePower', true);
        
      
       %% Channel Model 
      
       H = sqrt(0.5).*(randn(Nr, Nt) + randn(Nr, Nt).*1i);       % channel
       nf = sum(sum(H.*conj(H))) /Nt /Nr ;                       % normalization factor
       H = H./sqrt(nf);
       
       noise = sqrt(nv).*(randn(Nr, 1) + randn(Nr, 1).*1i);
      
       y = H(:,i1).*real(s) + H(:,i2).*imag(s).*1i + noise;     % recieved signal
       
       
      %% maximum likelihood estimation for QSM
      
      minV = inf;
      for mm = 1:M
          sr = real(Sym(mm));                                    
          si = imag(Sym(mm));
          for nn = 1:Nt
              hsr = H(:,nn)*sr;
              for kk = 1:Nt
                  hsi = H(:,kk)*si.*1i;
                 
                  temp = (y - (hsr+hsi));
                  d = temp' * temp;
                  
                  if d < minV
                      minV = d; 
                      EstSym = mm;
                      EstAnt1 = nn;
                      EstAnt2 = kk;
                  end
              end
          end 
      end
      
      
      
      %% QSM Demapping 
      
      DemAnt1 = de2bi(EstAnt1-1, N, 'left-msb');
      DemAnt2 = de2bi(EstAnt2-1, N, 'left-msb');
      DemSym = qamdemod(Sym(EstSym), M, 'gray','OutputType','bit', 'UnitAveragePower', true);
      
      BitsOut = [ DemAnt1 DemAnt2 DemSym.'];
      
      
      %% BER calculation 
      
      ber(ii) = ber(ii) + sum(abs(BitsIn - BitsOut));
      
      if ber(ii) > 2000
          run_loops = jj;
          
          txt = ['snr = ', num2str(snr(ii)), ' ber = ', num2str(ber(ii)./run_loops./TotalBits)];
          disp(txt);
          break;
      end
    end
    ber(ii) = ber(ii)./run_loops./TotalBits;
end

% figure
semilogy(snr, ber, 'k+--', 'linewidth', 1);  
    
    grid on
    
    xlabel ('SNR (dB)');
    ylabel ('Average BER');
    title  ('Quadrature Spatial Modulation');
