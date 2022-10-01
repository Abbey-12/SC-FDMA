%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Communication system %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SC-FDMA performance analaysis project%%%%%%%%%%%%%
%%%%%%%%%%%%% Bayleyegn Abebu Ademe & Wassie Solomon Fekadie %%%%%%%%%%%%%%%
%%%%%%%%%%%%%% feburary 17,2022                              %%%%%%%%%%%%%%
%%%%%%%%%%%%%% Unveristy of Trento, Italy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear;
close all;
%% Parameters 

N = 128;% The size of the transmitter IFFT and the receiver FFT (128/512)
M=8;
Q=N/M;% The bandwidth spreading factor
Q_D=Q-1;  % The bandwidth expnasion of distributed should be less than that of the interleaved 
cp =32 ; % CP length.
SP.subband =10;
% choose modulation type QAM or PSK
modulation_type='PSK';
order=2;

% Choose equalization type MMSE  (minimum mean squared-error equalizer)or ZFEQ (zero forcing equalizer)
equalizer_type='ZFEQ';
SNR = 0:1:30; % Simulated SNR range is from 0 dB to 30 dB.
N_simu = 10^5; % the number of simulation 

% Channels based on 3GPP TS 25.104.
channel = [1 10^(-9.7/20) 10^(-22.8/20)];
 channel = channel/sqrt(sum(channel.^2)); % Normalize the channel.


% Frequency expression of the channel response.
H = fft(channel,N );
%% 
% symbol error initialization
Ifdma_SER=zeros(length(SNR),1);
Dfdma_SER=zeros(length(SNR),1);
Lfdma_SER=zeros(length(SNR),1);
    
for n = 1:length(SNR)
% Initialize the error count.
Ifdma_error_number = 0;
Dfdma_error_number = 0;
Lfdma_error_number = 0;
for k = 1:N_simu
% Generate random data block.
% Input Generation
input_symbol=randi([0, (order-1)] ,1,M);
if modulation_type=='PSK'
input_signal = pskmod(input_symbol,order);

elseif  modulation_type=='QAM'
input_signal = qammod(input_symbol,order);
end 
%Apply FFT operation
input_signal_fft=fft(input_signal,M);

% Initialize subcarrier mappings 
Ifdma_mapping = zeros(1,N);
Dfdma_mapping = zeros(1,N);
Lfdma_mapping = zeros(1,N);

%Applying Mapping
Ifdma_mapping(1+SP.subband:Q:N) = input_signal_fft;
if Q ==1
Dfdma_mapping=Ifdma_mapping;
else 
Dfdma_mapping(1+SP.subband:Q_D:Q_D*M) = input_signal_fft; 

end 
Lfdma_mapping([1:M]+M*SP.subband) = input_signal_fft;

% Apply IFFT operation

Ifdma_IFFT = ifft(Ifdma_mapping,N);
Dfdma_IFFT = ifft(Dfdma_mapping,N);
Lfdma_IFFT = ifft(Lfdma_mapping,N);

% Cyclic Prefix Addition
Ifdma_cyclic = [Ifdma_IFFT(N-cp+1:N) Ifdma_IFFT ];
Dfdma_cyclic = [Dfdma_IFFT(N-cp+1:N) Dfdma_IFFT ];
Lfdma_cyclic = [Lfdma_IFFT(N-cp+1:N) Lfdma_IFFT ];
%% CHANNEL 

    % multi-path channel
    
    Ifdma_signal=filter(channel, 1, Ifdma_cyclic);
    Dfdma_signal=filter(channel, 1, Dfdma_cyclic);
    Lfdma_signal=filter(channel, 1, Lfdma_cyclic);
    % Generate AWGN 
    Noise  = (randn(1, N+cp)+1i*randn(1, N+cp))/sqrt(2); % N+cp by considering the added cyclic symbols 
    noisePower = 10^(-SNR(n)/10);
    % Add AWGN to the transmitted signal.
    Ifdma_signal =Ifdma_signal+sqrt(noisePower/M)*Noise;
    Dfdma_signal =Dfdma_signal+sqrt(noisePower/M)*Noise;
    Lfdma_signal =Lfdma_signal +sqrt(noisePower/M)*Noise;
    %% Receiver 
    % Removing cyclic prefix

    Ifdma_received =  Ifdma_signal (cp+1:N+cp);
    Dfdma_received = Dfdma_signal(cp+1:N+cp);
    Lfdma_received = Lfdma_signal (cp+1:N+cp);

   
    % Applying FFT , from time to frequency domain
    Ifdma_received = fft(Ifdma_received,N);
    Dfdma_received = fft(Dfdma_received,N);
    Lfdma_received = fft(Lfdma_received,N);
    
    % Applying demaping
    Ifdma_received = Ifdma_received(1+SP.subband:Q:N);
    if Q ==1
   Dfdma_received=Ifdma_received ;
    else 
    Dfdma_received = Dfdma_received(1+SP.subband:Q_D:Q_D*M);
      
     end 
   
    Lfdma_received = Lfdma_received ([1:M]+M*SP.subband);
    
    %  channel response of the subcarriers
    H_Ifdma=H(1+SP.subband:Q:N);
    if Q ==1
   H_Dfdma=  H_Ifdma ;
    else 
   H_Dfdma=H(1+SP.subband:Q_D:Q_D*M);
     end 
     
    H_Lfdma=H([1:M]+M*SP.subband);
   
%%    Perform frequency equalization
    % Zero Forcing 
    if equalizer_type =='ZFEQ'
    Ifdma_received = Ifdma_received./H_Ifdma;
    Dfdma_received = Dfdma_received./H_Dfdma;
    Lfdma_received = Lfdma_received./H_Lfdma;
    
    elseif equalizer_type=='MMSE'
    % MMSE
    C_Ifdma = conj(H_Ifdma)./(conj(H_Ifdma).*H_Ifdma + 10^(-SNR(n)/10));
    C_Dfdma = conj(H_Dfdma)./(conj(H_Dfdma).*H_Dfdma + 10^(-SNR(n)/10));
    C_Lfdma = conj(H_Lfdma)./(conj(H_Lfdma).*H_Lfdma + 10^(-SNR(n)/10));
    Ifdma_received  = Ifdma_received .* C_Ifdma;
    Dfdma_received  = Dfdma_received .* C_Dfdma;
    Lfdma_received  = Lfdma_received .* C_Lfdma;
    end 
    
    %% Applying  IFFT
    Ifdma_received = ifft(Ifdma_received);
    Dfdma_received = ifft(Dfdma_received);
    Lfdma_received = ifft(Lfdma_received);
    
    %%  Symbol detection 
    
    % QPSK Demodulation
     if modulation_type =='PSK'
    Ifdma_symbol = pskdemod( Ifdma_received , order);
    Dfdma_symbol = pskdemod( Dfdma_received , order);
    Lfdma_symbol = pskdemod( Lfdma_received , order);
    
    % QAM Demodulation
     elseif  modulation_type =='QAM'
    Ifdma_symbol = qamdemod( Ifdma_received , order);
    Dfdma_symbol = qamdemod( Dfdma_received , order);
    Lfdma_symbol = qamdemod( Lfdma_received , order);
    end 
    

    
    %% Error Calculation 
    
   % Number of correctly received symbol 
   Ifdma_correct_symbol = find((input_symbol-Ifdma_symbol) == 0);
   Dfdma_correct_symbol = find((input_symbol-Dfdma_symbol) == 0);
   Lfdma_correct_symbol = find((input_symbol-Lfdma_symbol) == 0);
   
   % the number of errors 
   Ifdma_error_number = Ifdma_error_number +(M-length(Ifdma_correct_symbol));
   Dfdma_error_number = Dfdma_error_number +(M-length(Dfdma_correct_symbol));
   Lfdma_error_number = Lfdma_error_number +(M-length(Lfdma_correct_symbol));
   
end 
   % Calculating symbol error rate 
    Ifdma_SER(n,:) = Ifdma_error_number/ (M*N_simu);
    Dfdma_SER(n,:) = Dfdma_error_number/ (M*N_simu);
    Lfdma_SER(n,:) = Lfdma_error_number/ (M*N_simu);
end


%% plotting

% figure  
 semilogy(SNR,Ifdma_SER, '-*',SNR, Dfdma_SER, '-o',SNR, Lfdma_SER, '-d')

first_zero1=find(Ifdma_SER==0,1,'first');

Ifdma_SER=Ifdma_SER(1:first_zero1);

first_zero2=find(Dfdma_SER==0,1,'first');
Dfdma_SER=Dfdma_SER(1:first_zero2);

first_zero3=find(Lfdma_SER==0,1,'first');
Lfdma_SER=Lfdma_SER(1:first_zero3);

Y_min=10^-8;

semilogy(SNR(1:first_zero1),max(Y_min,Ifdma_SER),'LineWidth', 2)
hold on, 

semilogy(SNR(1:first_zero2), max(Y_min, Dfdma_SER),'LineWidth', 2)
hold on, 

semilogy(SNR(1:first_zero3), max(Y_min,Lfdma_SER),'LineWidth', 2)
hold off

legend(sprintf('Interleaved user %K', Q),sprintf('Distributed user %K', Q),sprintf('Localized  user  %r', Q),'FontSize',10);

if modulation_type=='QAM'
title(['Symbol error rate diagram for  ',num2str(order),' - ' modulation_type, ' and input block size ',num2str(M) ]);
elseif modulation_type=='PSK' 
    if order==2
title(['Symbol error rate diagram for BPSK ',' and input block size ',num2str(M)]);
    else
 title(['Symbol error rate diagram for QPSK  ',' and input block size ',num2str(M)]);
    end 
end 

xlabel('Signal to Noise Ratio in [dB]', 'FontSize',10 )
ylabel('Symbol Error rate', 'FontSize',10 )

set(gca,'FontSize',12)
axis( [0 30 Y_min 1])



