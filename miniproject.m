% Ha Noi, 18/1/2024
% By Bui Van Canh

close all;clear ;clc;

N = 100000 ;    % number of bits or symbols

ip = rand(N,1)>0.5;  % generating 0,1 with equal probability

% Choosing the modulation techniques
for j = 1:3
 
if j == 1
    %FSK modualtion
    [mod,demod] = fsk();
    type = ' FSK ';
elseif j == 2
    %BPSK modulation
    [mod,demod] = bpsk();
    type = ' BPSK ';
else
    %QPSK modulation
    [mod,demod] = qpsk();
    type = ' QPSK ';
end

s = step(mod,ip);

Eb_N0_dB = -3:35; % multiple Eb/N0 values 

%Noise addition , equalization and error calculation
for i = 1:length(Eb_N0_dB)
   
   n = 1/sqrt(2)*(randn(length(s),1) + j*randn(length(s),1)); % white gaussian noise, 0dB variance 
   h = 1/sqrt(2)*(randn(length(s),1) + j*randn(length(s),1)); % Rayleigh channel
   % h = 1/sqrt(2)*(randn(1,N)+1i*randn(1,N));
   % Channel and noise Noise addition
   y = h.*s + 10^(-Eb_N0_dB(i)/20)*n;
   %y = h.*s;

   % no using equalization
        yHat1 = ones(size(h)) .* y;

   % using Zero forcing  equalization
        yHat2 = y./h;
      
   % using MMSE qualization
        yHat3 = (conj(h)./((abs(h)).^2+n)).*y;
       
% Bit in output:

%output with no using equalization 
   op1 = step(demod,yHat1);

%output with  Zero Forcing equalization 
   op2 = step(demod,yHat2);

%outpt with MMSE equalization
   op3 = step(demod,yHat3);

% Counting the errors:

%Number of bit errors with no equalization 
   nErr1(i) = size(find((ip-op1)),1);

%Number of bit errors with Zero Forcing equalization
   nErr2(i) = size(find((ip-op2)),1);

%Number of bit errors with MMSE equalization
   nErr3(i) = size(find((ip-op3)),1);

end

% BER CALCULATION:

% simulated ber with no using equalization
simBer1 = nErr1/N; 

% simulated ber with Zero Forcing equalization
simBer2 = nErr2/N; 

% simulated ber with MMSE equalization
simBer3 = nErr3/N; 

 % theoretical ber
theoryBerAWGN = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10)));

EbN0Lin = 10.^(Eb_N0_dB/10);

% theoryBer in Rayleigh Channel
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1))); 

% Figures
% plot
% Graph of no using equalization
figure
semilogy(Eb_N0_dB,theoryBerAWGN,'b-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,theoryBer,'r-','LineWidth',2);
semilogy(Eb_N0_dB,simBer1,'k-','LineWidth',2);
axis([-3 35 10^-5 0.5])
grid on
legend('AWGN','Rayleigh-Theory', 'Rayleigh-Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
head = strcat('BER for',type,'modulation in Rayleigh channel no using equalization' );
title(head);

% Graph of using Zero Forcing equalization
figure 
semilogy(Eb_N0_dB,theoryBerAWGN,'b-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,theoryBer,'r-','LineWidth',2);
semilogy(Eb_N0_dB,simBer2,'k-','LineWidth',2);
axis([-3 35 10^-5 0.5])
grid on
legend('AWGN','Rayleigh-Theory', 'Rayleigh-Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
head = strcat('BER for',type,'modulation in Rayleigh channel  using Zero Forcing' );
title(head);

% Graph of using MMSE equalization
figure 
semilogy(Eb_N0_dB,theoryBerAWGN,'b-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,theoryBer,'r-','LineWidth',2);
semilogy(Eb_N0_dB,simBer3,'k-','LineWidth',2);
axis([-3 35 10^-5 0.5])
grid on
legend('AWGN','Rayleigh-Theory', 'Rayleigh-Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
head = strcat('BER for',type,'modulation in Rayleigh channel  using MMSE' );
title(head);
 
end


% Function setup
% BPSK modulation and demodulation
function [mod,demod] = bpsk()
mod = comm.BPSKModulator();
demod = comm.BPSKDemodulator();
end

% FSK modulation and demodulation
function [mod,demod] = fsk()
mod = comm.FSKModulator();
demod = comm.FSKDemodulator();
end

% QPSK modulation and demodulation
function [mod,demod] = qpsk()
mod = comm.QPSKModulator();
mod.BitInput = true;
demod = comm.QPSKDemodulator();
demod.BitOutput = true;
end
