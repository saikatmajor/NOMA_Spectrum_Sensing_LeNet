clear;
clc;
close all;
SNRdB = -14;
SNR = 10^(SNRdB/10);

M = 28;                 % Number of secondary users
N = 100;                % Time slots in sensing interval
K = 20000;              % Training data length
MaxTrials = 100;        % Number of trials of Monte Carlo simulations with different channels
Omega = [4,1];          % Power ratios of PU
Omega = Omega/sum(Omega);
Q = 1;                  % Number of NOMA primary users (1 for downlink)
noiseVar = ones(1,M);

Theta = rand(K,1)>0.5;  % State of the PU
H = (randn(Q,M) + 1i*randn(Q,M))/sqrt(2); % Rayleigh channel between PU and SU
CC = zeros(K,2,M,M);

for k = 1:K
    theta = Theta(k,:);
    Tx = sqrt(Omega).*(randn(N,Q) + 1i*randn(N,Q))/sqrt(2);   % PU signal 
    S = sum(Tx,2);      % Power domain multiplexing (added).

    signalUser = S*H;
    noiseUser = sqrt(noiseVar).*(randn(N,M) + 1i*randn(N,M))/sqrt(2);

    Ps = trace(signalUser'*signalUser);     % Estimate signal power
    Pn = trace(noiseUser'*noiseUser);       % Estimate noise power
    signalPower = SNR/(Ps/Pn);

    rxUser = sqrt(signalPower)*(Theta(k).*S)*H + noiseUser;
    
    %signal_power = sum(abs(sqrt(signalPower)*signalUser).^2);
    %noise_power = sum(abs(noiseUser).^2);
    %snr_actual = 10*log10(sum(signal_power./noise_power)/M);
    
    Ry = rxUser'*rxUser/N;
    CC(k,1,:,:) = real(Ry);
    CC(k,2,:,:) = imag(Ry);
end
theta = Theta;
save('noma_downlink_2state_data.mat','CC','theta');
