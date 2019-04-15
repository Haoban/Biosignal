%% Task 1 Raw SSEP signals and ensemble averaging
load('521282S_datasample.mat');

raw_data = data_samples.raw_data;
Fs=data_samples.Fs;
amplitude = data_samples.amplitude_unit;
V_Scale = raw_data * amplitude* 1000000;

M = size(raw_data,1);
N = size(raw_data,2);
for n = 1:1:N
    sa(n)=sum(V_Scale(:,n))/M;
end
t=linspace(1,M,N);
figure
plot(t,V_Scale()','LineWidth',1);

hold on
plot(t, sa, '-r','LineWidth',5);
%% Task 2 SNR estimation
% Variance
for n = 1:1:N
        sigma_v_2(n) = sum((V_Scale(:,n)-sa(n)).^2)/M;
end

Sigma_v2 = sum(sigma_v_2)/N;

% mean squared signal energy
Es = 0;
for n = 1:1:N
    Es = Es+sa(n)^2;
end
Es = Es/N;

% SNR
SNR=10*log10(Es/Sigma_v2);

% Estimate averaging SNR dB
ave_SNR_dB = 20*log10(sqrt(M));

% New estimated SNR
New_estimated_SNR=ave_SNR_dB+SNR;
%% Task 3 Chirp modeling

newfeatures=chirp_features(data_samples);
chirp_signal=newfeatures.fitted_chirp;
hold on
plot(t,chirp_signal(:,1),'-g','LineWidth',3);

hold on
plot([newfeatures.fit_delay*1000 newfeatures.fit_delay*1000],[-100 100],'--k','LineWidth',1);
hold on
plot([newfeatures.fit_first_peak_location*1000 newfeatures.fit_first_peak_location*1000],[-100 100],'--k','LineWidth',1);
hold on
plot([newfeatures.fit_second_peak_location*1000 newfeatures.fit_second_peak_location*1000],[-100 100],'--k','LineWidth',1);
hold on
plot([0 100],[-newfeatures.fit_second_peak_amplitude -newfeatures.fit_second_peak_amplitude],'--k','LineWidth',1);
hold on
plot([0 100],[newfeatures.fit_first_peak_amplitude newfeatures.fit_first_peak_amplitude],'--k','LineWidth',1);




