load('521282S_data_2.mat')
% delta
Fs = 200;  % Sampling Frequency

N    = 800;      % Order
Fc1  = 1;        % First Cutoff Frequency
Fc2  = 4;        % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Hd = dfilt.dffir(b);
delta_signal = filter(Hd,signal) ;

% alpha
Fs = 200;  % Sampling Frequency

N    = 800;      % Order
Fc1  = 8;        % First Cutoff Frequency
Fc2  = 12;        % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Hd2 = dfilt.dffir(b);
alpha_signal = filter(Hd2,signal) ;

% theta
Fs = 200;  % Sampling Frequency

N    = 800;      % Order
Fc1  = 8;        % First Cutoff Frequency
Fc2  = 12;        % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Hd3 = dfilt.dffir(b);
theta_signal = filter(Hd3,signal) ;

% beta
Fs = 200;  % Sampling Frequency

N    = 800;      % Order
Fc1  = 12;        % First Cutoff Frequency
Fc2  = 25;        % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Hd4 = dfilt.dffir(b);
beta_signal = filter(Hd4,signal) ;

%% plot
figure('position',[246 99 560 990]);
ylims=[-0.03 0.03];
subplot(5,1,1);hold on;
plot(t, signal);
set(gca,'ylim',ylims);

subplot(5,1,2);hold on;
plot(t, delta_signal);
set(gca,'ylim',ylims);

subplot(5,1,3);hold on;
plot(t, theta_signal);
set(gca,'ylim',ylims);

subplot(5,1,4);hold on;
plot(t, alpha_signal);
set(gca,'ylim',ylims);

subplot(5,1,5);hold on;
plot(t, beta_signal);
set(gca,'ylim',ylims);

[s_upper,s_lower] = envelope(signal,30*Fs,'rms');
[d_upper,d_lower] = envelope(delta_signal, 30*Fs,'rms');
[t_upper,t_lower] = envelope(theta_signal, 30*Fs,'rms');
[a_upper,a_lower] = envelope(alpha_signal, 30*Fs,'rms');
[b_upper,b_lower] = envelope(beta_signal, 30*Fs,'rms');




subplot(5,1,1);
plot(t, s_upper,'r');
plot(t, s_lower,'r');

subplot(5,1,2);
plot(t, d_upper,'r');
plot(t, d_lower,'r');


subplot(5,1,3);
plot(t, t_upper,'r');
plot(t, t_lower,'r');


subplot(5,1,4);
plot(t, a_upper,'r');
plot(t, a_lower,'r');


subplot(5,1,5);
plot(t, b_upper,'r');
plot(t, b_lower,'r');


%% spectrogram


F = 0.1:0.1:32;
[S,F,T,P] = spectrogram(signal,Fs*30,Fs*29,F,Fs);
figure
subplot(3,1,1);hold on;
Image1 = imagesc(T/60,F,log10(P),[-7 -3]);
axis xy;

ylims=[0 1];
power_delta = (sum(P(10:40,:))./sum(P));
power_theta = (sum(P(40:80,:))./sum(P));
power_alpha = (sum(P(80:120,:))./sum(P));
power_beta = (sum(P(120:250,:))./sum(P));
subplot(3,1,2);hold on;

plot(T,power_delta,'r');
set(gca,'ylim',ylims);
hold on
plot(T,power_theta,'b');
set(gca,'ylim',ylims);
hold on
plot(T,power_alpha,'y');
set(gca,'ylim',ylims);
hold on
plot(T,power_beta,'g');
set(gca,'ylim',ylims);

legend('power_delta','power_theta','power_alpha','power_beta');

%% spectral
f1 = 0.1;
fh = 32;
n = fh/f1;
for k = 1:1:422
    SE(k) = 0;
    for i = 1:1:320
    SE(k) = SE(k) + P(i,k)*log(P(i,k));
    end
    SE(k) = (-SE(k)/log(Fs))/n;
end


subplot(3,1,3);
plot(T,SE);
