figure
hold on
plot(t, signal+spread,'-')
hold off

Fs = 200;  % Sampling Frequency

N    = 200;      % Order
Fc   = 90;       % Cutoff Frequency
flag = 'scale';  % Sampling Flag

win = hamming(N+1);

b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
LP = dfilt.dffir(b);

LP_signal = filter(LP,signal) 

%%%%%%%

Fs = 200;  % Sampling Frequency

Fnotch = 50;  % Notch Frequency
Q      = 35;  % Q-factor
Apass  = 3;   % Bandwidth Attenuation

BW = Fnotch/Q;

[b, a] = iirnotch(Fnotch/(Fs/2), BW/(Fs/2), Apass);
% Hd     = dfilt.df2(b, a);

Not_signal = filtfilt(b,a,LP_signal);

%%%%%%%

Fs = 200;  % Sampling Frequency

N     = 500;  % Order
Fstop = 0.1;  % Stopband Frequency
Fpass = 0.5;  % Passband Frequency
Wstop = 1;    % Stopband Weight
Wpass = 1;    % Passband Weight

% Calculate the coefficients using the FIRLS function.
b  = firls(N, [0 Fstop Fpass Fs/2]/(Fs/2), [0 0 1 1], [Wstop Wpass]);
Hp = dfilt.dffir(b);

HP_signal = filter(Hp,Not_signal);

figure
subplot(3,1,1);
plot(t, LP_signal+spread,'-');

subplot(3,1,2);
plot(t, Not_signal+spread,'-');

subplot(3,1,3);
plot(t, HP_signal+spread,'-');

LMSfit = dsp.LMSFilter('Length',11, 'Method' ,'LMS', 'StepSize', 0.6, 'WeightsOutputPort', 'false' )

EOG1 = 14
EOG2 = 15

channels = 15;

% eog_channel_in(:,1) = signal(:,EOG1)
% eog_channel_in(:,2) = signal(:,EOG2)


% LP_lms = LMSfilter(:,1:channels)


for j = 1:15
    for i = 1:10
        [ out1 , out2] = step(LMSfit, signal(:,EOG1) , signal(:,j));
    end
    newchannelEOG1(:,j) = out1;
    newchannelEEG1(:,j) = out2;
end



figure
hold on
plot(t, newchannelEEG1 + spread,'-');
hold off

for j = 1:15
    for i = 1:10
        [ out1 , out2] = step(LMSfit, newchannelEEG1(:,EOG2) , newchannelEEG1(:,j));
    end
    newchannelEOG2(:,j) = out1;
    newchannelEEG2(:,j) = out2;
end


figure
hold on
plot(t, newchannelEEG2+spread,'-')
hold off

[IC, A, W] = fastica(HP_signal');
plot(t,IC'/100+spread);

IC_fil=IC;
IC_fil(1:6,:)=0;
T=A*IC_fil;
plot(t,T'+spread);

