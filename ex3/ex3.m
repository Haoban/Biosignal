%% task 1
load('521282S_data_3.mat');

% help topoplot
figure
Vector = [];
topoplot(Vector,'channelLocations.locs','electrodes','ptslabels');

%% task 2
% spectrogram
figure

F = 0.1:0.1:32;
for i = 1:13
    [S,F,T,P] = spectrogram(signal(i,:),Fs*30,Fs*29,F,Fs);
    alpha_P = P(10:40,:);
    delta_P = P(80:120,:);
    power(i,:) = sum(P,1);
%     Sig_total_start(i) = log10(sum(P(:,1)));
%     Sig_total_2min(i) = log10(sum(P(:,120)));
%     Sig_total_5min(i) = log10(sum(P(:,300)));
%     Sig_total_7min(i) = log10(sum(P(:,420)));
end

    subplot(2,2,1);hold on;
    topoplot(power(:,1),'channelLocations.locs','electrodes','on','maplimits','maxmin');
    subplot(2,2,2);hold on;
    topoplot(power(:,120),'channelLocations.locs','electrodes','on','maplimits','maxmin');   
    subplot(2,2,3);hold on;
    topoplot(power(:,300),'channelLocations.locs','electrodes','on','maplimits','maxmin');
    subplot(2,2,4);hold on;
    topoplot(power(:,420),'channelLocations.locs','electrodes','on','maplimits','maxmin');

% 


% four channels
F8 = power(1,:);
Fpz = power(3,:);
F7 = power(4,:);
Fz = power(10,:);
frontal_electrodes = F8 + Fpz + F7 + Fz;
P3 = power(9,:);
P4 = power(7,:);
Pz = power(12,:);
Oz = power(13,:);
rear_electrodes = P3 + Pz + P4 + Oz;
figure 
plot(T/60, frontal_electrodes , 'r', T/60, rear_electrodes , 'b');

%% task 3
F = 0.1:0.1:32;
sum_P = 0;
for i = 1:13
    [S,F,T,P] = spectrogram(signal(i,:),Fs*30,Fs*29,F,Fs);
    alpha_Power(i,:) = sum(P(10:40,:),1);
    delta_Power(i,:) = sum(P(80:120,:),1);
    power(i,:) = sum(P,1);
    sum_P = sum_P + P;
end

power_delta_start = (sum(sum_P(10:40,1))./sum(P(:,1)));
power_alpha_start = (sum(sum_P(80:120,1))./sum(P(:,1)));
power_delta_2min = (sum(sum_P(10:40,120))./sum(P(:,120)));
power_alpha_2min = (sum(sum_P(80:120,120))./sum(P(:,120)));
power_delta_5min = (sum(sum_P(10:40,300))./sum(P(:,300)));
power_alpha_5min = (sum(sum_P(80:120,300))./sum(P(:,300)));
power_delta_7min = (sum(sum_P(10:40,420))./sum(P(:,420)));
power_alpha_7min = (sum(sum_P(80:120,420))./sum(P(:,420)));

    figure
    subplot(4,2,1);hold on;
    topoplot(alpha_Power(:,1),'channelLocations.locs','electrodes','on','maplimits','maxmin');
    subplot(4,2,2);hold on;
    topoplot(delta_Power(:,120),'channelLocations.locs','electrodes','on','maplimits','maxmin');   
    subplot(4,2,3);hold on;
    topoplot(alpha_Power(:,300),'channelLocations.locs','electrodes','on','maplimits','maxmin');
    subplot(4,2,4);hold on;
    topoplot(delta_Power(:,420),'channelLocations.locs','electrodes','on','maplimits','maxmin');
    subplot(4,2,5);hold on;
    topoplot(alpha_Power(:,1),'channelLocations.locs','electrodes','on','maplimits','maxmin');
    subplot(4,2,6);hold on;
    topoplot(delta_Power(:,120),'channelLocations.locs','electrodes','on','maplimits','maxmin');   
    subplot(4,2,7);hold on;
    topoplot(alpha_Power(:,300),'channelLocations.locs','electrodes','on','maplimits','maxmin');
    subplot(4,2,8);hold on;
    topoplot(delta_Power(:,420),'channelLocations.locs','electrodes','on','maplimits','maxmin');