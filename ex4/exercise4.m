%% task 1

load('521282S_eeg_data_4.mat');
F = 0.1:0.1:32;
[S,F,T,P] = spectrogram(signal,Fs*30,Fs*29,F,Fs);

figure
subplot(3,1,1);
plot(t, signal);
ylim([-1000 1000])

subplot(3,1,2);
Image1 = imagesc(T/60,F,log10(P),[-4 6]);
axis xy;

%% spectral entropy

f1 = 0.1;
fh = 32;
n = fh/f1;
for k = 1:1:991
    SE(k) = 0;
    for i = 1:1:320
    SE(k) = SE(k) + P(i,k)*log(P(i,k));
    end
    SE(k) = (-SE(k)/log(Fs))/n;
end

subplot(3,1,3);
plot(T,SE);
hold off

%% task 2
figure

t_start = find(t == 208/60);
t_end = find(t == 440/60);
plot(t2, signal(t_start:t_end));
ylim([-1000 1000]);



delta = 67;
t_ictal = find(t2 == 100);
t_pre = find(t2 == 90);
t_post = find(t2 == 132);
t_post_end = find(t2 == 142);

t_1 = find(t2 == 89);
t_2 = find(t2 == 160);

figure
plot([-1500 1500],[-1500 1500], '.');
hold on
% From 100 s to 132 s ------- seizure state
plot(signal(t_ictal+t_start : t_post+t_start),signal(t_ictal + delta+t_start : t_post + delta+t_start), '.','MarkerEdgeColor','c');
hold on
% From 1 s to 100-11 s ------- seizure state
plot(signal(t_start : t_1),signal(t_start + delta : t_1 + delta), '.','MarkerEdgeColor','r');
hold on
% From 132+11 s to 160 s ------- seizure state
plot(signal(t_post_end : t_2),signal(t_post_end + delta : t_2 + delta), '.','MarkerEdgeColor','r');
hold on
% From 90 s to 100 s ------- pre ictal state
plot(signal(t_pre+t_start : t_ictal+t_start),signal(t_pre + delta+t_start : t_ictal + delta+t_start), '.','MarkerEdgeColor','k');
hold on
% From 132 s to 142 s ------- post ictal state
plot(signal(t_post+t_start : t_post_end+t_start),signal(t_post + delta+t_start : t_post_end + delta+t_start), '.','MarkerEdgeColor','k');
hold on

