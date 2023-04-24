clear;
close all;
clc;
%% 1.1 Data Generation

[v_m,fs] = audioread("in-the-air.wav");
% sound(v_m,fs);

T_s = 1/fs;
N = length(v_m);
t = 0:T_s:(N-1)*T_s;
f = linspace(-fs/2,fs/2,N);
V_m = fftshift(fft(v_m))/sqrt(N);

% figure;
% subplot(2,1,1);
% plot(t,v_m);
% grid on;
% ylabel("v_m(t)");
% xlabel("t");
% 
% subplot(2,1,2);
% plot(f,V_m);
% grid on;
% ylabel("V_m(f)");
% xlabel("f");

% evaluation of the bandwidth

bw_v =find( V_m > 0.1*max(V_m) );
bw = max(bw_v / (length(f)/2))*10^4;

% bw = 0.7*10^4;
%% 1.2 Modulator
fc = 15*10^3;  % carrier frequency
k_AM = 0.02;  % modulation index

v_AM = ammod(v_m,fc,fs,0,k_AM);
V_AM = fftshift(fft(v_AM))/sqrt(N);

% figure;
% plot(f,V_AM);
% grid on;
% ylabel("V_AM(f)");
% xlabel("f");

% sound(v_AM,fs); 

% %% 1.3 Channel
% N_0 = 8*10^-4;
% z = (sqrt(N_0/2)*randn(1,N)).';
% x_r = v_AM + z;
% X_r = fftshift(fft(x_r))/sqrt(N);
% 
% figure;
% subplot(2,1,1);
% plot(t,x_r);
% grid on;
% ylabel("x_r(t)");
% xlabel("t");
% 
% subplot(2,1,2);
% plot(f,X_r);
% grid on;
% ylabel("X_r(f)");
% xlabel("f");
% 
% 
% %% 1.4 Demodulator
% 
% % 1.4.1
% x_L = bandpass(x_r,[fc-bw fc+bw],fs);
% X_L = fftshift(fft(x_L))/sqrt(N);
% 
% figure;
% subplot(2,1,1);
% plot(t,x_L);
% grid on;
% ylabel("x_L(t)");
% xlabel("t");
% 
% subplot(2,1,2);
% plot(f,X_L);
% grid on;
% ylabel("X_L(f)");
% xlabel("f");
% 
% % 1.4.2
% x_d = amdemod(x_L,fc,fs,0,k_AM);
% X_d = fftshift(fft(lowpass(x_d,bw,fs)))/sqrt(N);
% 
% % X_d1 = fftshift(fft(x_d))/sqrt(N);
% % X_d = lowpass(X_d1,bw,fs);
% 
% figure;
% subplot(2,1,1);
% plot(t,x_d);
% hold on
% plot(t,v_m);
% ylabel("v_m(t) and x_d(t)");
% xlabel("t");
% legend('x_d(t)','v_m(t)');
% 
% subplot(2,1,2);
% plot(f,X_d);
% hold on 
% plot(f,V_m);
% grid on;
% ylabel("V_m(f) and X_d(f)");
% xlabel("f");
% legend('X_d(f)','V_m(f)');
% 
% % sound(x_d,fs); %1.4.3
% 
% correlation = xcorr(x_d,v_m,0,'coeff'); %1.4.4
% 
% %% 1.5
% % 1.5.1
% 
% % 1.3 Channel
% N_0 = 0.02;
% z = (sqrt(N_0/2)*randn(1,N)).';
% x_r = v_AM + z;
% X_r = fftshift(fft(x_r))/sqrt(N);
% 
% figure;
% subplot(2,1,1);
% plot(t,x_r);
% grid on;
% ylabel("x_r(t)");
% xlabel("t");
% 
% subplot(2,1,2);
% plot(f,X_r);
% grid on;
% ylabel("X_r(f)");
% xlabel("f");
% 
% 
% % 1.4 Demodulator
% 
% % 1.4.1
% x_L = bandpass(x_r,[fc-bw fc+bw],fs);
% X_L = fftshift(fft(x_L))/sqrt(N);
% 
% figure;
% subplot(2,1,1);
% plot(t,x_L);
% grid on;
% ylabel("x_L(t)");
% xlabel("t");
% 
% subplot(2,1,2);
% plot(f,X_L);
% grid on;
% ylabel("X_L(f)");
% xlabel("f");
% 
% % 1.4.2
% x_d = amdemod(x_L,fc,fs,0,k_AM);
% X_d = fftshift(fft(lowpass(x_d,bw,fs)))/sqrt(N);
% 
% % X_d1 = fftshift(fft(x_d))/sqrt(N);
% % X_d = lowpass(X_d1,bw,fs);
% 
% figure;
% subplot(2,1,1);
% plot(t,x_d);
% hold on
% plot(t,v_m);
% ylabel("v_m(t) and x_d(t)");
% xlabel("t");
% legend('x_d(t)','v_m(t)');
% 
% subplot(2,1,2);
% plot(f,X_d);
% hold on 
% plot(f,V_m);
% grid on;
% ylabel("V_m(f) and X_d(f)");
% xlabel("f");
% legend('X_d(f)','V_m(f)');
% 
% sound(x_d,fs); %1.4.3
% 
% correlation(2) = xcorr(x_d,v_m,0,'coeff'); %1.4.4

%% 1.5.2
fd = 10^4;
v_FM = fmmod(v_m,fc,fs,fd);
V_FM = fftshift(fft(v_FM))/sqrt(N).';

%% 1.3 Channel
N_0 = 8*10^-4;
z = (sqrt(N_0/2)*randn(1,N)).';
x_r = v_FM + z;
X_r = fftshift(fft(x_r))/sqrt(N);

figure;
subplot(2,1,1);
plot(t,x_r);
grid on;
ylabel("x_r(t)");
xlabel("t");

subplot(2,1,2);
plot(f,X_r);
grid on;
ylabel("X_r(f)");
xlabel("f");


%% 1.4 Demodulator

% 1.4.1
x_L = bandpass(x_r,[fc-bw fc+bw],fs);
X_L = fftshift(fft(x_L))/sqrt(N);

figure;
subplot(2,1,1);
plot(t,x_L);
grid on;
ylabel("x_L(t)");
xlabel("t");

subplot(2,1,2);
plot(f,X_L);
grid on;
ylabel("X_L(f)");
xlabel("f");

% 1.4.2
x_d = fmdemod(x_L,fc,fs,fd);
X_d = fftshift(fft(lowpass(x_d,bw,fs)))/sqrt(N);

% X_d1 = fftshift(fft(x_d))/sqrt(N);
% X_d = lowpass(X_d1,bw,fs);

figure;
subplot(2,1,1);
plot(t,x_d);
hold on
plot(t,v_m);
ylabel("v_m(t) and x_d(t)");
xlabel("t");
legend('x_d(t)','v_m(t)');

subplot(2,1,2);
plot(f,X_d);
hold on 
plot(f,V_m);
grid on;
ylabel("V_m(f) and X_d(f)");
xlabel("f");
legend('X_d(f)','V_m(f)');

% sound(x_d,fs); %1.4.3

correlation(3) = xcorr(x_d,v_m,0,'coeff'); %1.4.4

%% 1.5
% 1.5.1

%% 1.3 Channel
N_0 = 0.02;
z = (sqrt(N_0/2)*randn(1,N)).';
x_r = v_FM + z;
X_r = fftshift(fft(x_r))/sqrt(N);

figure;
subplot(2,1,1);
plot(t,x_r);
grid on;
ylabel("x_r(t)");
xlabel("t");

subplot(2,1,2);
plot(f,X_r);
grid on;
ylabel("X_r(f)");
xlabel("f");


%% 1.4 Demodulator

% 1.4.1
x_L = bandpass(x_r,[fc-bw fc+bw],fs);
X_L = fftshift(fft(x_L))/sqrt(N);

figure;
subplot(2,1,1);
plot(t,x_L);
grid on;
ylabel("x_L(t)");
xlabel("t");

subplot(2,1,2);
plot(f,X_L);
grid on;
ylabel("X_L(f)");
xlabel("f");

% 1.4.2
x_d = fmdemod(x_L,fc,fs,fd);
X_d = fftshift(fft(lowpass(x_d,bw,fs)))/sqrt(N);

% X_d1 = fftshift(fft(x_d))/sqrt(N);
% X_d = lowpass(X_d1,bw,fs);

figure;
subplot(2,1,1);
plot(t,x_d);
hold on
plot(t,v_m);
ylabel("v_m(t) and x_d(t)");
xlabel("t");
legend('x_d(t)','v_m(t)');

subplot(2,1,2);
plot(f,X_d);
hold on 
plot(f,V_m);
grid on;
ylabel("V_m(f) and X_d(f)");
xlabel("f");
legend('X_d(f)','V_m(f)');

% sound(x_d,fs); %1.4.3

correlation(4) = xcorr(x_d,v_m,0,'coeff'); %1.4.4