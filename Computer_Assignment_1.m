clear;
close all;
%% 1.1 Data Generation

[v_m,fs] = audioread("in-the-air.wav");
%sound(v_m,fs);

T_s = 1/fs;
N = length(v_m);
t = 0:T_s:(N-1)*T_s;
f = linspace(-fs/2,fs/2,N);
V_m = fftshift(fft(v_m))/sqrt(N);

figure;
subplot(2,1,1);
plot(t,v_m);
grid on;
ylabel("v_m");
xlabel("t");

subplot(2,1,2);
plot(f,abs(V_m));
grid on;
ylabel("|V_m|");
xlabel("f");

% evaluation of the bandwidth
%bw = V_m>0.01*max(V_m);
%BW = 1.1khz

%% 1.2 Modulator
fc = 15*10^3;  % carrier frequency
k_AM = 0.02;  % modulation index

v_AM = ammod(v_m,fc, fs);

V_AM = fftshift(fft(v_AM))/sqrt(N);

figure;
plot(f,abs(V_AM));
grid on;
ylabel("V_AM");
xlabel("f");

sound(v_AM,fs); 