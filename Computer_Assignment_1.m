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
plot(f,V_m);
grid on;
ylabel("V_m");
xlabel("f");

