clear;
close all;
%% 1.1 Data Generation

[v_m,fs] = audioread("in-the-air.wav");
sound(v_m,fs);

T_s = 1/fs;
N = length(v_m);
t = 0:T_s:(N-1)*T_s;
f = linspace(-fs/2,fs/2,N);
V_m = fftshift(fft(v_m))/sqrt(N);

figure;
subplot(2,1,1);
plot(t,v_m);
grid on;
ylabel("v_m(t)");
xlabel("t");

subplot(2,1,2);
plot(f,abs(V_m));
grid on;
ylabel("|V_m(f)|");
xlabel("f");

% evaluation of the bandwidth
bw_v =find( V_m > 0.1*max(V_m) );
bw = max(bw_v / (length(f)/2))*10^4;

%% 1.2 Modulator
fc = 15*10^3;  % carrier frequency
k_AM = 0.02;  % modulation index

v_AM = ammod(v_m,fc,fs);

V_AM = fftshift(fft(v_AM))/sqrt(N);

figure;
plot(f,abs(V_AM));
grid on;
ylabel("V_AM(f)");
xlabel("f");

sound(v_AM,fs); 

%% 1.3 Channel
N_0 = 8*10^-4;
z = (sqrt(N_0/2)*randn(1,N)).';
x_r = v_AM + z;
X_r = fftshift(fft(x_r))/sqrt(N);

figure;
subplot(2,1,1);
plot(t,x_r);
grid on;
ylabel("x_r(t)");
xlabel("t");

subplot(2,1,2);
plot(f,abs(X_r));
grid on;
ylabel("|X_r(f)|");
xlabel("f");


%% 1.4 Demodulator

x_L = bandpass(x_r,[fc-bw/2,fc+bw/2],fs);
X_L = fftshift(fft(x_L))/sqrt(N);

figure;
subplot(2,1,1);
plot(t,x_L);
grid on;
ylabel("x_L(t)");
xlabel("t");

subplot(2,1,2);
plot(f,abs(X_L));
grid on;
ylabel("|X_L(f)|");
xlabel("f");

x_d = amdemod(x_L,fc,fs);
X_d1 = fftshift(fft(x_d))/sqrt(N);
X_d = lowpass(X_d1,bw,fs);

figure;
plot(f,abs(V_m));
hold on 
plot(f,abs(X_d),'--');
grid on;
ylabel("Amp");
xlabel("f");
legend('V_m(f)','X_d(f)');

sound(x_d,fs); %1.4.3

xcorr(x_d,v_m,0,'coeff'); %1.4.4









