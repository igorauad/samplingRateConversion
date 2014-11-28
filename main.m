clear all;
close all;
clc

set(0,'DefaultAxesFontSize',11)
set(0,'DefaultTextInterpreter','LaTex')

fs = 8000;          % sampling frequency
ts = 1 / fs;        % sampling period

%% Example signal

Nfft = 64;
nSymbols = 1e3;
x = generateOFDMSymbol( Nfft, nSymbols, 0.1, 0.95 );
psdFromFFT(x, Nfft, fs);

%% Upsampling

M = 4;
x_up = upsample(x, M);

% Freq Domain:
psdFromFFT(x_up, M*Nfft, M*fs);
title('Upsampled Signal PSD')

% Time Domain:
figure
subplot(211)
plot((0:(Nfft-1))*ts, x((Nfft+1):2*Nfft))
xlabel('Time')
ylabel('Amplitude')
title('Comparison in Time Domain')
subplot(212)
plot((0:(M*Nfft - 1))*(ts/M), x_up((M*Nfft + 1):2*M*Nfft), 'r')
xlabel('Time')
ylabel('Amplitude')
title('Comparison in Time Domain')

% NOTE: observe the gain reduced in upsampling

%% Interpolation

M = 4;
[x_interp b] = interp(x, M);

% Freq Domain:
psdFromFFT(x_interp, M*Nfft, M*fs);
title('Interpolated Signal PSD')

% Time Domain:
figure
subplot(211)
plot((0:(Nfft-1))*ts, x((Nfft+1):2*Nfft))
xlabel('Time')
ylabel('Amplitude')
title('Comparison in Time Domain')
subplot(212)
plot((0:(M*Nfft - 1))*(ts/M), x_interp((M*Nfft + 1):2*M*Nfft), 'r')
xlabel('Time')
ylabel('Amplitude')

% Filter
figure
freqz(b)
title('Anti-imaging filter')


%% Downsampling with aliasing

M = 2;
x_down = downsample(x, M);

% Freq Domain:
psdFromFFT(x_down, Nfft/M, fs/M);
title('Downsampled Signal PSD')

% Time Domain:
figure
subplot(211)
plot((0:(Nfft-1))*ts, x((Nfft+1):2*Nfft))
xlabel('Time')
ylabel('Amplitude')
title('Comparison in Time Domain')
subplot(212)
plot((0:((Nfft/M) - 1))*(M*ts), x_down((Nfft/M + 1):2*(Nfft/M)), 'r')
xlabel('Time')
ylabel('Amplitude')

% Note: observe the gain

%% Downsampling without aliasing

% Signal with shorter bandwidth:
x2 = generateOFDMSymbol( Nfft, nSymbols, 0.1, 0.4 );
psdFromFFT(x2, Nfft, fs);
title('Signal with narrower bandwidth')

M = 2;
x_down = downsample(x2, M);

% Freq Domain:
psdFromFFT(x_down, Nfft/M, fs/M);
title('Downsampled Signal PSD')

% Time Domain:
figure
subplot(211)
plot((0:(Nfft-1))*ts, x2((Nfft+1):2*Nfft))
xlabel('Time')
ylabel('Amplitude')
title('Comparison in Time Domain')
subplot(212)
plot((0:((Nfft/M) - 1))*(M*ts), x_down((Nfft/M + 1):2*(Nfft/M)), 'r')
xlabel('Time')
ylabel('Amplitude')
title('Comparison in Time Domain')

%% Decimation

M = 2;
x_decimated = decimate(x, M);

% Freq Domain:
psdFromFFT(x_decimated, Nfft/M, fs/M);
title('Decimated Signal PSD')

% Time Domain:
figure
subplot(211)
plot((0:(Nfft-1))*ts, x((Nfft+1):2*Nfft))
xlabel('Time')
ylabel('Amplitude')
title('Comparison in Time Domain')
subplot(212)
plot((0:((Nfft/M) - 1))*(M*ts), x_decimated((Nfft/M + 1):2*(Nfft/M)), 'r')
xlabel('Time')
ylabel('Amplitude')


%% Example FIR Filter

fc = 400;           % cuttof frequency

% Low-pass FIR Filter:
fir_order = 15;
fir_num = fir1(fir_order, fc / (fs/2));
fir_den = 1;

% Frequency response:
figure
freqz(fir_num, fir_den)
title('FIR Filter Freq. Response')

% Interpolated filter:
M = 4;              % upsampling factor
fir_num_up = interp(fir_num, M);

figure
freqz(fir_num_up, fir_den)
title('Interpolated FIR Filter Freq. Response')


%% Filtered signal vs interpolated filtered signal

% Filtered signal
y = filter(fir_num, fir_den, x);

% PSD:
ss
title('Filtered Signal PSD')

% Interpolated filtered signal
y_interp = interp(y, M);

% PSD:
psdFromFFT(y_interp, M*Nfft, M*fs);
title('Interpolated Filtered Signal PSD')

% Compare the filtered signal and its interpolated version
figure
subplot(211)
plot((0:(Nfft-1))*ts, y((Nfft+1):2*Nfft))
xlabel('Time')
ylabel('Amplitude')
title('Comparison in Time Domain')
subplot(212)
plot((0:(M*Nfft - 1))*(ts/M), y_interp((M*Nfft + 1):2*M*Nfft), 'r')
xlabel('Time')
ylabel('Amplitude')

%% Filter interpolated input signal and compared to original filtered signal

M = 4;
x_interp = interp(x, M);

y_oversampled = filter(fir_num_up, 1, x_interp);

% PSD:
psdFromFFT(y_oversampled, M*Nfft, M*fs);
title('Interpolated Filtered Signal PSD')

figure
subplot(211)
plot((0:(Nfft-1))*ts, y((Nfft+1):2*Nfft))
xlabel('Time')
ylabel('Amplitude')
subplot(212)
plot((0:(M*Nfft - 1))*(ts/M), y_oversampled((M*Nfft + 1):2*M*Nfft), 'r')
xlabel('Time')
ylabel('Amplitude')


%% Example from Analog filters

% Analog LP Butterworth filter:
order = 8;
Wn = 2*pi*fc;
[NUM, DEN] = butter(order, Wn,'s');

% Convert Analog to Digital filter:
[NUMd,DENd] = bilinear(NUM, DEN, fs);
figure
freqz(NUMd, DENd)

% Digital Filter with a higher sampling frequency:
[NUMd_over, DENd_over] = bilinear(NUM, DEN, M*fs);
figure
freqz(NUMd_over, DENd_over)

%% Compare

% 1) Regular filtering of the original signal:
y = filter(NUMd, DENd, x);

% 2) Interpolate the "regularly" filtered signal
M = 4;
y_interp = resample(y, M, 1);

% 3) Filter interpolated input signal using the interpolated digital filter
x_interp = interp(x, M);
y_interp2 = filter(NUMd_over, DENd_over, x_interp);


psdFromFFT(y, Nfft, fs);
psdFromFFT(y_interp, M*Nfft, M*fs);
psdFromFFT(y_interp2, M*Nfft, M*fs);

figure
subplot(311)
plot((0:(Nfft-1))*ts, y((Nfft+1):2*Nfft))
xlabel('Time')
ylabel('Amplitude')
title('Regular filtering')
subplot(312)
plot((0:(M*Nfft - 1))*(ts/M), y_interp((M*Nfft + 1):2*M*Nfft), 'r')
xlabel('Time')
ylabel('Amplitude')
title('Interpolated version of the above waveform')
subplot(313)
plot((0:(M*Nfft - 1))*(ts/M), y_interp2((M*Nfft + 1):2*M*Nfft), 'g')
xlabel('Time')
ylabel('Amplitude')
title('Filtering the interpolated input signal using oversampled filter')


