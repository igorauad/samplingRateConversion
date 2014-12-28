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
f_lo = 200;
f_hi = 3500;
x = generateOFDMSymbol( Nfft, nSymbols, f_lo, f_hi, fs);
psdFromFFT(x, Nfft, fs);

%% Upsampling

M = 4;
x_up = upsample(x, M);

% Freq Domain:
psdFromFFT(x_up, M*Nfft, M*fs);
title('Upsampled Signal PSD')

% Time Domain comparison:
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

% NOTE: observe the psd attenuation by a factor of M in the upsampled 
% signal. e.g. for M = 4, the psd reduces by 20*log10(4) ~= 12dB

%% Interpolation

M = 4;
[x_interp, b] = interp(x, M); % b is the vector with the taps from the
                              % anti-imaging filter

% Freq Domain:
psdFromFFT(x_interp, M*Nfft, M*fs);
title('Interpolated Signal PSD')

% Time Domain comparison:
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

% Anti-imaging filter frequency response
figure
freqz(b)
title('Anti-imaging filter')


%% Downsampling with aliasing

% Original signal to be downsampled:
psdFromFFT(x, Nfft, fs);
% Check whether its bandwidth goes beyond (fs/2)/M, case in which aliasing
% will occur when downsampling.
% From another point of view, the new sampling frequency fs/M (after
% downsampling) still has to be greater than twice the signal bandwith.
% Hence, aliasing will only be avoided with the original signal was
% oversampled, case in which there is room for some sampling frequency
% reduction.

M = 2;  % Downsampling factor
x_down = downsample(x, M);

% Freq Domain:
psdFromFFT(x_down, Nfft/M, fs/M);
title('Downsampled Signal PSD')

% Time Domain comparison:
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

%% Downsampling without aliasing

% Downsample a signal 'x2', whose bandwidth is narrower:
f_lo = 200;
f_hi = 1000;
x2 = generateOFDMSymbol( Nfft, nSymbols, f_lo, f_hi, fs);
psdFromFFT(x2, Nfft, fs);
title('Signal with narrower bandwidth')

M = 2;  % Downsampling factor
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

M = 2;  % Downsampling factor
x_decimated = decimate(x2, M);

% Freq Domain:
psdFromFFT(x_decimated, Nfft/M, fs/M);
title('Decimated Signal PSD')

% Time Domain:
figure
subplot(211)
plot((0:(Nfft-1))*ts, x2((3*Nfft+1):4*Nfft))
xlabel('Time')
ylabel('Amplitude')
title('Comparison in Time Domain')
subplot(212)
plot((0:((Nfft/M) - 1))*(M*ts), x_decimated((3*Nfft/M + 1):4*(Nfft/M)), 'r')
xlabel('Time')
ylabel('Amplitude')


