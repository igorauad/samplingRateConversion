function [ x_serial ] = generateOFDMSymbol(Nfft, nSymbols, f_lo, f_hi, fs)
% Generate a random Hermitian Symmetric OFDM symbol
% ------------------------------------------------
%   [ x_serial ] = generateOFDMSymbol( Nfft, nSymbols, fl, fup )
%
%   
%   Inputs:
%       Nfft        -   FFT Size
%       nSymbols    -   Number of Symbols to Generate
%       f_lo        -   Lowest normalized frequency to be used 
%       f_hi        -   Highest normalized frequency to be used
%
%   Output:
%       x_serial    -   Serialized time domain OFDM symbols


M = 16;         % QAM constellation

% Check if frequencies are valid:
if (f_lo > f_hi)
   error('Lower freq. can not be higher the Highest freq.'); 
end

if (f_lo > (fs/2) || f_hi > (fs/2))
   error('Frequency higher can not be higher than Nyquist.'); 
end

% Normalize Upper and Lower bounds:
f_lo = f_lo / (fs/2);
f_hi = f_hi / (fs/2);
% Note: this converts the frequencies to the range (0, 1), where 0 means DC
% and 1 means Nyquist.

% OFDM Modulation

% Modulator Object:
mod = modem.qammod('M', M, 'SymbolOrder', 'Gray');

% Preallocate Positive side of the spectrum
X_pos = zeros(Nfft/2 +1, nSymbols);

% Use tones within a limited bandwidth:
usedTones = (floor(f_lo*(Nfft/2 + 1)) + 1) : floor(f_hi*(Nfft/2 + 1));

% Modulate:
X_pos(usedTones, :) = mod.modulate(randi(M, length(usedTones), nSymbols) - 1);

% Hermitian symmetry:
X_hermitian = [X_pos; flipud(conj(X_pos(2:end-1, :)))];

% Force DC and Nyquist to be real valued:
X_hermitian(1) = abs(X_hermitian(1));
X_hermitian(Nfft/2 + 1) = abs(X_hermitian(Nfft/2 + 1));

% Ifft:
x = ifft(X_hermitian, Nfft);

% P/S conversion:
x_serial = x(:);

end

