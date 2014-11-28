function [ x_serial ] = generateOFDMSymbol( Nfft, nSymbols, fl, fup )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

M = 16;         % QAM constellation

% OFDM Modulation

% Modulator Object:
mod     =   modem.qammod('M', M, 'SymbolOrder', 'Gray');

% Preallocate Positive side of the spectrum
X_pos = zeros(Nfft/2 +1, nSymbols);

% Use tones within a limited bandwidth:
usedTones = floor(fl*(Nfft/2)) : floor(fup*(Nfft/2));

% Modulate:
X_pos(usedTones, :) = mod.modulate(randi(M, length(usedTones), nSymbols) - 1);

% Hermitian symmetry:
X_hermitian = [X_pos; flipud(conj(X_pos(2:end-1, :)))];

% Force DC and Nyquist to be real valued:
X_hermitian(1) = abs(X_hermitian(1));
X_hermitian(Nfft/2 + 1) = abs(X_hermitian(Nfft/2 + 1));

% Ifft:
x = ifft(X_hermitian);

% P/S conversion:
x_serial = x(:);

end

