function [ Sk_dbmHz, f, plotHnd ] = psdFromFFT( y, Nfft, fs)
% Simple PSD estimation using the FFT
% ------------------------------------
%   [ Sk_dbmHz, f, plotHnd ] = psdFromFFT( y, Nfft, fs)
%
%   This function implements a naive PSD estimation using the analytical
%   expressions with the DFT. Overlapping and windowing are not used. This
%   PSD estimation is very good when the signals being analyzed contain
%   exactly integer number of periods, which is the case for OFDM symbols.
%
%
%   Input:
%      y       ->  time-domain signal for which the PSD is to be estimated
%      Nfft    ->  FFT length
%      fs      ->  Sampling frequency 
% Note: 'fs' is used only for scaling the result appropriately
%
%   Output:
%      Sk_dbmHz ->  Estimated PSD in dBm/Hz
%      f        ->  frequency vector (axis)
%      plotHnd  ->  Plot Handle

deltaf = fs / Nfft;                     % Tone Spacing
nBlocks = floor(length(y) / Nfft);      % Number of FFT blocks

% S/P conversion:
y_splitted = reshape(y(1:(nBlocks*Nfft)), Nfft, nBlocks);

% FFT:
Y_fft = fft(y_splitted, Nfft);

% PSD estimator (Watts/tone):
Sk = mean( abs(Y_fft).^2, 2) / Nfft;

% Convert PSD estimator to Watts/Hz:
Sk_wattsHz = Sk / fs; 

% Convert PSD estimator to dB/Hz:
Sk_dbHz = 10*log10(Sk_wattsHz);     

% Convert PSD estimator to dBm/Hz:
Sk_dbmHz = Sk_dbHz + 30;          

% Analog frequency vector
f = (0:(Nfft-1))*deltaf;

if (nargout == 0 || nargout == 3)
    
    if(nargout == 3)
        figure('visible','off');
        plotHnd = gcf;
    else
        figure
    end
    
    % Check if better to display in Hz, kHz or MHz
    f_max = f(end);
    if (f_max >= 1e9)
        freq_vector = f / 1e9;
        x_label = 'Analog Frequency (GHz)';
        axis_f_max = fs/2 / 1e9;
        axis_deltaf = fs/8/1e9;
        
    elseif (f_max >= 1e6 && f_max < 1e9)
        freq_vector = f / 1e6;
        x_label = 'Analog Frequency (MHz)';
        axis_f_max = fs/2 / 1e6;
        axis_deltaf = fs/8/1e6;
        
    elseif (f_max >= 1e3 && f_max < 1e6)
        freq_vector = f / 1e3;
        x_label = 'Analog Frequency (kHz)';
        axis_f_max = fs/2 / 1e3;
        axis_deltaf = fs/8/1e3;
        
    else
        freq_vector = f;
        x_label = 'Analog Frequency (Hz)';
        axis_f_max = fs/2;
        axis_deltaf = fs/8;
        
    end     
    plot(freq_vector, Sk_dbmHz)    
    xlabel(x_label)
    ylabel('PSD (dbm/Hz)')
    % Limit freq axis       
    set(gca, 'XLim', [0 axis_f_max]);
    set(gca, 'XTick', 0:axis_deltaf:axis_f_max);    
    grid on
end

end

