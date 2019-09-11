% Check that convolution in time domain performed via the convolution
% theorem is time-sensitive (uncomment line 10 to see that the result actually changes for the inverse signal)
%% define signal
srate = 1000;
time = -2:1/srate:2;

f  = [4 10];
ff = linspace(f(1),f(2),length(time));
s1 = sin(2*pi*ff.*time);
% s1 = s1(end:-1:1);

% freq = 5;
% s1 = sin(2*pi*freq*time);
% coef = 1/400*(1:length(time));
% s1 = coef.*s1;

figure(1)
plot(time, s1)

%% define wavelet

wavelet_freq  = 10; % frequency of wavelet, in Hz

% create complex sine wave
sine_wave = exp( 1i*2*pi*wavelet_freq*time );

% create Gaussian window
s = 7 / (2*pi*wavelet_freq); % this is the standard deviation of the gaussian
gaus_win  = exp( (-time.^2) ./ (2*s^2) );

% now create Morlet wavelet
cmw  = sine_wave .* gaus_win;

%% define convolution parameters

nData = length(s1);
nKern = length(cmw);
nConv = nData + nKern - 1;
half_wav = floor( length(cmw)/2 )+1;

%% get fourier components

% FFT of wavelet, and amplitude-normalize in the frequency domain
cmwX = fft(cmw,nConv);
cmwX = cmwX ./ max(cmwX);

% FFT of data
dataX = fft(s1,nConv);

% now for convolution...
conv_res = dataX.*cmwX;
as = ifft(conv_res);

% cut 1/2 of the length of the wavelet from the beginning and from the end
as = as(half_wav-1:end-half_wav);

% compute hz for plotting
hz = linspace(0,srate/2,floor(nConv/2)+1);

%% some plots...

figure(2), clf

% plot power spectrum of data
subplot(311)
plot(hz,2*abs(dataX(1:length(hz))/length(s1)))
title('Power spectrum of data')
xlabel('Frequency'), ylabel('Amplitude')
set(gca,'xlim',[0 30])

% plot power spectrum of wavelet
subplot(312)
plot(hz,abs(cmwX(1:length(hz))))
title('Power spectrum of wavelet')
xlabel('Frequency'), ylabel('Amplitude')
set(gca,'xlim',[0 30])

% plot power spectrum of convolution result
subplot(313)
plot(hz,2*abs(conv_res(1:length(hz))/length(s1)))
title('Power spectrum of convolution result')
xlabel('Frequency'), ylabel('Amplitude')
set(gca,'xlim',[0 30])

%% Plot convolution in time domain

figure(3), clf
% plot the filtered signal (projection onto real axis)
subplot(311)
plot(time,real(as))
xlabel('Time (s)'), ylabel('Amplitude (\muV)')
set(gca,'xlim',[-2 2])


% plot power (squared magnitude from origin to dot-product location in
% complex space)
subplot(312)
plot(time,abs(as).^2)
xlabel('Time (s)'), ylabel('Power \muV^2')
set(gca,'xlim',[-2 2])


% plot phase (angle of vector to dot-product, relative to positive real
% axis)
subplot(313)
plot(time,angle(as))
xlabel('Time (ms)'), ylabel('Phase (rad.)')
set(gca,'xlim',[-2 2])