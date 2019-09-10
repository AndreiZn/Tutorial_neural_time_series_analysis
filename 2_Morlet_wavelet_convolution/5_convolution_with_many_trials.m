% mikexcohen.com
% Process many trials (combine into one long trial and then perform convolution)
%% Load data, define variables

load sampleEEGdata.mat;

% frequency parameters
min_freq = 2;
max_freq = 30;
num_frex = 40;
frex = linspace(min_freq,max_freq,num_frex);

% which channel to plot
channel2use = 'o1';

% other wavelet parameters
range_cycles = [4 10];
s = logspace(log10(range_cycles(1)), log10(range_cycles(end)), num_frex) ./ (2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;

%% Naive solution

% Naive solution is to iterate over channels (only one channel here),
% frequencies and trials. However, this approach can be optimized:
% 1) Combine data trials into one long "super trial", process it and then
% split it back into trials
% 2) Perform FFT of EEG.data for specific channel once, i.e. outside of the
% loop over frequencies

%% Process many trials optimally

% FFT parameters
nWave = length(wavtime);
nData = EEG.pnts * EEG.trials;
nConv = nWave + nData - 1;

% initialize output time-frequency data
tf = zeros(length(frex),EEG.pnts);

% now compute the FFT of all trials concatenated
alldata = reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,[]);
dataX = fft(alldata, nConv);

% loop over frequencies 
for fi=1:length(frex)
    % create wavelet and get its FFT
    % the wavelet doesn't change on each trial...
    wavelet = exp(1i*2*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX);
    
    % now run convolution in one step
    as = ifft(waveletX .* dataX);
    as = as(half_wave+1:end-half_wave);
    
    % and reshape back to time X trials
    as = reshape(as, EEG.pnts, EEG.trials);
    
    % compute power and average over trials
    tf(fi,:) = mean(abs(as).^2, 2); 
end

%% Plots

figure(1), clf
contourf(EEG.times,frex,tf,40,'linecolor','none')
set(gca, 'CLim', [0,5], 'ydir','normal','xlim',[-300,1000])

figure(2), clf
plot(EEG.times,tf(5,:))
set(gca,'xlim',[-300,1200])

title(['Power at ', num2str(round(10*frex(5))/10), ' Hz'])
xlabel('Time (ms)'), ylabel('Power (\muV^2)')
