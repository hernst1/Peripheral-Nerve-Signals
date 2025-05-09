close all

%% Import Data
load('data.mat');
%{
%% Vectorize Structs
stim = [VF Pinch Flex];
stim_names = ["VF" "Pinch" "Flex"];

%% Raw signal data visualization
figure('units','normalized','Position',[0.1,0.1,0.7,0.4])
    
for i = 1:length(stim)
    % Assign fields and triggers
    signal = stim(i).signal;
    labels = stim(i).trigger;
    TRIG = gettrigger(labels,0.5);
    TRIGend = gettrigger(-labels,-0.5);

    % Plot raw signals
    subplot(length(stim), 1, i)
    plot((1:length(signal))./fs,zscore(signal));
    hold on;
    plot((1:length(signal))./fs,zscore(labels),'y');
    stem(TRIG./fs,ones(length(TRIG),1)*max(zscore(labels)),'Color','g');
    stem(TRIGend./fs,ones(length(TRIG),1)*max(zscore(labels)),'Color','r');
    grid on; grid minor;
    xlim([0,length(signal)./fs]), ylim([-7.5,7.5])
    xlabel('Time (s)'), ylabel('Amplitude (uV)')
    title('Raw ' + stim_names(i) + ' signal with labels for stimulation periods')
end

%% PSD Estimates
figure('units','normalized','Position',[0.1,0.1,0.5,0.5])
colors = ['b','r','g']';
hold on;
for i = 1:length(stim)
    % Assign fields
    signal = stim(i).signal;
    labels = stim(i).trigger;

    % Identify interesting signal
    [rows_act,cols_act,values_act] = find(labels > 0);
    [rows_rest1,cols_rest,values_rest] = find(labels == 0);
    notOfInterest = signal(rows_rest1);
    signalOfInterest = signal(rows_act);

    % Calculate and plot PSD
    h = spectrum.welch; % creates the Welch spectrum estimator
    SOIf=psd(h,signalOfInterest,'Fs',fs); % calculates and plot the one sided PSD
    plot(SOIf); % Plot the one-sided PSD
    temp = get(gca);
    temp.Children(i).Color = colors(i);
    hold on
end
legend(stim_names,'Location','best')

%% Bandpass Filtering

fc1 = 800; % first cutoff frequency in Hz 
fc2 = 2200; % second cutoff frequency in Hz 

% normalize the frequencies
Wp = [fc1 fc2]*2/fs;

% Build a Butterworth bandpass filter of 4th order
% check the "butter" function in matlab

[b,a] = butter(4,Wp,"bandpass");

% Filter data of both classes with a non-causal filter
% Hint: use "filtfilt" function in MATLAB

figure('units','normalized','Position',[0.1,0.1,0.7,0.4])
for i = 1:length(stim)
    % Assign fields and triggers
    signal = stim(i).signal;
    labels = stim(i).trigger;
    TRIG = gettrigger(labels,0.5);
    TRIGend = gettrigger(-labels,-0.5);

    filteredSignal = filtfilt(b,a,zscore(signal));
    
    subplot(length(stim),1,i)
    plot((1:length(signal))./fs,filteredSignal);
    hold on;
    plot((1:length(signal))./fs,zscore(labels),'y');
    stem(TRIG./fs,ones(length(TRIG),1)*max(zscore(labels)),'Color','g');
    stem(TRIGend./fs,ones(length(TRIG),1)*max(zscore(labels)),'Color','r');
    grid on; grid minor;
    xlim([0,length(signal)./fs]),ylim([-7.5,7.5])
    xlabel('Time (s)'),ylabel('Amplitude (uV)')
    title('BP Filtered '+ stim_names(i) +' signal with labels for stimulation periods')
end

%% Filtered PSD

figure('units','normalized','Position',[0.1,0.1,0.5,0.5])
%colors = ['b','r','g'];
hold on;
for i = 1:length(stim)
    % Assign fields
    signal = stim(i).signal;
    labels = stim(i).trigger;

    % Identify interesting signal
    notOfInterest = filteredSignal(rows_rest1);
    signalOfInterest = filteredSignal(rows_act);

    % Calculate and plot PSD
    h = spectrum.welch; % creates the Welch spectrum estimator
    SOIf=psd(h,signalOfInterest,'Fs',fs); % calculates and plot the one sided PSD
    plot(SOIf); % Plot the one-sided PSD.
    %temp = get(gca);
    %temp.Children(i).Color = colors(i);
end
legend(stim_names,'Location','best')

%% 
%}

%% Execution Code
fs = 20000;
fc1 = 800;
fc2 = 2200;
order = 4;

% Raw signals
VF_signal = VF.signal;
VF_labels = VF.trigger;
subplot(3,1,1)
plotSignal(gca, VF_signal, VF_labels, fs, 'VF', 'Raw ');
%hold on
%plotSignal(gca, bpFilter(VF_signal,fc1,fc2,order,fs), VF_labels, fs, 'VF')

Pinch_signal = Pinch.signal;
Pinch_labels = Pinch.trigger;
subplot(3,1,2)
plotSignal(gca, Pinch_signal, Pinch_labels, fs, 'Pinch', 'Raw ')
%hold on
%plotSignal(gca, bpFilter(Pinch_signal,fc1,fc2,order,fs), Pinch_labels, fs, 'Pinch')

Flex_signal = Flex.signal;
Flex_labels = Flex.trigger;
subplot(3,1,3)
plotSignal(gca, Flex_signal, Flex_labels, fs, 'Flex', 'Raw ')
%hold on
%plotSignal(gca, bpFilter(Flex_signal,fc1,fc2,order,fs), Flex_labels, fs, 'Flex')
set(gcf,'units','normalized','Position',[0.1,0.1,0.7,0.6])

% PSD
    % PSD stim
figure('units','normalized','Position',[0.1,0.1,0.5,0.5])
psdSignal(VF_signal, VF_labels, fs, 'r');
hold on;
psdSignal(Pinch_signal, Pinch_labels, fs, 'b');
hold on;
psdSignal(Flex_signal, Flex_labels, fs, 'c');
hold on;

    % Attempt at PSD rest - why are the rest psd different? should i make the
    % assumption that they are the same? should i average them?
[rows_act,~,~] = find(VF_labels>0);
[rows_rest1,~,~] = find(VF_labels==0);
VF_rest = VF_signal(rows_rest1);
VF_act = VF_signal(rows_act);
h = spectrum.welch; % creates the Welch spectrum estimator
SOIf=psd(h,VF_rest,'Fs',fs); % calculates and plot the one sided PSD
%f = SOIf.Frequencies * 1000; % why does this not work for converting
%kHz to Hz
plot(SOIf); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'k';
xlabel('Frequency (kHz)');

    % PSD filtered VF
hold on;
psdSignal(filtered_VF_signal, VF_labels, fs, 'm');

legend({'VF', 'Pinch', 'Flex', 'VF Rest', 'VF Filtered'})

% Filtering
filtered_VF_signal = bpFilter(VF_signal,fc1,fc2,order,fs);
filtered_Pinch_signal = bpFilter(Pinch_signal,fc1,fc2,order,fs);
filtered_Flex_signal = bpFilter(Flex_signal,fc1,fc2,order,fs);

    % BP filtered signal plots without overlay
figure
subplot(3,1,1)
plotSignal(gca,filtered_VF_signal, VF_labels, fs, 'VF', 'BP Filtered ')
subplot(3,1,2)
plotSignal(gca, filtered_Pinch_signal, Pinch_labels, fs, 'Pinch', 'BP Filtered ')
subplot(3,1,3)
plotSignal(gca, filtered_Flex_signal, Flex_labels, fs, 'Flex', 'BP Filtered ')
set(gcf,'units','normalized','Position',[0.1,0.1,0.7,0.6])

    % BP filtered VF signal plot with overlay
figure
plotSignal(gca, VF_signal, VF_labels, fs, 'VF', 'Raw ');
hold on
plotSignal(gca, bpFilter(VF_signal,fc1,fc2,order,fs), VF_labels, fs, 'VF', 'BP')
legend({'Raw VF signal','Filtered VF signal'},'Location','best')

%% Plotting signals
function fig = plotSignal(ax, signal, labels, fs, stim_name, signal_type)
    TRIG = gettrigger(labels,0.5);
    TRIGend = gettrigger(-labels,-0.5);

    axes(ax);
    figure('units','normalized','Position',[0.1,0.1,0.7,0.4]);
    plot((1:length(signal))./fs,zscore(signal));
    hold on;
    plot((1:length(signal))./fs,zscore(labels),'y');
    stem(TRIG./fs,ones(length(TRIG),1) .* max(zscore(labels)),'Color','g');
    stem(TRIGend./fs,ones(length(TRIG),1) .* max(zscore(labels)),'Color','r');
    grid on; grid minor;
    xlim([0,length(signal)./fs])
    %ylim([-7.5,7.5])
    xlabel('Time (s)')
    ylabel('Amplitude (uV)')
    %title([signal_type stim_name ' signal with labels for stimulation periods'])
    % title for overlaying BP over raw
    title(['BP ' stim_name ' signal overlayed on raw ' stim_name ' signal'])
    fig = [];
end
%% PSD estimates
function psd_plot = psdSignal(signal, labels, fs, color)
    %figure('units','normalized','Position',[0.1,0.1,0.5,0.5])
    [rows_act,cols_act,values_act] = find(labels>0);
    [rows_rest1,cols_rest,values_rest] = find(labels==0);
    notOfInterest = signal(rows_rest1);
    signalOfInterest=signal(rows_act);
    h = spectrum.welch; % creates the Welch spectrum estimator
    SOIf=psd(h,signalOfInterest,'Fs',fs); % calculates and plot the one sided PSD
    %f = SOIf.Frequencies * 1000; % why does this not work for converting
    %kHz to Hz
    plot(SOIf); % Plot the one-sided PSD. 
    temp =get(gca);
    temp.Children(1).Color = color;
    title('Power Spectral Density Estimate for [VF, Pinch, Flex Stimuli], Rest, Filtered VF')
    xlabel('Frequency (kHz)');

    psd_plot = [];
end
%% Bandpass filtering
function bpSignal = bpFilter(signal, fc1, fc2, order, fs)
    Wp = [fc1 fc2]*2/fs;
    [b,a] = butter(order,Wp,"bandpass");
    bpSignal = filtfilt(b,a,signal);
end

%{
%% Example: Plot the raw flex signal
flex_signal=Flex.signal;
flex_labels=Flex.trigger;
flex_TRIG = gettrigger(flex_labels,0.5);
flex_TRIGend = gettrigger(-flex_labels,-0.5);

figure('units','normalized','Position',[0.1,0.1,0.7,0.4])
plot((1:length(flex_signal))./fs,zscore(flex_signal));
hold on;
plot((1:length(flex_signal))./fs,zscore(flex_labels),'y');
stem(flex_TRIG./fs,ones(length(flex_TRIG),1)*max(zscore(flex_labels)),'Color','g');
stem(flex_TRIGend./fs,ones(length(flex_TRIG),1)*max(zscore(flex_labels)),'Color','r');
grid on; grid minor;
xlim([0,length(flex_signal)./fs])
xlabel('Time (s)')
ylabel('Amplitude (uV)')
title('Raw Flex signal with labels for stimulation periods')
%}
%{
%% Example: PSD estimates
figure('units','normalized','Position',[0.1,0.1,0.5,0.5])
[rows_act,cols_act,values_act] = find(labels>0);
[rows_rest1,cols_rest,values_rest] = find(labels==0);
notOfInterest = signal(rows_rest1);
signalOfInterest=signal(rows_act);
h = spectrum.welch; % creates the Welch spectrum estimator
SOIf=psd(h,signalOfInterest,'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'b';
%}

