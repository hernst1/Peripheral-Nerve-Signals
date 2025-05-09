close all

%{
%% Extracting Features over overlapping windows
WSize = [0.05 0.1 0.3]; % window size in s
Olap = [0.0 0.25 0.75]; % overlap percentage

MAV_figs = figure('Visible','off');
VAR_figs = figure('Visible','off');
fig = figure('visible','off');

for i = 1:length(stim)

    for k = 1:length(WSize)
        window = floor(WSize(k)*fs);	% length of each data frame, 30ms
        nOlap = floor(Olap(i)*window);  % overlap of successive frames, half of WSize
        hop = window-nOlap;	    % amount to advance for next data frame

        signal = stim(i).signal;
        labels = stim(i).trigger; % labels of stimulus locations
    
        filteredSignal = filtfilt(b,a,zscore(signal));
    
        nx = length(signal);	            % length of input vector
        len = fix((nx - (window-hop))/hop);	%length of output vector = total frames
        
        % preallocate outputs for speed
        [MAV_feature, VAR_feature, featureLabels] = deal(zeros(1,len));

        Rise1 = gettrigger(labels,0.5); % gets the starting points of stimulations
        Fall1 = gettrigger(-labels,-0.5); % gets the ending points of stimulations
    
        for j = 1:len
            start_window = ((j-1)*hop+1);
            end_window = ((j-1)*hop+window);
    
            segment = filteredSignal(start_window:end_window);
            MAV_feature(j) = mean(abs(segment));
            VAR_feature(j) = mean((segment-mean(segment)).^2);
            
            % re-build the label vector to match it with the feature vector
            featureLabels(j) = sum(arrayfun(@(t) (start_window) >= Rise1(t) && (end_window) <= Fall1(t), 1:length(Rise1)));
        end

        %% Plotting the features
        % Note: when plotting the features, scale the featureLabels to the max of
        % the feature values for proper visualization
        sub_idx = (i-1)*length(WSize)+k;
        % MAV feature
        figure(MAV_figs)
        subplot(length(stim),length(WSize),sub_idx)
        plot(1:len,featureLabels.*MAV_feature,'r');
        hold on;
        plot(1:len,featureLabels.*max(MAV_feature),'m'); % can we assume that the max MAV feature is the max amplitude of stimulation, rather than just the max amplitude of response
        sgtitle('MAV feature extraction for filtered ' + stim_names(i) + ' signal')
        subtitle('WSize = '  + WSize(k) + ' Olap = ' + Olap(i))
        xlabel('Sample #'), ylabel('Amplitude (uV)')
        % VAR feature
        figure(VAR_figs)
        subplot(length(stim),length(WSize),sub_idx)
        plot(1:len,featureLabels.*VAR_feature,'b');
        hold on
        plot(1:len,featureLabels.*max(VAR_feature),'g');
        sgtitle('VAR feature extraction for filtered ' + stim_names(i) + ' signal')
        subtitle('WSize = '  + WSize(k) + ' Olap = ' + Olap(i))
        xlabel('Sample #'), ylabel('Amplitude (uV)')
       

    %% Feature Selection
        snr = 20*log10(mean(MAV_feature)/mean(MAV_rest));
    
    
    end
end

set(MAV_figs,'Visible','on')
set(VAR_figs,'Visible','on')
%}
%% Execution Code
%{
WSize = [0.05 0.1 0.3];
Olap = [0.0 0.25 0.75];
fs = 20000;

[VF_MAV_feature, VF_VAR_feature, VF_featureLabels, ~] = featureExtraction(filtered_VF_signal, VF_labels, Olap, WSize, fs);
featAvgStore = plotFeatures(VF_MAV_feature, VF_VAR_feature, VF_featureLabels, Olap, WSize, 'VF');
VF_SNR = featureSelection(featAvgStore);
disp(VF_SNR);

[Pinch_MAV_feature, Pinch_VAR_feature, Pinch_featureLabels, ~] = featureExtraction(filtered_Pinch_signal, Pinch_labels, Olap, WSize, fs);
featAvgStore = plotFeatures(Pinch_MAV_feature, Pinch_VAR_feature, Pinch_featureLabels, Olap, WSize, 'Pinch');
Pinch_SNR = featureSelection(featAvgStore);
disp(Pinch_SNR);

[Flex_MAV_feature, Flex_VAR_feature, Flex_featureLabels, ~] = featureExtraction(filtered_Flex_signal, Flex_labels, Olap, WSize, fs);
featAvgStore = plotFeatures(Flex_MAV_feature, Flex_VAR_feature, Flex_featureLabels, Olap, WSize, 'Flex');
Flex_SNR = featureSelection(featAvgStore);
disp(Flex_SNR);
%}

%% Execution code: Using structs based on Deland's code

WSize = [0.05 0.1 0.3];
Olap = [0.0 0.25 0.75];
fs = 20000;

stim_names = {'VF', 'Pinch', 'Flex'};
Filtered_signal_stor = {filtered_VF_signal, filtered_Pinch_signal, filtered_Flex_signal};
label_storage = {VF_labels, Pinch_labels, Flex_labels};

% Initialize Feat_av_stor if it doesn't exist
if ~exist('Feat_av_stor', 'var')
    Feat_av_stor = struct();
end

% Loop through each stimulus name
for i = 1:length(stim_names)
    % Create structure for the current stimulus name if it doesn't exist
    if ~isfield(Feat_av_stor, stim_names{i})
        Feat_av_stor.(stim_names{i}) = struct('MAV', struct(), 'VAR', struct(), 'lab', struct());
    end
    
    % Create structure for MAV feature if it doesn't exist
    if ~isfield(Feat_av_stor.(stim_names{i}), 'MAV')
        Feat_av_stor.(stim_names{i}).MAV = struct('stimuli', struct(), 'rest', struct(), 'feat', struct());
    end
    
    % Create structure for VAR feature if it doesn't exist
    if ~isfield(Feat_av_stor.(stim_names{i}), 'VAR')
        Feat_av_stor.(stim_names{i}).VAR = struct('stimuli', struct(), 'rest', struct(), 'feat', struct());
    end
end



for i = 1:length(stim_names)
    figure('units','normalized','Position',[0.1,0.1,0.7,0.6]);
    p = 1;
    for j = 1:length(Olap)
        for k = 1:length(WSize)
            subplot(3,3,p)
            [MAV_feature, VAR_feature, featureLabels] = feat_extract(Filtered_signal_stor{i}, label_storage{i}, WSize(k), Olap(j), fs);

            if (j == 2) && (k == 2)
                Feat_av_stor.(stim_names{i}).MAV.stimuli = mean(MAV_feature(featureLabels == 1));        
                Feat_av_stor.(stim_names{i}).MAV.rest = mean(MAV_feature(featureLabels == 0));
                Feat_av_stor.(stim_names{i}).MAV.feat = MAV_feature;
                Feat_av_stor.(stim_names{i}).lab = featureLabels;

                Feat_av_stor.(stim_names{i}).VAR.stimuli = mean(VAR_feature(featureLabels == 1));
                Feat_av_stor.(stim_names{i}).VAR.rest = mean(VAR_feature(featureLabels == 0));
                Feat_av_stor.(stim_names{i}).VAR.feat = VAR_feature;
                Feat_av_stor.(stim_names{i}).lab = featureLabels;
            end
            
            stem(find(featureLabels == 1), ones (1, length(find(featureLabels == 1) )).*max(max(MAV_feature), max(VAR_feature)), 'Color', 'r', 'Linewidth', 0.2, 'DisplayName', 'Stimulation Labels')
            hold on;
            stem(find(featureLabels == 0), ones (1, length(find(featureLabels == 0) )).*max(max(MAV_feature), max(VAR_feature)), 'Color', 'c', 'Linewidth', 0.2, 'DisplayName', 'Rest Labels')
            hold on;
            plot(1:length(MAV_feature), MAV_feature, 'Linewidth', 2, 'Color', 'k', 'DisplayName', 'MAV Feature')
            hold on;
            plot(1:length(VAR_feature), VAR_feature, 'Linewidth', 2, 'Color', 'b', 'DisplayName', 'VAR Feature')
            grid on; grid minor;
            set(gca, 'YScale', 'log')
            xlabel('Frame Count'), ylabel('Amplitude (uV)')
            sgtitle(['MAV and VAR Features: ' stim_names{i}])
            sbt = subtitle(['WSize: ' num2str(WSize(k)) ', Olap: ' num2str(Olap(j))]);
            sbt.FontSize = 12;
            lgd = legend('Location','east');
            %lgd.FontSize = 11;
            p = p + 1;
            %set(gca, 'FontSize', 15)
        end
    end
end

[SNR_MAV, SNR_VAR] = feat_select(Feat_av_stor, stim_names);
Feature = {'MAV', 'VAR'};

% Display MAV feature table
fprintf('MAV Feature:\n');
fprintf('| Stimulus Type | MAV SNR |\n');
fprintf('|---------------|---------|\n');
for i = 1:length(stim_names)
    fprintf('| %-13s | %-7.3f |\n', stim_names{i}, SNR_MAV(i));
end
fprintf('\n');

% Display VAR feature table
fprintf('VAR Feature:\n');
fprintf('| Stimulus Type | VAR SNR |\n');
fprintf('|---------------|---------|\n');
for i = 1:length(stim_names)
    fprintf('| %-13s | %-7.3f |\n', stim_names{i}, SNR_VAR(i));
end
fprintf('\n');


%% Copy of Deland's (TA) feat_extract function
function [MAV_feature, VAR_feature, featureLabels] = feat_extract(filteredSignal, labels, WSize, Olap, fs)
    nx = length(filteredSignal); % Calculate the length of input signal
    WSize = floor(WSize * fs); % Compute the length of each data frame
    nOlap = floor(Olap * WSize); % Compute the overlap of successive frames
    hop = WSize - nOlap; % Amount to advance for the next data frame
    len = fix((nx - (WSize - hop)) / hop); % Length of output vector = total frames
    
    % Initialize output arrays and preallocate for speed
    MAV_feature = zeros(1, len);
    VAR_feature = zeros(1, len);
    featureLabels = zeros(1, len);
    
    % Get triggers for labels
    Rise1 = gettrigger(labels, 0.5); % Starting points of stimulations
    Fall1 = gettrigger(-labels, -0.5); % Ending points of stimulations
    
    % Extract features
    for k = 1:len
        start_idx = (k - 1) * hop + 1;
        end_idx = (k - 1) * hop + WSize;

        segment = filteredSignal(start_idx:end_idx);
        MAV_feature(k) = mean(abs(segment));
        VAR_feature(k) = mean((segment - mean(segment)).^2);
        
        % Re-build the label vector to match it with the feature vector
        featureLabels(k) = sum(arrayfun(@(t) start_idx >= Rise1(t) && end_idx <= Fall1(t), 1:length(Rise1)));
    end
end
%% Feature Selection via SNR - my attempt
function [SNR_MAV, SNR_VAR] = feat_select(feature_avg_struct, stim_names)
    SNR_MAV = ones(length(stim_names), 1);
    SNR_VAR = ones(length(stim_names), 1);

    for i = 1:length(stim_names)
        SNR_MAV(i) = 20*log10(feature_avg_struct.(stim_names{i}).MAV.stimuli/feature_avg_struct.(stim_names{i}).MAV.rest);
        SNR_VAR(i) = 20*log10(feature_avg_struct.(stim_names{i}).VAR.stimuli/feature_avg_struct.(stim_names{i}).VAR.rest);
    end
end

% - beginning of my code that worked except feature selection
%{

%% Extracting Features over overlapping windows - my code
function [MAV_feature, VAR_feature, featureLabels, len] = featureExtraction(filteredSignal, labels, Olap, WSizes, fs)
    
    numOlap = length(Olap);
    numWSizes = length(WSizes);

    len = 1;
    
    MAV_feature = zeros(numOlap, numWSizes, len);
    VAR_feature = zeros(numOlap, numWSizes, len);
    featureLabels = zeros(numOlap, numWSizes, len);

    for i = 1:numOlap
        for j = 1:numWSizes
            WSize = floor(WSizes(j)*fs);	    % length of each data frame
            nOlap = floor(Olap(i)*WSize);  % overlap of successive frames, half of WSize
            hop = WSize-nOlap;	    % amount to advance for next data frame
            nx = length(filteredSignal);	            % length of input vector
            len = fix((nx - (WSize-hop))/hop);	%length of output vector = total frames
            
            Rise1 = gettrigger(labels,0.5); % gets the starting points of stimulations
            Fall1 = gettrigger(-labels,-0.5); % gets the ending points of stimulations

            for k = 1:len
                start_idx = (k-1)*hop+1;
                end_idx = (k-1)*hop+WSize;

                segment = filteredSignal(start_idx:end_idx);
                MAV_feature(i,j,k) = mean(abs(segment));
                VAR_feature(i,j,k) = mean((segment-mean(segment)).^2);
                
                % re-build the label vector to match it with the feature vector
                featureLabels(i,j,k) = sum(arrayfun(@(t) start_idx >= Rise1(t) && end_idx <= Fall1(t), 1:length(Rise1)));
            end
        end
    end
end

%% Plotting Features for visual feature selection
function featAvgStore = plotFeatures(MAV_feature, VAR_feature, featureLabels, Olap, WSize, stim_name)
    % Plot MAV and VAR features for each combination of Olap and WSize
    figure('units','normalized','Position',[0.1,0.1,0.7,0.4]);
    featAvgStore = zeros(1,4);

    for i = 1:length(Olap)
        for j = 1:length(WSize)

            sub_idx = (i - 1) * length(WSize) + j;
            subplot(length(Olap), length(WSize), sub_idx);
            
            % Extract features for the current combination of Olap and WSize
            current_MAV_feature = squeeze(MAV_feature(i, j, :));
            current_VAR_feature = squeeze(VAR_feature(i, j, :));
            current_featureLabels = squeeze(featureLabels(i, j, :));

            if (i == 2) && (j == 2)
                format long
                featAvgStore(1,1) = mean(current_MAV_feature(find(current_featureLabels == 1)));
                featAvgStore(1,2) = mean(current_MAV_feature(find(current_featureLabels == 0)));
                featAvgStore(1,3) = mean(current_VAR_feature(find(current_featureLabels == 1)));
                featAvgStore(1,4) = mean(current_VAR_feature(find(current_featureLabels == 0)));
                %disp(featAvgStore)
            end

            % Plot stimulation stems
            stem(find(current_featureLabels == 1), ones(1, length(find(current_featureLabels == 1))).*max(max(current_MAV_feature), max(current_VAR_feature)), 'Color', 'r', 'Linewidth', 0.2, 'DisplayName', 'Stimulation Labels');
            hold on;

            % Plot MAV and VAR features
            plot(1:length(current_MAV_feature), current_MAV_feature, 'Linewidth', 2, 'Color', 'k', 'DisplayName', 'MAV Feature');
            hold on;
            plot(1:length(current_VAR_feature), current_VAR_feature, 'Linewidth', 2, 'Color', 'b','DisplayName', 'VAR Feature');
           
            % Set title, axes, and scale
            sgtitle(['MAV and VAR Feature Extraction - ' stim_name])
            subtitle(['Olap = ' num2str(Olap(i)) ', WSize = ' num2str(WSize(j))]);
            xlabel('Frame Count'), ylabel('Amplitude (uV)')
            set(gca, 'YScale', 'log')
        end
    end
    legend('Location','best');
end

%% Feature Selection via SNR
function SNR = featureSelection(featAvgStore)
    SNR = zeros(2,1);
    for i = length(featAvgStore)/2
        format long
        SNR(i) = 20*log10(featAvgStore(i*2-1)/featAvgStore(i*2));
    end
end
% - this is end of my final code that worked except for feature selection
%}
%% Extracting Features over overlapping windows
%{
function [MAV_feature, VAR_feature, featureLabels, len] = featExt(filteredSignal, labels, Olap, WSizes, fs)
    % Determine the dimensions of the output arrays
    num_Olap = length(Olap);
    num_WSizes = length(WSizes);
    
    % Compute len based on the maximum WSize and minimum Olap
    max_WSize = max(WSizes);
    min_Olap = min(Olap);
    len = floor((length(filteredSignal) - max_WSize * fs * (1 - min_Olap)) / (max_WSize * fs * (1 - min_Olap)));

    % Preallocate output arrays
    MAV_feature = zeros(num_Olap, num_WSizes, len);
    VAR_feature = zeros(num_Olap, num_WSizes, len);
    featureLabels = zeros(num_Olap, num_WSizes, len);

    for i = 1:num_Olap
        for j = 1:num_WSizes
            % Compute window parameters
            WSize = floor(WSizes(j) * fs);      % Length of each data frame
            nOlap = floor(Olap(i) * WSize);     % Overlap of successive frames, half of WSize
            hop = WSize - nOlap;                % Amount to advance for next data frame
            
            Rise1 = gettrigger(labels, 0.5);    % Gets the starting points of stimulations
            Fall1 = gettrigger(-labels, -0.5);  % Gets the ending points of stimulations
            
            for k = 1:len
                start_idx = (k - 1) * hop + 1;
                end_idx = (k - 1) * hop + WSize;

                segment = filteredSignal(start_idx:end_idx);
                MAV_feature(i, j, k) = mean(abs(segment));
                VAR_feature(i, j, k) = mean((segment - mean(segment)).^2);
                
                % Rebuild the label vector to match it with the feature vector
                featureLabels(i, j, k) = sum(arrayfun(@(t) start_idx >= Rise1(t) && end_idx <= Fall1(t), 1:length(Rise1)));
            end
        end
    end
end
%}
%******
%{
function [MAV_feature, VAR_feature, featureLabels, len] = featExt(filteredSignal, labels, Olap, WSizes, fs)
    
    MAV_feature = cell(length(Olap), length(WSizes));
    VAR_feature = cell(length(Olap), length(WSizes));
    featureLabels = cell(length(Olap), length(WSizes));

    for i = 1:length(Olap)
        for j = 1:length(WSizes)
            WSize = floor(WSizes(j)*fs);	    % length of each data frame
            nOlap = floor(Olap(i)*WSize);  % overlap of successive frames, half of WSize
            hop = WSize-nOlap;	    % amount to advance for next data frame
            nx = length(filteredSignal);	            % length of input vector
            len = fix((nx - (WSize-hop))/hop);	%length of output vector = total frames
            
            % preallocate outputs for speed
            [MAV_feature{i,j}, VAR_feature{i,j}, featureLabels{i,j}] = deal(zeros(1,len));
            
            Rise1 = gettrigger(labels,0.5); % gets the starting points of stimulations
            Fall1 = gettrigger(-labels,-0.5); % gets the ending points of stimulations
            
            for k = 1:len
                start_idx = (k-1)*hop+1;
                end_idx = (k-1)*hop+WSize;

                segment = filteredSignal(start_idx:end_idx);
                MAV_feature{i,j}(k) = mean(abs(segment));
                VAR_feature{i,j}(k) = mean((segment-mean(segment)).^2);
                
                % re-build the label vector to match it with the feature vector
                featureLabels{i,j}(k) = sum(arrayfun(@(t) start_idx >= Rise1(t) && end_idx <= Fall1(t), 1:length(Rise1)));
            end
        end
    end
end
%}
%******

%% Plotting features
%{
function [MAV_fig, VAR_fig] = plotFeatures(MAV_feature, VAR_feature, featureLabels, Olap, WSize, len, stim_name)
    % Note: when plotting the features, scale the featureLabels to the max of
    % the feature values for proper visualization
    for i = 1:length(Olap)
        for j = 1:length(WSize)

            sub_idx = (i-1)*length(Olap)+j;

            % MAV feature
            MAV_fig = figure('Visible','off');
            figure(MAV_fig)
            subplot(length(Olap),length(WSize),sub_idx)
            stem(find(featureLabels == 1), ones(1, max(featureLabels == 1)) .* max(MAV_feature), 'r');
            hold on;
            stem(find(featureLabels == 0), ones(1, max(featureLabels == 0)) .* max(MAV_feature), 'c');
            hold on;
            plot(1:length(MAV_feature), MAV_feature, 'k');
            %plot(1:len,featureLabels.*MAV_feature,'k');
            hold on;
            %plot(1:len,featureLabels.*max(MAV_feature),'r'); % can we assume that the max MAV feature is the max amplitude of stimulation, rather than just the max amplitude of response
            sgtitle(['MAV feature extraction for filtered ' stim_name ' signal'])
            subtitle(['WSize = ' num2str(WSize(j)) ' Olap = ' num2str(Olap(i))])
            xlabel('Sample #'), ylabel('Amplitude (uV)')

            % VAR feature
            VAR_fig = figure('Visible','off');
            figure(VAR_fig)
            subplot(length(Olap),length(WSize),sub_idx)
            plot(1:len,featureLabels.*VAR_feature,'k');
            hold on
            plot(1:len,featureLabels.*max(VAR_feature),'b');
            sgtitle(['VAR feature extraction for filtered ' stim_name ' signal'])
            subtitle(['WSize = ' num2str(WSize(j)) ' Olap = ' num2str(Olap(i))])
            xlabel('Sample #'), ylabel('Amplitude (uV)')
        end
    end
end
%}
%{
function plotFeatures(MAV_feature, VAR_feature, featureLabels, Olap, WSizes, stim_name)
    % Plot MAV and VAR features
    figure('units','normalized','Position',[0.1,0.1,0.7,0.4]);
    
    num_Olap = size(MAV_feature, 1);
    num_WSizes = size(MAV_feature, 2);
    
    for i = 1:num_Olap
        for j = 1:num_WSizes
            sub_idx = (i - 1) * num_WSizes + j;
            subplot(num_Olap, num_WSizes, sub_idx);
            
            current_MAV_feature = squeeze(MAV_feature(i, j, :));
            current_VAR_feature = squeeze(VAR_feature(i, j, :));
            current_featureLabels = squeeze(featureLabels(i, j, :));
            
            stem(find(current_featureLabels == 1), ones(1, max(current_featureLabels == 1)) .* max(max(current_MAV_feature), max(current_VAR_feature)), 'r');
            hold on;
            stem(find(current_featureLabels == 0), ones(1, max(current_featureLabels == 0)) .* max(max(current_MAV_feature), max(current_VAR_feature)), 'c');
            hold on;
            plot(1:length(current_MAV_feature), current_MAV_feature, 'k');
            hold on;
            plot(1:length(current_VAR_feature), current_VAR_feature);
            sgtitle(['MAV and VAR Feature Extraction - ' stim_name])
            subtitle(['Olap = ' num2str(Olap(i)) ', WSize = ' num2str(WSizes(j))]);
            xlabel('Sample #'), ylabel('Amplitude (uV)')
            grid on;
        end
    end
end
%}

%{
function plotFeaturess(MAV_feature, VAR_feature, featureLabels, Olap, WSize, stim_name)
    % Plot MAV features
    figure('units','normalized','Position',[0.1,0.1,0.7,0.4]);
    
    for i = 1:length(Olap)
        for j = 1:length(WSize)
            sub_idx = (i - 1) * length(WSize) + j;
            subplot(length(Olap), length(WSize), sub_idx);
            
            % Extract features for the current combination of Olap and WSize
            current_MAV_feature = MAV_feature{i, j};
            current_VAR_feature = VAR_feature{i, j};
            current_featureLabels = featureLabels{i, j};

            % Plot MAV feature
            %stem(find(current_featureLabels == 1), ones(1, length(find(current_featureLabels == 1)).*max(max(current_MAV_feature), max(current_VAR_feature))))
            %hold on;
            %stem(find(current_featureLabels == 0), ones(1, length(find(current_featureLabels == 0)).*max(max(current_MAV_feature), max(current_VAR_feature))))
            %hold on;
            plot(1:length(current_MAV_feature), current_featureLabels .* max(current_MAV_feature));
            hold on;
            plot(1:length(current_MAV_feature), current_MAV_feature, 'k');
            hold on;
            plot(1:length(current_VAR_feature), current_featureLabels .* max(current_VAR_feature));
            hold on
            plot(1:length(current_VAR_feature), current_VAR_feature, 'r');
            sgtitle(['MAV and VAR Feature Extraction - ' stim_name])
            subtitle(['Olap = ' num2str(Olap(i)) ', WSize = ' num2str(WSize(j))]);
            set(gca, 'YScale', 'log')
            xlabel('Frame Count'), ylabel('Amplitude (uV)')
            %grid on;
            
        end
    end
end
%}
%******
%{
function plotFeatures(MAV_feature, VAR_feature, featureLabels, Olap, WSize, stim_name)
    % Plot MAV features
    figure('units','normalized','Position',[0.1,0.1,0.7,0.4]);

    for i = 1:length(Olap)
        for j = 1:length(WSize)
            sub_idx = (i - 1) * length(WSize) + j;
            subplot(length(Olap), length(WSize), sub_idx);
            
            % Extract features for the current combination of Olap and WSize
            current_MAV_feature = MAV_feature{i, j};
            current_VAR_feature = VAR_feature{i, j};
            current_featureLabels = featureLabels{i, j};

            % Iterate over elements of current_featureLabels
            for k = 1:length(current_featureLabels)
                % If the label is 1, plot a red stem
                if current_featureLabels(k) == 1
                    stem(k, max(max(current_MAV_feature), max(current_VAR_feature)), 'r', 'Linewidth', 0.2, 'DisplayName', 'Labels for stimulation');
                    hold on;
                % If the label is 0, plot a cyan stem
                elseif current_featureLabels(k) == 0
                    stem(k, max(max(current_MAV_feature), max(current_VAR_feature)), 'c', 'Linewidth', 0.2, 'DisplayName', 'Labels for rest');
                    hold on;
                end
            end

            % Plot MAV and VAR features
            plot(1:length(current_MAV_feature), current_MAV_feature, 'Linewidth', 2, 'Color', 'k', 'DisplayName', 'MAV Feature');
            hold on;
            plot(1:length(current_VAR_feature), current_VAR_feature, 'Linewidth', 2, 'Color', 'b','DisplayName', 'VAR Feature');
            hold on;
            
            % Set title, axes, and scale
            sgtitle(['MAV and VAR Feature Extraction - ' stim_name])
            subtitle(['Olap = ' num2str(Olap(i)) ', WSize = ' num2str(WSize(j))]);
            xlabel('Frame Count'), ylabel('Amplitude (uV)')
            set(gca, 'YScale', 'log')
        end
    end
    legend('Stimulation','Rest','MAV Feature','VAR Feature','Location','best')
end
%}
%******
%{
    % Plot VAR features
    figure('units','normalized','Position',[0.1,0.1,0.7,0.4]);
    
    for i = 1:length(Olap)
        for j = 1:length(WSize)
            sub_idx = (i - 1) * length(WSize) + j;
            subplot(length(Olap), length(WSize), sub_idx);
            
            % Extract features for the current combination of Olap and WSize
            current_VAR_feature = VAR_feature{i, j};
            current_featureLabels = featureLabels{i, j};
            
            % Plot VAR feature
            %plot(1:length(current_VAR_feature), current_featureLabels .* current_VAR_feature, 'k');
            plot(1:length(current_VAR_feature), current_VAR_feature, 'k');
            hold on;
            plot(1:length(current_VAR_feature), current_featureLabels .* max(current_VAR_feature), 'b');
            sgtitle(['VAR Feature Extraction - ' stim_name])
            subtitle(['Olap = ' num2str(Olap(i)) ', WSize = ' num2str(WSize(j))]);
            xlabel('Sample #'), ylabel('Amplitude (uV)')
            grid on;
        end
    end
%}
%{
function [MAV_feature, VAR_feature, featureLabels, len] = featureExtraction(filteredSignal, labels, Olap, WSizes, fs)
    
    numOlap = length(Olap);
    numWSizes = length(WSizes);

    len = 1;
    
    MAV_feature = zeros(numOlap, numWSizes, len);
    VAR_feature = zeros(numOlap, numWSizes, len);
    featureLabels = zeros(numOlap, numWSizes, len);

    for i = 1:numOlap
        for j = 1:numWSizes
            WSize = floor(WSizes(j)*fs);	    % length of each data frame
            nOlap = floor(Olap(i)*WSize);  % overlap of successive frames, half of WSize
            hop = WSize-nOlap;	    % amount to advance for next data frame
            nx = length(filteredSignal);	            % length of input vector
            len = fix((nx - (WSize-hop))/hop);	%length of output vector = total frames
            
            Rise1 = gettrigger(labels,0.5); % gets the starting points of stimulations
            Fall1 = gettrigger(-labels,-0.5); % gets the ending points of stimulations

            for k = 1:len
                start_idx = (k-1)*hop+1;
                end_idx = (k-1)*hop+WSize;

                segment = filteredSignal(start_idx:end_idx);
                MAV_feature(i,j,k) = mean(abs(segment));
                VAR_feature(i,j,k) = mean((segment-mean(segment)).^2);
                
                % re-build the label vector to match it with the feature vector
                featureLabels(i,j,k) = sum(arrayfun(@(t) start_idx >= Rise1(t) && end_idx <= Fall1(t), 1:length(Rise1)));
            end
        end
    end
end

function featAvgStore = plotFeatures(MAV_feature, VAR_feature, featureLabels, Olap, WSize, stim_name)
    % Plot MAV and VAR features for each combination of Olap and WSize
    figure('units','normalized','Position',[0.1,0.1,0.7,0.4]);
    featAvgStore = zeros(1,4);

    for i = 1:length(Olap)
        for j = 1:length(WSize)

            sub_idx = (i - 1) * length(WSize) + j;
            subplot(length(Olap), length(WSize), sub_idx);
            
            % Extract features for the current combination of Olap and WSize
            current_MAV_feature = squeeze(MAV_feature(i, j, :));
            current_VAR_feature = squeeze(VAR_feature(i, j, :));
            current_featureLabels = squeeze(featureLabels(i, j, :));

            if (i == 2) && (j == 2)
                featAvgStore(1,1) = mean(current_MAV_feature(find(current_featureLabels == 1)));
                featAvgStore(1,2) = mean(current_MAV_feature(find(current_featureLabels ~= 1)));
                featAvgStore(1,3) = mean(current_VAR_feature(find(current_featureLabels == 1)));
                featAvgStore(1,4) = mean(current_VAR_feature(find(current_featureLabels ~= 1)));
            end

            % Plot stimulation stems
            stem(find(current_featureLabels == 1), ones(1, length(find(current_featureLabels == 1))).*max(max(current_MAV_feature), max(current_VAR_feature)), 'Color', 'r', 'Linewidth', 0.2, 'DisplayName', 'Stimulation Labels');
            hold on;

            % Plot MAV and VAR features
            plot(1:length(current_MAV_feature), current_MAV_feature, 'Linewidth', 2, 'Color', 'k', 'DisplayName', 'MAV Feature');
            hold on;
            plot(1:length(current_VAR_feature), current_VAR_feature, 'Linewidth', 2, 'Color', 'b','DisplayName', 'VAR Feature');
           
            % Set title, axes, and scale
            sgtitle(['MAV and VAR Feature Extraction - ' stim_name])
            subtitle(['Olap = ' num2str(Olap(i)) ', WSize = ' num2str(WSize(j))]);
            xlabel('Frame Count'), ylabel('Amplitude (uV)')
            set(gca, 'YScale', 'log')
        end
    end
    legend('Location','best');
end

function SNR = featureSelection(featAvgStore)
    SNR = zeros(2,1);
    for j = length(featAvgStore)/2
        SNR(j) = 20*log10(featAvgStore(j*2-1)/featAvgStore(j*2));
    end
end
%}
