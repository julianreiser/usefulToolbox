%% readiness_potential_quickndirty_modular.m
%
% Quick preprocessing for RP analysis (modular version)
% This script processes merged EEG set files with names of the form:
%   [experiment][condition][repetition]_[subject]_merged.set
% Example: readinessPotentialSit1_002_merged.set
%
% Author: Julian Elias Reiser, 2025-04-10

clear all; close all; clc;  % clear workspace and command window

%% Define Paths

pathVars.inLocation = '/Users/julianreiser/Downloads/mbt_workshop/data/merged/';
pathVars.outLocation = '/Users/julianreiser/Downloads/mbt_workshop/data/prepro/';
pathVars.eegLab = '/Users/julianreiser/owncloud/projects/functions/eeglab2023.1';

%% Modular Parameters
% Define the experiment name, conditions, repetitions, and subject IDs.
experimentName = 'readinessPotential';
conditions     = {'Sit', 'Walk'};   % adjust conditions as needed
repetitions    = [1, 2];             % for example, repetition numbers
subjects       = {'001', '002', '003', '004'};  % subject IDs (as strings)

% Preprocessing parameters
HPF  = 1;    % high-pass filter cutoff (Hz)
LPF  = 30;   % low-pass filter cutoff (Hz)
REJ1 = 3;    % threshold for ICA artifact rejection (if needed)
REJ2 = 5;    % threshold for epoch artifact rejection
FROM = -0.5; % epoch start time (s)
TO   = 0.8;  % epoch end time (s)

% List channels to ignore (e.g., accelerometers, gyros, etc.)
CHANS = 'EEG';

%% Start EEGLAB

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%% Preallocation for ERP structure (optional)
ERPs = struct;

%% Loop over Subjects, Conditions, and Repetitions

for s = 1:length(subjects)
    for c = 1:length(conditions)
        for r = 1:length(repetitions)
            
            % Construct the file name using the modular naming convention:
            % [experimentName][Condition][Repetition]_[Subject]_merged.set
            fileName = sprintf('%s%s%d_%s_merged.set', experimentName, conditions{c}, repetitions(r), subjects{s});

            fprintf('Processing file: %s\n', fileName);
            
            % Load the merged EEG set file
            EEG = pop_loadset('filename', fileName, 'filepath', pathVars.inLocation);
            
            % Channel editing: add channel location info if missing and remove unwanted channels
            besaFile = dir(fullfile(pathVars.eegLab, 'plugins', 'dipfit*' , '**', '*BESA*','standard-10-5-cap385.elp'));
            besaFile = fullfile(besaFile.folder,besaFile.name);
            EEG = pop_chanedit(EEG, 'look up', besaFile);
            
            % find channels that are EEG-related
            chansIdx = find(strcmpi({EEG.chanlocs.type},CHANS));
            
            % copy the data before deletion
            OLD = EEG;
            EEG = pop_select(EEG, 'channel', chansIdx);
            
            % Filtering
            EEG = pop_eegfiltnew(EEG, 'locutoff', HPF);
            EEG = pop_eegfiltnew(EEG, 'hicutoff', LPF);
            
            % Store this processed dataset in ALLEEG
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            %% Epoching Around Specified Triggers
            % For this example we assume the triggers of interest are '11' and '22'
            EVENTS = {'11', '22'};
            % Find reference channels (e.g., those containing 'TP' for mastoid channels)
            refchans = find(~cellfun(@isempty, strfind({EEG.chanlocs.labels}, 'TP')));
            
            EEG = pop_epoch(EEG, EVENTS, [FROM TO], 'epochinfo', 'yes');
            
            % Artifact rejection on epochs (using joint probabilities and kurtosis)
            EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, REJ2, REJ2, 0, 1, 0, [], 0);
            EEG = pop_rejkurt(EEG, 1, 1:EEG.nbchan, REJ2, REJ2, 0, 1, 0, [], 0);
            
            % Re-reference the data (e.g., to TP channels if appropriate)
            EEG = pop_reref(EEG, refchans);
            
            % Update dataset name for clarity
            procName = sprintf('%s%s%d_%s_preproc-epo', experimentName, conditions{c}, repetitions(r), subjects{s});
            EEG.setname = procName;
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            %% (Optional) Run ICA and Remove Artifactual Components
            
            % Run ICA
            EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'interrupt', 'on');
            EEG = pop_iclabel(EEG, 'default');
            
            % Identify bad ICs using ICLabel classifications:
            badICs = find( EEG.etc.ic_classification.ICLabel.classifications(:,2) >= 0.8 | ... % muscle
                           EEG.etc.ic_classification.ICLabel.classifications(:,3) >= 0.8 | ... % eye
                           EEG.etc.ic_classification.ICLabel.classifications(:,6) >= 0.8 );     % other bad components
            
            if ~isempty(badICs)
                EEG = pop_subcomp(EEG, badICs, 0);
                fprintf('Removed ICs: %s\n', mat2str(badICs));
            end
            
            EEG.setname = [procName, '_ICAcorr'];
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            %% Save the Processed Set
            outFileName = sprintf('%s%s%d_%s_prep.set', experimentName, conditions{c}, repetitions(r), subjects{s});
            EEG = pop_saveset(EEG, 'filename', outFileName, 'filepath', pathVars.outLocation);
            
            %% (Optional) Store ERP Data for Further Analysis
            % Store selected fields (e.g., data, times, srate) in a structure.
            erpField = sprintf('sub%s_%s%d', subjects{s}, conditions{c}, repetitions(r));
            ERPs.(erpField).data  = EEG.data;
            ERPs.(erpField).times = EEG.times;
            ERPs.(erpField).srate = EEG.srate;
            ERPs.(erpField).event = EEG.event;
            ERPs.(erpField).chanlocs = EEG.chanlocs;
            
        end % repetition loop
    end % condition loop
end % subject loop

%% Optionally save the ERPs structure for later analysis
save(fullfile(pathVars.outLocation, 'ERPs.mat'), 'ERPs');

load(fullfile(pathVars.outLocation, 'ERPs.mat'), 'ERPs');

% Define experimental grouping parameters
conditionsToAverage = {'Sit', 'Walk'};    % experimental conditions
responseFactors = {'11', '22'};           % response codes
leftElectrodes =  {'C3', 'FC5', 'CP5'};
rightElectrodes = {'C4', 'FC6', 'CP6'};


% Define motor cortex electrode indices

leftIndices = findChans(ERPs,leftElectrodes);   % e.g., C3, FC5, CP5
rightIndices = findChans(ERPs,rightElectrodes);  % e.g., C4, FC6, CP6

% Initialize containers
leftGrand = struct(); rightGrand = struct();
erpFields = fieldnames(ERPs);

% Loop through posture and response combinations
for c = 1:length(conditionsToAverage)
    cond = conditionsToAverage{c};
    condFields = erpFields(contains(lower(erpFields), lower(cond)));

    for r = 1:length(responseFactors)
        resp = responseFactors{r};
        keyName = sprintf('%s%s', cond, resp);

        subjLeft = [];
        subjRight = [];

        for i = 1:length(condFields)
            erpStruct = ERPs.(condFields{i});
            
            if isfield(erpStruct, 'event')  % if behavioral labels are stored
                [~, firstIdx] = unique([erpStruct.event.epoch],'first');
                epochIdx = find(strcmp({erpStruct.event(firstIdx).type}, resp));
            else
                warning('No "resp" field found. Using all epochs.');
                epochIdx = 1:size(erpStruct.data,3);
            end
            if isempty(epochIdx), continue; end

            subjERP = mean(erpStruct.data(:,:,epochIdx), 3);
            subjLeftERP = mean(subjERP(leftIndices, :), 1);
            subjRightERP = mean(subjERP(rightIndices, :), 1);

            subjLeft(end+1, :) = subjLeftERP;
            subjRight(end+1, :) = subjRightERP;
        end

        leftGrand.(keyName) = mean(subjLeft, 1);
        rightGrand.(keyName) = mean(subjRight, 1);
    end
end

% Plotting
timeVec = ERPs.(erpFields{1}).times;
figure;

% Left Motor Cortex
subplot(1,2,1); hold on;
labelsL = fieldnames(leftGrand);
for i = 1:length(labelsL)
    plot(timeVec, leftGrand.(labelsL{i}), 'LineWidth', 1);
end
xlabel('Time (s)'); ylabel('Amplitude (µV)');
title('Left Motor Cortex: Grand Avg ERP'); legend(labelsL, 'Location', 'Best'); grid on;

% Right Motor Cortex
subplot(1,2,2); hold on;
labelsR = fieldnames(rightGrand);
for i = 1:length(labelsR)
    plot(timeVec, rightGrand.(labelsR{i}), 'LineWidth', 1);
end
xlabel('Time (s)'); ylabel('Amplitude (µV)');
title('Right Motor Cortex: Grand Avg ERP'); legend(labelsR, 'Location', 'Best'); grid on;

function chanIdx = findChans(erpStruct,chanNames)
    chanIdx = [];
    erpFields = fieldnames(erpStruct);
    chanlocs = erpStruct.(erpFields{1}).chanlocs;
    for chans = 1:length(chanNames)
        chanIdx(chans) = find(strcmpi({chanlocs.labels},chanNames{chans}))
    end % for chans
end % function ChanIdx
