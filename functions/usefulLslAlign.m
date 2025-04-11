%% Merging script
% This script is going to merge all files that you recorded for all
% subjects within a specified folder. The only file that we need for sure
% is the .xdf file.
% Please make sure that you include the .xdf (labrecorder), the .edf (faros),
% and the bdf (mbt amplifier). In a folder within the subject folder,
% please include the exported pupil labs folder for the subject.
% NOTE: if you want to include pupil labs, make sure that you run the
% python script to create a mergable file before.
%
% Author: Julian Elias Reiser, 2025-04-10

%% Define paths
clear all
close all

% indicate folder with subjects
pathVars.inLocation = '/Users/julianreiser/Downloads/mbt_workshop/data/melanieData/';
pathVars.outLocation = '/Users/julianreiser/Downloads/mbt_workshop/data/merged/';

% indicate how your experiment is named
exp = 'readinessPotentialWalk2';

% indicate if we need to skip any subjects
subjects2skip = {};

% this determines subjects on its own
subjects        = struct();
subjects.XDF = dir(fullfile(pathVars.inLocation, [exp '_*']));
subjects.XDF = regexp({subjects.XDF.name}, '_(.*)$', 'tokens');
subjects.XDF = cellfun(@(x) x{1}{1}, subjects.XDF, 'UniformOutput', false);

subjects.COR = dir(fullfile(pathVars.outLocation, [exp '_*merged.set']));
subjects.COR = cellfun(@(x) x{1}{1}, arrayfun(@(x) regexp(x.name, '_([0-9]+)_', 'tokens'), subjects.COR, 'UniformOutput', false),'UniformOutput', false);

%% iterate through the participants
for dataset = 1:length(subjects.XDF)
    
    % determine the subject
    subject = subjects.XDF{dataset};
    disp(['... Computing subject ' subject(end-1:end) ' ...'])
    subfolder = [exp '_' subject];

    % check if already computed
    if any(strcmpi(subjects2skip,subject)) 
        disp('... skipping subject as indicated. Continuing ...')
        continue
    elseif any(strcmpi(subjects.COR,subject))
        disp('... xdf data have already been merged. Continuing ...')
        continue
    else
        %% load xdffile
        % define path to subject files
        subPath = fullfile(pathVars.inLocation,subfolder);
        
        % find names of all necessary files
        xdfFile = dir(fullfile(subPath,[exp '_' subject '.xdf']));
        xdfFileName = xdfFile.name;
        
        % check if edfFile exists
        if ~isempty(dir(fullfile(subPath,[exp '*.edf'])))
            edfFile = dir(fullfile(subPath,[exp '_*.edf']));
            edfFileName = edfFile.name;
            disp('... .edf file present found in the subject folder...')
            ecgData = edfread(fullfile(subPath, edfFileName));
            ecgInfo = edfinfo(fullfile(subPath, edfFileName));
        else disp('... no .edf file present in the subject folder...')
        end
        
        % check if bdfFile exists
        if ~isempty(dir(fullfile(subPath,[exp '*.bdf'])))
            bdfFile = dir(fullfile(subPath,[exp '_*.bdf']));
            bdfFileName = bdfFile.name;
            disp('... .bdf file present found in the subject folder...')
        else disp('... no .bdf file present in the subject folder...')
            bdfFileName = [];
        end

        % load xdf file
        if isempty(xdfFileName)
            disp('... no EEG data present in the subject folder. Continuing ...')
            continue
        elseif ~isempty(bdfFileName)
            disp(['... Loading EEG files for subject ' subject ' from amlifier data...'])
            EEG = pop_biosig(fullfile(subPath,bdfFileName));
            streams = load_xdf(fullfile(subPath,xdfFileName));
            xdfEEG = pop_loadxdf(fullfile(subPath,xdfFileName));
        elseif ~isempty(xdfFileName) & isempty(bdfFileName)
            disp(['... Loading EEG files for subject ' subject ' from xdf...'])
            EEG = pop_loadxdf(fullfile(subPath,xdfFileName));
            streams = load_xdf(fullfile(subPath,xdfFileName));
        end % if isempty
        
        %% find names of streams
        for stream = 1:length(streams)
            lslnames{stream} = streams{1,stream}.info.name;
        end % for stream

        % if no eeg stream is found, just skip the subject
        if ~any(contains([lslnames{:}], {'Counter - ','Sync'}, 'IgnoreCase', true)) & ~any(contains([lslnames{:}], 'EEG', 'IgnoreCase', true))
            disp('... no counter stream or EEG stream present in the xdf streams. Continuing ...')
            continue
        end % if any

        %% adding the lsl_times to the EEG struct
        % if the EEG was written directly into the xdf; no counter
        if isempty(bdfFileName) & any(contains([lslnames{:}], 'EEG', 'IgnoreCase', true))

            % find indices in lsl_times to add to EEG struct
            for x = 1:length(streams)
                if isfield(streams{x}.info,'name') && contains(streams{x}.info.name, 'EEG')
                    disp(['... Found EEG in streams'])
                    xdfEEG = streams{x};
                    break;
                end
            end % for x = 1:length(streams)

            % check if lsl_times are already written
            if isfield(EEG, 'lsl_times') & length(xdfEEG.time_stamps) == length(EEG.times) % EEG.times are in steps of 2, xdf_counter.time_series in steps of 1 (thus /2) the first indices of EEG.lsl_times will be 0 due to no matching value in the xdf (thus those need to be substracted)
                disp('... lsl_times have already been added to EEG struct. Continuing ...')
            % otherwise write it
            else
                % check if the field is already present, otherwise initialize
                if ~isfield(EEG, 'lsl_times')
                    EEG.lsl_times = [];
                else disp('EEG.lsl_times array is incomplete. Overwriting existing values.')
                end
            end % if ~isfield(EEG, 'lsl_times')
            
            % now copy the time_stamps from xdf to EEG-struct
            EEG.lsl_times = xdfEEG.time_stamps;
        
        % if there is a counter and a bdf file
        elseif ~isempty(bdfFileName) & any(contains([lslnames{:}], {'Counter','Sync'}, 'IgnoreCase', true))
            
            % load counter data and append lsl_times to EEG struct
            if exist('xdfCounter','var')
                disp('... xdfCounter variable already exists. Continuing ...')
            else
                % look for the counter stream and copy information to xdfCounter
                for x = 1:length(streams)
                    if isfield(streams{x}.info,'name') && contains([lslnames{:}], {'Counter','Sync'}, 'IgnoreCase', true)
                        disp(['... Found eeg in streams'])
                        xdfCounter = streams{x};
                        break;
                    end
                end % fzor x = 1:length(streams)
            end % if exist('xdfCounter','var')

            % 
            if isfield(EEG, 'lsl_times') && length(xdfCounter.time_stamps)/(1000/xdfCounter.info.nominal_srate) == length(EEG.lsl_times) - sum(EEG.lsl_times == 0) % EEG.times are in steps of 2, xdf_counter.time_series in steps of 1 (thus /2) the first indices of EEG.lsl_times will be 0 due to no matching value in the xdf (thus those need to be substracted)
                disp('... lsl_times have already been added to EEG struct. Continuing ...')
            else
    
                if ~isfield(EEG, 'lsl_times')
                    EEG.lsl_times = [];
                else
                    disp('EEG.lsl_times array is incomplete. Overwriting existing values.')
                end % if ~isfield(EEG, 'lsl_times')
    
                % find indices in lsl_times to add to EEG struct
                [valuetmp,eegTimeIdx,lslTimeIdx] = intersect(EEG.data(find(strcmpi({EEG.chanlocs.labels},'cnt')),:),xdfCounter.time_series);
                EEG.lsl_times(1:eegTimeIdx(1)-1) = NaN;
                EEG.lsl_times(eegTimeIdx) = xdfCounter.time_stamps(lslTimeIdx);
                EEG.lsl_times(eegTimeIdx(end):length(EEG.data)) = NaN;

            end % if isfield(EEG,'lsl_times'...

            %% [Counter Branch] Merge events from xdfEEG.event using counter numbers
            if isfield(xdfEEG, 'event') && ~isempty(xdfEEG.event)
                disp('... copying events from xdfEEG (using counter mapping) to EEG event structure ...');
                
                % Ensure EEG.event exists
                if ~isfield(EEG, 'event') || isempty(EEG.event)
                    EEG.event = [];
                end
                
                % Loop over each event in xdfEEG.event
                for iEvent = 1:length(xdfEEG.event)
                    % Get the event’s latency index (into xdfEEG.data)
                    eventLatency = xdfEEG.event(iEvent).latency;
                    
                    % Check that the latency index is within the available columns
                    if round(eventLatency) > size(xdfEEG.data,2)
                        warning('Event latency index %d exceeds available data columns in xdfEEG.data.', round(eventLatency));
                        continue;
                    end
                    
                    % Retrieve the counter value from the first channel of xdfEEG.data.
                    % These counter values match xdfCounter.time_series.
                    counterValue = xdfEEG.data(1, round(eventLatency));
                    
                    % Find the corresponding index in xdfCounter.time_series where the counter value matches.
                    counterIdx = find(xdfCounter.time_series == counterValue, 1, 'first');
                    if isempty(counterIdx)
                        warning('Counter value %.3f not found in xdfCounter.time_series for event %d.', counterValue, iEvent);
                        continue;
                    end
                    
                    % Get the event’s timestamp (in seconds) from xdfCounter.time_stamps
                    evtTime = xdfCounter.time_stamps(counterIdx);
                    
                    % Map the event timestamp to a sample index in EEG.lsl_times using interpolation.
                    eventSample = find(EEG.lsl_times == evtTime);
                    
                    % Build the new EEG event structure.
                    newEvent.latency = eventSample;
                    newEvent.type    = xdfEEG.event(iEvent).type;  % using the counter event type from xdfEEG
                    
                    % Optionally, include event duration if available.
                    if isfield(xdfEEG.event(iEvent), 'duration')
                        newEvent.duration = xdfEEG.event(iEvent).duration;
                    end
                    
                    % Append the new event into EEG.event
                    EEG.event(end+1).type = newEvent.type;
                    EEG.event(end).latency = newEvent.latency;
                    EEG.event(end).duration = newEvent.duration;
                end
                disp('... Successfully merged counter-mapped events from xdfEEG.event into EEG.event.');
            else
                disp('... No events found in xdfEEG to merge using counter mapping.');
            end
        else disp('... no counter, no EEG stream, and no bdf file found. Continuing...')
            continue
        end % if isempty(bdfFileName)

        %% find indices of interesting streams
        disp('Detected streams:')
        for streami = 1:length(streams)
            if contains(streams{1,streami}.info.type,'audio')
                xdfAudio = streams{1,streami};
                disp('... cellphone audio stream found ...')
            elseif contains(streams{1,streami}.info.type,'Gaze')
                xdfGaze = streams{1,streami};
                disp('... neon gaze stream found ...')
            elseif contains(streams{1,streami}.info.type,'Event') & contains(streams{1,streami}.info.name,'Neon')
                xdfGazeEvent = streams{1,streami};
                disp('... neon event stream found ...')
            elseif contains(streams{1,streami}.info.type,'ECG') & contains(streams{1,streami}.info.name,'faros')
                xdfEcg = streams{1,streami};
                disp('... faros ECG found ...')
            elseif contains(streams{1,streami}.info.type,'Marker') & contains(streams{1,streami}.info.name,'Audio classifier')
                xdfAudioClassifier = streams{1,streami};
                disp('... Audio Classifier found ...')
            end
        end
        
        %% Merge gaze files or fallback to XDF gaze if no aligned CSV
        csvGaze = [];
        csvFix = [];
        pupil_events_local = [];
        loadXdfGaze = false;
        
        if ~isempty(dir(fullfile(subPath, '*_pupil')))
            pupil_folder = dir(fullfile(subPath, '*_pupil'));
            pupil_folder = pupil_folder.name;
            pupil_subfolder = dir(fullfile(subPath, pupil_folder));
            pupil_subfolder = pupil_subfolder(3).name;
        
            disp(['... ' subject ' pupil folder found. Processing.'])
            pupil_path = fullfile(subPath, pupil_folder, pupil_subfolder);
        
            % Check for time-aligned gaze file
            if isfile(fullfile(pupil_path, 'time_aligned_gaze.csv'))
                disp(['... time_aligned_gaze.csv found for ' subject])
                csvGaze = readtable(fullfile(pupil_path, 'time_aligned_gaze.csv'));
            else
                disp(['... No time_aligned_gaze.csv found for ' subject '. Using XDF gaze instead.'])
                loadXdfGaze = true;  % Only use xdf if no aligned gaze CSV
            end
        
            % Optionally load fixations if available
            if isfile(fullfile(pupil_path, 'fixations.csv'))
                disp(['... fixations.csv found for ' subject])
                csvFix = readtable(fullfile(pupil_path, 'fixations.csv'));
            else
                disp(['... No fixations.csv found for ' subject])
            end
        
            % Optionally load events if available
            if isfile(fullfile(pupil_path, 'events.csv'))
                disp(['... Found events.csv for ' subject])
                pupil_events_local = readtable(fullfile(pupil_path, 'events.csv'));
                pupil_events_local.lsl_times = str2double(cell(height(pupil_events_local), 1));
            else
                disp(['... No events.csv file for ' subject])
            end
        else
            % No pupil folder at all → fallback to xdf
            loadXdfGaze = true;
        end
        
        %% Add gaze from xdfGaze if present and fallback was triggered
        if loadXdfGaze && exist("xdfGaze", "var")
            % 1. Remove NaNs from EEG timestamps (to find the valid portion)
            validEEG = ~isnan(EEG.lsl_times);
            eegTimes = EEG.lsl_times(validEEG);
        
            % 2. Gaze data and timestamps
            gazeData = xdfGaze.time_series;
            gazeTimes = xdfGaze.time_stamps;
        
            % 3. Interpolate each gaze channel to EEG timestamps
            interpGaze = nan(size(gazeData,1), length(EEG.lsl_times));
        
            for ch = 1:size(gazeData,1)
                interpVals = interp1(gazeTimes, gazeData(ch,:), eegTimes, 'linear');
                fullInBounds = false(1, length(EEG.lsl_times));
                fullInBounds(validEEG) = eegTimes >= gazeTimes(1) & eegTimes <= gazeTimes(end);
                temp = nan(1, length(EEG.lsl_times));
                temp(validEEG) = interpVals;
                temp(~fullInBounds) = NaN;
                interpGaze(ch,:) = temp;
        
                % Add channel label and type if available
                if isfield(xdfGaze.info.desc, 'channels') && ...
                   isfield(xdfGaze.info.desc.channels, 'channel') && ...
                   length(xdfGaze.info.desc.channels.channel) >= ch
                    label = xdfGaze.info.desc.channels.channel{ch}.label;
                    type = xdfGaze.info.desc.channels.channel{ch}.type;
                else
                    label = ['Gaze_' num2str(ch)];
                    type = 'Gaze';
                end
        
                EEG.chanlocs(end+1).labels = label;
                EEG.chanlocs(end).type = ['Gaze_' type];
            end
        
            % 4. Append interpolated gaze data to EEG
            EEG.data(end+1:end+size(interpGaze,1),:) = interpGaze;
        
            disp('... Imported gaze channels from xdfGaze into EEG.data.')
        end

        %% add ECG from xdfEcg if present
        % check if faros file exists
        if exist("xdfEcg", "var") & ~exist('ecgData','var')
            % 1. Remove NaNs from EEG timestamps (to find the valid portion)
            validEEG = ~isnan(EEG.lsl_times);
            eegTimes = EEG.lsl_times(validEEG);
        
            % 2. ECG data and timestamps
            ecgTimes = xdfEcg.time_stamps;
            ecgData = double(xdfEcg.time_series);  % convert int16 to double for interpolation
        
            % 3. Interpolate ECG to EEG timestamps
            interpECG = nan(1, length(EEG.lsl_times));
        
            interpVals = interp1(ecgTimes, ecgData, eegTimes, 'linear');
            fullInBounds = false(1, length(EEG.lsl_times));
            fullInBounds(validEEG) = eegTimes >= ecgTimes(1) & eegTimes <= ecgTimes(end);
        
            temp = nan(1, length(EEG.lsl_times));
            temp(validEEG) = interpVals;
            temp(~fullInBounds) = NaN;
            interpECG(:) = temp;
        
            % 4. Add channel to EEG
            if isfield(xdfEcg.info.desc, 'channels') && ...
               isfield(xdfEcg.info.desc.channels, 'channel') && ...
               isfield(xdfEcg.info.desc.channels.channel{1}, 'label')
                label = xdfEcg.info.desc.channels.channel{1}.label;
            else
                label = 'ECG';
            end
        
            EEG.chanlocs(end+1).labels = label;
            EEG.chanlocs(end).type = 'ECG';
            EEG.data(end+1,:) = interpECG;
        
            disp('... imported xdf ECG into EEG.data. Continuing ...')
        elseif exist("xdfEcg", "var") & exist('ecgData','var')
            %% Replace LSL ECG with high-quality offline ECG from EDF
            % Assumes: ecgData (from edfread), ecgInfo (from edfinfo), xdfEcg (LSL ECG), EEG struct exists
            
            % --- Extract and format ECG signal ---
            newECG = ecgData.ECG;
            
            % Ensure it's numeric and a row vector
            if iscell(newECG)
                newECG = cell2mat(newECG);
            end
            if iscolumn(newECG)
                newECG = newECG';
            end
            
            fs_edf = ecgInfo.NumSamples(1);  % Sampling rate from EDF header
            
            % --- Get LSL ECG and interpolate to EDF rate ---
            lsl_ecg_ts = xdfEcg.time_stamps;
            lsl_ecg_data = double(xdfEcg.time_series);  % Convert to double
            
            interpTime = lsl_ecg_ts(1):1/fs_edf:lsl_ecg_ts(end);
            interpECG = interp1(lsl_ecg_ts, lsl_ecg_data, interpTime, 'pchip');
            
            % --- Cross-correlate for alignment ---
            tmpECG = newECG;
            cutoff = round(0.2 * length(newECG));  % Ignore edges to avoid noise
            tmpECG([1:cutoff, end-cutoff:end]) = 0;
            
            [r, lags] = xcorr(interpECG, tmpECG);
            [~, lagIdx] = max(r);
            offset = lags(lagIdx);
            
            % --- Apply offset ---
            if offset > 0
                newECG(1:offset) = [];
            elseif offset < 0
                newECG = [zeros(1, abs(offset)), newECG];
            end
            
            % --- Trim to match length of interpolated LSL ECG ---
            minLen = min(length(interpECG), length(newECG));
            newECG = newECG(1:minLen);
            newTime = (0:minLen-1) / fs_edf + lsl_ecg_ts(1);
            
            % --- Replace LSL ECG in xdfEcg with aligned EDF ECG ---
            xdfEcg.time_series = newECG;
            xdfEcg.time_stamps = newTime;
            
            disp('... Replaced LSL ECG with aligned EDF ECG in xdfEcg.');
            
            %% Append aligned ECG to EEG
            validEEG = ~isnan(EEG.lsl_times);
            eegTimes = EEG.lsl_times(validEEG);
            
            % Interpolate to EEG timebase
            interpECG = nan(1, length(EEG.lsl_times));
            interpVals = interp1(newTime, newECG, eegTimes, 'linear');
            fullInBounds = false(1, length(EEG.lsl_times));
            fullInBounds(validEEG) = eegTimes >= newTime(1) & eegTimes <= newTime(end);
            
            temp = nan(1, length(EEG.lsl_times));
            temp(validEEG) = interpVals;
            temp(~fullInBounds) = NaN;
            interpECG(:) = temp;
            
            % Append to EEG structure
            EEG.data(end+1,:) = interpECG;
            EEG.chanlocs(end+1).labels = 'ECG_EDF';
            EEG.chanlocs(end).type = 'ECG';
            
            disp('... appended aligned ECG to EEG.data as channel "ECG_EDF".');
        end

        %% add audio from xdfAudio if present
        if exist("xdfAudio", "var")
            % 1. Remove NaNs from EEG timestamps (to find the valid portion)
            validEEG = ~isnan(EEG.lsl_times);
            eegTimes = EEG.lsl_times(validEEG);
            
            % 2. Audio data and timestamps
            audioTimes = xdfAudio.time_stamps;
            audioData = xdfAudio.time_series;
        
            % 3. Interpolate each audio channel to EEG timestamps
            interpAudio = nan(size(audioData,1), length(EEG.lsl_times));
        
            for ch = 1:size(audioData,1)
                interpVals = interp1(audioTimes, audioData(ch,:), eegTimes, 'linear');
                fullInBounds = false(1, length(EEG.lsl_times));
                fullInBounds(validEEG) = eegTimes >= audioTimes(1) & eegTimes <= audioTimes(end);
                temp = nan(1, length(EEG.lsl_times));
                temp(validEEG) = interpVals;
                temp(~fullInBounds) = NaN;
                interpAudio(ch,:) = temp;
        
                % Try to get channel label/type info if available
                if isfield(xdfAudio.info.desc, 'channels') && ...
                   isfield(xdfAudio.info.desc.channels, 'channel') && ...
                   length(xdfAudio.info.desc.channels.channel) >= ch
                    label = xdfAudio.info.desc.channels.channel{ch}.label;
                    type = xdfAudio.info.desc.channels.channel{ch}.type;
                else
                    label = ['Audio_' num2str(ch)];
                    type = 'Audio';
                end
        
                EEG.chanlocs(end+1).labels = label;
                EEG.chanlocs(end).type = ['Audio_' type];
            end
        
            % 4. Append interpolated audio to EEG.data
            EEG.data(end+1:end+size(interpAudio,1),:) = interpAudio;
        
            disp('... imported all data from xdf audio to EEG.data. Continuing ...')
        end

        %% Merge Audio Classifier stream as EEG.events
        if exist('xdfAudioClassifier', 'var') && ~isempty(xdfAudioClassifier)
            disp('... merging Audio Classifier events into EEG.event ...');
        
            % Initialize EEG.event if it does not exist or is empty
            if ~isfield(EEG, 'event') || isempty(EEG.event)
                EEG.event = [];
            end
            
            % Check for required fields in the Audio Classifier stream
            if isfield(xdfAudioClassifier, 'time_stamps') && isfield(xdfAudioClassifier, 'time_series')
                % Loop over each event in the Audio Classifier stream
                for iEvent = 1:length(xdfAudioClassifier.time_series)
                    % Get the event type from the cell array (e.g., 'Speech', 'Silence', etc.)
                    evType = xdfAudioClassifier.time_series{iEvent};
                    
                    % Retrieve the event's corresponding time stamp (in seconds)
                    evtTime = xdfAudioClassifier.time_stamps(iEvent);
                    
                    % Map the event's time stamp to a sample index in EEG.lsl_times
                    eventSample = interp1(EEG.lsl_times, 1:length(EEG.lsl_times), evtTime, 'nearest', NaN);
                    if isnan(eventSample)
                        warning('Audio Classifier event %d with time %.3f s could not be aligned to EEG.lsl_times.', iEvent, evtTime);
                        continue;
                    end
                    
                    % Build the new event structure
                    newEvent.latency = eventSample;
                    newEvent.type = ['AudioClassifier_' evType];  % prepend a tag for clarity
        
                    % Append the event to EEG.event
                    EEG.event(end+1).latency = newEvent.latency;
                    EEG.event(end).duration = 1;
                    EEG.event(end).type = newEvent.type;
                end
                disp('... Successfully merged Audio Classifier events into EEG.event.');
            else
                disp('... Audio Classifier stream does not have the required time_stamps or time_series fields.');
            end
        end
        %% Finish the merger
        % check if everything is ok with the dataset
        EEG = eeg_checkset(EEG);
        % check if especially the events are ok
        EEG = eeg_checkset(EEG,'eventconsistency');
    
        % now save
        disp(['... saving subject ' subject ' to merged folder. Continuing ...'])
        EEG = pop_saveset(EEG,'filename',[exp '_' subject '_merged.set' ],'filepath',[pathVars.outLocation]);

    end % if any
end % for subject
