% help:
% Function to convert voltage to force Fz (normal to tapping box) based on
% measurements performed with a Biosemi AIB and a ATI mini 40 - FT4421
% 
% provide file name of the .mat file containing the time series of the
% 6 components recorded from the sensor, data should be filtered for 50 Hz 
% and baseline corrected (5 to 10 s) for offests between channels, keep 
% only the 6 channels corresponding to the force components then, send data
% to workspace and save .mat file
% provide the sampling rate (sr)
% varargin = filname.mat,sr 
%
% save data using import matrix from matlab : Fz
%
%
% Julien Lambert & Cédric Lenoir 27/01/2023, RnB-Lab
%
function conversion_volt_to_force(participant,sr)

%% Preamble
% Load parameters file
params = motor_get_params(); 

% Import version 6 of Letswave
import_lw(6)

% Define paths
raw_dir = fullfile(params.path_data, ...
    sprintf('sub-1%03d', participant)); 

% Compute trial duration
trial_dur = params.tempi/1000 * ...
            params.n_events_per_pattern * ...
            params.n_patterns_per_trial;

% Indicate calibration coefficients
calib_coeff_fz = [-10.7662351676953 -0.14181979074703 ...
    -10.2649711848043 0.17631695409008 -10.4154237298717 ...
    0.4954806017287];


for session = (1:2:3) % for the pre- and post-training sessions (1 & 3)

    %% Load .bdf data
    file_name = sprintf(['ses-%03d/eeg/' ...
        'sub-1%03d_task-motor_session-%01d.bdf'], ...
        session, participant, session);
    [header, data] = RLW_import_BDF(fullfile(raw_dir, file_name)); 
    
    %% Extract force data
    % Indicate relevant channels
    chan_names = {'Ana1', 'Ana2', 'Ana3', 'Ana4', 'Ana5', 'Ana6'}; 

    % Extract only force data
    [header_force, data_force] = ...
        RLW_arrange_channels(header, data, chan_names); 

    % Segment data ot keep tapping trials only (with a 5-s baseline)
    [header_force, data_force] = RLW_segmentation(...
            header_force, data_force, {'2'}, ...
            'x_start', -params.dur_epoch_buffer, ...
            'x_duration', trial_dur + 2 * params.dur_epoch_buffer);

    % Resize data to have a 3-D matrix
    data_force = squeeze(data_force);

    % Create an empty matrix to store the data
    matrix_session = zeros(size(data_force, 3), ...
            size(data_force, 1) + 1);

    for trial = 1:5

        % Extract data for the current trial
        data_current_trial = squeeze(data_force(trial, :, :));

        % Create an empty matrix to store the data
        forces = zeros(size(data_current_trial, 1), ...
            size(data_current_trial, 2));

        for component = 1:6

            %% Perform baseline substraction
            % Compute mean value on the first 5 s
            mean_baseline = mean(data_current_trial(component, 1:(sr*5)));

            for data_point = 1:size(data_current_trial, 2)

                % Extract baseline value from each data point
                data_current_trial(component, data_point) = ...
                    data_current_trial(component, data_point) - ...
                    mean_baseline;

            end

            %% Compute force
            % Weight the component
            forces(component, :) = calib_coeff_fz(component) .* ...
                data_current_trial(component, :);


        end

        % Recompose the signal
        fz = sum(forces);

        % From µV to V
        fz = -fz/1000000;

        %% Plot
        % Build time vector
        x = linspace(-5, size(data_current_trial, 2)/sr - 5, ...
            size(data_current_trial, 2));

        % Plot the data
        figure('Position', [0,0,2000,500])
        plot(x, fz)
        ax = get(gca);
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'Force (N)';
        hold on
        xline(0, '--r') % beginning of trial
        xline(trial_dur, '--r') % end of trial
        hold off
        title('Normal force during tapping')

        %% Store data
        matrix_session(: , trial + 1) = fz;

    end

    %% Save data
    % Add time
    matrix_session(: , 1) = x;

    % Add column names
    col_labels = {'time', 'trial_1', 'trial_2', 'trial_3', 'trial_4', ...
        'trial_5'};
    final_table = array2table(matrix_session, 'VariableNames', col_labels);

    % Write table
    export_name = sprintf('sub-1%03d_ses-%03d_force.csv', ...
        participant, session);
    export_path = fullfile(params.path_output, ...
        'data/segmented/force', ...
        sprintf('sub-1%03d', participant));

    % Create participant folder if needed
    if ~isfolder(export_path)
        mkdir(export_path)
    end

    % Write table
    writetable(final_table, fullfile(export_path, export_name));

end

