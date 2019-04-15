function data_with_features = chirp_features(data_samples)
% Copyright (c) 2015, Eero Väyrynen, version 1.01
% CHIRP_FEATURES.m function implements method presented in (PAPER HERE)
% The function uses an optimization method (Particle Swarm Optimization Toolbox)
% downloadable from http://www.mathworks.com/matlabcentral/fileexchange/7506-particle-swarm-optimization-toolbox
% A chirp model (chirpmodel.m) is fitted to the averaged data using an objective function (objFun.m).
% Optional preprocessing using interquartile mean for outlier removal and a
% zero-phase filtering to remove low frequency drifts. Optimization of
% parameters is performed in a set variable range varRange) and goodness-of-fit 
% measures are computed. Amplitude/location based basic SEP features are produced 
% using the fitted chirp model.
%
% INPUTS
%
% data_samples(i)                  	A structure containing raw SEP data and recording parameters for each i sample
% Fields:
%
% .Fs                               Sampling frequency
% .raw_data                         Matrix containing raw SEP data where each N sample is stored in the rows of the matrix
% .amplitude_unit                   Scale factor for signal amplitudes (conversion to volts)
% 
% OUTPUTS
%
% data_with_features(i)             A structure containing all the input structure fields (above) and fields produced in the fitting 
% Fields:
%
% .filtered                         A flag denoting if prefiltering was used
% .IQM                              A flag denoting if interquartile mean was used to eliminate outlier data
% .peakdirection                    Stores the original first peak direction for peak flipping 
% .F_model                          Chirp model used in fitting (lin, poly3, poly4, exp, linexp)
% .mean_signal                      Mean signal waveform calculated from the raw data (filtered if filtering is enabled)
% .fitted_chirp                     Fitted chirp signal waveform
%
% .fit_delay                        Chirp delay (tau)
% .fit_linweight_start              Linear weight start (a)
% .fit_linweight_stop               Linear weight stop (b)
% .fit_nonlinear_weight_exp         Nonlinear weight (lambda_w)
% .fit_start_F                      Start frequency (f0)
% .fit_end_F                        End frequency (f1)
% .fit_F_linearity                  Nonlinear frequency parameter (lambda_f)
% .fit_Gbest                        Best objective function value
%
% .fit_NSS_residual                 Normalised sum of squares residual of the fitting
%
% .fit_first_peak_amplitude         First SEP peak amplitude
% .fit_second_peak_amplitude        Second SEP peak amplitude
% .fit_first_peak_location          First SEP peak delay
% .fit_second_peak_location         Second SEP peak delay
% .fit_peak2peak_amplitude          SEP first and second peak Peak-to-peak amplitude 
% .fit_peak2peak_delay              SEP first and second peak Peak-to-peak delay
%
% .fit_RMSE                         Root Mean Squared Error of the fit

% Version 1.01 changes
% - New BP filter applied [10 280] Hz
% - New PSO parameter defaults for speed improvements (see code, pso_Trelea_vectorized call)

global fs
global data_points
global origData
global chirp_length
global model

model = 'lin';       % set chirp model (lin, poly3, poly4, exp, linexp)
prefiltering = 0;       % use filtering
filteroutliers = 0;     % use InterQuartile Mean to filter outliers

fs = data_samples(1).Fs;
data_points = size(data_samples(1).raw_data,2);
min_delay = 0.006;
max_delay = 0.02;
chirp_length = floor(data_points-fs*max_delay); % set chirp_length to data length minus max delay
data_with_features = data_samples; % include original data and parameters
[data_with_features.filtered] = deal(0); % set default filtered flag
[data_with_features.IQM] = deal(0); % set default outlier filtering (interquartile mean) flag

%% preprocessing FIR filter
if prefiltering,
    Fc = [10 280]; % IIR BP filter cutoff at 10Hz and 280Hz
    N = 4;
    d = designfilt('bandpassiir','DesignMethod','butter', 'FilterOrder',N,'HalfPowerFrequency1',Fc(1),'HalfPowerFrequency2',Fc(2),'SampleRate',fs);
    [data_with_features.filtered] = deal(1); % set filtered flag
   
end

%% Main loop
for i = 1:1:length(data_samples)
    
    % Calcualate and select mean signal for chirp optimisation
    if filteroutliers,% filter outliers using trimmean
        origData = trimmean(data_with_features(i).raw_data,50,1);
        [data_with_features.IQM] = deal(1); % set filtered flag
    else
        origData = mean(data_samples(i).raw_data,1); % mean from raw data        
    end
    
    %% preprocessing
    if prefiltering
        offset = mean(origData(1:20)); % Get bias offset before filtering to avoid step in mirroring
        filtData = filtfilt(d,[-1*origData(end:-1:1)+offset origData-offset]); % mirror, flip, tune offset, and zerophase filter
        origData = filtData((data_points+1):end)+offset; % cut, return offset and replace original
    end
    
    % Scale the signal to microvolts
    origData = origData.*data_with_features(i).amplitude_unit.*10^6;
    data_with_features(i).mean_signal = origData; % store the mean signal
    
    % flip the signal so that maximum value is positive (emphasize beginning)
    tempSignal = abs(origData.*linspace(1,0.01,data_points)); % linearly weighted abs signal
    [~, maxlocation] = max(tempSignal(floor(fs*min_delay)+1:end)); % get max absolute value after min delay
    maxlocation = maxlocation + floor(fs*min_delay); % adjust maxlocation to original data
    data_with_features(i).peakdirection = sign(origData(maxlocation)); % store sign of the peak amplitude
    origData = data_with_features(i).peakdirection*origData; % flip first peak to positive side for optimisation


    
    %% Optimisation
    % Tries to find parameters delay, a, b, lambda_w, f1, f2, lambda_f with for
    % maximal fitting between original data and optimised chirpfunction.
    
    varRange = [min_delay max_delay ;... % tau
        0.5*max(abs(origData)) 1.5*max(abs(origData));... % a (assumes signal first peak is positive to reduce computations and increase robustness)
        0 0.5*max(abs(origData));... % b
        1/5 1/0.002;...% attenuation nonlinearity lambda_w
        20 140;... % start frequency f0
        0.01 10 ;... % end target frequency f1
        0 1]; % frequency nonlinearity lambda_f (0 only exponential, 1 only linear)
    dim = size(varRange,1); % Dimension corresponds to number of optimised parameters
    minmax = 0; % 0 = minimized, 1 = objFcn function is maximized
    maxvel = 4; % maximum velocity of particles
    
    % PSO
    [optOUT,~,~] = pso_Trelea_vectorized('objFcn',dim,maxvel,varRange,minmax,[0 2000 50 2 2 0.9 0.4 1500 1e-25 250 NaN 1 0]);
    
    %% store optimized parameter outputs
    data_with_features(i).F_model = model; % chirp model
    
    data_with_features(i).fit_delay = optOUT(1); % chirp delay (tau)
    data_with_features(i).fit_linweight_start = optOUT(2); % linear weight start (a)
    data_with_features(i).fit_linweight_stop = optOUT(3); % linear weight stop (b)
    data_with_features(i).fit_nonlinear_weight_exp = optOUT(4); % nonlinear weight (lambda_w)
    data_with_features(i).fit_start_F = optOUT(5); % start frequency (f0)
    data_with_features(i).fit_end_F = optOUT(6); % end frequency (f1)
    if strcmp(model,'lin'),
        data_with_features(i).fit_F_linearity = 1; % nonlinear frequency parameter (lambda_f) linear case       
    elseif strcmp(model,'exp'),
        data_with_features(i).fit_F_linearity = 0; % nonlinear frequency parameter (lambda_f) exp case       
    else
        data_with_features(i).fit_F_linearity = optOUT(7); % nonlinear frequency parameter (lambda_f)
    end
    data_with_features(i).fit_Gbest = optOUT(8); % best objective function value
    
    
    %% Generate the optimised chirp function
    
    t = linspace(0,data_points/fs,data_points)'; % linear time
    x = chirpmodel(t, data_with_features(i).fit_delay, data_with_features(i).fit_linweight_start, data_with_features(i).fit_linweight_stop, data_with_features(i).fit_nonlinear_weight_exp, data_with_features(i).fit_start_F, data_with_features(i).fit_end_F, data_with_features(i).fit_F_linearity, data_with_features(i).F_model);
    data_with_features(i).fitted_chirp = data_with_features(i).peakdirection*x;

    
    %% Calculate normalized sum of squares residual of the fitted chirp
    
    data_with_features(i).fit_NSS_residual = mean((data_with_features(i).mean_signal-data_with_features(i).fitted_chirp').^2)/(mean(data_with_features(i).fitted_chirp.^2));
    
    %% calculate model peaks and locations(delays), peak2peak amplitude/delay values
    
    [data_with_features(i).fit_first_peak_amplitude, data_with_features(i).fit_first_peak_location] = max(data_with_features(i).peakdirection*data_with_features(i).fitted_chirp);
    
    [data_with_features(i).fit_second_peak_amplitude, data_with_features(i).fit_second_peak_location] = min(data_with_features(i).peakdirection*data_with_features(i).fitted_chirp);
    data_with_features(i).fit_second_peak_amplitude = abs(data_with_features(i).fit_second_peak_amplitude);
    
    data_with_features(i).fit_first_peak_location = data_with_features(i).fit_first_peak_location/fs;
    
    data_with_features(i).fit_second_peak_location = data_with_features(i).fit_second_peak_location/fs;
    
    data_with_features(i).fit_peak2peak_amplitude = data_with_features(i).fit_first_peak_amplitude + data_with_features(i).fit_second_peak_amplitude;
    data_with_features(i).fit_peak2peak_delay = data_with_features(i).fit_second_peak_location - data_with_features(i).fit_first_peak_location;

    %% fit statistics
    data_with_features(i).fit_RMSE = sqrt(mean((data_with_features(i).mean_signal-data_with_features(i).fitted_chirp').^2));
    
    %% debug plots and outputs
%     figure
%     plot(t, data_with_features(i).mean_signal)
%     hold on
%     plot(t, data_with_features(i).fitted_chirp)
%     data_with_features(i)

end

end