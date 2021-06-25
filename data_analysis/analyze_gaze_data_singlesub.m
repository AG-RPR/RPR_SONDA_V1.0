function [x_result, y_result] = analyze_gaze_data_singlesub(varargin)
% All arguments are "optional" (depending on how the filenames are saved), order of input
% is not relevant. The separate input arguments are useful for cycling
% through multiple conditions, to quickly analyze one subject the syntax
% load_data('subj', <full_file_name>); will work fine.
%
% INPUT ARGUMENTS:
% subj: string, subject name or code
% vel: integer, [1:4] velocity level
% vfd: integer, [0:3] simulated vfd: 0=control; 1=central defect; 2=peripheral defect; 3=hemifield defect
% contrast: integer, [0,1] 0=low contrast; 1=high contrast
% pursuit: integer, [0,1] 0=smooth pursuit; 1=saccadic pursuit
% ui: ignore the input arguments and load data by User Interface: in this case it extracts the file name and analyzes all the trials with for the selected subjectXconditions
[open_ui, filename] = load_data(varargin{:}); % get the filename depending on the string chosen by the user
if open_ui == 1 % if the user asked the UI mode, get the filename from the UI
    [filename, filepath] = uigetfile('matlab');
end
for trial_number = 1:6 % load the data from each trial
    if iscell(filename);  filename = filename{:}; end % if filename comes as a cell, make it a string
    if exist('filepath','var')
        load([filepath, filename(1:end-5), num2str(trial_number)]); % filename(1:end-5) is to remove file extension and trial number, which will be added procedurally
    else
        load([pwd,'/', filename, num2str(trial_number)]);
    end
    
    % time series trimming to remove possible artifacts from start/end of recording
    [x_stim, x_resp] = trim_timeseries(x_stim, x_resp);
    [y_stim, y_resp] = trim_timeseries(y_stim, y_resp);
    
    % convert time series from pixels to degrees of visual field
    x_stim = pix2deg(x_stim, xScreen, monitWidth, viewDist) - pix2deg(xScreen/2, xScreen, monitWidth, viewDist);
    x_resp = pix2deg(x_resp, xScreen, monitWidth, viewDist) - pix2deg(xScreen/2, xScreen, monitWidth, viewDist);
    y_stim = pix2deg(y_stim, xScreen, monitWidth, viewDist) - pix2deg(xScreen/2, xScreen, monitWidth, viewDist);
    y_resp = pix2deg(y_resp, xScreen, monitWidth, viewDist) - pix2deg(xScreen/2, xScreen, monitWidth, viewDist);
    
    % remove possible outliers (values outside screen boundaries) NB the threshold value depends on your setup!
    x_resp(abs(x_resp)>30) = NaN;
    y_resp(abs(y_resp)>30) = NaN;
    
    % remove blinks
    blink_thresh = 300; % threshold (deg/sec) to separate autentic saccades from blink artifacts
    [x_resp, ~, xnanflag] = blink_filter(x_resp,blink_thresh,0);
    [y_resp, ~, ynanflag] = blink_filter(y_resp,blink_thresh,0);
    [x_stim, ~, xnanflag] = blink_filter(x_stim,blink_thresh,0);
    [y_stim, ~, ynanflag] = blink_filter(y_stim,blink_thresh,0);
    
    % low-pass filtering
    lpf = designfilt('lowpassiir','FilterOrder',1,'HalfPowerFrequency',0.5,'DesignMethod','butter'); % low pass filter, adjust filter power depending on jitter in the data. For SONDA data 0.5 is recommended
    x_resp = filtfilt(lpf, x_resp);
    y_resp = filtfilt(lpf, y_resp);
    
    % calculate cosine similarity
    xCS(trial_number) = cosine_similarity(x_stim,x_resp);
    yCS(trial_number) = cosine_similarity(y_stim,y_resp);
    
    % convert time series from position to velocity: calculate the
    % derivative of the x and y positions over time, both for stimulus and
    % gaze, horizontal and vertical
    x_stim_vel = diff(x_stim);
    x_resp_vel = diff(x_resp);
    y_stim_vel = diff(y_stim);
    y_resp_vel = diff(y_resp);
    
    one_second = round(1/ifi); % how long is one second in frames. it will be used to trim the central part of the crosscorrelogram
    % "ccg" does cross correlation, normalization and trimming
    [x_axis, x_vel_cor] = ccg(x_resp_vel, x_stim_vel,one_second*2);
    [~, y_vel_cor] = ccg(y_resp_vel, y_stim_vel,one_second*2);
    x_axis = x_axis*ifi;
    
    % stack the cross-correlograms for the averaging later
    x_vel_tot(:,trial_number)=x_vel_cor;
    y_vel_tot(:,trial_number)=y_vel_cor;
      
    ecc_bin = -20:20; % array which represent eccentricities from -20 deg to + 20 deg. it's going to be used for the error distribution
    % calculate the difference (in position) between gaze and stimulus, both horizontal and vertical
    x_pos_err{trial_number} = x_resp - x_stim;
    y_pos_err{trial_number} = y_resp - y_stim;
end
% averaging cross-correlograms
x_vel_avg = nanmean(x_vel_tot,2);
y_vel_avg = nanmean(y_vel_tot,2);

% compute error distributions
x_pos_err_tot = [x_pos_err{:}];
hx=hist(x_pos_err_tot, ecc_bin);
y_pos_err_tot = [y_pos_err{:}];
hy=hist(y_pos_err_tot, ecc_bin);

% gaussian fits to error cross-correlograms and distributions
[x_gmodel, x_good]=fit(x_axis',x_vel_avg,'gauss1');
[y_gmodel, y_good]=fit(x_axis',y_vel_avg,'gauss1');
[hx_gauss, hx_good]=fit(ecc_bin',(hx/sum(hx))','gauss1');
[hy_gauss, hy_good]=fit(ecc_bin',(hy/sum(hy))','gauss1');

% spatio-temporal features
x_result.amp = x_gmodel.a1;
x_result.lag = x_gmodel.b1;
x_result.width = x_gmodel.c1;
x_result.adjR2 = x_good.adjrsquare;
x_result.erramp = hx_gauss.a1;
x_result.errmean = hx_gauss.b1;
x_result.errstd = hx_gauss.c1;
x_result.errR2 = hx_good.adjrsquare;
x_result.CoSi = nanmean(xCS);

y_result.amp = y_gmodel.a1;
y_result.lag = y_gmodel.b1;
y_result.width = y_gmodel.c1;
y_result.adjR2 = y_good.adjrsquare;
y_result.erramp = hy_gauss.a1;
y_result.errmean = hy_gauss.b1;
y_result.errstd = hy_gauss.c1;
y_result.errR2 = hy_good.adjrsquare;
y_result.CoSi = nanmean(yCS);

%% Plotting
clf;
subplot(121)
plot(x_axis(1:4:end), x_vel_avg(1:4:end),'k.', 'markersize',15);
hold on;
plot(x_axis, x_gmodel(x_axis), 'k', 'LineWidth',2);
ccgylim = [-0.05 0.2];
ccgxlim = [-0.2 1];
ylim(ccgylim); xlim(ccgxlim);  plot([0 0],[-1 1],'k--','LineWidth',0.5);
xlabel('time (s)'); ylabel('xcorr'); title('cross-correlogram');

subplot(122);
bar(ecc_bin,hx/sum(hx),'k','LineWidth',1); xlim([-20 20]);
hold on;
plot(ecc_bin, hx_gauss(ecc_bin),'k');
ylim([0 0.5]);
xlim([-20 20]);
alpha(0.5);
xlabel('positional error (deg)');
ylabel('probability (%)');
title('error distribution');
shg;
end