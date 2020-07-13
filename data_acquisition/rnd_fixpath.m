% function rnd_fixpath(duration_frames, xScreen, yScreen)
%
% create a set of brownian motion trajectories with different velocity
% levels given the amount of frames (duration of the tracking session) and
% screen resolution
function [x_pos, y_pos]=rnd_fixpath(duration_frames, xScreen, yScreen)
iter = 0;
clc;
for i = 1:6
    disp(['Computing trajectory #' num2str(i)]);
    % 8 is the optimal value for 120hz
    vel = 8; % with 10 it takes forever but it still finishes (procedurally, it's more difficult to fit a a faster trajectory within the screen boundaries)
    while 1 % generate random trajectories until the criteria are satisfied (must not exceed screen boundaries)
        x_vel = randn(1,duration_frames); % allocate an array of random velocity values for x and y components, as long as the number of frames of the duration of the gaze acquisition session
        y_vel = randn(1,duration_frames);
        
        x_velgain = vel*4; % scale the velocity gain for horizontal and vertical components so that x = 2y. These values are arbitrary, but work pretty well on common 21-24" monitors
        y_velgain = vel*2;
        
        x_vel(1) = 0; % start from 0 velocity (after the integration to obtain the position it will help in setting the starting position of the stimulus arbitrarily)
        y_vel(1) = 0;
        x_vel=x_vel.*x_velgain; % apply the gain to the array of random velocity values
        y_vel=y_vel.*y_velgain;
        
        % low-pass filter to avoid excessive jitter
        cutoff = 10; % cutoff frequency (Hz)
        sigma = cutoff/sqrt(2*log(2)); % standard deviation of the gaussian filter
        sz = duration_frames;    % length of gaussian filter vector
        x = linspace(-sz / 2, sz / 2, sz);
        g_filter = exp(-x .^ 2 / (2 * sigma ^ 2)); % construction of the gaussian filter
        g_filter = g_filter / sum (g_filter); % normalization of the gaussian filter
        
        lp_x_vel = conv(x_vel, g_filter, 'same'); % apply the gaussian filter to low-pass horizontal and vertical velocity arrays
        lp_y_vel = conv(y_vel, g_filter,'same');
        
        x_pos = cumsum(lp_x_vel); % integration of velocity over time to obtain the array of positions (x and y coordinates of the stimulus)
        y_pos = cumsum(lp_y_vel);
        
        x_pos = x_pos + xScreen/2; % since it had 0 velocity, after the integration the starting position is [0, 0]. Adding half the screen resolution will set the starting position of the stimulus to the center of the screen
        y_pos = y_pos + yScreen/2;
        
        if max(x_pos) > xScreen || min(x_pos) < 0 || max(y_pos) > yScreen || min(y_pos) < 0 % check if the final position array violates the screen boundaries
            iter = iter + 1; % count iterations (for debugging purposes)
        else
            break; % if no boundaries are violated, stop generating trajectories and move on
        end
    end
    
    fsave = ['v0_path' num2str(i) '.mat']; % store the trajectories in a .mat file. v0 is the baseline velocity
    save(fsave, 'x_pos', 'y_pos');
    disp('Done!');
end

%%
disp('First set of trajectory created, starting to compute velocity levels');

for j = 2:5 % 4 velocity levels, determined as ratios of the baseline velocity (v0)
    disp(['Computing velocity #' num2str(j-1)]);
    for i = 1:6 % trial number
        fname = ['v0_path' num2str(i) '.mat']; % load the baseline velocity
        load(fname);
        
        % linear scaling different of different velocity levels, by resampling the baseline with different
        % ratios and trimming the length afterward
        x_pos = resample(x_pos,4,j); % j is the resampling numerator, with a denominator = 4. So the resulting velocities levels are v1 = 2/4, v2 = 3/4 , v3 = 4/4 and v4 = 5/4 of the baseline velocity (v3 effectively is equal to v0)
        y_pos = resample(y_pos,4,j); 
        x_pos = x_pos(20:end); % remove the first 20 samples, as the resampling often introduces artifacts at the very beginning and end of the time series
        y_pos = y_pos(20:end);
        
        % remove the jitter introduced by the resampling using a moving
        % average with a temporal window of 3 frames, both for x and y
        % coordinates
        tmp = x_pos;
        tmp(2:2:end-1) = (tmp(1:2:end-2)+tmp(3:2:end))/2;
        x_pos = tmp;
        
        tmp = y_pos;
        tmp(2:2:end-1) = (tmp(1:2:end-2)+tmp(3:2:end))/2;
        y_pos = tmp;
               
        fsave = ['v' num2str(j-1) '_path' num2str(i) '.mat']; % save the trajectory with the appropriate name
        save(fsave,'x_pos','y_pos');
        
    end
    disp('Done!');
end
disp('Everything done!');
