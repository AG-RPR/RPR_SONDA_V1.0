% Juxtapose epochs of 2 seconds each from the already exisiting trajectories to create the "saccadic pursuit" trajectories
% hz is the refresh rate of the screen and duration_frames is the duration of the tracking session in frames, both are automatically detected from acquire_gaze_gata.m
function make_saccadic_path(hz, duration_frames)
epoch_length = hz*2; % length of the epoch in frames
try
    for v = 1:4 % for each velocity level
        for n = 1:6 % for each trial
            loadname = ['v' num2str(v) '_path' num2str(n) '.mat']; % load one trajectory at a time
            load(loadname);
            xsmooth(n,1:duration_frames,v) = x_pos(1:duration_frames); % reshape smooth trajectories into a TRIAL by LENGTH by VELOCITY matrix
            ysmooth(n,1:duration_frames,v) = y_pos(1:duration_frames);
        end
    end
catch
    disp('Path files not found, check that you are in the correct folder'); % if cannot load the trajectories, print a warning message
    return;
end
xsacc=nan(6,size(xsmooth,2),4); % allocate a matrix with the same shape for saccadic trajectories
ysacc=nan(6,size(xsmooth,2),4);
for v = 1:4 % for each velocity
    for n = 1:6 % for each trial
        xsacc(n,1:hz*2,v) = xsmooth(n,1:hz*2,v); % the first epoch of each saccadic trajectory is the same as the smooth one
        ysacc(n,1:hz*2,v) = ysmooth(n,1:hz*2,v);
        t=2*hz+1; % start index of the epoch (it's length is fixed and given by epoch_length)
        while t+epoch_length <= size(xsmooth,2)+1 % check that the current epoch doesn't end beyond the end of the duration of the current trial
            tmp=Shuffle(setdiff(1:6,n)); % Shuffle is a PTB function (/Psychtoolbox/PsychProbability/Shuffle.m)
            rn = tmp(1); % random n
            if t+epoch_length == size(xsmooth,2)+1 % check if we are in the last epoch (needs a special indexing)
                xsacc(n,t:epoch_length,v) = xsmooth(rn,epoch_length+1-2*hz:epoch_length,v); % reminder for me in the future: this is one of your usual fucked-up indexing, just trust yourself this time, it works.
                ysacc(n,t:epoch_length,v) = ysmooth(rn,epoch_length+1-2*hz:epoch_length,v);
            else
                rt=randi([1 epoch_length-2*hz]); % if we are not in the last epoch, assign to the current epoch of saccadic trajectory a random epoch of smooth trajectory
                xsacc(n,t:t+2*hz,v) = xsmooth(rn,rt:rt+2*hz,v);
                ysacc(n,t:t+2*hz,v) = ysmooth(rn,t:t+2*hz,v);
            end
            t=t+epoch_length; % increment the index to the next epoch
        end
    end   
end
xsmooth = xsmooth(:,1:duration_frames,:); % trim the trajectories to the approriate length
ysmooth = ysmooth(:,1:duration_frames,:);
xsacc = xsacc(:,1:duration_frames,:);
ysacc = ysacc(:,1:duration_frames,:);
save('matrix_paths', 'xsmooth', 'ysmooth', 'xsacc', 'ysacc'); % save the trajectory matrix
end