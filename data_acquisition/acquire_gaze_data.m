% function acquire_gaze_data(varargin)
%
% 'acquire_gaze_data' is the main function to collect gaze data from an EyeLink eye-tracker. The data will be already formatted to be used with the
% SONDA analysis pipeline.

function acquire_gaze_data(varargin)
%% cleaning workspace & default variables
close all; % close all windows (not screens!)
fclose('all'); % close all files
clc; % clear command window
commandwindow;

%% default variables
myPath = [pwd '/'];
cd(myPath); % go to path location
dummymode = 1; % flag to switch on / off dummymode (trial without eye-tracker)
subj ='demo1'; % string variable for subject name
session_duration = 20; % duration of each tracking session
nTrial = 6; % number of trials per each condition
cont_cond = [0]; %contrast
pursuit_cond = [0 1]; % 0 = smooth pursuit, 1 = saccadic pursuit
vel_cond = [3]; % velocity level (1 - 4). These levels are pre-computed in rnd_fixpath.m
vfd_cond = [0]; % gaze contingent simulated visual field defect condition (0 = no defect, 1 = central defect, 2 = peripheral defect, 3 = hemifield defect)
contrast = [0.10 0.50]; % contrast levels
monitWidth = 540; % physical width of the display monitor (mm)
viewDist = 600; % viewing distance of the observer (mm)
skip_trial = 0; % flag to allow the user to skip trial
go_out = 0; % flag to allow the user to quit a session any time

%% parse input arguments
for nk=1:2:nargin % for each pair of input argument
    switch varargin{nk}
        case {'dummymode', 'dummy' ,'dm'}
            dummymode = varargin{nk+1};
        case {'subj','Subj','subject','Subject','name'}
            subj=varargin{nk+1};
        case {'ntrial','ntrials'}
            nTrial=varargin{nk+1};
        case {'duration'}
            session_duration=varargin{nk+1};
            if session_duration > 20
                disp('The duration of each trial cannot exceed 20 seconds');
                return
            end
        case {'pursuit','Pursuit','pat'}
            pursuit_cond = varargin{nk+1};
        case {'contrast','Contrast', 'lum'}
            cont_cond =varargin{nk+1};
        case {'vel','Vel','velocity','Velocity'}
            vel_cond=varargin{nk+1};
        case{'vfd', 'cond'}
            vfd_cond = varargin{nk+1};
        otherwise
            error(sprintf('invalid parameter %s\n',varargin{nk}));
    end
end

%% Create list of conditions for the data acquisition (see https://doi.org/10.1101/2020.04.20.20072603 for details)
list = rnd_list(subj, myPath, pursuit_cond,  cont_cond, vel_cond, vfd_cond);

%% Initialze PTB (PsychToolBox) screen & miscellaneous graphic options
Screen('Preference', 'SkipSyncTests', 1); % skip PTB sync-test
Screen('Preference', 'SkipSyncTests', 2);
Screen('Preference', 'VisualDebuglevel', 3); % skip PTB warning
screenNumber = max(Screen('Screens'));  % get the pointer to the screen (select secondary screen in case of multiple displays)
white = WhiteIndex(screenNumber); % define white (usually 255)
grey=GrayIndex(screenNumber); % define gray (usually 127)
[w, windowRect] = Screen('OpenWindow', screenNumber, grey); % open the window object and get a pointer to access that window, initialize the window with a grey background
Screen('Flip', w); % initial flip to clean the screen
ifi = Screen('GetFlipInterval', w); % get interframe interval
hz = Screen('NominalFrameRate', w,1); % get refresh rate
MaxPriority(w); % set top priority
[xCenter, yCenter] = RectCenter(windowRect); % get central coordinates (resolution-dependent)
[xScreen, yScreen]= Screen('WindowSize', w); % get screen coordinates (resolution-dependent)
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % enable alpha blending to smooth the edge of the stimulus and the VFD conditions
Screen('TextSize', w, 30); % set default text size to 30
AssertOpenGL; % ensure that OpenGL is running properly before running the data acquisition

%% Define some experimental constants and pre-allocate arrays for storing results
duration_frames = round(hz) * session_duration; % duration of each tracking session in frames
diameter_vfd = 15; % pre-defined diameter of the central VFD (in degrees of visual field)
diameter = deg2pix(diameter_vfd,xScreen,monitWidth,viewDist); % convert the size of the central VFD from degrees of visual to pixels
radius_vfd = diameter/2; % get the VFD radius
x_resp = nan(1,duration_frames); % stores the x coordinates of the tracked eye
y_resp = nan(1,duration_frames); % stores the y coordinates of the tracked eye
x_stim = nan(1,duration_frames); % stores the x coordinates of the moving stimulus
y_stim = nan(1,duration_frames); % stores the y coordinates of the moving stimulus

%% Experimental variables & misc stuff
try
    load([myPath 'stim_trajectories/matrix_paths.mat']); % try to access the matrix containing the pre-generated trajectories of the stimulus, if it fails, it creates a new one
catch
    for ddd=1:10; disp('.'); end
    create_new_paths = input('Data containing the trajectory not found. Do you want to create it? [1: YES | 0: NO] '); % ask confirmation to create the trajectories
    if create_new_paths == 1 % if YES, create the trajectories in the subfolder paths
        mkdir(myPath,'stim_trajectories'); % create new folder "stim_trajectories"
        cd([myPath '/stim_trajectories']); % move to the new folder
        rnd_fixpath(duration_frames*2,xScreen,yScreen); % create trajectories (smooth pursuit) based on system's specs (refresh rate and resolution)
        make_saccadic_path(hz); % "scramble" the smooth pursuit trajectories to obtain the saccadic ones (see Grillini et al. 2020, Supplementary Materials for more details)
    else
        disp('Check that you are in the correct folder.'); % if NOT, print an error message and quit
        disp('End of testing.');
        sca;
        return;
    end
end

load([myPath 'stim_trajectories/matrix_paths.mat']); % pre-load the trajectories of the stimulus
% the matrices have the following structure:
% TRIAL X FRAME X VELOCITY
% EXAMPLES:
% xsmooth(4,:,3); horizontal component, smooth pursuit, 4th trial, 3rd velocity level
% ysacc(2,:,4); vertical component, saccadic pursuit, 2nd trial, 4th  velocity level

%% Eyelink initialization & calibration

el=EyelinkInitDefaults(w); % create a structure with graphic environment and tracker infos
EyelinkInit(dummymode,1) % initialization, if no eye-tracker is connect start in dummy mode
[~, vs]=Eyelink('GetTrackerVersion'); % get eye tracker version
switch dummymode
    case 0
        fprintf('Running the test on a ''%s'' tracker.\n', vs ); % if not in dummy mode, print the actual version of the eye-tracker (sometimes it is useful to debug issues)
    case 1
        disp('Running the test in dummy mode [no eyetracker connected]'); % if in dummy mode, warn the user
end
Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA'); % make sure that we get gaze data from the Eyelink (this is some kind of data formatting)

EyelinkDoTrackerSetup(el); % calibrate the eye-tracker (standard 9-point calibration)
EyelinkDoDriftCorrection(el); % do a final calibration check using drift correction
Screen('FillRect',w,grey); % clean the screen and fill it with a grey background
Screen('Flip', w);
eye_used = -1; % -1 because we don't know yet which eye is available. The EyeLink will change this automatically once it picks-up either eye

%% Data acquisition session
for ii = 1:size(list,1) % for each item on the list of conditions (the one created with rnd_list.m)
    % Format strings to be used in the filename, depending on the current conditions
    if go_out; ListenChar; break; end
    pat = list(ii,1); % pattern = type of trajectory (smooth = 0; saccadic = 1)
    lum = list(ii,2); % luminance = contrast level (low = 0; high = 1)
    vel = list(ii,3); % velocity level (1 to 4, see rnd_fixpath.m)
    vfd = list(ii,4); % type of visual field defect simulation (1 to 4)
    disp(['Trial #' num2str(ii)]);
    switch pat
        case 0
            disp('Pursuit modality: SMOOTH');
            pat_string = 'sm'; % smooth pursuit
        case 1
            disp('Pursuit modality: SACCADIC');
            pat_string = 'sc'; % saccadic pursuit
    end
    switch lum
        case 0
            disp('Contrast: LOW');
            lum_string = 'l'; % low luminance
        case 1
            disp('Contrast: HIGH');
            lum_string = 'h'; % high luminance
    end
    disp(['Velocity level: ' num2str(vel)]);
    switch vfd
        case 0
            disp('Simulated Visual Field Defect: NONE');
            vfd_string = 'c'; % control, no defect
        case 1
            disp('Simulated Visual Field Defect: CENTRAL SCOTOMA');
            vfd_string = 'm'; % pseudo-macular degeneration
        case 2
            disp('Simulated Visual Field Defect: PERIPHERAL SCOTOMA');
            vfd_string = 'g'; % pseudo-glaucoma
        case 3
            disp('Simulated Visual Field Defect: HEMIANOPIA');
            vfd_string = 'e'; % pseudo-hemianopia
    end
    switch vel
        case 1
            vel_string = 'v1';
        case 2
            vel_string = 'v2';
        case 3
            vel_string = 'v3'; % v3 is usually the baseline (see rnd_fixpath.m)
        case 4
            vel_string = 'v4';
    end
    
    %% Stimulus initialization (Gaussian luminance blob)
    size_blob_deg = 0.86; % diameter of the stimulus in degrees. The resulting size of the stimulus correspond to a Goldmann III stimulus = 0.43 deg radius
    size_blob=round(deg2pix(size_blob_deg,xScreen,monitWidth, viewDist)); % size of the stimulus in pixels
    [xg,yg]=meshgrid(-size_blob/2:size_blob/2-1,-size_blob/2:size_blob/2-1); % create the meshgrid of coordinates to be used to create the PTB texture of the stimulus
    gauss_sigma = size(xg)/4; gauss_sigma=gauss_sigma(1); % standard deviation of the gaussian aperture to mask the stimulus and giving it smooth edges
    gaussian = exp(-((xg/gauss_sigma).^2)-((yg/gauss_sigma).^2)); % gaussian aperture
    gauss_blob(:, :, 2)=  gaussian * white; % alpha layer of the texture containing the gaussian aperture. 255 = transparent, 0 = opaque
    gauss_blob(:, :, 1) = ones(length(gauss_blob(:,1,2)),length(gauss_blob(:,1,2)))* white; % background layer of the texture. 255 = white, 0 = black (for colored stimuli, create 4 layers, one for alpha and 3 for RGB)
    if dummymode % if in dummy mode load a little icon of an eye as a place-holder for the mouse cursor
        immat = imread('/Users/Ale/Google_Drive/REPERIO/EMC/Source_code/eye-icon.png');
        immat(immat == 0) = grey;
    end
    x_mouse = xCenter; % position mouse coordinates at the center of the screen (useful only in dummy mode)
    y_mouse = yCenter;
    if dummymode
        texim = Screen('MakeTexture', w, immat); % if in dummy mode, create a PTB texture out of the eye icon
    end
    rectim = [0 0 75 75]; % scale the eye icon to 75x75 pix
    gauss_blob(:, :, 1) = ones(length(gauss_blob(:,1,2)),length(gauss_blob(:,1,2)))* round(white*contrast(lum+1)) + grey; % adjust stimulus contrast depending on the current condition. "contrast" is a 1x2 vector
    tex_blob=Screen('MakeTexture', w, gauss_blob); % create a PTB texture for the stimulus
    rect_blob = [0 0 size_blob size_blob]; % scale the stimulus to the appropriate size in pixel
    %% Simulating visual field defects
    [xs,ys]=meshgrid(-xScreen:xScreen, -xScreen:xScreen); % meshgrid coordinates to create the VFD texture to be dipsplayed on screen
    transLayer=2; % transparency layer
    maskblob=ones(2*xScreen+1, 2*xScreen+1, transLayer) * grey;
    d=sqrt(xs.^2+ys.^2); % circular area based on the mashgrid coordinates
    pp = 4; % exponent of the gaussian mask. The higher, the sharper the edge
    dm = 0; % center of the gaussian mask. Default at 0
    dsd=radius_vfd; % size. When equal to radius_vfd it looks okay
    switch vfd
        case 0
        case 1 % central loss
            maskblob(:,:,transLayer)=round(exp(-((d-dm)/dsd).^pp) * white);  % create a transparency mask with a gaussian occlusion at the center
            texScot=Screen('MakeTexture', w, maskblob); % create a PTB texture
            dstRectScot=Screen('Rect', texScot);
        case 2  % peripheral loss
            maskblob(:,:,transLayer)=round(-exp(-((d-dm)/dsd).^pp) * white)+white; % create a transparency mask with a gaussian opening at the center
            texScot=Screen('MakeTexture', w, maskblob); % create a PTB texture
            dstRectScot=Screen('Rect', texScot);
        case 3 % hemifield loss
            maskHemi = ones(yScreen, xScreen*2)*grey;
            maskHemi(:,:,2)=ones(yScreen, xScreen*2)*white;
            for c=1:yScreen % this little loop creates a vertical smoothed edge
                a(c,:)=5:5:255;
            end
            maskHemi(:,1:51,2)=a;
            texHemi=Screen('MakeTexture', w, maskHemi); % create a PTB texture
            rectHemi = [0 0 xScreen*2 yScreen]; % resize the hemifield defect to fit the screen, doubled x dimension to be ensure that the defect is always present even if the gaze falls outside the screen
    end
    %% session start (6 trials) (STOPPED HERE) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    for k = 1:nTrial
        switch pat
            case 0 % smooth
                x_pos = xsmooth(k,:,vel);
                y_pos = ysmooth(k,:,vel);
            case 1 % saccadic
                x_pos = xsacc(k,:,vel);
                y_pos = ysacc(k,:,vel);
        end
        x_pos = x_pos(1:duration_frames); % shouldn't be necessary but does no harm. It ensures that the arrays to store the results have proper length
        y_pos = y_pos(1:duration_frames);
        cd(myPath);
        % if the edfFile name gets longer than 7 digits it will give some trouble with
        % the old Eyelink (it will be saved but not retrieved)....Doesn't
        % matter, since we keep the edf files only as a backup
        edfFile=[subj vel_string vfd_string lum_string pat_string num2str(k) '.edf']; %open file with subject name to record data to
        Eyelink('Openfile', edfFile);
        if go_out; ListenChar; break; end  % in case the user at any time requested to quit the trial
        %% Option screen with user controls
        % creates a minimal interface to show the user what condition is
        % coming next, and to allow trial skipping or quitting
        % (with/without saving data)
        while 1
            Screen('FillRect',w,grey);
            ygap = 50;
            low_border = 50;
            yspacing = low_border:ygap:low_border+ygap*5;
            DrawFormattedText(w, '[SPACEBAR]: START', 'center', 'center', white);
            DrawFormattedText(w, '[S]: Skip Next Trial', 20, yspacing(3) , white);
            DrawFormattedText(w, '[C]: Calibration', 20, yspacing(2) , white);
            DrawFormattedText(w, '[ESC]: Exit', 20, yspacing(1), white);
            Screen('TextSize', w, 20);
            ygap = 30;
            yspacing = low_border:ygap:low_border+ygap*5;
            DrawFormattedText(w, ['Trial #' num2str(k)], 20, yScreen-yspacing(6), white);
            switch pat
                case 0
                    DrawFormattedText(w,'Pursuit modality: SMOOTH', 20, yScreen-yspacing(5), white);
                case 1
                    DrawFormattedText(w,'Pursuit modality: SACCADIC', 20, yScreen-yspacing(5), white);
            end
            switch lum
                case 0
                    DrawFormattedText(w,'Contrast: LOW', 20, yScreen-yspacing(4), white);
                case 1
                    DrawFormattedText(w,'Contrast: HIGH', 20, yScreen-yspacing(4), white);
            end
            DrawFormattedText(w, ['Velocity Level: ' num2str(vel)], 20, yScreen-yspacing(3), white);
            switch vfd
                case 0
                    DrawFormattedText(w,'Simulated Visual Field Defect: NONE (CONTROL)', 20, yScreen-yspacing(2), white);
                case 1
                    DrawFormattedText(w,'Simulated Visual Field Defect: CENTRAL', 20, yScreen-yspacing(2), white);
                case 2
                    DrawFormattedText(w,'Simulated Visual Field Defect: PERIPHERAL', 20, yScreen-yspacing(2), white);
                case 3
                    DrawFormattedText(w,'Simulated Visual Field Defect: HEMIFIELD', 20, yScreen-yspacing(2), white);
            end
            Screen('TextSize', w, 30);
            Screen('Flip',w);
            WaitSecs(0.2);
            [~, keyCode, ~] = KbWait(-3); % wait for user input. -3 listens to all keyboard devices connected
            ansKey=find(keyCode); % NB: these keycodes are based on MacOS and change depending on the OS
            if ansKey == 41 % ESC (quit & save)
                ListenChar;
                ShowCursor;
                if ~dummymode
                    Eyelink('StopRecording');
                    Eyelink('CloseFile');
                    commandwindow;
                    saveExit = input('Save file? [1 = YES | 0 = NO] '); % ask the user if they want to save the data after quitting
                    if saveExit
                        try
                            fprintf('Receiving data file ''%s''\n', edfFile);
                            status=Eyelink('ReceiveFile'); % try to retrieve edf file (not necessary for the analysis, but useful as a backup)
                            if status > 0
                                fprintf('ReceiveFile status %d\n', status);
                            end
                            if 2==exist(edfFile, 'file')
                                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd); % if successful, indicate where the edf can be found
                            end
                        catch
                            fprintf('Problem receiving data file ''%s''\n', edfFile ); % if fail, print error
                        end
                    end
                end
                fclose('all');
                sca; % close the PTB screen
                return
            elseif ansKey == 6 % C (calibration)
                EyelinkrackerSetup(el); % re-calibrate the eye-tracker
                EyelinkDoDriftCorrection(el);
            elseif ansKey == 44 % SPACEBAR
                break; % continue the data acquisition normally
            elseif ansKey == 22 % S (skip next trial)
                skip_trial = 1; % flag to skip next trial
                break;
            else
                DrawFormattedText(w, '[WARNING] INVALID INPUT','center','center', white); % if none of the keys pressed by the user is valid, print a warning and wait for another input
                Screen('Flip', w);
                WaitSecs(1);
            end
        end
        
        %% Display moving stimulus on screen
        vbl=Screen(w,'Flip'); % vbl is a time-stamp associated with the "flip" of the PTB screen
        rect_blob = CenterRectOnPoint(rect_blob,xCenter,yCenter);
        Eyelink('StartRecording'); % start recording eye position
        WaitSecs(0.1); % record a few samples before we actually start displaying
        Eyelink('Message', ['START_TRIAL_' num2str(k) '_' vel_string '_' vfd_string '_' lum_string '_' pat_string]); % mark zero-plot time in data file
        SetMouse(xCenter,yCenter,w);
        HideCursor;
        for i = 1:duration_frames
            if dummymode
                rectim = CenterRectOnPoint(rectim, x_mouse, y_mouse);
                Screen('DrawTexture', w, texim, [], rectim);
            end
            if skip_trial
                skip_trial = 0;
                break;
            end
            [~,~,keyCode,~]=KbCheck;
            keyPressed = find(keyCode);
            if keyPressed == 41 % ESC during the trial for quitting
                sca;
                go_out = 1;
                break;
            end
            if keyPressed == 22 % S skip trial
                break;
            end
            rect_blob=CenterRectOnPoint(rect_blob,x_pos(i),y_pos(i));
            Screen('DrawTexture', w, tex_blob,[],rect_blob);
            if dummymode % if in dummy mode, treat the mouse coordinates as if they were the gaze position
                [x_mouse, y_mouse]=GetMouse(w);
                x_stim(i)=x_pos(i);
                y_stim(i)=y_pos(i);
                x_resp(i)=x_mouse;
                y_resp(i)=y_mouse;
                switch vfd % if the condition requires a simulated VFD, display on screen the appropriate texture that we created before
                    case 0
                    case 1
                        dstRectScot=CenterRectOnPoint(dstRectScot,x_mouse,y_mouse); % the position of the VFD will be controlled by mouse
                        Screen('DrawTexture', w, texScot,[],dstRectScot); % draw the VFD on the screen
                    case 2
                        dstRectScot=CenterRectOnPoint(dstRectScot,x_mouse,y_mouse); % the position of the VFD will be controlled by mouse
                        Screen('DrawTexture', w, texScot,[],dstRectScot); % draw the VFD on the screen
                    case 3
                        rectHemi = CenterRectOnPoint(rectHemi,x_mouse+xScreen-25, yCenter); % the position of the VFD will be controlled by mouse
                        Screen('DrawTexture', w, texHemi, [], rectHemi)  % draw the VFD on the screen
                end
            else % if the EyeLink is connected
                if Eyelink( 'NewFloatSampleAvailable') > 0  % check for presence of a new sample update
                    evt = Eyelink( 'NewestFloatSample'); % get the sample in the form of an event structure
                    if eye_used ~= -1 % do we know which eye to use yet?
                        x = evt.gx(eye_used+1); % if we do, get current gaze position from sample
                        y = evt.gy(eye_used+1); % +1 as we're accessing MATLAB array
                        if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0  % do we have valid data and is the pupil visible?
                            switch vfd % if yes, check if we need to display a VFD or not
                                case 0 % if no VFD is required, just skip to storing coordinates
                                case 1
                                    dstRectScot=CenterRectOnPoint(dstRectScot,x,y); % the position of the VFD will be controlled by gaze
                                    Screen('DrawTexture', el.window, texScot,[],dstRectScot); % draw the VFD on the screen
                                case 2
                                    dstRectScot=CenterRectOnPoint(dstRectScot,x,y); % the position of the VFD will be controlled by gaze
                                    Screen('DrawTexture', el.window, texScot,[],dstRectScot); % draw the VFD on the screen
                                case 3
                                    rectHemi = CenterRectOnPoint(rectHemi,x+xScreen-25, yCenter); % the position of the VFD will be controlled by gaze
                                    Screen('DrawTexture', el.window, texHemi, [], rectHemi); % draw the VFD on the screen
                            end
                            x_stim(i)=x_pos(i); % store stimulus coordinates
                            y_stim(i)=y_pos(i);
                            x_resp(i)=x; % store gaze coordinates
                            y_resp(i)=y;
                        else % if data is invalid (e.g. during a blink) draw scotoma in the last valid position
                            if i>1
                                x_stim(i)=x_pos(i);
                                y_stim(i)=y_pos(i);
                                x_resp(i) = x_resp(i-1);
                                y_resp(i) = y_resp(i-1);
                            end
                        end
                    else % if we don't know what eye is being tracked, find out
                        eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
                        if eye_used == el.BINOCULAR % if both eyes are tracked
                            eye_used = el.LEFT_EYE; % use left eye
                        end
                    end
                end
            end
            vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi); % flip the screen, updating the position of the stimulus and the vbl timestamp
        end
        
        %% Store gaze and stimulus coordinates data
        Eyelink('Message', ['END_TRIAL_' num2str(k) '_' vel_string '_' vfd_string '_' lum_string '_' pat_string]); % print a message in the .edf file to mark the end of the data acquisition
        Eyelink('StopRecording');
        if ~go_out
            DrawFormattedText(w, 'Saving data...','center','center', [255 255 255]);
            Screen('Flip',w);
        end
        try % try to fetch edf file
            fprintf('Receiving data file ''%s''\n', edfFile);
            status=Eyelink('ReceiveFile'); % try to retrieve edf file (not necessary for the analysis, but useful as a backup)
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(edfFile, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd); % if successful, indicate where the edf can be found
            end
        catch
            fprintf('Problem receiving data file ''%s''\n', edfFile ); % if fail, print error
        end
        x_resp = x_resp(1:length(x_stim));
        y_resp = y_resp(1:length(y_stim));
        filename = [myPath 'correlograms/' subj vel_string vfd_string lum_string pat_string num2str(k)]; % concatenate all the condition strings into a single filename
        try
            save(filename, 'ifi','hz', 'vel', 'pat', 'vfd','lum','monitWidth', 'viewDist', 'myPath', 'subj', 'x_stim', 'x_resp', 'y_stim', 'y_resp', 'xScreen', 'yScreen', 'xCenter', 'yCenter');
        catch
            mkdir(myPath, 'correlograms'); % if the folder "correlogram" doesn't exist yet, create it
            save(filename, 'ifi','hz', 'vel', 'pat', 'vfd','lum','monitWidth', 'viewDist', 'myPath', 'subj', 'x_stim', 'x_resp', 'y_stim', 'y_resp', 'xScreen', 'yScreen', 'xCenter', 'yCenter');
        end
    end
end
ShowCursor;
ListenChar;
cd(myPath);
sca;

