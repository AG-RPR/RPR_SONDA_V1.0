function [open_ui, filename] = load_data(varargin)

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
% ui: ignore the input arguments and load data by User Interface

for k=1:2:nargin
    switch varargin{k}
        case {'subj'}
            subj=varargin{k+1};
        case {'vel'}
            vel_f=varargin{k+1};
            vel_str = ['v' num2str(vel_f)];
            info_data.velocity = vel_str;
        case {'vfd'}
            cond_f=varargin{k+1};
        case {'contrast'}
            lum_f=varargin{k+1};
        case {'pursuit'}
            pat_f=varargin{k+1};
        case {'plot', 'Plot'}
            do_plot=varargin{k+1};
        case {'ui'}
            open_ui = 1;
        otherwise
            error(sprintf('invalid parameter %s\n',varargin{k}));
    end
end

if exist('open_ui','var')
    filename=0;
    return;
else 
    open_ui = 0;
end

if exist('subj','var')
    loadmat = 1; % we're going to open a .mat file
    if subj < 10 % if the input to subj is numeric and less than 10, pad it with a 0 and convert it to string
        subj_str = ['0' num2str(subj)];
        info_data.subjID = subj_str;
    else
        subj_str = num2str(subj);
        info_data.subjID = subj_str;
    end
end

if exist('cond_f','var')
    switch cond_f
        case 0
            cond_str = 'c'; % "Control
        case 1
            cond_str = 'm'; % "Macular" (central defect)
        case 2
            cond_str = 'g'; % "Glaucoma" (peripheral defect)
        case 3
            cond_str = 'h'; % "Hemianopia" (hemifield defect)
    end
    info_data.condition = cond_str;
end

if exist('pat_f','var')
    switch pat_f
        case 0
            pat_str= {'sm'};
        case 1
            pat_str= {'sc'};
    end
    info_data.pursuit = pat_str;
end

if exist('lum_f','var')
    switch lum_f
        case 0
            lum_str = 'l';
        case 1
            lum_str = 'h';
    end
    info_data.contrast = lum_str;
end

filename = strcat(subj_str,vel_str,cond_str, lum_str, pat_str);

end