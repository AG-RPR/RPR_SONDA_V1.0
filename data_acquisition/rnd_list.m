% function list  = rnd_list(subj, myPath, [pat], [lum], [vel], [vfd])
% 
% subj: string containing subject ID
% myPath: string containing the path where the list is created
% Creates a list with all the possible (unique) combinations of conditions. 
% Default conditions are:
%
% pat = [0 1];      % smooth and saccadic pursuits (default: both)
% lum = [0 1];      % low and high contrast (default: high)
% vel = [1 2 3 4];  % velocity levels. The default velocity is level 3
% vfd = [0 1 2 3];  % simulated gaze-contingent visual field defect. 0: No defect | 1: central defect | 2: peripheral defect | 3: hemifield defect (default: 0)

function rnd_list_out = rnd_list(subj, myPath, varargin)
%% initialize default parameters
pat = [0 1];
lum = 1;
vel = 3;
vfd = 0;

%% parse input arguments
if nargin >= 3
    switch numel(varargin)
        case 1
            if ~isempty(varargin{1}); pat = varargin{1}; end
        case 2
            if ~isempty(varargin{1}); pat = varargin{1}; end
            if ~isempty(varargin{2}); lum = varargin{2}; end
        case 3
            if ~isempty(varargin{1}); pat = varargin{1}; end
            if ~isempty(varargin{2}); lum = varargin{2}; end
            if ~isempty(varargin{3}); vel = varargin{3}; end
        case 4
            if ~isempty(varargin{1}); pat = varargin{1}; end
            if ~isempty(varargin{2}); lum = varargin{2}; end
            if ~isempty(varargin{3}); vel = varargin{3}; end
            if ~isempty(varargin{4}); vfd = varargin{4}; end            
    end
end
%% create list of conditions
list = combvec(pat, lum, vel, vfd)';
if size(list,1) > 1
    rnd_list_out = Shuffle(list,2); % Shuffle is a Psychtoolbox function
else
    rnd_list_out = list;
end
try
    fl=fopen([myPath 'lists/' subj '_list.txt'],'w+');
    fprintf(fl,'%s\t%s\t%s\t%s\n', 'pat', 'lum', 'vel', 'vfd');
catch
    mkdir(myPath,'lists');
    fl=fopen([myPath 'lists/' subj '_list.txt'],'w+');
    fprintf(fl,'%s\t%s\t%s\t%s\n', 'pat', 'lum', 'vel', 'vfd');
end
for ii = 1 : size(rnd_list_out,1)
    fprintf(fl,'%d\t%d\t%d\t%d\n',rnd_list_out(ii,1), rnd_list_out(ii,2), rnd_list_out(ii,3), rnd_list_out(ii,4));
end

fclose(fl);

end