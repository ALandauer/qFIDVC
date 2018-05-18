function [u, cc, dm, m, tSwitch] = funIDVC(varargin)
% u = funIDVC(filename, sSize, runMode) is the main function that performs
% IDVC on a time increment of volumetric images.  
% 
% INPUTS
% -------------------------------------------------------------------------
%   filename: string for the filename prefix for the volumetric images in 
%             the current directory.  
%             Input options:
%             --- If image is not within a cell) ---
%             1) 'filename*.mat' or 'filename*' 
% 
%             --- If image is within a cell that contains multichannels ---
%             2) filename{1} = 'filename*.mat' or 'filename*' and
%                filename{2} = channel number containing images you want to
%                              run IDVC on.
%                (if the channel is not provided, i.e. length(filename) = 1
%                , then channel = 1
%                 
%   sSize: interrogation window (subset) size for the first iterations.  
%          Must be, 32,64,96, or 128 voxels and a three column
%          array (one for each dimenision) or scalar (equal for all
%          dimensions).
%   runMode: string that defines the method of running IDVC. Options: 
%             cumulative (time0 -> time1, time0 -> time2, ...)
%             (Allowable inputs: 'c','cum','cumulative')
%             or 
%             incremental (time0 -> time1, time1 -> time2, ...)
%             (Allowable inputs: 'i','inc','incremental')
%             
% OUTPUTS
% -------------------------------------------------------------------------
%   u: displacement field vector defined at every meshgrid point with 
%      spacing dm. Format: cell array, each containing a 3D matrix for each
%      time point
%         (components in x,y,z)
%         u{time}{1} = displacement in x-direction
%         u{time}{2} = displacement in y-direction
%         u{time}{3} = displacement in z-direction
%         u{time}{4} = magnitude
%   cc: peak values of the cross-correlation for each interrogation point
%   dm: meshgrid spacing (8 by default)
%   m: mesh grid used for the in the displacement grid
%
% NOTES
% -------------------------------------------------------------------------
% none
% 
% For more information please see 
%

%% ---- Opening & Reading the First Image into CPU Memory ----
[fileInfo, sSize0, sSizeMin, runMode, u_] = parseInputs(varargin{:});

I{1} = loadFile(fileInfo,1);

%% ---- Opening and Reading Subsequent Images ---
numImages = length(fileInfo.filename);
u = cell(numImages-1,1); cc = cell(numImages-1,1);
cc = cell(numImages-1,1);
                  
for i = 2:numImages % Reads Volumes Starting on the Second Volumes
    tStart = tic;
    I{2} = loadFile(fileInfo,i);
    
    %Start DVC
    disp(['Current file: ' fileInfo.filename{i}])
    [u_, cc{i-1}, dm, m, tSwitch(i-1)] = IDVC(I,sSize0,sSizeMin,u_,runMode);
    
    % Saving iterations of the DVC
    u{i-1}{1} = -u_{1};  u{i-1}{2} = -u_{2};  u{i-1}{3} = -u_{3}; u{i-1}{4} = u_{4};
    
    if strcmpi(runMode(1),'i')
        I{1} = I{2};
        u_ = num2cell(zeros(1,3));

    end
    
    if strcmpi(runMode(1),'h') == 1 % If hybrid mode
        if tSwitch(end) == 1
            disp('Updating reference image due to poor correlation.')
            disp('Rerunning DIC on the time step')
            I{1} = loadFile(fileInfo,i-1); %Update reference image
            u_ = num2cell(zeros(1,4)); % No predictor displacement in next time step
                  
            % Run DIC again on updated time step
            [u_,cc{i-1},dm,m,tSwitch(i-1)] = IDVC(I,sSize0,sSizeMin,u_);
            tSwitch(i-1) = tSwitch(i-1) + 1;
            
            %Save iterations of DIC
            u{i-1}{1} = -u_{1};  u{i-1}{2} = -u_{2}; u{i-1}{3} = -u_{3}; u{i-1}{4} = u_{4};
        end
    else
        u_ = u_(1:3);
    end
                  
    disp(['Elapsed time for all iterations: ',num2str(toc(tStart))]);
end
                  if strcmpi(runMode(1),'h') == 1
                  % Update displacement for to cumulative
                  option = 'spline';
                  tStart = tic;
                  disp('Update displacement for Hybrid mode');
                  [u] = FIDVCinc2cum(u,tSwitch,dm,m,option);
                  disp(['Elapsed Time for updating displacements: ',num2str(toc(tStart))]);
                  end
                  
end

function I = loadFile(fileInfo,idx)
I = load(fileInfo.filename{idx});
fieldName = fieldnames(I);
I = getfield(I,fieldName{1});
if iscell(I)
    if numel(I), I = I{1};
    else
        I = I{fileInfo.dataChannel};
    end
end
end

function varargout = parseInputs(varargin)
%  = parseInputs(filename, sSize, sSizeMin, runMode)


% Parse filenames
filename = varargin{1};
if iscell(filename)
    if length(filename) == 1, fileInfo.datachannel = 1;
    else fileInfo.datachannel = filename{2};
    end
    filename = filename{1};
end


[~,filename,~] = fileparts(filename);
filename = dir([filename,'.mat']);
fileInfo.filename = {filename.name};

if isempty(fileInfo), error('File name doesn''t exist'); end

% Ensure dimensionality of the subset size
sSize = varargin{2};
if numel(sSize) == 1
    sSize = sSize*[1 1 1];
elseif numel(sSize) ~=3
    error('Subset size must be a scalar or a three column array');
end

%get minium subset size
sSizeMin = varargin{3};
                  
% Ensure range of subset size
if min(sSize) < 16 || max(sSize > 128)
   error('Subset size must be within 16 and 128 pixels');
end

% Ensure even subset size
% if sum(mod(sSize,4)) > 0
%     error('Subset size must be even');
% end

if sum(mod(sSize,32)) ~= 0
    error('Subset size must be 16, 32, 64, 96, or 128 voxels in each dimension');
end

% Check run method input
runMode  = varargin{4};

if ~(strcmpi(runMode(1),'c') || strcmpi(runMode(1),'i') || strcmpi(runMode(1),'h'))
    error('Run method must be incremental or cumulative or hybrid');
end

% Initial guess of displacement field = [0 0 0];
u0 = num2cell(zeros(1,3));

% Outputs
varargout{      1} = fileInfo;
varargout{end + 1} = sSize;
varargout{end + 1} = sSizeMin;
varargout{end + 1} = runMode;
varargout{end + 1} = u0;

end
