function [u, cc, dm, mFinal, decorrFlag] = IDVC(varargin)
% [u, cc] = IDVC(I,sSize,u0,className);
% I = filterDisplacements(I0,filterSize,z) applies a low-pass convolution
% filter to the displacement field to mitigate divergence based on
%
% F. F. J. Schrijer and F. Scarano. Effect of predictor corrector filtering
% on the stability and spatial resolution of iterative PIV interrogation.
% Exp. Fluids, 45(5):927{941, May 2008. doi: 10.1007/s00348-008-0511-7
%
% INPUTS
% -------------------------------------------------------------------------
%   I0: cell containing the undeformed, I0{1}, and deformed, I0{2} 3-D
%       images
%   sSize: interrogation window (subset) size
%   u0: pre-estimated displacement field (typically zeros)
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: displacement field vector defined at every meshgrid point with
%      spacing dm. Format: cell array, each containing a 3D matrix
%         (components in x,y,z)
%         u{1} = displacement in x-direction
%         u{2} = displacement in y-direction
%         u{3} = displacement in z-direction
%         u{4} = magnitude
%   cc: peak values of the cross-correlation for each interrogation
%   dm: meshgrid spacing (8 by default)
%
% NOTES
% -------------------------------------------------------------------------
% To run you need a compatible C compiler. Please see
% (http://www.mathworks.com/support/compilers/R2014a/index.html)
%
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast
% iterative digital volume correlation algorithm for large deformations.
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2


% PRESET CONSTANTS
maxIterations = 10; % maximum number of iterations
dm = 8; % desired output mesh spacing
convergenceCrit = [0.25, 0.5, 0.0625]; % convergence criteria
ccThreshold = 0.75; % bad cross-correlation threshold (mean - ccThreshold*stdev for q-factor distribution)

sizeThresh = 1331;  %Threshold for maximum size of bad correlation regions
percentThresh = 25; %Threshold for max % of total measurement pts failing q-factor testing
stDevThresh = 0.15; %Threshold for max st. dev. of the fitted gaussian to the second peak

cc = cell(1);
cc{1} = struct('max',[],'maxIdx',[],'qfactors',[],'q_thresh',[],...
    'qfactors_accept',[]);


[I0, sSize, sSizeMin, sSpacing, padSize, DVCPadSize, u] = parseInputs(varargin{:});

% START ITERATING
i = 2; converged01 = 0; SSE = []; I = I0;

t0 = tic;
while ~converged01 && i - 1 < maxIterations
    ti = tic;
    
    % Check for convergence
    [converged01, SSE(i-1) , sSize(i,:), sSpacing(i,:)] = ...
        checkConvergenceSSD(I,SSE,sSize,sSizeMin,sSpacing,convergenceCrit);
    
    if ~converged01
        finalSize = sSize(i,:);
        
        [I, m] = parseImages(I,sSize(i,:),sSpacing(i,:));
        %warp images with displacement guess if a cumulative step
        if i == 2 %only on the first iteration
            if numel(u{1}) == 1 %on the first image the disp guess is zero
                u{1} = zeros(length(m{1}),length(m{2}),length(m{3}));
                u{2} = zeros(length(m{1}),length(m{2}),length(m{3}));
                u{3} = zeros(length(m{1}),length(m{2}),length(m{3}));
            end
        else
            u{1} = inpaint_nans3(u{1}); %inpaint nans from last time's edge pts
            u{2} = inpaint_nans3(u{2});
            u{3} = inpaint_nans3(u{3});
            I = volumeMapping(I,m,u); %otherwise map the images w/ the initial guess
            [I, m] = parseImages(I,sSize(i,:),sSpacing(i,:));
            %             u{1} = 0; %reset disp to zero
            %             u{2} = 0;
        end
        
        % run cross-correlation to get an estimate of the displacements
        [du, cc{i-1}] = DVC(I,sSize(i,:),sSpacing(i,:),DVCPadSize,ccThreshold);
        
        % add the displacements from previous iteration to current
        [u, ~, cc{i-1}, mFinal] = addDisplacements(u,du,cc{i-1},m,dm,i-1);
        
        % filter the  displacements using a predictor filter
        u = filterDisplacements(u,sSize(i,:)/dm);
        
        %flag all outliers (q-factor and fluctuation-based)
        [cc{i-1},~] = flagOutliers(u,cc{i-1},1.25,0.1);
        
        % remove outliers in displacement field
        [u,cc{i-1}.alpha_mask,cc{i-1}.nan_mask,cc{i-1}.edge_pts] = replaceOutliers(u,cc{i-1});
        %         u = removeOutliers(u);
        
        % mesh and pad images based on new subset size and spacing
        [I, mFinal] = parseImages(I0,sSize(i,:),sSpacing(i,:));
        
        % map volumes based on displacment field
        I = volumeMapping(I,mFinal,u);
        
        disp(['Elapsed time (iteration ',num2str(i-1),'): ',num2str(toc(ti))]);
        i = i + 1;
    end
    
end

%prep outputs
decorrFlag = decorrelationCheck(cc,sizeThresh,percentThresh,stDevThresh);
u{1} = cc{end}.edge_pts.*u{1};
u{2} = cc{end}.edge_pts.*u{2};
[u,cc,mFinal] = parseOutputs(u,cc,finalSize,padSize,mFinal);

disp(['Convergence at iteration ',num2str(i)]);
disp(['Total time: ',num2str(toc(t0))]);
end



%% ========================================================================
function varargout = parseImages(varargin)
% pads images and creates meshgrid

I{1} = single(varargin{1}{1});
I{2} = single(varargin{1}{2});
sSize = varargin{2};
sSpacing = varargin{3};


prePad = sSize/2;
postPad = sSize/2;

sizeI = size(I{1});
I{1} = padarray(I{1},prePad,'replicate','pre');
I{1} = padarray(I{1},postPad,'replicate','post');

I{2} = padarray(I{2},prePad,'replicate','pre');
I{2} = padarray(I{2},postPad,'replicate','post');


idx = cell(1,3);
for i = 1:3, idx{i} = (1:sSpacing(i):(sizeI(i) + 1)) + sSize(i)/2; end

% [m{1},m{2},m{3}] = meshgrid(idx{:});

varargout{    1} = I;
varargout{end+1} = idx;

end

%% ========================================================================
function varargout = parseInputs(varargin)
% parses inputs and pads images so that there is an divisable meshgrid number.

I0{1} = single(varargin{1}{1});
I0{2} = single(varargin{1}{2});

% I0{1} = permute(I0{1},[2 1 3]);
% I0{2} = permute(I0{2},[2 1 3]);

sSize = varargin{2};
sSize = [sSize(2), sSize(1), sSize(3)];

sSizeMin = varargin{3};

sSpacing = sSize/2;
u0 = varargin{4};

DVCPadSize = sSpacing/2;

sizeI0 = size(I0{1});
sizeI = ceil(sizeI0./sSpacing).*sSpacing;
prePad = ceil((sizeI - sizeI0)/2);
postPad = floor((sizeI - sizeI0)/2);

I{1} = padarray(I0{1},prePad,'replicate','pre');
I{1} = padarray(I{1},postPad,'replicate','post');

I{2} = padarray(I0{2},prePad,'replicate','pre');
I{2} = padarray(I{2},postPad,'replicate','post');

varargout{    1} = I;
varargout{end+1} = sSize;
varargout{end+1} = sSizeMin;
varargout{end+1} = sSpacing;
varargout{end+1} = [prePad; postPad];
varargout{end+1} = DVCPadSize;
varargout{end+1} = u0;
end


function [u,cc,m] = parseOutputs(u,cc,finalsSize,padSize,m)
% parses outputs and unpads the displacment field and cc.

% % % % unpadSize(1,:) = ceil(padSize(1,:)/filterSpacing);
% % % % unpadSize(2,:) = floor((padSize(2,:)+1)/filterSpacing);
% % % % % +1 from the extra meshgrid point during the meshing of the DVC algorithm. EBK (10-23-2013)
% % % %
% % % % for i = 1:3
% % % %     u{i} = u{i}(1+unpadSize(1,1):end-unpadSize(2,1),...
% % % %         1+unpadSize(1,2):end-unpadSize(2,2),...
% % % %         1+unpadSize(1,3):end-unpadSize(2,3));
% % % % end
% % % % u{4} = sqrt(u{1}.^2 + u{2}.^2 + u{3}.^2);
% % % %
% % % % cc = cc(1+unpadSize(1,1):end-unpadSize(2,1),...
% % % %     1+unpadSize(1,2):end-unpadSize(2,2),...
% % % %     1+unpadSize(1,3):end-unpadSize(2,3));
filterSpacing = finalsSize/2;
u{4} = sqrt(u{1}.^2 + u{2}.^2 + u{3}.^2);
for i = 1:3
    m{i} = m{i}-padSize(1,i)-filterSpacing(i);
end

%remove the memory-intesive and generally unneeded part of the cc struct
for jj = 1:length(cc)
    cc{jj}.max = [];
end

end



function decorrFlag = decorrelationCheck(cc,sizeThresh,percentThresh,stDevThresh)
%This function takes in a cc structure and checks the q-factor maps for
%decorrelation indications.

%Initilize
decorrFlag = 0;

if nargin < 4
    stDevThresh = 0.135;
elseif nargin < 3
    percentThresh = 15;
elseif nargin < 2
    sizeThresh = 1300;
end

trim_size = ceil(cc{end}.sSize/2)+1;

%get the q-factors maps, exclude boundary pixels since they are
%reliably poor.
qf1 = cc{end}.qfactors_accept{1}(1+trim_size(1):end-trim_size(1),...
    1+trim_size(2):end-trim_size(2),1+trim_size(3):end-trim_size(3));

qf2 = cc{end}.qfactors_accept{2}(1+trim_size(1):end-trim_size(1),...
    1+trim_size(2):end-trim_size(2),1+trim_size(3):end-trim_size(3));

%compute the % of bad correlation pts
num_pts = numel(qf1)+numel(qf2);
num_fail = sum(isnan(qf1(:))+isnan(qf2(:)));

if num_fail/num_pts*100 >= percentThresh
    decorrFlag = 1;
    disp('percent decorr')
end

%binarize the qf maps
qf1_bin = isnan(qf1);
qf2_bin = isnan(qf2);

%get connect regions
qf1_connectivity = bwconncomp(qf1_bin,8);
qf2_connectivity = bwconncomp(qf2_bin,8);

%find the number of px in each connect region
qf1_numPixels = cellfun(@numel,qf1_connectivity.PixelIdxList);
qf2_numPixels = cellfun(@numel,qf2_connectivity.PixelIdxList);

%find the largest area of the connected regions
[~,idx1] = max(qf1_numPixels);
[~,idx2] = max(qf2_numPixels);

%compute the region properties for the largest region, normalized by is
%eccentricity (since its harder to correctly predict far away from known data
%with inpaint nans)
if ~isempty(idx1)
    S1 = regionprops(qf1_connectivity,'Area');
    area(1) = S1(idx1).Area;
    norm_area(1) = area(1);
else
    norm_area(1) = 0;
end

if ~isempty(idx2)
    S2 = regionprops(qf2_connectivity,'Area');
    area(2) = S2(idx2).Area;
    norm_area(2) = area(2);
else
    norm_area(2) = 0;
end
% qf_holeSize(idx1) = 0;
% [qf_second,idx2] = max([qf1_numPixels,qf2_numPixels]);
% norm_area
if max(norm_area(:)) > sizeThresh
    decorrFlag = 1;
    disp('size decorr')
end

% do the st dev based thresholding
complete_qfactors = cc{end}.qfactors(1:2,:);
for ii = 1:size(complete_qfactors,1)
    
    x = sort(complete_qfactors(ii,:));
    
    %do the parameter estimation
    options = statset('MaxIter',2000, 'MaxFunEvals',4000);
    % options.FunValCheck = 'off';
    
    %     disp('Parameters estimated for single peak Gaussian')
    paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
    paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
        paramEsts(2)];
    
    qfactor_fit_means(ii) = paramEsts(3);
    qfactor_fit_stds(ii) = paramEsts(5);
    
end

if cc{end}.sSize(1) == 16
    stDevThresh = stDevThresh + 0.01;
elseif cc{end}.sSize(1) == 64
    stDevThresh = stDevThresh - 0.01;
end

if qfactor_fit_stds(1) >= stDevThresh || qfactor_fit_stds(2) >= stDevThresh
    decorrFlag = 1;
    disp('stdev decorr')
end


end