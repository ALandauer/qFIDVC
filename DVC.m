function [u, cc] = DVC(varargin)
% [du, cc] = DVC(I,sSize,sSpacing,ccThreshold) estimates
% displacements between two volumetric images through digital volume
% correlation.
%
% INPUTS
% -------------------------------------------------------------------------
%   I: cell containing the undeformed, I{1}, and deformed, I{2} 3-D images
%   sSize: interrogation window (subset) size
%   sSpacing: interrogation window (subset) spacing.  Determines window
%             overlap factor
%   ccThreshold: threshold value that defines a bad cross-correlation
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the displacement field (u{1:3} = {u_x, u_y, u_z})
%   cc: peak values of the cross-correlation for each interrogation
%
% NOTES
% -------------------------------------------------------------------------
% all functions are self contained
%
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast
% iterative digital volume correlation algorithm for large deformations.
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

% Parse inputs and create meshgrid
[I,m,mSize,sSize,sSizeFull,sSpacing,MTF,M,ccThreshold] = parseInputs(varargin{:});

% Initialize variables
mSize_ = prod(mSize);
u123 = zeros(mSize_,3);
% cc = zeros(mSize_,1);
cc = struct('max',[],'sSpacing',[],'sSize',[],'maxIdx',[],'qfactors',[],'q_thresh',[],...
    'qfactors_accept',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%inject small noise to avoid flat subsets
noise1 = rand(size(I{1}))/10000;
noise2 = rand(size(I{2}))/10000;
I{1} = I{1}+noise1;
I{2} = I{2}+noise2;

% A = cell(1,mSize_);

for k = 1:mSize_
    
    
    tStart = tic; % begin timer
    %--------------------------------------------------------------------------
    % grab the moving subset from the images
    subst = I{1}(m{1}(k,:),m{2}(k,:),m{3}(k,:));
    B = I{2}(m{1}(k,:),m{2}(k,:),m{3}(k,:));
    
    % multiply by the modular transfer function to alter frequency content
    subst = MTF.*subst; B = MTF.*B;
    
    % run cross-correlation
    cc.A{k} = xCorr3(subst,B,sSize);

    % find maximum index of the cross-correlaiton
    [max_val(k), maxIdx{k}] = max(cc.A{k}(:));
    
    % compute voxel resolution displacements
    [u1, u2, u3] = ind2sub(sSize,maxIdx{k});
    
    % gather the 3x3x3 voxel neighborhood around the peak
    try xCorrPeak = reshape(cc.A{k}(u1 + (-1:1), u2 + (-1:1), u3 + (-1:1)),27,1);
        % last squares fitting of the peak to calculate sub-voxel displacements
        du123 = lsqPolyFit3(xCorrPeak, M{1}, M{2});
        u123(k,:) = [u1 u2 u3] + du123' - sSize/2 - 1;
        %--------------------------------------------------------------------------
    catch
        u123(k,:) = nan;
    end
    %     xCorrPeak = reshape(A(u1 + (-1:1), u2 + (-1:1), u3 + (-1:1)),27,1);
    %
    %     % least squares fitting of the peak to calculate sub-voxel displacements
    %     du123 = lsqPolyFit3(xCorrPeak, M{1}, M{2});
    %     u123(k,:) = [u1 u2 u3] + du123' - (sSize/2) - 1;
    %--------------------------------------------------------------------------
    
    
end

cc.maxIdx = maxIdx;
cc.max = max_val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reshape displacements and set bad correlations to zero

if size(sSpacing,1) == 1
    spacingChange = 1;
elseif sSpacing(end,1) == sSpacing(end-1,1)
    spacingChange = 0;
else
    spacingChange = 1;
end

cc.max = reshape(double(cc.max),mSize);
[cc, ccMask_] = ...
    removeBadCorrelations(cc,ccThreshold,spacingChange,mSize);
cc.sSpacing = sSpacing(end,:);
cc.sSize = sSize;

% % cc = reshape(double(cc),mSize);
% cc = permute(cc,[2 1 3]);
% % [cc, ccMask] = removeBadCorrelations(I,cc,ccThreshold);

for ii = 1:2
    ccMask{ii} = reshape(double(ccMask_(:,ii)),mSize);
end

% u = cell(1,3);
% for i = 1:3
%     u{i} = reshape(double(u123(:,i)),mSize).*ccMask;
%     u{i} = permute(u{i},[2 1 3]);
% end
u{1} = reshape(double(u123(:,2)),mSize);%.*ccMask;
u{2} = reshape(double(u123(:,1)),mSize);%.*ccMask;
u{3} = reshape(double(u123(:,3)),mSize);%.*ccMask;

end

%% ========================================================================
function varargout = parseInputs(varargin)
% Parse inputs and create meshgrid

I{1} = varargin{1}{1};
I{2} = varargin{1}{2};
sSizeFull = varargin{2};
sSpacing = varargin{3};
padSize = varargin{4};
ccThreshold = varargin{5};

% pad images with zeros so that we don't grab any subset outside of the image
% domain. This would produce an error
I{1} = padarray(I{1},padSize,'replicate','both');
I{2} = padarray(I{2},padSize,'replicate','both');
sizeV = size(I{1});

sSize = sSizeFull(end,:);

% Initialize Mesh Variables
idx = cell(1,3);
for i = 1:3
    idx{i} = (1+padSize(i)) : sSpacing(i) : (sizeV(i)-sSize(i)-padSize(i)+1);
end
[m{1},m{2},m{3}] = ndgrid(idx{:});


% sSize = [sSize(2) sSize(1) sSize(3)];
mSize = size(m{1});
mSize_ = prod(mSize);

m_ = cell(1,3);
for k = 1:3
    m_{k} = zeros([mSize_,sSize(k)],'uint16');
    repmat_ = repmat((1:sSize(k))-1,mSize_,1);
    m_{k} = bsxfun(@plus, repmat_,m{k}(:));
end

% Initialize quadratic least squares fitting coefficients
[mx, my, mz] = meshgrid((-1:1),(-1:1),(-1:1));
m = [mx(:), my(:), mz(:)];

M{1} = zeros(size(m,1),10);
for i = 1:size(m,1)
    x = m(i,1); y = m(i,2); z = m(i,3);
    M{1}(i,:) = [1,x,y,z,x^2,x*y,x*z,y^2,y*z,z^2];
end

M{2} = M{1}'*M{1};

% Generate Moduluar transfer function (see eq. 3)
[~,~,MTF] = generateMTF(sSize);

%% Parse outputs

varargout{    1} = I;
varargout{end+1} = m_;
varargout{end+1} = mSize;
varargout{end+1} = sSize;
varargout{end+1} = sSizeFull;
varargout{end+1} = sSpacing;
varargout{end+1} = MTF;
varargout{end+1} = M;
varargout{end+1} = ccThreshold;

end

%% ========================================================================
function A = xCorr3(A,B,sSize)
% performs fft based cross correlation of A and B (see equation 2)

A = fftn(A,sSize);
B = fftn(B,sSize);
B = conj(B);
A = A.*B;
A = ifftn(A);
A = real(A);
A = fftshift(A);
end

%% ========================================================================
function    duvw = lsqPolyFit3(b, M, trMM)
% LeastSqPoly performs a 3D polynomial fit in the least squares sense
% Solves M*x = b,
% trMM = transpose(M)*M
% trMb = tranpose(M)*b
%
% If you need to generate the coefficients then uncomment the following
% [mx, my, mz] = meshgrid(-1:1,-1:1,-1:1);
% m = [mx(:), my(:), mz(:)];
%
% for i = 1:size(m,1)
%    x = m(i,1); y = m(i,2); z = m(i,3);
%    M1(i,:) = [1,x,y,z,x^2,x*y,x*z,y^2,y*z,z^2];
% end
%
% trMM1 = M'*M;

% b = log(b);
trMb = sum(bsxfun(@times, M, b));

x = trMM\trMb'; %solve for unknown coefficients

A = [x(6), 2*x(5), x(7);
    2*x(8),  x(6) x(9)
    x(9),    x(7), 2*x(10)];

duvw = (A\(-x([2 3 4])));
end

%% ========================================================================
function varargout = generateMTF(sSize)
% MTF functions taken from
% J. Nogueira, A Lecuona, P. A. Rodriguez, J. A. Alfaro, and A. Acosta.
% Limits on the resolution of correlation PIV iterative methods. Practical
% implementation and design of weighting functions. Exp. Fluids,
% 39(2):314{321, July 2005. doi: 10.1007/s00348-005-1017-1

%% equation 4

if prod(single(sSize == 32)) || prod(single(sSize == 64)) || prod(single(sSize == 16))
    sSize = sSize(1);
    
    x = cell(1,3);
    [x{1}, x{2}, x{3}] = meshgrid(1:sSize,1:sSize,1:sSize);
    
    nu{1} = 1;
    for i = 1:3
        x{i} = x{i} - sSize/2 - 0.5;
        x{i} = abs(x{i}/sSize);
        nu{1} = nu{1}.*(3*(4*x{i}.^2-4*x{i}+1));
    end
    
    %% equation 5
    [x{1}, x{2}, x{3}] = meshgrid(1:sSize,1:sSize,1:sSize);
    
    for i = 1:3, x{i} = x{i} - sSize/2 - 0.5; end
    
    r = abs(sqrt(x{1}.^2 + x{2}.^2 + x{3}.^2)/sSize);
    nu{2}  = zeros(size(r));
    nu{2}(r < 0.5) = 24/pi*(4*r(r < 0.5).^2-4*r(r < 0.5)+1);
    
    %% equation 6
    [x{1}, x{2}, x{3}] = meshgrid(1:sSize,1:sSize,1:sSize);
    
    nu{3} = 1;
    for i = 1:3
        x{i} = x{i} - sSize/2 - 0.5;
        x{i} = (x{i}/sSize);
        
        nu{3} = nu{3}.*(12*abs(x{i}).^2 - 12*abs(x{i}) + 3 + ...
            0.15*cos(4*pi*x{i}) + 0.20*cos(6*pi*x{i}) + ...
            0.10*cos(8*pi*x{i}) + 0.05*cos(10*pi*x{i}));
        
    end
    nu{3}(nu{3} < 0) = 0;
    
else
    
    nu{1} = ones(sSize(1),sSize(2),sSize(3));
    nu{2} = nu{1};
    nu{3} = nu{1};
    
end

nu = cellfun(@(x) x/sum(x(:)), nu, 'UniformOutput',0);
nu = cellfun(@sqrt, nu, 'UniformOutput',0);

varargout = nu;

end

%% ========================================================================
function [cc, ccMask] = ...
    removeBadCorrelations(cc,ccThreshold,sizeChange,mSize)

%get peak locations and cc_min maps (i.e. cc - cc(min))
[peak,cc_min] = cellfun(@(x) cc_max_find(x),cc.A,'UniformOutput',0);

%compute two primary quality metrics, as given in "Xue Z, Particle Image
% Velocimetry Correlation Signal-to-noise Metrics, Particle Image
% Pattern Mutual Information and Measurement uncertainty Quantification.
% MS Thesis, Virginia Tech, 2014.

%peak to corr. energy ratio
pce = cellfun(@(x,y) (y^2)/(1/numel(x)*(sum(abs(x(:).^2)))),cc_min,peak,'UniformOutput',0);
%min value -> 1 (worst case)

%peak to entropy ratio
ppe = cellfun(@(x) q_entropy(x),cc_min,'UniformOutput',0);%peak to cc (information) entropy
%min value -> 0 (worst case)

qfactors_ = cell2mat(...
    cellfun(@(x,y) [x(:);y(:)], pce,ppe,'UniformOutput',0))';

for k = 1:2
    qf_ = (qfactors_(:,k)-min(qfactors_(:,k)));
    cc.qfactors(:,k) = qf_/max(qf_);
end

if sizeChange == 1
    %recompute threshold, only use pce & ppe since these give the best
    %results emprically.
    for ii = 1:2
        
        [qf_para{ii},single_distro] = bimodal_gauss_fit(cc.qfactors(ii,:));
        
        if single_distro == 0%(qf_para{ii}(2) + 2*qf_para{ii}(4)) < (qf_para{ii}(3) - 2*qf_para{ii}(5))
            cc.q_thresh{ii} = qf_para{ii}(3) - ccThreshold*qf_para{ii}(5);
        elseif single_distro == 1
            cc.q_thresh{ii} = qf_para{ii}(3) - ccThreshold*qf_para{ii}(5);
        else
            cc.q_thresh{ii} = qf_para{ii}(3) - ccThreshold*qf_para{ii}(5);
        end
    end
    q_trim = [cc.q_thresh{1};cc.q_thresh{2}];
else
    q_trim = [cc.q_thresh{1};cc.q_thresh{2}];
end

%NaN the qfactor values that are below the threshold
temp = bsxfun(@le,cc.qfactors(:,1:2),q_trim');
qfactors_accept = cc.qfactors(:,1:2);
qfactors_accept(temp) = NaN;

for ii = 1:2
    cc.qfactors_accept{ii} = reshape(double(qfactors_accept(:,ii)),mSize);
end

ccMask = ones(size(qfactors_accept)) + ...
    zeros(size(qfactors_accept)).*qfactors_accept;

end

%% ========================================================================
function [paramEsts,single_distro] = bimodal_gauss_fit(x)
%This function takes a dataset and fits a bimodal Gaussian distro to it.

x = sort(x);

%set function for bimodal Gaussian
pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
    p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
pdf_single = @(x,mu1,sigma1) ...
    normpdf(x,mu1,sigma1);

%starting params, biased mixture toward "good" values,
%centered at quartiles, equal std dev.
pStart = 0.25;
muStart = quantile(x,[.10 .75]);
sigmaStart(1) = sqrt(var(x(1:round(length(x)/5))));
%- 0.25*diff(quantile(x,[0.01 0.25])).^2);
sigmaStart(2) = sqrt(var(x(ceil(length(x)/10):ceil(3*length(x)/4))));
%... - 0.25*diff(quantile(x,[0.25 0.75])).^2);%1:round(length(x)/2)
start = [pStart muStart sigmaStart];

%set lower and upper bounds
lb = [0 -inf -inf 0.00001 0.00001];
ub = [1 inf inf inf inf];

%do the parameter estimation
options = statset('MaxIter',1800, 'MaxFunEvals',3600);
% options.FunValCheck = 'off';
try
    single_distro = 0;
    paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start, ...
        'lower',lb, 'upper',ub, 'options',options);%,'optimfun','fmincon'
    
    if paramEsts(2)-paramEsts(4) >= paramEsts(3)+paramEsts(5) || ...
            paramEsts(2)+paramEsts(4) <= paramEsts(3)-paramEsts(5)
        
        single_distro = 1;
        %     disp('Parameters estimated for single peak Gaussian')
        paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
        paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
            paramEsts(2)];
        
    end
    
catch
    single_distro = 1;
    %     disp('Parameters estimated for single peak Gaussian')
    paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
    paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
        paramEsts(2)];
end

% % %show the result
% % figure
% % % [~, bins] =
% % histogram(x,100);
% % % bins = -2.5:.5:7.5;
% % % h = bar(bins,histc(x,bins)/(length(x)*0.5),'histc');
% % % histogram(x,100)
% % % h.FaceColor = [0.9 0.9 0.9];
% % xgrid = linspace(1.1*min(x),1.1*max(x),200);
% % pdfgrid = pdf_normmixture(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),...
% %     paramEsts(4),paramEsts(5));
% % hold on
% % plot((paramEsts(3) - 2*paramEsts(5)),pdfgrid,'or')
% % plot((paramEsts(2) + 2*paramEsts(4)),pdfgrid,'*r')
% % plot(xgrid,pdfgrid,'-b')
% % hold off
% % xlabel('x')
% % ylabel('Probability Density')

end

function [peak,cc_min] = cc_max_find(cc)
%find the peak and zero-adjusted cc map

cc_min = cc - min(cc(:));%zero-adjust
% cc_filt = imgaussfilt3(cc_min); %filter to remove noise from peak value

[peak,~] = max(cc_min(:)); %get the index of the peak


end

function [ppe] = q_entropy(cc_min)
%compute entropy q-factor for a given cc map

[cc_hist,~] = histcounts(cc_min,150); %get histogram values

entropy = 0;
p = cc_hist/sum(cc_hist); %compute probablities
for i = 1:length(p)%compute entropy
    if p(i) == 0
        entropy = entropy+p(i);
    else
        entropy = entropy+p(i)*log(1/p(i));
    end
end

ppe = 1/entropy; %peak to cc (information) entropy
%min value -> 0 (worst case)


end

