%% Example run file for the IDVC - LDTFM
% The two images that are going to be run are 'vol00.mat' and 'vol01.mat.'
% the deformation are defined as a four-pole Gaussian prescribed
% displacement field.  See below for distriub
%
% Central z (x_3) plane. Number denotes width of gaussian in voxels.  See
% section 3.1 and 3.2 in Bar-Kochba et al. (2014)
% --------------------------
% |                        |
% |     32           64    |
% |                        |
% |                        |
% |                        |
% |     96           128   |
% |                        |
% |                        |
% --------------------------
%
% VARIABLES OPTIONS
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
%          Must be 32,64,96, or 128 voxels and a three column
%          array (one for each dimenision) or scalar (equal for all
%          dimensions).
%   
%   sSize: interrogation window (subset) size minimum value.
%
%   runMode: string that defines the method of running IDVC. Options:
%             cumulative (time0 -> time1, time0 -> time2, ...)
%             (Allowable inputs: 'c','cum','cumulative')
%             or
%             incremental (time0 -> time1, time1 -> time2, ...)
%             (Allowable inputs: 'i','inc','incremental')
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u:  displacement field vector calculated from FIDVC. Format: cell array,
%      which is a 3D vector (components in x,y,z)  per each time point
%      (units are in voxels)
%         u{time}{1} = displacement in x-direction at t=time of size MxNxP
%         u{time}{2} = displacement in y-direction at t=time of size MxNxP
%         u{time}{3} = displacement in z-direction at t=time of size MxNxP
%   cc: peak values of the cross-correlation for each interrogation
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

clear; close all;
%%
sSize = [64 64 64];
sSizeMin = 32;
runMode = 'i';
filename = 'vol_series_e1*.mat';
% filename = 'Crop*.mat';
% filename = 'test*.mat';

% Estimate displacements via qIDVC
[u, ~, dm, m] = funIDVC(filename, sSize, sSizeMin, runMode);

save('results_qDVC_uniform_e1_inc.mat','u','sSize','sSizeMin','dm', 'm', 'runMode','-v7.3');
% save('resultsFIDVCnewsSize16.mat','u','cc','dm', 'm');
% save('data_20180510_1449','-v7.3')

%%
% figure
% subplot(1,3,1)
% hist(u{1}{1}(:),100)
% subplot(1,3,2)
% hist(u{1}{2}(:),100)
% subplot(1,3,3)
% hist(u{1}{3}(:),100)




