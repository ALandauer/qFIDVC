function himage = imshow3(I,varargin)
%IMSHOW3 Display 3D image stack.
%   IMSHOW3(I) generates a figure for viewing slices of a 3D image. A
%   slider at the bottom of the figure allows users to switch between
%   frames.
%
%   IMSHOW3(___,Name,Value) uses additional parameter name-value 
%   pairs. Valid parameters include:
%
%       'Parent'        Handle of the parent figure to use. Can be a
%                       numeric scalar or the string 'new'.
%
%                       Default: sum(mfilename+0)
%
%       'Scale'         Scalar indicating how to resize the image. Passed
%                       as an input to IMRESIZE.
%
%                       Default: 1 (no resizing)
%
%   Example
%	-------
%   View a 3D MRI sample:
%
%       filename = 'mri.tif';
%       info = imfinfo(filename);
%       for ii=1:length(info)
%           I(:,:,ii) = imread(filename,'Index',ii,'Info',info);
%       end
%       imshow3(I);
%
%   See also IMSHOW.

% Copyright 2016 Matthew R. Eicholtz

% Default parameter values
default = struct(...
    'Parent',sum(mfilename+0),...
    'Scale',1);

% Parse inputs
[h,scale] = parseinputs(default,varargin{:});
I = imresize(I,scale,'nearest');
[row,col,layers] = size(I);

% Setup figure
if ishandle(h); close(h); end %close the figure if it exists
h = figure(h);
h.Visible = 'on';
h.Color = [1 1 1];
h.MenuBar = 'none';
h.ToolBar = 'none';
h.Resize = 'off';

marginx = 12;
marginheader = 28;
marginfooter = 12;
if layers>1; marginfooter = 40; end
pos = h.Position; %initial figure position
pos(2) = min(pos(2),200);
h.Position = [pos(1:2), col+2*marginx, row+marginheader+marginfooter];
axes('Units','pixels','Position',[marginx,marginfooter,col,row]);
if layers>1
    uicontrol(...
        'Style','slider',...
        'Callback',@updatefigure,...
        'Max',layers,...
        'Min',1,...
        'Units','pixels',...
        'Position',[0 0 h.Position(3)+2 28],...
        'SliderStep',[1/(layers-1) 1/(layers-1)],...
        'Tag','slider',...
        'Value',1);
end

% Set application-defined data
setappdata(h,'I',I);

% Render the first image
colormap gray
himage = image(I(:,:,1),'Tag','image');
himage.CDataMapping = 'scaled';

hold on;
text(col-4,-24,sprintf('%d/%d',1,layers),'FontSize',12,'FontWeight','bold','Color','k',...
    'Interpreter','none','HorizontalAlignment','right','VerticalAlignment','top');
hold off;
axis off;

end

%% Helper functions
function varargout = parseinputs(default,varargin)
%PARSEINPUTS Custom input parsing function.
    p = inputParser;
    
    p.addParameter('parent',default.Parent,...
        @(x) isnumeric(x)|strcmp(x,'new'));
    p.addParameter('scale',default.Scale,...
        @(x) validateattributes(x,{'numeric'},{'scalar','nonempty','nonsparse','>',0,'finite'}));
    
    p.parse(varargin{:});
    
    [h,scale] = struct2vars(p.Results);
    
    if strcmp(h,'new')
        h = length(findobj('type','figure'))+1;
    end
    
    varargout = {h,scale};
end
function updatefigure(obj,~)
%UPDATEFIGURE Update the figure when the slider value changes.
    ind = round(obj.Value);
    obj.Value = ind;
    
    h = obj.Parent;
    I = getappdata(h,'I');
    
    himage = findobj(h,'Tag','image');
    htext = findobj(h,'Type','text');
    
    himage.CData = I(:,:,ind);
    htext.String = sprintf('%d/%d',ind,size(I,3));
end

