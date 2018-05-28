%Another function that loops through a series of 3d synth images and grabs a
%centerline slice.
%A Landauer, 3.7.17

clearvars

folder = 'high_density2';
ext = '.mat';
files = dir(strcat('./',folder,'/*',ext));

s = length(files);

X = [1,512];
Y = [1,512];
Z = 1:16;

%Make a new output folder if none exists
if exist(['./',folder,'_slices/'],'dir') ~= 7
    mkdir(['./',folder,'_slices/']);
end

cnt = 0;
for ii = s:-1:1
    cnt = cnt+1
    
    load(['./',folder,'/',files(ii).name]);
    
    I0_series{ii} = sum(I0(X(1):X(2),Y(1):Y(2),Z),3);
    I0_series{ii} = uint8(255*I0_series{ii}/max(max(I0_series{ii})));
    
    cur_img = I0_series{ii};
    
    imwrite(cur_img,['./',folder,'_slices/',num2str(cnt+999),'_',files(ii).name(1:end-4)],...
        'tif','Compression','none');
    
end

save(['img_stack_2d_',files(1).name(1:20)],'I0_series')
    
    
    
    
    
    
    
