%Another function the loops through a series of 3d synth images and grabs a
%centerline slice.
%A Landauer, 3.7.17

folder = 'high_density';
ext = '.mat';
files = dir(strcat('./',folder,'/*',ext));

s = length(files);

X = [1,512];
Y = [1,512];
Z = 8;

cnt = 0;
for ii = s:-1:1
    cnt = cnt+1;
    
    load(files(ii).name);
    
    I0_series{ii} = uint8(I0(X(1):X(2),Y(1):Y(2),Z));
    
    cur_img = I0_series{ii};
    
    imwrite(cur_img,['./',folder,'_slices/',num2str(cnt+999),'_',files(ii).name(1:end-4)],...
        'tif','Compression','none');
    
end

save(['img_stack_2d_',files(1).name(1:20)],'I0_series')
    
    
    
    
    
    
    
