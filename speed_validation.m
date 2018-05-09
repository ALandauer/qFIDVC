clear

sz = 1024;

img = rand(sz,sz,sz/4);

[x1,x2,x3] = ndgrid(1:sz,1:sz,1:sz/4);

[subpx_pts1,subpx_pts2,subpx_pts3] = ndgrid(1:0.98:sz,1:0.93:sz,1:0.79:sz/4);

img_interp = zeros(size(img));
t2 = tic;
img_interp = mirt3D_mexinterp(img,subpx_pts1,subpx_pts2,subpx_pts3);
toc(t2);

img_interp = zeros(size(img));
t1 = tic;
img_interp = interpn(x1,x2,x3,img,subpx_pts1,subpx_pts2,subpx_pts3,'linear');
toc(t1);




