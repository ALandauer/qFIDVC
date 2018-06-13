%This function computes the convolution of the bead intensity and the point
%spread function of pinhole optics

clear all
close all
clc

%% CONVOLUTION
stretch = 0.5; %uniaxial stretch level

% for stretch = 1:-0.15:0.25

% set run parameters
L = 532*10^-9; %incident wavelength

% d^2/l
I_0 = 255; %reference intensity

r = 0.6; %bead radius

dz = 1.2; %spacing of image plane from pinhole

d = sqrt(0.61*L*dz);   %pinhole radius, based on Young, 1971 "pinhole optics"

sigma = [0.5,0.5,0.5]; %3d gaussian distros

%compute semi-axis lengths
a(1) = r*stretch;
a(2) = r/sqrt(stretch);
a(3) = r/sqrt(stretch);

%compute the convolution
st = -5.0;
sp = 0.2;
ed = 5.0;

pool = gcp('nocreate');

if isempty(pool)
    curCluster = parcluster('local');
    curCluster.NumWorkers = 16;
    saveProfile(curCluster);
    pool = parpool(16);
end

k = 500;

t0 = tic;
t1 = tic;
v = st:sp:ed;
I_field = zeros(4,length(st:sp:ed)*length(st:sp:ed)*length(st:sp:ed));
I_mat = zeros(length(v),length(v),length(v));

cnt = 0;
for ii = 1:length(v) %st:sp:ed
    m = v(ii);
    fprintf('\n %i of %i and loop time %f',cnt,length(I_field),toc(t1));
    t1 = tic;
    for jj = 1:length(v) %st:sp:ed
        n = v(jj);
%         t2 = tic;
        parfor kk = 1:length(v) %st:sp:ed
            
            p = v(kk)+dz;
%             FofI = @(x,y,z) I_0*exp(-((x.^2)/(2*sigma(1)^2) + y.^2/(2*sigma(2)^2) +z.^2./(2*sigma(3)^2)))./...
%      (1+exp(-2*k*(1-((m-x).^2/a(1)^2)-((n-y).^2/a(2)^2)-((p-z).^2/a(3)^2))));
%            FofI = @(x,y,z) bead_img_conv(x,y,z,m,n,p,a,d,L,I_0,sigma);
            
            FofI = @(x,y,z) I_0*exp(-((x.^2)/(2*sigma(1)^2) + y.^2/(2*sigma(2)^2) +z.^2./(2*sigma(3)^2)))./...
    (1+exp(-2*k*(1-((m-x).^2/a(1)^2)-((n-y).^2/a(2)^2)-((p-z).^2/a(3)^2))));

            I_cur = integral3(FofI,-1.5*sigma(1),1.5*sigma(1),-1.5*sigma(2),1.5*sigma(2),...
                -(1.5*sigma(3)+dz),(1.5*sigma(3)+dz));%,'AbsTol',5e-4,'RelTol',5e-3);

%             I_cur = triplequad(FofI,-1.5*sigma(1),1.5*sigma(1),-1.5*sigma(2),1.5*sigma(2),...
%                 -(1.5*sigma(3)+dz),(1.5*sigma(3)+dz));
            
%             for qq = 1:8
%                 I_cur_(qq) = integral3(FofI,o{qq});
%             end
%             I_cur = sum(I_cur_);

%             I_cur = triplequad(FofI,st,ed,st,ed,st+5,ed+5);
            
            
            bead(ii,jj,kk) = 1./(1+exp(-2*k*(1-((m).^2/a(1)^2)-((n).^2/a(2)^2)-((p).^2/a(3)^2))));

            I_field(:,cnt+kk) = [m,n,p,I_cur];
            
            I_mat(ii,jj,kk) = I_cur;
            
        end
%         toc(t2)
        cnt = cnt+length(v);
    end
end
fprintf('\nTotal integration time: %f\n',toc(t0))

clk = clock;
save(['bead_image_output_lambda_',num2str(stretch),'_',num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
    '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'])

% end
