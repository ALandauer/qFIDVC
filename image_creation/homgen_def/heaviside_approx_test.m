
ccc

particle_radius = 5;

a(1) = particle_radius*1.2;
a(2) = particle_radius/sqrt(1.2);
a(3) = particle_radius/sqrt(1.2);

[xd,yd,zd] = meshgrid(...
    linspace(-1.2*a(1),1.2*a(1),5*particle_radius),...
    linspace(-1.2*a(2),1.2*a(2),5*particle_radius),...
    linspace(-1.2*a(3),1.2*a(3),5*particle_radius));
% disp_pts = [xd,yd,zd];

%define ellipsoid distances from center, the surface is at el_dist == 1
%since x^2/a^2 + y^2/b^2 + y^2/c^2 = 1
el_dist = sqrt(xd.^2./(a(1).^2) + yd.^2./(a(2).^2) + zd.^2./(a(3).^2));

%
%     for

k = 10000;

ell_tophat = heavisideApprox(-(el_dist-1),k);


%%exp(-(x^2 + y^2)/(0.3528((l/d)(sqrt(x^2+y^2+z^2)))^2))
%%(1+exp(-2k(1-((m-x)^2/a^2)-((n-y)^2/b^2)-((p-z)^2/c^2)))

l = 532*10^-9;
d = 0.0001;
% d^2/l
k = 100000; %large number
I_0 = 255; %reference intensity

r = 0.0001;
stretch = 1;

a(1) = r*1.2;
a(2) = r/sqrt(1.2);
a(3) = r/sqrt(1.2);

I_func = @(x,y,z) I_0*exp(-(x.^2 + y.^2)./(0.3528*((l/d).*(sqrt(x.^2+y.^2+z.^2))).^2))./...
    (1+exp(-2*k*(1-((m-x).^2/a(1)^2)-((n-y).^2/a(2)^2)-((p-z).^2/a(3)^2))));



%%%  http://www.wolframalpha.com/input/?i=exp(-(x%5E2+%2B+y%5E2)%2F(0.3528((l%2Fd)(sqrt(x%5E2%2By%5E2%2Bz%5E2)))%5E2))%2F(1%2Bexp(-2k(1-((m-x)%5E2%2Fa%5E2)-((n-y)%5E2%2Fb%5E2)-((p-z)%5E2%2Fc%5E2))))

% 
% sizeOut = size(ell_tophat);
% hFigRotated = figure;
% hAxRotated  = axes;
% slice(double(ell_tophat),sizeOut(2)/2,sizeOut(1)/2,sizeOut(3)/2);
% grid on, shading interp, colormap gray, axis image
% 
% 
% sizeOut = size(el_dist);
% hFigRotated = figure;
% hAxRotated  = axes;
% slice(double(el_dist),sizeOut(2)/2,sizeOut(1)/2,sizeOut(3)/2);
% grid on, shading interp, colormap gray, axis image


