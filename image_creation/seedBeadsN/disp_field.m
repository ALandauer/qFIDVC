function [u] = disp_field(x_,y_,z_,step,sizeI)
%This function takes in the curent location, and returns the cartisean
%displacement at that location for the Lame problem.

nu = 0.3;
E = 1e6;

%change the coord system
x_ = x_ - sizeI(1)/2;
y_ = y_ - sizeI(2)/2;
z_ = z_ - sizeI(3)/2;

a = 25;
% % r = sqrt(x.^2 + y.^2);
% % t = atan2(y,x);
% % 
% % Sinf = step*50000;

eT = zeros(3);
eT(1,1) = 3*step;
% eT(1,2) = step/2;
% eT(2,1) = step/2;
eT(2,2) = step;
eT(3,3) = step;

%
% %the line-load in a half-space.
% u_r = -(2*(1 - nu^2)/(pi*E))*P*cos(t).*log(r)-(1+nu)*(1-2*nu)/(pi*E)*P.*t.*sin(t);
% u_t = (2*(1 - nu^2)/(pi*E))*P.*sin(t).*log(r)+(1+nu)/(pi*E)*P*sin(t)-...
%     2*(1+nu)*(1-2*nu)/(pi*E)*P.*t.*cos(t);
%
% u_x = u_r*cos(u_t);
% u_y = u_r*sin(u_t);
%
% u = [u_x,u_y];

%the lame problem
% % % % if r >= a
% % % %     Srr_int = Sinf/2*(r*t-2*a*t+(r+4*a-(9*a^4/r^3)+(4*a^2/r))*sin(t)*cos(t));
% % % %     Stt_int = Sinf/2*(r*t-a^2/r-(r-(9*a^4/r^3)-10*a)*sin(t)*cos(t));
% % % %     Srt_int = Sinf/4*(r+(9*a^4/r^3)-a^2/r-8*a)*cos(2*t);
% % % %     
% % % %     u_r = 1/E*(Srr_int-nu*Stt_int);
% % % %     u_t = 1/E*(Stt_int-nu*Srr_int);
% % % %     u_z = 1/E*(-nu*(Srr_int+Stt_int));
% % % % %     urt = 1/E*(1+nu)*Srt_int;
% % % % else
% % % %     u_r = 0;
% % % %     u_t = 0;
% % % %     u_z = 0;
% % % % %     urt = 0;
% % % % end
% % % % 
% % % % u_x = u_r*cos(u_t);
% % % % u_y = u_r*sin(u_t);
% % % % 
% % % % u = [u_x,u_y,u_z];

%eshelby soln for a sphere with transform strain eT
x = [x_,y_,z_];
R = sqrt(x(1).^2 + x(2).^2 + x(3).^2);
pT = (E/(1+nu))*(eT + nu*(eT(1,1)+eT(2,2)+eT(3,3))*eye(3)/(1-2*nu));

for i = 1:3
   u(i) = ((1+nu)*(a^3))/(2*(1-nu)*E) * ...
       (((2*pT(i,:)*x' + (pT(1,1)+pT(2,2)+pT(3,3))*x(i))/(15*R^5))*(3*a^2-5*R^2) + ...
       (x*pT*x'*x(i)/R^7)*(R^2 - a^2) + (4*(1-nu)*pT(i,:)*x')/(3*R^3));
end

end