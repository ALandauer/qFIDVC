function H = heavisideApprox(x,k)
%analytical approx to the Heaviside step function, larger K yeilds a better
%approximation.

H = 1./(1+exp(-2*k*x));

end