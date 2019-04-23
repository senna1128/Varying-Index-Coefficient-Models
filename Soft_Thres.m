% This function is soft-thresholding function for Lasso 
% Input
  % x: matrix
  % kappa: threshold
% Output
  % x: soft(x,kappa)

function [ x ] = Soft_Thres(x, kappa)

x(abs(x)<=kappa) = 0;
x(x<-kappa) = x(x<-kappa) + kappa;
x(x>kappa) = x(x>kappa) - kappa;

end

