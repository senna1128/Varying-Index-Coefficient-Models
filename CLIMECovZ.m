%% This function estimate the precision matrix of Z, i.e. SNPs

function [OmegaHat] = CLIMECovZ(Z, n, d2, gam)

SigmaHat = Z'*Z/n;

cvx_begin
   variable OmegaHat(d2,d2)
   minimize( sum(sum(abs(OmegaHat))) )
   subject to
      max(max(abs(SigmaHat*OmegaHat - eye(d2)))) <= gam
cvx_end

end

