function r=estcoh(x,n,p,prob,gam) 
%
% Author: Andrew Walden
% when used with fzero( )
% finds gamh=x  given n,p gam which gives cdf=prob
%
% n         no. complex degrees of freedom
% p         no. series
% x         defines P(gamh <x)
% prob      corresponding prob
% gam       actual mag sqd multiple coherence
%
%
r=prob-cdfcohfin(n,p,x,gam);
 