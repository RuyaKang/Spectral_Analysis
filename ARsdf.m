function S = ARsdf(p, phi, f)
% compute periodogram for a range of frequency
%     input:
%        p = order of the process
%        phi = column vector of parameters
%        f = column vector of frequencies
%     output: 
%        vector of sdf
vec = 1-exp(-j*2*pi*f.*(1:p))*phi;
S = 1./(vec.*conj(vec));
end