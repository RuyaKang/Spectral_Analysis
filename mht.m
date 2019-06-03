function RHH = mht(gamma, K, p, alpha)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

L = length(gamma);

% sort the estimates in descending order
R = sort(gamma, 'descend');

% test trough the sorted estimates, stop when a single one pass the test
l = 1;
reject = 0;

while l <= L
    C = fzero(@(x) estcoh(x,K,p,(1-alpha)^(1/(L-l+1)),0), [0.000001, 0.999999]);
    if R(l) < C
        break
    else
        reject = reject + 1;
    end
    
    l = l+1;

end

RHH = reject/L;

end

