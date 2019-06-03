function [N, X, E, M, E_source, lambda, norm_spectra, norm_cross_spectra, coherences] = MvHawkesSimulation(v,A,beta,lambda0,varargin)
% HAWKES_Mv Simulates mutually exciting Hawkes process, and computes the
% theoretical spectrum
% 
% SYNOPSIS: [N, X, E] = hawkes_Mv(alpha,beta, mu)
%
% INPUT v - reversion level (p vector)
%       A - jump scale (a pxp matrix)
%       beta - exponential decay (a pxp matrix)
%       lambda0 - initial intensity (p vector)
% Optional:
%		Delta - bin width of summarisation
%       Nmax - number of events 
%       Tmax - length of time-interval 0->T
%       Jump - "fixed" or "exp" fixed or random jumps?
%       Shift - vector of shift parameters
%
% OUTPUT N - discrete counting process
%        X - binned process
%        E - set of events for process realisation
%        M - a set of marks/class correponding to event
%        lambda - list of lambdas updated after each event
%        norm_spectra - the theoretical spectra of each of the processes
%        norm_cross_spectra - the theoretical cross spectra
%        coherences - theoretical coherences        
%
% REMARKS - This is Algorithm 1 from 
% `Simulation and Calibration of a Fully Bayesian Marked Multidimensional 
% Hawkes Process with Dissimilar Decays'
% Adapted and improved version of the multivariate version of the Hawkes
% simulation algorthim in Dassios 2013.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input parsing
pa = inputParser;
defaultDelta = 1;
defaultNmax = 10000;
defaultTmax = 50;

addOptional(pa,'Delta',defaultDelta,@isnumeric);
addOptional(pa,'Nmax',defaultNmax,@isnumeric);
addOptional(pa,'Tmax',defaultTmax,@isnumeric);

parse(pa,varargin{:});  
Delta = pa.Results.Delta;   % Defualt is 1
Tmax = pa.Results.Tmax;
Nmax = pa.Results.Nmax;

% Start main routine
% Get dimension from length of vector
p = length(lambda0);    

% Rates are recorded marginally for stream
% lambdaTm=diag(lambda0);

lambdaTp(:,:,1) = ones(p);
lambdaTp(logical(eye(size(lambdaTp)))) = lambda0;

% Initialise necessary lists and arrays
j=2;    % Due to matlab indexing k=2 corresponds to first event
E=0;  % Init to zero time
M=0;   % Init marks
E_source = 0;   % Init log of source of time point

while(j<Nmax && E(end)<Tmax)
    % Generate an array of interarrival times, simulated for each 
    % process and due to each of the p processes (i.e. pxp array)
    
    % Simulate immigrant processes ~ Exp(v)
    recip_v = 1./v;
    % Being explicit with recip_v(:,1) for Python interpreter
    rng shuffle % ensures that a different seed is used each time for rand
    imm_times = -recip_v(:,1).*log(rand([p,1])); % px1 array of simulated immigrants
    
    % Matrix used to simulate the interarrival times 
    % D_arr = (beta.*log(ones(p)-rand(p))) ./ (lambdaTp(:,:,j-1)); % this will be negative
    D_arr = (beta.*log(rand(p))) ./ (lambdaTp(:,:,j-1)); 
    
    % Enforce condition
    cond_D_arr = ones(p) + D_arr;
    % cond_D_arr = ones(p) - D_arr;
    accepted_indices = cond_D_arr > 0; % Gives array of logical true/false
    
    % Provided condition holds, set array of times (log works element-wise)
    % Values not holding the condition are set to zero using failed_indices
    % interarrival_j = (-1./beta).*log(ones(p)-D_arr).*accepted_indices; %
    % seems an error in the paper, below ensures +ve interarrival times
    interarrival_j = (-1./beta).*log(ones(p)+D_arr).*accepted_indices;
    
    % Set zero values to infinity as these don't satisfy the condition
    interarrival_j(interarrival_j==0) = inf;
    % Find the shortest inter-arrival time from pxp array and immigrants
    interarrivals = [imm_times,interarrival_j]; % now p*(p+1)
    minimum_interarrival = min(min(interarrivals));
    
    % Assign the minimum value to the array of event times
    E = [E;E(j-1)+minimum_interarrival];
    
    % Get the indices of the minimum time and update array of marks
    % min_index_x: between 1 and p, denotes the process the time relates to
    % min_index_y - 1: denotes the source process of the new event time
    [min_index_x, min_index_y] = find(interarrivals == minimum_interarrival);
    
    M = [M;min_index_x];
    % Set the event sources: if 0, then an immigrant, if E_source(i)=M(i),
    % then self-excitation, otherwise, cross-exciation
    E_source = [E_source;min_index_y - 1]; 
    % Update the conditional intensity function matrix, add
    % A(:,min_index_x) to min_index_x th column of lambdaTp
    add_intensity = zeros(p);
    add_intensity(:,min_index_x) = A(:,min_index_x);
    lambdaTp(:,:,j) = lambdaTp(:,:,j-1).*exp(-beta.*(E(end)-E(end-1))) + add_intensity;
     
    % lambdaTm = (lambdaTp(:,:,j-1)-diag(v)).*exp(-beta.*(E(end)-E(end-1)))+diag(v);
    % lambdaTp(:,:,j) = lambdaTm + diag(A(:,min_index_x));
    j=j+1;  % Increase the index counter by 1
end

% Adding in option to add a ping process
% E = [E; Ping];
% E = sort(E);
lambda = lambdaTp;
% Convert to binned counts
K=floor(Tmax/Delta); 
N=zeros([p,K]);
Nt=zeros([1,p]);
for j=1:K
    for m=1:p
%     Xtemp = length(E(E>=(k-1)*Delta & E<=k*Delta));
        Xtemp(m) = sum(squeeze(M(E>=(j-1)*Delta & E<=j*Delta)==m));
        Nt(m)= Nt(m)+Xtemp(m);
        N(m,j) = Nt(m);
        X(m,j) = Xtemp(m);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the theoretical spectrum for Delta the bin width of summarisation
% of the data

m = ceil(Tmax/Delta); % m is the number of bins
mid_point = floor(m/2); % Equivalently: floor(T/(2*delta));

% Represents the Nyquist frequencies to plot the spectrum against  
x = linspace(0,1/(2*Delta),mid_point + 1);  

counter = 1;
theoretical_spectral_matrix = zeros(length(v),length(v),length(x)); % length(A) gives the dimension of the Mv HP
testing_lambda = diag((eye(length(A)) - (A./beta))\v);
for iter_s = x
    iter_s_val = 2*pi*sqrt(-1)*iter_s;
    theoretical_spectral_matrix(:,:,counter) = (eye(length(A)) - A./(beta + iter_s_val)) \ testing_lambda / (eye(length(A)) - transpose(A./(beta - iter_s_val)));
    counter = counter + 1;
end

pairs = nchoosek(1:length(v), 2);

norm_cross_spectra = zeros(length(x),size(pairs,1));
for cross_edge = 1:size(pairs,1)
    edge_pair = pairs(cross_edge,:);
    norm_cross_spectra(:,cross_edge) = abs(squeeze(theoretical_spectral_matrix(edge_pair(1),edge_pair(2),:)));
end

norm_spectra = zeros(length(x), length(v));
for self_edge = 1:length(v)
    norm_spectra(:,self_edge) = abs(squeeze(theoretical_spectral_matrix(self_edge,self_edge,:))); % this is real but has +/- 0i terms so abs makes real
end

coherences = zeros(length(x), size(pairs,1));
for coh_pairs = 1:size(pairs,1)
    edge_pair = pairs(coh_pairs,:);
    coherences(:,coh_pairs) = (norm_cross_spectra(:,coh_pairs)).^2./(norm_spectra(:,edge_pair(1)).*norm_spectra(:,edge_pair(2)));
end

end
