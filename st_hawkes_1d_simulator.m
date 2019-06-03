function [proc1_freq, E, lambda, theoretical_spectrum, S11] = st_hawkes_1d_simulator(v,alpha,beta,lambda0,varargin)
% HAWKES_1D Simulates self exciting Hawkes process with exponential
% excitation function
% 
%
% INPUT v - reversion level
%       alpha - jump scale
%       beta - exponential decay
%       lambda0 - initial intensity (if negative draw from stationary)
% Optional:
%		Delta - bin width of summarisation
%       Nmax - number of events 
%       Tmax - length of time-interval 0->T
%       Jump - "fixed" or "exp" fixed or random jumps?
%
% OUTPUT proc1_freq - binned process
%        E - set of events for process realisation
%        lambda - list of lambdas updated after each event
%        theoretical_spectrum - the spectrum of the simulated HP
%        S11 - periodogram estimate of the spectrum        
%
% REMARKS - This is an implementation of the univariate
%       simulation algoirhtm (Alg 1) of "Exact simulation of Hawkes process
%       with exponentially decaying intensity", Dassios (2013)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing
pa = inputParser;
defaultDelta = 1;
defaultNmax = 10000;
defaultTmax = 1000;

addOptional(pa,'Delta',defaultDelta,@isnumeric);
addOptional(pa,'Nmax',defaultNmax,@isnumeric);
addOptional(pa,'Tmax',defaultTmax,@isnumeric);

parse(pa,varargin{:});  
Delta = pa.Results.Delta;   % Defualt is 1
Tmax = pa.Results.Tmax;
Nmax = pa.Results.Nmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(lambda0 < 0)
    % Sample from stationary distribution
    % p7 Dassios (2013) with a equiv v, delta equiv beta, beta equiv alpha
    lambdaTp(1) = v + gamrnd((v/beta),(beta*alpha-1)/beta);
    lambdaTm(1) = lambdaTp(1);
else
    lambdaTp(1)=lambda0;
    lambdaTm(1)=lambda0;
end

k=2; % Due to matlab indexing k=2 corresponds to first event
E=[0];  % Init to zero time

while(k<Nmax && E(end)<Tmax)
    % Draw random decision variable
    D = 1 + ( beta*log(rand(1)) / (lambdaTp(k-1)-v) );
    S1 = -(1/beta)*log(D); 
    S2 = -(1/v)*log(rand(1));
    if(D>0)
       S = min(S1,S2); 
    elseif(D<0)
       S = S2; 
    end
    E = [E;E(k-1)+S];  % add a new event
    
    lambdaTm(k) = (lambdaTp(k-1)-v)*exp(-beta*(E(end)-E(end-1)))+v;
    lambdaTp(k) = lambdaTm(k) + alpha;
    
    k=k+1;  % Update count
   
end

E = sort(E);

lambda = lambdaTp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated spectrum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of events per bin for series of event times 
proc1_freq = transpose(histcounts(E,0:Delta:Tmax));
proc1_freq = proc1_freq - mean(proc1_freq); % Deduct the mean of the counts
% Compute periodogram estimator 
N = length(proc1_freq); 
H = 1/sqrt(N) * ones(N, 1);
% S11 stores the estimated periodogram estimate (note that it is symmetric)
S11 = abs(fftshift(fft(H.*proc1_freq))).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theoretical values of the spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stationary lambda, as given in Hawkes' original paper (1971a)
stat_lambda = lambda0*beta/(beta-alpha);
% Theoretical spectrum, from Hawkes' paper. Note that instead of angular
% frequence (omega), frequency is used (omega=2*pi*f). 
x = linspace(0,1/(2*Delta),N/2 + 1);
theoretical_spectrum = transpose(stat_lambda*((beta.^2 + (x*2*pi).^2)./((beta - alpha)^2 + (x*2*pi).^2)));

end


