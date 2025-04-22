function [jmptout,sttout] = simCTMC(Q,T,nsim,instt)
%
%[jmptout,sttout] = simCTMC(Q,T,nsim,[instt])
% 
% This function simulates jump times and associated states for a
% continuous-time Markov chain with associated Q-matrix Q (see Norris
% (1997), Markov Chains, Cambridge University Press).
%
% Inputs:
%   Q: square Q-matrix for Markov chain
%   T: simulation over time range [0,T); set T=Inf if you want to simulate
%       until an absorbing state is reached (if T=Inf but there is no
%       absorbing state, or if there is, but it's possible that the
%       absorbing state will never be reached, then setting T=Inf returns
%       an error).
%   nsim: number of simulations to run (only used if instt is not passed in)
%   instt: optional vector of initial states; if passed in, nsim = size of
%       instt; otherwise, nsim draws are made from the stationary
%       distribution of the Markov chain (if there are multiple stationary
%       distributions, an error is returned).
%
% Outputs:
%   jmptout: if nsim=1, row vector containing the jump times, where the
%       first element is always 0, and the last element is the last jump
%       prior to period T; if nsim>1, jmptout is a cell vector whose j-th
%       element is the row vector of jump times for the j-th simulation.
%   sttout: if nsim=1, row vector whose k-th entry contains the state of 
%       the system immediately following the jump at time jmptout(k); the
%       first element is always the initial state; if nsim>1, sttout is a
%       cell vector whose j-th element is the row vector of states for the
%       j-th simulation.
%
% To simulate a time series for the states, first obtain jmptout and sttout
% using this function, then pass the results to the function simCTMC_t.

%% initial checks

dQ = diag(Q);

if any(abs(sum(Q,2))>1e-12*abs(dQ))
    error('Rows of Q must sum to 0.')
end

nQ = size(Q,1);
if ~ismatrix(Q) || (size(Q,2) ~= nQ)
    error('Q must be a square matrix.')
end

if ~isscalar(T) || (T <= 0)
    error('T must be a positive scalar.')
end

%% Get jump transition matrix Pi

Pi = Q - diag(dQ);
Pi = Pi + diag(dQ==0);
Pi = Pi./sum(Pi,2);

cuPi = cumsum(Pi,2); % cumulate rows of Pi; to be used to get random draws

%% Exponential distribution parameters for jumps

lam = abs(dQ);
mu = 1./lam;

%% Compute stationary distribution

stdst = null(Q')';
stdst = stdst./sum(stdst,2);

if isinf(T)
    if any(sum(stdst>0,2) > 1)
        error(['T=Inf specified but system has a stationary distribution '...
            'that doesn''t have an absorbing state.'])
    end
end

%% Draw initial state from stationary distribution if necessary

if exist('instt','var') % if initial state is passed in
    instt = instt(:);       % convert to column vector
    % check that initial states are valid (integers between 1 and nQ)
    if any(mod(instt,1) ~= 0) || (max(instt) > nQ) || (min(instt) < 1)
        error(['Initial state(s) not valid. Must be an integer between 1 '...
            'and the dimension of Q.'])
    end
    nsim = numel(instt);    % number of simluations is number of initial states
    curstt = instt;         % vector that holds current state value
else                    % otherwise draw from stationary distribution
    
    if mod(nsim,1) ~= 0 || (nsim <= 0)
        error('nsim must be a positive integer.')
    end
    
    if size(stdst,1) > 1
        error(['Multiple stationary distributions. Please specify initial '...
            'states using ''instt''.'])
    end
    % draw from stationary distribution
    custdst = cumsum(stdst);
    p = rand(nsim,1);
    
    % the following finds the index of the first value of custdst that
    % exceeds the random draw p; this is the index of the inital state
    [~,curstt] = max(p<=custdst,[],2);
end

%% Allocate memory and initialize arrays

mxaddcols = 1000;

mxlam = max(lam);
njin = min(ceil(T*mxlam),mxaddcols); % guess of upper limit for total number of jumps 
                % each simulation will make; just used for pre-allocating
                % memory

% matrix to hold times of jumps
jmpts = zeros(nsim,njin);

% matrix to hold states
stts = zeros(nsim,njin);
stts(:,1) = curstt;

%% Simulation

curt = zeros(nsim,1); % vector to keep track of current states
fl = true(nsim,1);  % flag where true = need to keep simulating

j = 1;  % loop counter
nj = njin;  % variable to track current # of columns of jmpts, stts

while any(fl)   % as long as at least one simulation hasn't reached T periods
    
    j = j+1;
    
    % add more columns to jmpts and stts if necessary
    if j > nj
        njadd = min(ceil((T-min(curt))*mxlam),mxaddcols);
        jmpts = [jmpts, zeros(nsim,njadd)];
        stts = [stts, zeros(nsim,njadd)];
        nj = nj + njadd;
    end
    
    % current states and times for simulations we still need draws for
    curstt2 = curstt(fl); 
    curt2 = curt(fl);
    ndr = nnz(fl);
    
    % draw exponential random variables for periods until next jump
    jmp = -mu(curstt2) .* log(rand(ndr,1));
    
    % jmpfl = true if we need to draw a new state, which is the case if the
    % next jump happens before T periods; otherwise the current state
    % remains until the end of the simulation period
    jmpfl = (curt2+jmp < T);
    njmpfl = nnz(jmpfl);
    
    % draw new states where necessary
    p = rand(njmpfl,1);
    cuPijmp = cuPi(curstt2(jmpfl),:);
    % this uses the same trick as drawing from stationary distribution above
    [~,newstt] = max(p<=cuPijmp,[],2);
    
    % assign new states to appropriate elements of curstt2, and set all
    % other elements to 0
    curstt2(jmpfl) = newstt;
    curstt2(~jmpfl) = 0;
    
    curstt(fl) = curstt2; % assign to curstt for use on the next iteration
    stts(fl,j) = curstt2; % save state
    
    
    curt2 = curt2+jmp; % assign time of next jumps to curt2
    curt(fl) = curt2;   % assign to curt for use on next iteration
    jmpts(fl,j) = curt2; % save jump times
    
    fl = (curt<T);  % update flag for whether we have enough draws
end

%% Prepare output

if nsim == 1
    lstjmp = find(stts==0,1)-1;
    jmptout = jmpts(:,1:lstjmp);
    sttout = stts(:,1:lstjmp);
else
    jmptout = cell(nsim,1);
    sttout = cell(nsim,1);
    for j = 1:nsim
        jmptsj = jmpts(j,:);
        sttsj = stts(j,:);
        lstjmp = find(sttsj==0,1)-1;
        jmptout{j} = jmptsj(1,1:lstjmp);
        sttout{j} = sttsj(1,1:lstjmp);
    end
end

