function [X,Y]=StepSimulatorC(M,N)
% Given the number of time points N and the model M, simulate a staircase
% signal created according to M.P0 and M.C
%
% This 'C' version simulates a continuous-time Markov process.
% Exponentially-distributed dwell times are computed, and the
% continuous-time function is sampled according to the duty cycle to
% produce the output Y.  X is sampled effectively with duty cycle = 0.
%
% Output:  X noiseless signal
%          Y signal with noise of s.d. M.Sigma.
%
% M.DutyCycle is a number between 0 and 1 representing the integration time
% of the camera relative to the sampling interval.  Zero means there are
% no intermediate steps; 0.5 means half of the time there will be
% intermediates.
%
% 1 Apr 08 fs

[nu ns ns1]=size(M.C);
X=zeros(N,1);  % Noiseless, shows instantaneous transitions
Y=X;           % Starts out noiseless, but includes intermediate points.

% Construct the intermediate step switch
dc=M.DutyCycle;    % e.g. 0.9
dc1=1-M.DutyCycle; % e.g. 0.1

c0=M.C;
% Make the target state matrix and the a_ii vector.
if ns==1  % special case: the no-transition probability is C(1,1,1).
    a0=c0(1);
    c0(1)=0;
else
    for i=1:ns
        a0(i)=sum(c0(:,i,i));
        c0(:,i,i)=0;
    end;
end;

% Pick the starting state according to P0
v=sum(M.P0(:,:));  % sum over the first index, which is the step size.
i=PickJ(v);

% % Find the amplitude of the initial value
% s=PickJ(M.P0(:,i))-1;
% if s>nu/2  % handle negative values.
%     s=s-nu;
% end;

% No, force the starting amplitude to zero
s=0;

t0=1;   % starting time
t=t0;
X(t:N)=s;
Y(t:N)=s;

while t0<N-1
    t=t0-1/(1-a0(i))*log(rand);   % increment by an exponentially-distributed dwell time
    t=min(t,N);

    it=floor(t);  % Integer and fractional parts.
    ft=t-it;

    % Find the new molecular state
    v=sum(c0(:,i,:));  % sum over the first index, which is the step size.
    j=PickJ(v);

    % find the step size
    s=PickJ(c0(:,i,j))-1;
    if s>nu/2
        s=s-nu;
    end;

    % Allow for partial steps.  ps is the partial step fraction.
    if ft<dc1  % When the duty cycle is low, take no step most of the time.
        ps=0;
    elseif dc>0 % If there is some integration
        ps=(ft-dc1)/dc;  % The rest of the time, take a value between 0 and 1
    else
        ps=1;   % Case that shouldn't occur.
    end;

    Y(it)=Y(it)+ps*s;  % Increment the present point.  Y includes the partial step
    if it < N   % If there are any more points left
        Y(it+1:N)=Y(it+1:N)+s;
        X(it+1:N)=X(it+1:N)+s; % X skips the partial step.
    end;
    i=j;
    t0=t;
end;

Y=Y+M.Sigma*randn(N,1);  % M.Sigma is already in y-units
Y=round(Y);


% meant=sumd/sumn
% nps=nps/sumn

function j=PickJ(v)
% Given a normalized vector v (that is, sum(v)==1) find a random index j that is weighted by the elements
% of the vector.

q=rand*sum(v);
j=find((cumsum(v) >= q),1,'first');

% ns=numel(v);
% q=rand*sum(v);
% s=v(1);
% j=1;
% while (s<q)&&(j<ns)
%     j=j+1;
%     s=s+v(j);
% end;
