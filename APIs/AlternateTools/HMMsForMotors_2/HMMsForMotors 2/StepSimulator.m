function [X,Y]=StepSimulator(M,N)
% function [X,Y]=StepSimulator(M,N)
% Discrete-time step simulator.
% Given the number of time points N and the model M, simulate a staircase
% signal created according to M.P0 and M.C
% Output:  X noiseless signal
%          Y signal with noise of s.d. M.Sigma.
% M.DutyCycle is a number between 0 and 1 representing the integration time
% of the camera relative to the sampling interval.  Zero means there are
% no intermediate steps; 0.5 means half of the time there will be
% intermediates.
%
% Modified 19 Apr 07 to include M.Dutycycle instead of having this be an
% argument. -fs
% Modified 12 Sep 07 to allow negative steps. -fs


% Here is how one would initialize the random number generators to a
% standard state:
% randn('state',0);
% rand('state',0);

[nu ns ns1]=size(M.C);
X=zeros(N,1);

% Construct the intermediate step switch
partsw=[M.DutyCycle 1-M.DutyCycle]';

% Pick the starting state
    v=sum(M.P0(:,:));  % sum over the first index, which is the step size.
    i=PickJ(v);
    % Find the corresponding step
    s=PickJ(M.P0(:,i)/v(i))-1;
    if s > nu/2
        s = s-nu;  % allow a negative value
    end;
    X(1)=s;

    for t=2:N
    % Find the new molecular state
    v=sum(M.C(:,i,:));  % sum over the first index, which is the step size.
    j=PickJ(v);

    % find the step size
    delta=PickJ(M.C(:,i,j)/v(j))-1;
    if delta > nu/2
        delta = delta-nu;  % negative steps 
    end;
    % Allow for partial steps
    factor=1;
    if PickJ(partsw)==1  % See if there is to be a partial step.
        factor=rand;
    end;
    X(t)=s+round(delta*factor);  % This point may show a partial step
    s=s+delta;  % The next step will be full size.
    i=j;
end;

Y=X+M.Sigma*randn(N,1);  % M.Sigma is already in y-units
Y=floor(Y);


function j=PickJ(v)
% Given a vector v, find a random index j that is weighted by the elements
% of the vector.
ns=numel(v);
q=rand;
s=v(1);
j=1;
while (s<q)&&(j<ns)
    j=j+1;
    s=s+v(j);
end;

