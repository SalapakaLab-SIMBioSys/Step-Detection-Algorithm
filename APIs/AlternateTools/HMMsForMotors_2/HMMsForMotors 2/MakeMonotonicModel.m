function M=MakeMonotonicModel(nu, yQuantum, noiseSigma, transProb, stepSizes, stepSigma)
% function M=MakeMonotonicModel(nu, yQuantum, noiseSigma, transProb, stepSizes, stepSigma)
% Create the  model data structure M, based on simple parameters.
%   We assume the monotonic sequence of states is 1->2, 2->3,..., n->1.
% nu is the number of position states (e.g. 128).
% yQuantum is the number of nm per unit (e.g. 2 -> 2nm per y unit).
% noiseSigma is the noise standard deviation, in nm.
% transProb is a vector of transition probabilities a12, a23, etc.
% stepSizes is a vector of the position increments, in nm; element i gives the
%   increment between molecular state i and i+1 in nm.  These values are
%   the means of Gaussian distributions that are assigned to the step
%   size distribution.  The SD of each of these distributions is given by
%   the stepSigma vector (which can alternatively be a scalar) also in
%   nm.  If the jth element of stepSizes is zero, this means that no
%   steps will be taken.
%   Negative stepSigma values are taken to mean that the step size
%   distribution consists of delta functions at stepSizes(i)+stepSigma and
%   stepSizes(i)-stepSigma.  stepSigma is optional and may be a scalar.
%   Default value: 0.
% Fixed step size bug for ns=1 fs 27Dec08
%
% The returned model structure M contains
%    M.C the big TPM (matrix of vectors) nu x ns x ns, in y-units
%    M.P0 initial probabilities nu x ns
%    M.Sigma noise s.d. scalar, in y units
%    M.SigmaOld (set to zero)
%    M.YQuantum (copied from yQuantum)
%    M.DutyCycle = 0 (default)
%    M.Epsi = 1 (default)
%
% Example parameters:
% nu=32;
% yQuantum=2;
% noiseSigma=8;
% transProb=[0.05 0.1 0.05 0.1];
% stepSizes=[20 0 -20 0];
% stepSigma=128;  % very broad, for initializing with a nearly flat distribution.
% dutycycle=1;
if nargin < 6
    stepSigma=0;
end;

ns=numel(transProb);
if numel(stepSigma)<ns
    stepSigma=repmat(stepSigma(1),ns,1);
end;

% Scale all displacements by yQuantum
stepSigma=stepSigma/yQuantum;
StepP=circshift(diag(stepSizes/yQuantum),[0 1]);
% stepSizes=stepSizes/yQuantum;

% Make an A matrix allowing only transitions to the right.
I=diag(1-transProb);
Is=circshift(diag(transProb),[0 1]); % upper diagonal
A=(I+Is);

if ns < 2  % Special case of no molecular transitions
    A=transProb;
    StepP=stepSizes(1)/yQuantum; % There have to be steps.
end;
% ---Create model from parameters---

%  Make the basic coordinates
x=[-nu/2:nu/2-1]';
x=fftshift(x);

% Create the C array
[ns ns1]=size(A);
halves=[0.5 0.5]';
M.C=zeros(nu,ns,ns);
for i=1:ns
    for j=1:ns
        if StepP(i,j)~=0  % We make a distribution of step sizes
            if stepSigma(i)>0  % distribution is Gaussian
                c=exp(-(x).^2/(2*stepSigma(i)^2));
                c=c/sum(c);   % Normalize it.
                M.C(:,i,j)=A(i,j)*circshift(c,round(StepP(i,j)));
            else  % stepSigma is zero or negative.
                for pol=-1:2:1  % polarity of deviation
                    d=pol*stepSigma(i);  % pick + and - stepSigma
                    k=mod(round(StepP(i,j)+d)+1,nu);
                    M.C(k,i,j)=M.C(k,i,j)+0.5*A(i,j);
                end;
            end;
        end;
        M.C(1,i,j)=A(i,j)-sum(M.C(2:nu,i,j));  % force the column sum to be aij
    end;
end;

if ns==1  % Force the zero-step component to be 1-a12
    M.C(1,i,j)=1-A;
end;

% Initial probability vector: we assume all are equally probable.
M.P0=ones(nu,ns)/(nu*ns);
M.YQuantum=yQuantum;
M.Sigma=noiseSigma/yQuantum;  % Convert to y units.
M.SigmaOld=0;  % Mark the starting model.
M.DutyCycle=0;
M.Epsi=1;  % Essentially no prior.
