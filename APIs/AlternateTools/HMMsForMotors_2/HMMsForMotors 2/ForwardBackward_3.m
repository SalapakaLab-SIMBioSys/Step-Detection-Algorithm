function [LL newM gamma] = ForwardBackward_3(M,Y,I0,test,BigB)
% function [LL newM gamma] = ForwardBackward_3(M,Y,I0,test,BigB)
% Forward-backward algorithm for step-HMM
% Using the multi-state model structure M.  Only the first two arguments
% are required.  This Version 3 makes MAP estimation of transition
% probabilities
%
% The M data structure defines the model
%  M.C big transition probability matrix, TPM (matrix of vectors), nu x ns x ns
%  M.P0 initial probabilities nu x ns
%  M.Sigma noise s.d. (scalar).  Set to 0 to force an approximate guess.
%  M.DutyCycle relative integration time, between 0 and 1
%  M.SigmaOld previous sigma estimate.  Set to zero to initiate the
%  optimizer.
%  M.Epsi is the prior probability parameter for a_ii.  It should be set to
%  roughly 1/tmax, where tmax is the longest expected dwell time.  To turn
%  off MAP estimation set M.Epsi to a value larger than 1, e.g. 100.
%
% Y is the data vector.
% The function returns the updated model and the log-likelihood value.
% I0 is the intensity vector, one element for each time point, with the
% nominal value equal to 1.  The noise sigma at time point t is taken to be
% sigma/sqrt(I0(t)).
% If I0 is absent, or a scalar, then the same intensity (default=1) is assumed for all
% points.  Points where the position is indeterminant can be indicated by NaN
% values.
% If test=1, timing is printed out.

% Fred's version using the mex functions WtConvol and QuickShift.
% 10 Apr 2007
% Placed bounds on sigma search 31 Dec 07 fs.
% Fixed initial sigma estimation 5 Mar 10 fs.
% Inserted Fiona's corrected ComputeA2 function 23 Mar 10 fs. 

% Functions called: Makeb to obtain the b-function
%                   QuickShift as replacement for circshift (~2x speed)
%                   WtConvol as fast weighted convolution in obtaining Xi
%                   MinOfParab for sigma re-estimation
% this function also calls itself recursively in the sigma re-estimation.

ErfMode=1;  % Zero means to use the expanded Gaussian to make the b calculation.
            % 1 means the more exact erf calculation is used.
% BigB=1; % If I0 isn't given, still store a b matrix for each time point.
        % --this is faster but uses more memory, about 8*nt*nu^2 bytes.

LOnly=nargout<2;  % if no re-estimation is desired, compute only the loglikelihood.
newM=M;           % By default, copy the old parameters.

if nargin<3       % Default: noise sigma is constant
    I0=1;
end;

if nargin<4       % If test=1, timing information is printed out.
    test=0;
end;

if nargin<5
    BigB=1;      % Default is to make the large b-function lookup matrix 
                 % (nu x nu x nt in size).  Set this to zero to save memory
                 % for large datasets
end;

nt=numel(Y);       % number of time steps
[nu ns ns1]=size(M.C);  % nu and ns: numbers of states

if ~((numel(I0)==nt) || (numel(I0)==1))
    disp('I0 must be a scalar or the same size as Y.  Value forced to 1');
    I0=1;
end;

%  Handle the case that there is no sigma estimate.
if M.Sigma==0     % No value at all for sigma: get one
    d=zeros(nt,1);
    d(1:nt-1)=abs(Y(1:nt-1)-Y(2:nt));
    M.Sigma=median(d.*sqrt(I0));          % Fixed. Works for either vector or scalar I0
    M.SigmaOld=0;
end;


% Make the hint matrices
hint=zeros(ns,ns1); %this matrix looks at each distribution in C and finds if it is
                    %only a Delta function
for i=1:ns
    for j=1:ns1
        q=sum(M.C(:,i,j));
        hint(i,j)= (q==0)*0+(q==M.C(1,i,j))*1+(q>M.C(1,i,j))*2;
    end;
end;

% We calculate the C matrix here; The C matrix is the general transition
% prob matrix. Essentially, this part takes M.C(:,i,j) which is a vector as
% a function of w (step size) and converts it into a matrix which is a
% function of u and v. Here, w=v-u. The original vector from M.C is
% calculated for u=0. So, below we calculate additional columns that
% correspond to u=1,2,3... . This is done by vertically shifting the
% columns by 1 with respect to the last column. Lastly, the matrix is
% transposed so that the values corresponding to different u's run down the
% columns and rows shifted down by 1 (the last step may have to do with
% the previous fftshift that was used in generating M.C(x,*,*).
% (this is done intuitively by Fred [commented out] but done in 2 shorter lines
% and in weird way below by sheyum)
% Make the general TPM matrix Cmat(u,v,i,j).
%
Cmat=zeros(nu,nu,ns,ns1);
for i=1:ns
    for j=1:ns1
        c1=M.C(:,i,j); mat=c1;
        for q=1:nu-1
            mat=[mat,QuickShift(c1,q)];
        end;
        Cmat(:,:,i,j)=mat';
    end;
end;

if test, tic, end;

sigma=M.Sigma;
dutycycle=M.DutyCycle;

% construct the new noise or 'b' function. Here b is constructed for all
% Y(t) and I0(t) if I0 is a vector; otherwise it's computed only once.
%

if numel(I0) == nt % VarI is a flag for a valid I0 array.  Otherwise we ignore intensity.
    VarI=1;
    BigB=1;  % Of necessity we'll make the b matrix for each t.
else
    VarI=0;  % We can still use the BigB option if it's set above.
end;

% Make the default b matrices
bnull=ones(nu,nu)/nu; % flat pdf for the case I=0.
b0=Makeb(nu,sigma/sqrt(I0(1)),dutycycle,ErfMode); % Constant intensity: copy this for every datapoint

% If desired, make the big b matrix.
if BigB
    b=zeros(nu,nu,nt); %initialize the matrix
    for t=1:nt
        if VarI
            b0=Makeb(nu,sigma/sqrt(I0(t)),dutycycle,ErfMode);
        end;
        b(:,:,t)=QuickShift( b0, [Y(t),Y(t)] );
    end;
end;

if test
    TimeBMatrix=toc
    tic
end;


% *****CALCULATION OF THE FORWARD VARIABLES******
alpha=zeros(nu,ns,nt);
anorm=zeros(nt,1); %normalization

%Forward variable at t=1 is defined differently than the rest
if BigB
    bt=b(:,:,1);
elseif Y(1) ~= NaN
    bt=QuickShift(b0,[Y(1),Y(1)]);
else
    bt=bnull;
end;
for i=1:ns
    alpha(:,i,1)=M.P0(:,i).*diag(bt);
end;
anorm(1)=squeeze(sum(sum(alpha(:,:,1))));
alpha(:,:,1)=alpha(:,:,1)/anorm(1); %normalize the first alpha matrix

%Forward variables from t=2 to end
%compute alpha(v,j)_t=\sum_u,i (alpha(u,i)_{t-1}'*b.*Cmat(u,v,i,j))'
for t=2:nt
    if BigB
        bt=b(:,:,t);
    elseif ~isnan(Y(t))
        bt=QuickShift(b0,[Y(t),Y(t)]);
    else
        bt=bnull;
    end;
    atemp=zeros(nu,ns);
    for i=1:ns
        for j=1:ns
            if hint(i,j)==2 % general transition
                atemp(:,j)=atemp(:,j)+( alpha(:,i,t-1)'*(bt.*Cmat(:,:,i,j)) )'; %summing over i, for general case
            elseif hint(i,j)==1 % ij transition only; no step
                atemp(:,j)=atemp(:,j)+alpha(:,i,t-1).*diag(bt)*M.C(1,i,j); %summing over i,for this case
            end;
        end;
    end;
    anorm(t)=sum(atemp(:));  %summing over u
    alpha(:,:,t)=atemp/anorm(t); %alpha(v,j,t) normalized
end;
% *****END OF THE FORWARD VARIABLES******
if test
    TimeForward=toc
end;

LL=sum(log(anorm)); %calculate the log-likelihood

if LOnly
    return;   % Quit here if we are only computing the likelihood LL.
end;


if test, tic, end;

% *****CALCULATION OF THE BACKWARD VARIABLES******
beta=zeros(nu,ns,nt);
beta(:,:,nt)=ones(nu,ns)/anorm(nt);   %initializing first beta (at last time point)

%Backward variables from second last time point to beginning
%compute beta(u,i)_t=\sum_v,j (bt(:,:).*Cmat(:,:,i,j))*(beta(:,j,t+1)
%NOTICE: as opposed to the alpha's, the beta's are summed over j (not i) and
%v (not i). The matrix multiplication here looks a little different
%from what we have above in the FORWARD case.

for t=nt-1:-1:1
    if BigB
        bt1=b(:,:,t+1);
    elseif Y(t+1) ~= NaN
        bt1=QuickShift(b0,[Y(t+1),Y(t+1)]);
    else
        bt1=bnull;
    end;
    btemp=zeros(nu,ns);
    for i=1:ns
        for j=1:ns
            if hint(i,j)==2 % general transition
                btemp(:,i)=btemp(:,i)+(bt1.*Cmat(:,:,i,j))*beta(:,j,t+1); %summing over i, for general case
            elseif hint(i,j)==1 % ij transition only; no step
                btemp(:,i)=btemp(:,i)+(diag(bt1).*beta(:,j,t+1))*M.C(1,i,j); %summing over i,for this case
            end;
        end;
    end;
    beta(:,:,t)=btemp/anorm(t); %normalization
end;
% *****END OF THE BACKWARD VARIABLES******

if test
    TimeBackward=toc
    tic
end;

%*****CALCULATION OF GAMMA*****
% Note, we hardly need gamma unless we're outputting it.

gamma = alpha.*beta;  %eq (13) gamma, but  not normalized yet
for t=1:nt
    gamma(:,:,t)=gamma(:,:,t)/sum(sum(gamma(:,:,t))); % gamma normalized at each t value.
end;
gamma=max(0,gamma);  % Force gamma to be non-negative
%*****END OF GAMMA*****

%=>=>=>=>Re-estimation of the initial probability matrix
newM.P0=gamma(:,:,1);

%*****CALCULATION OF Xi; Eq (15)*****
% xi(w,i,j,t)=sum_u  {alpha_t(u,i)*c(w,i,j)*beta_{t+1}(u+w,j)*b(u+w,Y(t+1)).


Xi=zeros(nu,ns,ns,nt-1);  % local re-estimate variable Xi(w,i,j,t)
for t=1:nt-1
    if BigB
        bt1=b(:,:,t+1);
    else
        bt1=QuickShift(b0,[Y(t+1),Y(t+1)]);
    end;

    for i=1:ns
        for j=1:ns
            if hint(i,j)==2 % general transition
                Xi(:,i,j,t)=M.C(:,i,j).*...
                    WtConvol(alpha(:,i,t),beta(:,j,t+1),bt1);

            elseif hint(i,j)==1 % no-step transition
                Xi(1,i,j,t)=M.C(1,i,j)*(alpha(:,i,t)'*(diag(bt1).*beta(:,j,t+1)));

            end;
        end;
    end;
end;
sumXi_t=sum(Xi,4); %sum over t -> sumXi_t(w,i,j)
sumXi_tjw=sum(sum(sumXi_t,3)); % sum over w and j, leaving it as function of i only

% Get the diagonal TPM elements a, then normalize the rest to correspond.

% Standard code: c(w,i,j)=sum_t{Xi(w,i,j,t)}/sum_t,w',j'{Xi(w,i,j',t)}
% The following is the MAP code, assuming the prior 
% p(aii) ~ 1/(1 - aii + M.Epsi).  To turn off MAP estimation, make M.Epsi
% large (say = 1).

c=zeros(nu,ns,ns);
for qi=1:ns
    a=ComputeA(sumXi_t(1,qi,qi),sumXi_tjw(qi),M.Epsi);
    c(:,qi,:)=sumXi_t(:,qi,:)./sumXi_tjw(qi);    %eq. (18), re-estimating c_ij(w)
    c(1,qi,qi)=0;
    normFactor=(1-a)/sum(sum(c(:,qi,:),3));
    c(:,qi,:)=c(:,qi,:)*normFactor;
    c(1,qi,qi)=a;  % Now the sum of c(:,qi,:) should be unity.
end;
newM.C=c;

if test
    TimeReestC=toc
end;

% ***** Sigma Re-estimation *****
% We optimize sigma on the basis of the old model parameters.  It is to be
% able to look at the old parameters that we include this within the F-B
% routine.

% Parameters for re-estimation
minDiff=0.02;
maxFactor=0.6;

if test, tic, end;

% Check to see if the former sigma value has been set.
if M.SigmaOld==0
    M.SigmaOld = 0.8 * M.Sigma;  % Default if there's not past history for optimizer.
end;

si=M.Sigma;
diff=abs(si-M.SigmaOld);  % step size.
diff=min(si*maxFactor,max(diff,minDiff));  % diff is > MinDiff and < MaxFactor*si.

% The middle value of x is already given, and its corresponding L has been computed.
x(2)=log(si);
L(2)=LL;
for i=1:2:3 % Pick up the other two values.
    s=si+diff*(i-2);
    x(i)=log(s);   % L is a nicer parabolic function of log(sigma) than of sigma.
    M.Sigma=s;
    L(i)=ForwardBackward_3(M,Y,I0);  % Call ourselves recursively.
end;

newM.SigmaOld=si;
xmax=MinOfParab(L,x);
xmax=max(x(1)-0.5,min(xmax,x(3)+0.5));  % force the value to be within the range searched.
newM.Sigma=exp(xmax);

if test
    TimeReestSigma=toc
    L
    sigmas=exp(x)
end;


%*********Local functions*********

% function a=ComputeA(b,B,epsi)
% % Use quadratic formula to solve for diagonal tpm elements
% beta=1+epsi;
% q=beta.*B+b+1;
% a=(q - sqrt(q.^2-4*B.*(beta.*b+1)))./(2*B);

function a=ComputeA(b,B,epsi)  % Fiona's corrected code
% Use quadratic formula to solve for diagonal tpm elements
beta=1+epsi;
q=beta.*B+b-1;
a=(-q + sqrt(q.^2+4*b.*(1-B).*(beta)))./(2*(1-B));

function abscissa_x=MinOfParab(funcpnts,xpnts)
a=xpnts(1); b=xpnts(2); c=xpnts(3);
fa=funcpnts(1); fb=funcpnts(2); fc=funcpnts(3);

top=(b-a)^2*(fb-fc)-(b-c)^2*(fb-fa);
bottom=(b-a)*(fb-fc)-(b-c)*(fb-fa);
abscissa_x=b-(0.5*top/bottom);
