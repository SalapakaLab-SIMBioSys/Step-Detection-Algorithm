function [EstY EstI LP Wrap]=ViterbiRestoration(M,Y,I0)
%
% Viterbi algorithm for step-HMM
% Using the multi-state model structure M
% F. Sigworth 23 Nov 04
% revision 7 Dec 04
% rev 10 Apr 07 to call Makeb function.  -fs
%
% M structure
%  M.C big TPM (matrix of vectors) nu x ns x ns
%  M.P0 initial probabilities nu x ns
%  M.Sigma noise s.d. scalar
% returns the restoration EstY, the estimated sequence of molecular states
% EstI, the log probability of the optimum path LP,
% and a vector Wrap that can be added to Y or EstY to "wrap" them to
% the range 1...nu.

ErfMode=1;  % Use erf for b, instead of gaussian.
if nargin<3
    I0=1;
end;
nt=numel(Y);  % number of time steps
if numel(I0)<nt  % Handle the case where I0 is a scalar.
    I0=I0(1)*ones(nt,1);
end;

[nu ns ns1]=size(M.C);  % number of levels, states

% make log b
% x0=[-nu/2:nu/2-1]';
% x=fftshift(x0);
% for t=1:nt
%     %b(:,t)=exp(-x.^2/(2*(M.Sigma*I0(t)).^2));
%     b(:,t)=exp(-x.^2/(2*(M.Sigma/I0(t)).^2));  %*****Sigma divided by I*******
%                                                %****'b' is now 2D***
% end;

% ***make log b
Logb=zeros(nu,nu,nt); %initialize the matrix
for t=1:nt
    b0=Makeb(nu,M.Sigma/sqrt(I0(t)),M.DutyCycle,ErfMode);
    b0=fftshift(b0);
    Logb(:,:,t)=log(QuickShift( b0, [Y(t),Y(t)] ));
end;

% Old code:
% [u v]=ndgrid(-nu/2:nu/2-1,-nu/2:nu/2-1);
% u=fftshift(u); v=fftshift(v);
% Logb=zeros(nu,nu,nt); %initialize the matrix
% sigma=M.Sigma;
% for t=1:nt
%     totsigma=zeros(nu,nu);
%     totsigma=(sigma/I0(t))^2+(v-u).^2/12;
%     holder=1./sqrt(2*pi*totsigma).*exp(-(((v+u)/2).^2)./(2*totsigma)); %contruct b
%     top=circshift(holder,[Y(t),Y(t)]);
%     bottom=sum((sum(top)));
%     Logb(:,:,t)=log(top/bottom); %normalize the matrix at every time point
% end;

% Make the big A matrix
N=nu*ns;  % dimension of the whole problem.
A=zeros(N,N);
C=zeros(nu,nu);
for i=1:ns
    for j=1:ns
        %  The columns of the submatrix C are the log of the C vector.
        ct=log(max(M.C(:,i,j),eps));
        for k=1:nu
            C(k,:)=(circshift(ct,k-1))';
            %ct=circshift(ct,1);
        end;
        % insert the submatrices
        ii=nu*(i-1);
        jj=nu*(j-1);
        A(ii+1:ii+nu,jj+1:jj+nu)=C;
    end;
end;
d=zeros(N,nt); % for debugging

% t=1 point
P=reshape(M.P0,N,1);    % Convert to a single long vector
delta=max(log(max(P,eps))+repmat(diag(Logb(:,:,1)),ns,1)); % log maximum starting probability; eq. 21
delta=repmat(delta,N,1);
psi=zeros(N,nt);  % vectors that point back to the previous best state
%     d(:,1)=delta;


% Forward recursion;
for t=2:nt
    [q p]=max(A+repmat(delta,1,N)+repmat(Logb(:,:,t),ns,ns)); %q is max(phi*c) in eq 22
    psi(:,t)=p'; %arguments (indices) of max[delta_i*a_ij
    delta=q';   %delta is phi(t+1)in eq 22 and adding B since it is log
end;

% Compute the log probability of the best path
mstate=zeros(nt,1);  % Trace of metastates
[LP s]=max(delta);
mstate(nt)=s; %is reshaped matrix contin i's and u's

% Backward recursion
for t=nt:-1:2
    s=psi(s,t);
    mstate(t-1)=s; %contains pos ad state of max prob at each time
end;

EstU=zeros(nt,1);
EstI=zeros(nt,1);
EstY=zeros(nt,1);
Wrap=zeros(nt,1);

% Convert the metastate numbers to amplitudes and state indices
EstU=mod(mstate-1,nu);
EstI=floor((mstate-1)/nu)+1;

% Unwrap the amplitudes

% EstY=EstU+nu*round((Y-EstU)/nu);
Wrap=nu*round((EstU-Y)/nu);

% handle NaNs in the data
pt=find(isnan(Wrap));
for i=pt'
    Wrap(i)=Wrap(i-1);
end;

EstY=EstU-Wrap;



figure(3);
% % Here is an output overview
% % plot([Y+Wrap EstU X+Wrap 5*EstI-14]);
% subplot(1,2,1);
% xvals=(1:nt)';
% plot([Y+Wrap EstU]);
% subplot(1,2,2);
% plot([Y EstY]);

