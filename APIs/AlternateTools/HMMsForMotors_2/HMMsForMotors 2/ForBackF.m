function [LL newM gamma] = ForBackF(M,Y)
% Forward-backward algorithm for step-HMM
% Using the multi-state model structure M
% Uses FFT speed-up.
% F. Sigworth 22 Nov 04
% Last modification 7 Dec 04

% The M data structure defines the model
%  M.C is the big TPM (matrix of vectors) nu x ns x ns
%  M.P0 initial probabilities nu x ns
%  M.Sigma noise s.d. scalar

% Y is the data vector.  Missing points are marked as NaN values.
% the function returns the updated model and the log-likelihood value.

LOnly=nargout<2;  % if no re-estimation is desired, compute only the loglikelihood.
newM=M;               % By default, copy the old parameters.

nt=numel(Y);  % number of time steps
[nu ns ns1]=size(M.C);  % numbers of states

I0=1;  % This version doesn't handle I0 variations.

%  Handle the case that there is no sigma estimate.
if (M.Sigma<=0 || isnan(M.Sigma))    % No value at all for sigma: get one
    d=zeros(nt,1);
    d(1:nt-1)=abs(Y(1:nt-1)-Y(2:nt));
    M.Sigma=median(d./sqrt(I0));          % Works for either vector or scalar I0
    M.SigmaOld=0;
end;

% Make the hint matrices
hint=zeros(ns,ns);
for i=1:ns
    for j=1:ns
        q=sum(M.C(:,i,j));
        if q==0
            hint(i,j)=0;    % case 0: no transition
        elseif q==M.C(1,i,j)
            hint(i,j)=1;    % case 1: state transition but no step
        else
            hint(i,j)=2;    % case 2: general transition
        end;
    end;
end;

% hint=ones(ns,ns)*2;

% Make a wrapped-around vector x() like FFT points, such that
% x(1)=0; x(nu)=-1;
% x(nu/2)=nu/2-1; x(nu/2+1)=-nu/2.

x0=[-nu/2:nu/2-1]';
x=fftshift(x0);

%     tic;

% The b function is a Gaussian with s.d. sigma
% we assume it is an even function.
b=exp(-x.^2/(2*M.Sigma^2));
b=b/sum(b);


% disp('Forward');
alpha=zeros(nu,ns,nt);
anorm=zeros(nt,1);

for i=1:ns
    alpha(:,i,1)=M.P0(:,i).*circshift(b,Y(1));
end;
anorm(1)=squeeze(sum(sum(alpha(:,:,1))));
alpha(:,:,1)=alpha(:,:,1)/anorm(1);  % bug fixed: alpha(1) is now normalized.

% FFT version is faster for nu> 64, much faster for 256 or greater.
fC=fft(M.C);  % 1d FFT of the c vectors

for t=2:nt
    atemp=zeros(nu,ns);
    if isnan(Y(t))
        bt=ones(nu,1)/nu;  % no data available at this point
    else
        bt=circshift(b,Y(t));
    end;
    for i=1:ns
        for j=1:ns
            if hint(i,j)==2 % general transition
                fa=fft(alpha(:,i,t-1));
                atemp(:,j)=atemp(:,j)+bt .* real(ifft(fC(:,i,j).*fa));

            elseif hint(i,j)==1 % ij transition only
                atemp(:,j)=atemp(:,j)+M.C(1,i,j)*bt.*alpha(:,i,t-1);
            end;
        end;
    end;
    anorm(t)=sum(atemp(:));
    alpha(:,:,t)=atemp/anorm(t);
end;

LL=sum(log(anorm));

if ~LOnly

    % disp('Backward');
    beta=zeros(nu,ns,nt);
    beta(:,:,nt)=ones(nu,ns)/anorm(nt);

    fCt=conj(fC);
    for t=nt-1:-1:1
        btemp=zeros(nu,ns);
        if isnan(Y(t+1))
            bt1=ones(nu,1)/nu;  % no data available at this point
        else
            bt1=circshift(b,Y(t+1));
        end;
        for i=1:ns
            for j=1:ns
                if hint(i,j)==2  % general case
                    fb=fft(beta(:,j,t+1).*bt1);
                    btemp(:,i)=btemp(:,i)+real(ifft(fCt(:,i,j).*fb));

                elseif hint(i,j)==1
                    btemp(:,i)=btemp(:,i)+M.C(1,i,j)*bt1.*beta(:,j,t+1);
                end;
            end;
        end;
        beta(:,:,t)=btemp/anorm(t);
    end;

    % normalization of gamma, to match definition
    gamma = alpha.*beta;
    for t=1:nt    % gamma is normalized at each t value.
        gamma(:,:,t)=gamma(:,:,t)/sum(sum(gamma(:,:,t)));
    end;
   % Force gamma to be non-negative
    gamma=max(0,gamma);


    % disp('Re-estimation');

    % Re-estimate sigma
    ssum=0;
    s0=0;

    for t=1:nt
        if ~isnan(Y(t))  % i.e. valid data point
            x2=circshift(x,Y(t)).^2;
            s2=0;
            for i=1:ns
                s2=s2+gamma(:,i,t)'*x2;
            end;
            ssum=ssum+s2;
            s0=s0+1;
        end;
    end;
    sigmares=sqrt(ssum/s0);

    cre=zeros(nu,ns,ns);  % local re-estimate variable
    cres=zeros(nu,ns,ns);  % accumulator for reestimates of C
    cris=zeros(1,ns);  % accumulator for sum_jwt{cre(wijt)}, a function of i only.
    for t=1:nt-1
        if isnan(Y(t+1))
            bti1=ones(nu,ns)./(nu*ns);
        else
            bti1=repmat(circshift(b,Y(t+1)),1,ns);
        end;
        fb=fft(beta(:,:,t+1).*bti1);
        fa=conj(fft(alpha(:,:,t)));
        for i=1:ns
            for j=1:ns
                if hint(i,j)==2
                    cre(:,i,j)=M.C(:,i,j).*real(ifft(fa(:,i).*fb(:,j)));
                elseif hint(i,j)==1
                    %                         cre(1,i,j)=M.C(1,i,j)*real(fa(:,i)'*fb(:,j))/nu;
                    cre(1,i,j)=M.C(1,i,j)*alpha(:,i,t)'*(beta(:,j,t+1).*bti1(:,j));
                end; % if
            end; % j
        end; % i
        cre=cre/sum(cre(:));  % cre is now the xi variable, xi(w,i,j).
        cres=cres+cre;
        cris=cris+sum(sum(cre),3);  % sum over w and j
    end; % for t
    for i=1:ns
        cres(:,i,:)=cres(:,i,:)/cris(i);
    end;

    newM.C=cres;
    newM.Sigma=sigmares;
    newM.SigmaOld=M.Sigma;
    newM.P0=gamma(:,:,1);

    %         figure(2);
    %         DisplayModel(M,L,iter);
end; % if LOnly
