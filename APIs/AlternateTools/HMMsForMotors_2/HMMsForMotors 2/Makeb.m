function b=Makeb(nu,sigma,DutyCycle,ErfMode);
% ErfB.m
% Calculate b using error functions.
% If ErfMode = 0 then we compute b from a Gaussian function.

if nargin<3
    DutyCycle=0;
end;
if nargin<4
    ErfMode=1;
end;


% The b function is the distribution (H(u)-H(v))/(u-v), convolved with the
% Gaussian with variance sigma.  (H is the step function.)  After
% convolution this is 1/2*(erf(u/s2s)-erf(v/s2s)).  To handle the case u=v,
% in which case the distribution is a Gaussian, we can choose u and v to be
% very similar, by adding and subtracting a small value epsi.

if ErfMode

    % actual computation
    s2s=sqrt(2)*sigma;
    epsi=.001;  % delta for differentiating the error function to get Gaussian.
    p=-nu/2;
    q=nu/2-1;
    [u v]=ndgrid(p-epsi:q-epsi,p+epsi:q+epsi);
%     u=fftshift(u);
%     v=fftshift(v);

    if DutyCycle > 0
        eucol=0.5*erf(u(:,1)/s2s);
        eu=repmat(eucol,1,nu);
        evrow=0.5*erf(v(1,:)/s2s);
        ev=repmat(evrow,nu,1);
        b1=(eu-ev)./(u-v);
    else
        b1=0;
    end;
    % b1=0.5*(erf(u/s2s)-erf(v/s2s))./(u-v);  % Accomplishes the same thing
    % % but makes erf evaluations of each element.

    % In the case of DutyCycle<1 we add in a Gaussian dependent on v.
    if DutyCycle < 1
        rowg=exp(-(v(1,:).^2/(2*sigma^2)));
        rowg=rowg/sum(rowg);
        b0=repmat(rowg,nu,1);
        %     b0=exp(-(v.^2/(2*sigma^2)));
        %     b0=b0/sum(b0(1,:));

    else
        b0=0;
    end;
    % Final b function
    b=DutyCycle*b1+(1-DutyCycle)*b0;

else % not ErfMode

    % % Old gaussian b function code
    %
    [u v]=ndgrid(-nu/2:nu/2-1,-nu/2:nu/2-1);
%     u=fftshift(u); v=fftshift(v);

    totsigma=(sigma)^2+(v-u).^2/12;
    b=1./sqrt(2*pi*totsigma).*exp(-(((v+u)/2).^2)./(2*totsigma)); %contruct b
    % b=b/sum(sum(b)); %normalize the matrix at every time point

end;
