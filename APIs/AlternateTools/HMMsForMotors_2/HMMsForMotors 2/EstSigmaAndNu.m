function [sigma nu]=EstSigmaAndNu(Y,I0)
% Find reasonable values for these variables based on the step data vector
% Y.

if nargin<2
    I0=1;
end;
nt=numel(Y);
d=abs(diff(Y));
sigma=median(d./sqrt(I0));          % Works for either vector or scalar I0

dx1=max(d);
dx2=max(abs(Y(1:nt-2)-Y(3:nt)));

nu=2.4*max(dx1,dx2);  % 20% margin.


