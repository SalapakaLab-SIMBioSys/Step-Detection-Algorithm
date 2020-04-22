function xs=Step125(x)
% Find the value in a 1..2..5 x 10^n sequence that is >= x.
%
n=floor(log10(x));
x0=x/10^n;
testvals=[0.5 1 2 5 10];
index=find((x0<=testvals),1);
xs=testvals(index)*10^n;
