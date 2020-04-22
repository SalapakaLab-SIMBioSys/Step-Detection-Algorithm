function n=NextNiceNumber(x,f)
% function n=NextNiceNumber(x)
% Find the next even integer >= x whose largest
% prime factor is f.  This is used for finding a nice dimension of vectors
% for use with mixed-radix FFTs.

n=ceil(x);
% force x to be even
if mod(n,2)>0
    n=n+1;
end;
% Increment n until the largest factor is <= f.
factors=factor(n);
while max(factors)>f
    n=n+2;
    factors=factor(n);
end;
