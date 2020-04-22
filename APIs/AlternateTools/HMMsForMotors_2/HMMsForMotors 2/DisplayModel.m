function DisplayModel(M, L, iterations, string, BottomHalf)
% Display the HMM described by M, the log likelihood value,
% the iteration count, and an arbitrary string.
% If BottomHalf=1, we use only the bottom half of the figure.
% Modified for M.YQuantum field 8 Apr 07 -fs
% Modified to use cla, not clf when drawing 18 Jun 07 -fs

% Handle default values
if nargin<2
  L=0;
end;
if nargin<3
  iterations=0;
end;
if nargin<4
    string=[];
end;
if nargin<5
    BottomHalf=0;
end;

[nu ns ns1 ni]=size(M.C);
if mod(nu,2)==1 % odd nu
    x0=[-(nu-1)/2:(nu-1)/2]';
else  % even nu
    x0=[-nu/2:nu/2-1]';
end;
x1=x0*M.YQuantum;  % X values in nm

% Figure out which C entries are degenerate (that is, have all their mass
% in the first element of the C(:,i,j) vector.  We keep track of all the
% non-degenerate entries to make a graph of each one.
hint=zeros(ns,ns);
ixa=0;  % number of A-matrix elements to show
ixd=0;  % number of distributions to show

% if ns>1  % multi-state models
for i=1:ns
    for j=1:ns
        q=sum(M.C(:,i,j));
        if (q>0)&&((q==M.C(1,i,j)) || i==j) % Transition with no step, or only transition
            ixa=ixa+1;
            inai(ixa)=i;
            inaj(ixa)=j;
            vala(ixa)=M.C(1,i,j);
        end;
        if (q>0) && (q>M.C(1,i,j))     % there is a distribution
            ixd=ixd+1;
            indi(ixd)=i;
            indj(ixd)=j;
        end;
    end;
end;

% Compute the number of rows and columns in the figure
b=(BottomHalf>0);
nr=floor(sqrt(ixd+1))+b;
nc=ceil((ixd+1)/(nr-b));

% List the parameters in the first subplot area
subplot(nr,nc,1+b*nc);
cla;
axis off
dy=0.1;
x=0; y=0.95;

% Draw the string
if numel(string)>0
    %text(x,y,string); The line below was added to suppress TeX
    %interpretation of the string...it was giving error each time!
    text('Position',[x,y],'Interpreter','none','String',string)
    y=y-dy;
end;

% Draw the number of iterations
if iterations>0
  text(x,y,['Iterations: ' num2str(iterations)]);
  y=y-dy;
end;

% Draw the L value
text(x,y,['L = ' num2str(L)]);
y=y-dy;
text(x,y,['sigma = ' num2str(M.Sigma*M.YQuantum)]);
y=y-dy;

% Draw each transition probability
for i=1:ixa
  string=['a_{' num2str(inai(i)) num2str(inaj(i)) '} = ' num2str(vala(i))];
  if inai(i)==inaj(i)
      string=[string '  (dwell = ' num2str(1/(1-vala(i))) ')'];
  end;
  text(x,y,string);
  y=y-dy;
%   if inai(i)==inaj(i)
%       text(x,y,['  dwell = ' num2str(1/(1-vala(i)))]);
%       y=y-dy;
%   end;
end;

% Draw each distribution
for i=1:ixd
  subplot(nr,nc,i+1+b*nc);
  c=M.C(:,indi(i),indj(i));
  if indi(i)==indj(i)
    c(1)=0;
  end;
%   plot(x1,fftshift(c));
 PlotAndLabelPeaks(x1,fftshift(c));
 meanc=x1'*fftshift(c)/sum(c);
%    [mxval mxinds]=max(fftshift(c));
 % Show the mean and mode values below the graph
  xlabel(['Step size, nm   C_{' num2str(indi(i)) num2str(indj(i)) '}   mean = ' num2str(meanc)]);
  ylabel('Step probability');
end;

drawnow;
