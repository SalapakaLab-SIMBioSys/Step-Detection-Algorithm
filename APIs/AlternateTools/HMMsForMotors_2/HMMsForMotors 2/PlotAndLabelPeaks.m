function PlotAndLabelPeaks(x,y,minsize,xscale)
% function PlotAndLabelPeaks(x,y,minsize,xscale)
% Plot a histogram of y with the peaks labeled with the (scaled) x values.
% Function called by DisplayModel.
% minsize is the size of the smallest peak to be labeled, relative to the
% largest peak; default is .02.  xscale is the bin width; default is 1.

if nargin<4
    xscale=1;
end;
if nargin<3
    minsize=.02;
end;

border=0.1; % Extra space at top of plot for labels on peaks
yspace=.005; % extra space below label

n=numel(y);
dy=diff(y);
yp=y;

% Get the global maximum value.
[mxval mxpt]=max(yp);
mxval1=mxval;

% Get all other peak values.
i=1;
while mxval>minsize*mxval1
    mxv(i)=mxval;
    mxp(i)=mxpt;
    i=i+1;

    % march downhill from the peak, marking the first local minima found.
    p1=mxpt;
    while (p1<n) && (y(p1+1)<y(p1))
        p1=p1+1;
    end;
    p0=mxpt;
    while (p0>1) && (y(p0-1)<y(p0))
        p0=p0-1;
    end;
    % blank the points surrounding the present peak.
    yp(p0:p1)=0;
    % find the next peak
    [mxval mxpt]=max(yp);
end;
npts=i-1;

% Plot the histogram, leaving border space above it.
bar(x,y,'b');
absmax=max(y);
axis([-inf inf 0 absmax*(1+border)]);

% Label the peaks
for i=1:npts
    j=x(mxp(i));
    text(j,mxv(i)+yspace*absmax,num2str(j*xscale),'HorizontalAlignment','center','VerticalAlignment', 'baseline','fontsize',9,'color',[0 0 0]);
end;
