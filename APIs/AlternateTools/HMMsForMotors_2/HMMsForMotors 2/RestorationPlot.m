function RestorationPlot(Y, EstY, wrap, yquantum, EYbase, markerstr)
% function RestorationPlot(Y, EstY, wrap, yquantum, EYbase,markerstr) Fancy
% plotting of the original data Y and the Viterbi reconstruction EstY,
% labeling the steps.  Each step label shows the step size, corrected by
% yquantum. The argument wrap sets the span of the plot; a typical value is
% nu. The optional vector EYbase (same size as EstY) is the baseline
% estimate, added to EstY for plotting.  The argument markerstr is the
% description of the marker Matlab uses for plotting Y.  Default is 'r+',
% i.e. red crosses.

% Modified to use yquantum 8 Apr 07 -fs
% Added the EYbase argument 18 June 07 -fs
% Added TextXPadFactor etc. to prevent overlapping labels 11 Aug 07 -fs

if nargin<4
    yquantum=1;
end;
if nargin<5
    EYbase=0;
end;
if nargin<6
    markerstr='r+';
end;

Hyst=0.1;  % Hysteriesis (fraction of wrap) for wrapping.
border=.05; % border (fraction of wrap) above and below the plot
drifteps=0.9; % allowable linear drift not detected as steps.
TextXPadFactor=1; % min spacing of text labels to avoid overlap.  Set to zero to turn off spacing check.
TextXShift=0.5;   % Shift text over closer to lines.

nt=numel(Y);

% Find the steps in the estimated Y, taken to be always integers
Ysteps=round(EstY(2:nt)-EstY(1:nt-1));
Ysteps(nt)=0;

% Scale everything into physical units.
Ysteps=yquantum*Ysteps;
Y=yquantum*Y;
EstY=yquantum*EstY;
wrap=wrap*yquantum;
EYbase=EYbase*yquantum;


StepPtrs=find(Ysteps);  %vector of time pnts where nonzero step happens
%to get step sizes do hist(Ysteps(StepPtrs))

% Set up for a display that wraps around to make it not too tall
Shift=zeros(nt,1);  % Displacement to add to plotted data points
PShift=Shift;       % Same, but with NaNs to prevent connecting points at wraps
maxy=(1+Hyst)*wrap;
ceily=1.5*wrap; %max values of y in each plot; 1.5 there to allow noisy data to stay in one wrap
miny=-2*Hyst*wrap;
floory=-0.5*wrap;

for t=1:nt
    diff=EstY(t)-Shift(t);
    if diff > wrap
        Shift(t+1:nt)=EstY(t);  % Force it to zero; contains offset values for neigh plots
        PShift(t+1)=NaN;    % To avoid plotting a point, we assign it NaN
        Y(t+1)=NaN;  %%%
        maxy=min(ceily,max(maxy,diff));
    elseif diff < -Hyst*wrap
        Shift(t+1:nt)=EstY(t)-wrap;
        PShift(t+1)=NaN;
        miny=max(floory,min(miny,diff));
    end;
end;

PShift=PShift(1:nt);
PShift=PShift+Shift;

x=1:nt;
% PltY=EstY-Shift;

% Do the plotting here.
plot(x,Y-Shift,markerstr,x,EstY-PShift,'k-','markersize',3);
axis([0 inf miny-border*wrap maxy+border*wrap]);
qtextx=nt/10;
qtexty=maxy+2*border*wrap;
extent1=zeros(1,4);
xok=1;
yok=1;

for i=1:numel(StepPtrs)
    p=StepPtrs(i);             %x coord of step label.  Add 0.5 to get to center of step.
    stepVECTOR(i)=Ysteps(p);
    str=num2str(stepVECTOR(i));    % step label
    tx=p+TextXShift;
    ty=EstY(p)-Shift(p)+Ysteps(p)/2;  %y coord of step label
    ty=max(miny,(min(maxy,ty)));
    if TextXPadFactor>0
        % I guess we have to draw some text to see how big it is...
        h2=text(qtextx,qtexty,str,'HorizontalAlignment','left','fontsize',9,'color',[1 1 1],...
            'fontweight','normal');  % Draw invisible (white) text just to get its size.
        extent2=get(h2,'Extent');
        xok=tx>extent1(1)+TextXPadFactor*(extent1(3)+extent2(3));
        yok=ty<extent1(2)-extent1(4)/2 || ty>extent1(2)+3*extent1(4)/2;
    end;

    if xok || yok % not too close to previous: draw the label
        h1=text(tx,ty,str,'HorizontalAlignment','right','fontsize',9,'color',[.2 0 .5],...
            'fontweight','normal');
            extent1=get(h1,'Extent'); % extent=[left bottom width height]
    end;

end;
%stepVECTOR'  %print the steps to screen
xlabel('Time, frames');
ylabel('Position, nm');
