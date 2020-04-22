function ms=QuickShiftMex(m,shift)
% Faster replacement for circshift for time-critical applications.
% About twice as fast as circshift.
% Calls the mex routine QuickShiftMex, if available.
%
    ms=circshift(m,shift);    
