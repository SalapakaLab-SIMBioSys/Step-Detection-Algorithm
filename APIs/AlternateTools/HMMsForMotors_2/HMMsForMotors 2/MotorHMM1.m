% MotorHMM1.m
% Analyze noisy motor protein data, using the "one state" HMM and automatic
% assignment of parameters.

% The user needs to provide a data trace in the vector Y, or else use the
% following 6 lines of code to make simulated data.
% --alternating steps of 8 and 28 nm (+/- .3 nm), dwell time = 10

nt=500;  % no. of time points
simnoise=5;  % noise, in nm rms
M0=MakeMonotonicModel(100, 1, simnoise, [0.1 0.1], [8 28], 0.3);
rand('state',0); randn('state',0);  % Reset the random number generators (if desired)
M0.DutyCycle=1;  % Simulate the full CCD integration time
[X Y]=StepSimulatorC(M0,nt);  % Continuous-time simulator
% [X Y]=StepSimulator(M0,nt);   % Alternative, simple discrete-time simulator

% *** Analyze the data vector Y, and if given, the "noiseless" version X
% for comparison. ***

% Analysis parameters
SimpleFB = 1;  % 1: Use the simple but fast FFT-based calculation. No prior.
               % 0: use the full HMM.
quantFrac=0.5;  % Maximum size of the quantum, as fraction of noise sigma.
tol=1e-5;       % Step-to-step tolerance 
maxiters=200;

% If an X vector isn't provided, just copy Y for display purposes.
try ntx=numel(X); catch err; ntx=0; end;
if ntx ~= nt
    X=Y;
end;

% Pick values for the quantum and the number of position states nu.
[s n]=EstSigmaAndNu(Y);
yquantum=Step125(quantFrac*s/2.5); % yquantum is always smaller than the estimated sigma times quantFrac.
nu=NextNiceNumber(n/yquantum,7);  % FFT-friendly value.

% quantize the data and display it.
Y1=reshape(round(Y'/yquantum),nt,1);  % force it to be a column vector.
X1=reshape(round(X/yquantum),nt,1);
figure(1); clf;
plot([X1,Y1]);
xlabel('Time');
ylabel(['Position in quanta of ' num2str(yquantum) ' nm']);
drawnow;

% Starting model.  The wild guess is that the dwell time is about sqrt(nt).
M=MakeMonotonicModel(nu,yquantum,0,1/sqrt(nt),1,inf); % flat step-size distribution
M.DutyCycle=1;       % Only used when SimpleFB=0.
M.Epsi=1/sqrt(nt);   % For MAP prior.  Only used when SimpleFB=0.
Mold=M;
d=inf;
iter=1;
figure(2);  % This figure shows the progress as we iterate.

% Main loop
while (d>tol) && (iter<maxiters)
    if SimpleFB
        [L M gamma] = ForBackF(M,Y1);  % Fast, FFT-based
    else
        [L M gamma] = ForwardBackward_3(M,Y1);
    end;
    d=RelModelChange(M,Mold);
    Mold=M;
    str=sprintf('nu = %d  quant=%4.1g   diff = %8.4g',nu,M.YQuantum,d);
    DisplayModel(M, L, iter, str, 0);
    iter=iter+1;
end;

% Done optimizing the model M.  Now get the Viterbi restored timecourse
EstY=ViterbiRestoration(M,Y1);
% Make the expanded plot(s) of the restoration with the data
figure(3); clf;
if ntx>0  % We plot the X trace, detecting steps assuming it is noiseless.
    subplot(2,1,2);
    RestorationPlotM(EstY,X1,nu,M.YQuantum,0,'b+');
    title('Restored timecourse (dots) and X');
    subplot(2,1,1);
end;
RestorationPlotM(Y1,EstY,nu,M.YQuantum);
title('Y and restored timecourse');

% Re-plot the original data with the restoration
figure(1);
plot([X1 Y1 EstY]);
legend('X', 'Y', 'Restored', 'location', 'southeast');
xlabel('Time');
ylabel(['Position in quanta of ' num2str(yquantum) ' nm']);
