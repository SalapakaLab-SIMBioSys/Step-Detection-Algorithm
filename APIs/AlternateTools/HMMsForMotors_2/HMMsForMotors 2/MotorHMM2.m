% MotorHMM2.m
% Analyze noisy motor protein data, using a two-state HMM
% More complicated problem, uses MAP estimation.

% --alternating steps of 8 and 16 nm, dwell time = 5
nt=400;  % no. of time points
simnoise=10;  % noise, in nm rms
M0=MakeMonotonicModel(100, 1, simnoise, [0.2 0.2], [8 16], 0);

rand('state',0); randn('state',0);  % Reset the random number generators (if desired)
M0.DutyCycle=1;  % Simulate the full CCD integration time

[X Y]=StepSimulatorC(M0,nt);  % Continuous-time simulator

% Analysis parameters
SimpleFB = 0;   % 0: use the full HMM.
tol=1e-6;       % Step-to-step tolerance for convergence test
maxiters=200;

% Pick values for the quantum and the number of position states nu.
yquantum=4;
nu=64;

% quantize the data and display it.
X=X/yquantum;
Y=Y/yquantum;


figure(1); clf;
plot([X,Y]);
xlabel('Time');
ylabel(['Position in quanta of ' num2str(yquantum) ' nm']);
drawnow;

% Starting model.  The wild guess is that the dwell time is about sqrt(nt).
M=MakeMonotonicModel(nu,yquantum,0,[.4 .4],[-10 -5],200); % nearly flat step-size distributions
M.DutyCycle=1;       % Only used when SimpleFB=0.
M.Epsi=.01;          % MAP prior: maximum dwell ~100 time points.
Mold=M;
d=inf;
iter=1;
figure(2);  % This figure shows the progress as we iterate.

% Main loop
while (d>tol) && (iter<maxiters)
    if SimpleFB
        [L M gamma] = ForBackF(M,Y);  % Fast, FFT-based
    else        [L M gamma] = ForwardBackward_3(M,Y);
    end;
    d=RelModelChange(M,Mold);
    Mold=M;
    str=sprintf('nu = %d  quant=%4.1g   diff = %8.4g',nu,M.YQuantum,d);
    DisplayModel(M, L, iter, str, 0);
    iter=iter+1;
end;

% Done optimizing the model M.  Now get the Viterbi restored timecourse
EstY=ViterbiRestoration(M,Y);
% Make the expanded plot(s) of the restoration with the data
figure(3); clf;
% We plot the X trace, detecting steps assuming it is noiseless.
subplot(2,1,2);
RestorationPlot(EstY,X,nu,M.YQuantum,0,'b+');
title('Restored timecourse (dots) and X');
subplot(2,1,1);
RestorationPlot(Y,EstY,nu,M.YQuantum);
title('Y and restored timecourse');

% Re-plot the original data with the restoration
figure(1);
plot([X Y EstY]);
legend('X', 'Y', 'Restored', 'location', 'southeast');
xlabel('Time');
ylabel(['Position in quanta of ' num2str(yquantum) ' nm']);
