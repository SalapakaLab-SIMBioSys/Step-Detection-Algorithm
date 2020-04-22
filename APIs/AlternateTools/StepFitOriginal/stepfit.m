function est = stepfit(data,varargin)
%
%USAGE:
%x = stepfit(y)
%
%   x           Step Fit
%   y           Data on which to fit steps
%
%x = stepfit(y,[tuning parameter],value)
%
%   Optional Tuning Parameters:
%   Fs          Sampling time. Required if tau is being specified.
%   tau         Time constant of probe response
%   outputnoise Standard deviation of the noise measured at output, Estimates noise at the output if unspecified.
%               Noise estimation is not reliable under extreme dynamical distortion and high stepping speed.
%   measnoise   Standard deviation of the noise in the measurement system.
%               Assumes zero measurement noise if unspecifed
%   passes      Number of iterations in histogram refinement
%   verbose     verbose = 1: Display progress of histogram.  verbose = 2: Save the progress in fig files
%
%Author: Tanuj Aggarwal
%Algorithm Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3743561/pdf/nihms470920.pdf

global Fs tau outputnoise measnoise passes verbose

numarg = nargin-1;
if(numarg<0)
    help(mfilename)
    return
end

Fs = 1;
tau = 0;
outputnoise = noise_std(data)';
measnoise = 0;
passes = 3;
verbose = 0;
for i = 1:2:numarg
    evalc([varargin{i} ' = ' num2str(varargin{i+1})]);
end
% figure()
% plot([1,2,3],[tau, tau, tau],'linewidth',2);
% hold on
% plot([1,2,3],[Fs, Fs, Fs],'linewidth',2);
% hold off
% pause on
% pause(5)


fprintf(' Fs = %.2g \n tau = %.2g \n outputnoise = %.2g \n measnoise = %.2g \n passes = %.2g \n verbose = %.2g\n',Fs,tau,outputnoise,measnoise,passes,verbose);
est = stepfit1(data(:),Fs,tau,outputnoise,measnoise,passes,verbose);
% est=smoothfitsub(est);
est = reshape(est,size(data));
% LLR = stepfitLLR(est,data);
% fprintf('\nStepping Likelihood = %.2g\n',LLR);
if(nargout==0 && verbose==0)
    plot([data(:) est(:)])
end
disp('Done')


function sm1=smoothfitsub(fit)
fit=fit(:);
sm=fit;
sm1=fit;
steplocations=find(diff(fit));
if(isempty(steplocations))
    LLR=0;
else
    steplocations=[0 ; steplocations ;length(fit)];
    ind=zeros(length(steplocations)-1,1);
    for i=1:length(steplocations)-1
        ind(i)=round((steplocations(i)+1+steplocations(i+1))/2);
        sm(ind(i))=fit(ind(i));
    end
    sm1(ind(1):ind(end))=interp1(ind,sm(ind),ind(1):ind(end));
%     LLR=(-.5*(sum((y-fit).^2)-sum((sm1-y).^2))/var(y-fit))/length(y);
end

function est = stepfit1(data,varargin)
global noise resolution pass b a envelopeflag y endbounding
global Fs tau outputnoise measnoise passes verbose
y = data(:)';
try
    if(tau ~= 0)
        sysc = tf(1,[tau 1]);
        sysd = c2d(sysc,1/Fs,'zoh');
        [b,a] = tfdata(sysd,'v');
    else
        b = [1];
        a = [1 0];
    end
catch
    %     disp('Warning: Control toolbox required: Using unoptimized discretization for dynamics representation!');
    b = [ 1/Fs/tau];
    a = [1 1/Fs/tau-1];
end%find dynamics coefficients

if(verbose)
    disp(['b: [' num2str(b) ']'])
    disp(['a: [' num2str(a) ']'])
end

[h,w] = freqz(b,a,1000);
noiseamp = sqrt(sum((abs(h)).^2)/length(h));
noise = max(0,(outputnoise-measnoise)/noiseamp); % this is thermal noise SD
est = y;
oldest = y;
oldnorm = 0;
envelopeflag = 0;
previousenvelopeflag = 0;
endbounding = 0;
for pass = 1:passes
    est = vitpass(y,est);
    steps = diff(est);
    if(verbose)
        disp([', ' num2str(sum(logical(steps))) ' steps']);
    end
    newnorm = norm(oldest-est,2);
    if(newnorm == 0)
        errchck = 0;
    else
        errchck = abs(newnorm-oldnorm)/oldnorm;
    end
    if(resolution<outputnoise/50);errchck = 0;end
    %     disp([newnorm oldnorm errchck])
    if(sum([newnorm oldnorm errchck]) == 0 && envelopeflag)% || logical(previousenvelopeflag-envelopeflag))
        %         disp('Histogram Converged: END');
        %                 break;
        %         return;
    end
    oldnorm = newnorm;
    oldest = est;
    endbounding = max(endbounding,previousenvelopeflag-envelopeflag);
    previousenvelopeflag = envelopeflag;
    %     disp([envelopeflag endbounding])
end
% hgsave('histogramiter');
% est(2:end) = est(1:end-1);
% est = meanify(est,y);
if(nargout>1)
    stepfitconfidence(y,est);
end
%---------------------------------------------------------------------
function est = vitpass(y,estin)
global userargs outputnoise measnoise pass wt noise envelopeflag
global verbose resolution lowsave highsave trackcost b a endbounding
% partition data into blocks and run viterbi on each block
% global trackcost
%% windowsize
N = length(y);
imp = zeros(1,1000);
imp(1) = 1;
impf = filter(b,a,imp);
bw = find(impf>0.2*max(impf),1,'last');
window = max(10,min(length(y),bw));
%% noiseamp
if(pass >= 1)%noiseamp
    [h,w] = freqz(b,a,10000);
    noiseamp = sqrt(sum((abs(h)).^2)/length(h));
    noise = (outputnoise-measnoise)/noiseamp;
    %     maxamp = max(abs(h));
end
%% ncores
% try
%     ncores = matlabpool('size');
% catch
%     ncores = 1;
% end%ncores
%% Blocksize

steps = diff(estin);
ind = find(steps);
if(pass == 1)
    meandwell = max(1000,mean(floor(length(estin)/2)));
else
    meandwell = max(1000,round(mean(diff(ind))));
end
blocksize = max(10000,max(10*meandwell,ceil(N/4/ceil(N/10/meandwell))));
r = floor((N-1)/blocksize);
overlap = max(.1*blocksize,2*meandwell);
est = zeros(1,length(y));
miny = min(y);
% maxy = max(y);
%% noise reestimation
if(pass>1)%outputnoise
    %     steps = diff(estin);
    %     steplocations = find(steps);
    %     dwelltime = diff(steplocations);
    %     meandwelltime = mean(dwelltime);
    outputnoise = noise_std(y-filter(b,a,estin)); %re-estimate noise
end
%% envelope
if(pass == 1)%envelope
    [low high maxstep minstep miny maxy] = findenvelope2(estin,window,resolution);
    
    lowsave = low;
    highsave = high;
    resolution = mean(high-low)/100;
else
    %     [low high maxstep minstep miny maxy] = findenvelope(estin,y,window,resolution,outputnoise);
    [low high maxstep minstep miny maxy resolution ] = findenvelope3(estin,y,window,resolution,outputnoise,lowsave,highsave);
    % [low high maxstep minstep miny maxy res] = findenvelope4(estin,resolution,lowsave,6);
    lowsave = low;
    highsave = high;
    % resolution = min(res,resolution);
end

low = round((low-miny)/resolution)+1;
high = round((high-miny)/resolution)+1;
nstates = max(high-low+1);
% trackcost = inf*ones(nstates,length(low));
% if(endbounding && ~envelopeflag)
%     est = estin;
%     return;
% end
message = ['Pass:' num2str(pass) ', Resolution = ' sprintf('%4.3f',resolution) ];
if(verbose)
    fprintf([message]);
end
%% compute weight

if(pass <= 1 )
    range = minstep-1:resolution:maxstep;
    [val, in] = min(abs(range));
    range = range-range(in);
    wt.zeroindex = in;
%     wt.pdf = abs(sign(range));
    wt.pdf = (abs(range));
    wt.pdf = 1 + (abs(range));
%     wt.pdf = wt.pdf/max(range);
    wt.pdf(wt.zeroindex) = 0;
else
    wt = computeweight(estin,minstep,maxstep);
    % wt = computeweight2(estin,resolution,adapt,smooth,minstep,maxstep,trackcost
    % ,tracklow);
end%compute weight
start = zeros(1,r+1);
stop = zeros(1,r+1);
eststart = zeros(1,r+1);
eststop = zeros(1,r+1);

%% parallel operation
msg = char(zeros(r+1,1)+'|');
fprintf(msg);
for j = 1:r+1
    i = j-1;
    start(j) = max(blocksize*i+1-overlap,1);
    stop(j) = min(blocksize*(i+1)+overlap,N);
    eststart(j) = blocksize*i+1;
    eststop(j) = min(N,blocksize*(i+1));
    yme{j} = y(start(j):stop(j));
    lowtmp{j} = low(start(j):stop(j));
    hightmp{j} = high(start(j):stop(j));
    %     trackcosttmp{j} = trackcost(:,start(j):stop(j));
end %%% slice data for parallel operation
param.noise = noise;
param.measnoise = measnoise;
param.outputnoise = outputnoise;
param.wt = wt;
param.resolution = resolution;
param.b = b;
param.a = a;
param.pass = pass;
envelope = envelopeflag;
for j = 1:r+1  %*****************parfor**************
    [tmp costtmp] = viterbistepdetector(yme{j},miny,lowtmp{j},hightmp{j},param,envelope);
    myest{j} = tmp;
end

for j = 1:r+1
    tempest = myest{j};
    est(eststart(j):eststop(j)) = tempest(eststart(j)-start(j)+1:eststop(j)-start(j)+1);
end  %% merge sliced data

%% display progress
if(verbose)
    figure(4);clf;
    %         nr = floor(sqrt(passes));
    %         nc = ceil(passes/nr);
    subplot(2,1,1);
    v = 1:N;
    rv = N:-1:1;
    polyx = [v rv];
    polyenv = [(low-1)*resolution+miny (fliplr(high)-1)*resolution+miny];
    polyy = [y rv*0];
    py = fill(polyx,polyy,'k');hold on;set(py,'facealpha',0,'edgealpha',0.1);
    pf = fill(polyx,polyenv,'r');set(pf,'facealpha',0.2,'edgealpha',0);hold on;
    plot(est,'r','linewidth',2);
    xlabel('Samples');
    ylabel('Position');
    %     legend('Data','Envelope','Fit','location','best');
    axis tight;
    
    subplot(2,1,2);
    cla
    x = resolution*((1:length(wt.pdf))-wt.zeroindex);
    smoothedhistogram = exp(-wt.pdf);
    smoothedhistogram = smoothedhistogram-min(smoothedhistogram);
    smoothedhistogram(wt.zeroindex) = 0;
    %     Area(x-resolution/2,smoothedhistogram,'facecolor','r','linewidth',1);
    [bins nhist] = stephist(est,linspace(min(x),max(x),25),gcf,{'blue'},'normalized','bar');
    if(~any(nhist))
        return
    end
    if(max(nhist))
        bar(bins,nhist*max(smoothedhistogram)/max(nhist),'b');
    end
    hold on;
    plot(x-resolution/2,smoothedhistogram,'r','linewidth',1);
    xlabel('Step-size');
    ylabel('Probability');
    axis tight;
    drawnow();
    if(verbose == 2)
        str = sprintf('StepFitIterationProgress %d .fig',pass);
        hgsave(str);
    end
end

%% compute weight

% if(isequal(b,[0 1]) && isequal(a,[1 0]))
%     est = meanify(est,y);
% end
% est = est(1:N);
% wt = computeweight(est,minstep,maxstep);
%---------------------------------------------------------------------
function wt = computeweight(est,minstep,maxstep)
global resolution outputnoise noise pass b a y passes
% smooth = 0.9/(1-passes);
% smooth = max(2,round(1*noise/(resolution)*min(2,pass*smooth+1-smooth)));

redfactor = 1;%1-0.8*pass/passes;
smooth = max(3,round(redfactor*noise/resolution));
steps = diff(est);
s = minstep-2*smooth*resolution:resolution:maxstep+2*smooth*resolution;
[s0 zindex] = min(abs(s));
s = s-s(zindex);
% add penalty to size of step
scale = abs(s);
scale = -scale;
scale = scale - min(scale);
scale = scale/max(scale);
% scale = 1;
N = histc(steps,s-resolution/2).*scale;

if(smooth)
    t = N(zindex);
    Ntmp = N;
    Ntmp(zindex) = 0;
    Ntmp(zindex+1) = 0 ;
    Ntmp(zindex-1) = 0 ;
    d = smooth;
    H = gaussfir(1/d,d,1);
    % keyboard();
    % Ntmp = [Ntmp zeros(1,2*d+1)];
    % Ntmp = filter(H,1,Ntmp);
    % Ntmp(1:end-d) = Ntmp(d+1:end);
    % Ntmp = Ntmp(1:length(s));
    % Ntmp(zindex) = 0;
    % N = Ntmp;
    % N(zindex) = max(H)*t;
    % keyboard();
    Ntmpright = [Ntmp(zindex+1:end) zeros(1,2*d+1)];
    Ntmpleft = [Ntmp(1:zindex-1) zeros(1,2*d+1)];
    Ntmpright = filter(H,1,Ntmpright);
    Ntmpleft = filter(H,1,Ntmpleft);
    Ntmpright(1:end-d) = Ntmpright(d+1:end);
    Ntmpleft(1:end-d) = Ntmpleft(d+1:end);
    Ntmpleft = Ntmpleft(1:end-2*d-1);
    Ntmpright = Ntmpright(1:end-2*d-1);
    Ntmp = [Ntmpleft 0 Ntmpright];
    s(zindex) = 0;
    %% normalize with step size
    % ar = sum(Ntmp);
    % Ntmp = Ntmp./abs(s);
    % Ntmp = Ntmp/sum(Ntmp)*ar;
    %%
    Ntmp(zindex) = t*max(H);
    % Ntmp = max(N,Ntmp);
    % figure(3);N1 = N;N1(zindex) = 0;Ntmp1 = Ntmp;Ntmp1(zindex) = 0;plot([N1' Ntmp1']);pause(1);
    N = Ntmp;
end
N = N+1e-5;
N(zindex+1) = 0 ;
N(zindex-1) = 0 ;
N = N/sum(N);
N = -log(N);
wt.zeroindex = zindex;
wt.pdf = N;
wt.pdf(zindex) = 0;
%---------------------------------------------------------------------
function [est costtrack] = viterbistepdetector(y,miny,low,high,param,envelopeflag)
noise = param.noise;
% noise = sqrt(movingvar(y,50));
measnoise = param.measnoise;
outputnoise = param.outputnoise;
wt = param.wt;
resolution = param.resolution;
b = param.b;
a = param.a;
pass = param.pass;
N = length(y);
if(pass>1 ) %%sigmabar definition
    sigmabar = 2*(sum(b.^2)*(noise.^2)+sum(a.^2)*(measnoise^2));
%     sigmabar = 9*(noise.^2+measnoise^2+resolution^2);
else
    %     wt.pdf = wt.pdf/max(wt.pdf);
    %     wt.pdf(wt.zeroindex) = 0;
    sigmabar = 9*(noise.^2+measnoise^2+resolution^2);
end
% nstates = ceil((maxy-miny)/resolution)+1;
maxstep = max(high(2:end)-low(1:end-1));
minstep = min(low(2:end)-high(1:end-1));
maxband = max(high-low+1); % spread of envelope
S = (low(1):low(1)+maxband-1)'*ones(1,length(y)); % Initialize memory for surviving states.
cost = ((S(:,1)-1)*resolution+miny-y(1)).^2; % Initialize cost associated with each end state.
levelsup = 0:maxstep;
levelsdown = minstep:1:-1;
levels = max(1,[levelsdown levelsup]'+wt.zeroindex);
levels = min(levels,length(wt.pdf));
pdfcost1 = wt.pdf(levels);
pdfcostoffset1 = length(levelsdown)+1;
zhat = S*0+y(1); % Initialize memory for surviving output.
costtrack = 0;
for t = max(length(a)):length(y)%% Viterbi
    lastrow = (low(t-1):high(t-1))';
    llastrow = length(lastrow);
    row = (low(t):high(t));
    lrow = length(row);
    vr = ones(llastrow,1);
    vc = ones(1,lrow);
    dx = -lastrow*vc+vr*row+pdfcostoffset1;
    z = b(1)*vr*((row-1)*resolution+miny); %************************************
    lastrowtmp = lastrow;
    for i = 2:max(length(a),length(b))
        if(t-i>= 2)
            idx = lastrowtmp-low(t-i+1)+1;
            if(i <= numel(b)); z = z+b(i)*((lastrowtmp-1)*resolution+miny)*vc; end %********************
            if(i <= numel(a))
                if(pass == 1 ); z = z-a(i)*zhat(idx,t-i+1)*vc;else %dont erase
                    z = z-a(i)*y(t-i+1)*vr*vc ;
                end
            end
            lastrowtmp = low(t-i)-1+S(idx,t-i+1);
        end
    end
    if(min(size(dx)) == 1)
        noisecost = reshape(pdfcost1(dx),size(dx))*sigmabar;
    else
        noisecost = pdfcost1(dx)*sigmabar;
    end
    measurementcost = (y(t)-z).^2;
    tempcost = cost(1:llastrow)*vc+measurementcost+noisecost;
    [mincost,ind] = min(tempcost);     %select candidate with least cost
    zhat(1:lrow,t) = z(ind+(0:length(ind)-1)*llastrow)';
    cost(1:lrow) = mincost';
    S(1:lrow,t) = ind;
end
[mincost,ind] = min(cost(1:high(end)-low(end)+1));
est = reconstruction(S,low,ind,resolution,miny);
%---------------------------------------------------------------------
function [low high maxstep minstep miny maxy] = findenvelope2(signal,window,resolution)
% global trackcost
low = minfilter(signal,window);
high = maxfilter(signal,window);
mid = (low+high)/2;
low = 2*low-mid-2;
high = 2*high-mid+2;
maxstep = max(high-low);
minstep = min(low-high);
miny = min(low);
maxy = max(high);
%---------------------------------------------------------------------
function [low high maxstep minstep miny maxy resolution] = findenvelope3(signal,y,window,resolution,noise,low,high)
global envelopeflag tau Fs endbounding pass topthickness bottomthickness
if(pass <= 2)
    topthickness = max(5*resolution,mean(high-signal));
    bottomthickness = max(5*resolution,mean(signal-low));
end
if(~(min(abs(high-signal)) <= resolution || min(abs(low-signal)) <= resolution))
    envelopeflag = 1;
    topthickness = topthickness/2;
    bottomthickness = bottomthickness/2;
    low = signal-bottomthickness;
    high = signal+topthickness;
else
    envelopeflag = 0;
    %     disp('Envelope too narrow');
    topper = topthickness*(logical(abs(high-signal) <= 6*resolution));
    lower = bottomthickness*(logical(abs(signal-low) <= 6*resolution));
    low = low-lower;
    high = high+topper;
end
low = min(low,minfilter(signal-2*resolution,window));
high = max(high,maxfilter(signal+2*resolution,window));
finelevel = max(100,1);
resolution = min(resolution,.5*(min(high-low)/10 +max(high-low)/finelevel));

high = max(low+2*resolution,high);
maxstep = max(high(2:end)-low(1:end-1));
minstep = min(low(2:end)-high(1:end-1));
miny = min(low);
maxy = max(high);
function [low high maxstep minstep miny maxy resolution] = findenvelope4(signal,y,window,resolution,noise,low,high)
global envelopeflag tau Fs endbounding pass topthickness bottomthickness
if(pass <= 2)
    topthickness = min(10*resolution,mean(high-signal));
    bottomthickness = min(10*resolution,mean(signal-low));
end
if(~(min(abs(high-signal)) <= resolution || min(abs(low-signal)) <= resolution))
    envelopeflag = 1;
    if(1)%~endbounding)
        topthickness = topthickness/2;
        bottomthickness = bottomthickness/2;
        low = signal-bottomthickness;
        high = signal+topthickness;
    end
else
    envelopeflag = 0;
    %     disp('Envelope too narrow');
    if(~endbounding)
        topper = topthickness*(logical(abs(high-signal) <= 4*resolution));
        lower = bottomthickness*(logical(abs(signal-low) <= 4*resolution));
    else
        topper = topthickness+5*resolution;
        lower = bottomthickness+5*resolution;
    end
    low = low-lower;
    high = high+topper;
end
if(~(envelopeflag && endbounding))
    
    low = minfilter(low,window);
    high = maxfilter(high,window);
    finelevel = max(50,1);
    % if(envelopeflag == 1)
    resolution = min(resolution,.5*(min(high-low)/10 +max(high-low)/finelevel));
    % end
end
% if(endbounding);envelopeflag = 1;end
maxstep = max(high-low);
minstep = min(low-high);
miny = min(low);
maxy = max(high);
%---------------------------------------------------------------------
function y = maxfilter(x,window)
N = length(x);
y = x;
for i = 1:N
    start = max(1,i-window);
    stop = min(i+window,N);
    y(i) = max(x(start:stop));
end
%---------------------------------------------------------------------
function y = minfilter(x,window)
N = length(x);
y = x;
for i = 1:N
    start = max(1,i-window);
    stop = min(i+window,N);
    y(i) = min(x(start:stop));
end
%---------------------------------------------------------------------
function est = meanify(est,y)
% disp('meanifying');
steplocations = [0 find(diff(est)) length(est)];
for i = 1:length(steplocations)-2
    ind = steplocations(i)+1:steplocations(i+1);
    est(ind) = mean(y(ind));
end
%---------------------------------------------------------------------
function noiseD = noise_std(y)
p = abs(fft(y));
for j = 1:9
    ind = p > (mean(p) + 3*std(p));
    p(ind) = mean(p)+1*std(p);
end
noiseD = 1.1*sqrt(sum(p.^2)/length(p)^2);

%---------------------------------------------------------------------
function [v,s] = movingvar(x,m)
offset = floor(m/2);
x = x(:);
s = movingmean(x,m);
N = length(x);
npts = floor(N/m);
locations = round(linspace(1,N,npts));
xt = spline(locations,s(locations),1:N);
xt = xt(:);
x = x-xt;
x = x.^2;
b = ones(1,m)/m;
a = 1;
v = filter(b,a,x);
v(offset+1:end-m+offset) = v(m+1:end);
v(1:offset) = v(offset+1);
v(end-m+offset+1:end) = v(end-m+offset);
%---------------------------------------------------------------------
function xm = movingmean(x,m)
x = x(:);
b = ones(1,m)/m;
a = 1;
xm = filter(b,a,x);
%---------------------------------------------------------------------
function [y i] = lmax(x)
N = length(x);
s = zeros(1,N);
der = diff(x);
der = sign(der);
i = find(diff(der)<0);
i = i+1;
y = x(i);
%---------------------------------------------------------------------
function c = stepfitconfidence(y,fit)
sm = fit;
sm1 = fit;
steplocations = [0 find(diff(fit)) length(fit)];
for i = 1:length(steplocations)-1
    ind(i) = round((steplocations(i)+1+steplocations(i+1))/2);
    sm(ind(i)) = fit(ind(i));
end
sm1(ind(1):ind(end)) = interp1(ind,sm(ind),ind(1):ind(end));
c = (-.5*(sum((y-fit).^2)-sum((sm1-y).^2))/var(y-fit))/length(y);

function est = reconstruction(S,low,lastindex,resolution,miny)
N = size(S,2);
est(N) = (low(N)+lastindex-1-1)*resolution+miny;
ind = lastindex;
for k = N:-1:2
    ind = S(ind,k);
    est(k-1) = (ind+low(k-1)-1-1)*resolution+miny;
end

function LLR = stepfitLLR(y,fit)
y = y(:);
fit = fit(:);
sm = fit;
sm1 = fit;
steplocations = [0 ;find(diff(fit)) ;length(fit)];
for i = 1:length(steplocations)-1
    ind(i) = round((steplocations(i)+1+steplocations(i+1))/2);
    sm(ind(i)) = fit(ind(i));
end
sm1(ind(1):ind(end)) = interp1(ind,sm(ind),ind(1):ind(end));
LLR = (-.5*(sum((y-fit).^2)-sum((sm1-y).^2))/var(y-fit))/length(y);


function [bincenters Nhist] = stephist(stepsignal,bincenters,fignum,col,varargin)
s = size(stepsignal);
l = max(s);
m = min(s);
temp = reshape(stepsignal,l,m);
tempcol = col;
bw = bincenters(2)-bincenters(1);
figure(fignum);
for j = 1:m
    stepsignal = temp(:,j);
    steps = diff(stepsignal);
    steps = steps(logical(steps));
    Nhist(:,j) = histc(steps,bincenters-(bincenters(2)-bincenters(1))/2);
    
    if(nargin>4)
        if(strcmp(varargin{1},'normalized') || strcmp(varargin{1},'Normalized'))
            Nhist(:,j) = Nhist(:,j)/sum(Nhist(:,j));
        elseif(strcmp(varargin{1},'smooth') || strcmp(varargin{1},'Smooth'))
            Nhist(:,j) = gaussfir_unbiased(Nhist(:,j),5);
            
        end
    end
end
