function est=stepfit1(data,varargin)
%    
%USAGE:
%x=stepfit1(y,'parameter',value)
%
%x           Step Fit
%y           Data on which to fit steps
%Optional Parameters:
%Fs                 Sampling time. Required if tau is being specified.
%tau                Time constant of first order dynamics 
%outputnoise Standard deviation of the noise measured at output, Estimates noise at the output if unspecified.
%                      Noise estimation is not reliable under extreme dynamical distortion and high stepping speed.
%measnoise   Standard deviation of the noise in the measurement system.
%                       Assumes zero measurement noise if unspecifed
%passes          Number of iterations in histogram refinement
%verbose         verbose=1: Display progress of histogram.  verbose=2: Save the progress in fig files

numarg=nargin-1;
Fs=1;
tau=0;
evalc('outputnoise=noise_std(data)');
measnoise=0;
passes=10;
verbose=0;
for i=1:2:numarg
    msg=[varargin{i} '=' num2str(varargin{i+1})];
    evalc(msg);
%     disp(msg);
end
msg=sprintf('Fs=%.2e \n tau= %.2e \n outputnoise=%.2e \n measnoise=%.2e \n passes= %.2e \n verbose=%.2e',Fs,tau,outputnoise,measnoise,passes,verbose);
disp(msg);
est=stepfit_sub(data,Fs,tau,outputnoise,measnoise,passes,verbose);
LLR=stepfitLLR(data,est);
fprintf('Stepping LLR=%.3f',LLR)
function est=stepfit_sub(data,varargin)
global userargs Fs tau noise outputnoise measnoise resolution
global pass passes verbose b a envelopeflag y endbounding
y=data(:)';
for i=1:6
    userargs{i}=[];
end
for i=1:nargin-1
    userargs{i}=varargin{i};
end
Fs=1;
tau=0;
evalc('outputnoise=noise_std(y)');
measnoise=0;
resolution=.1;
passes=10;
verbose=0;
if(nargin<1);disp('Atleast 1 argument required');return;end   
if(nargin>1);if(~isempty(varargin{1}));Fs=varargin{1};end;end        
if(nargin>2);if(~isempty(varargin{2}));tau=varargin{2};end;end   
if(nargin>3);if(~isempty(varargin{3}));outputnoise=varargin{3};end;end
if(nargin>4);if(~isempty(varargin{4}));measnoise=varargin{4};end;end 
if(nargin>5);if(~isempty(varargin{5}));passes=varargin{5};end;end   
if(nargin>6);if(~isempty(varargin{6}));verbose=varargin{6};end;end   
try%check for parallel processing capability
    s=matlabpool('size');
    if(s<=1)
        disp('starting parallel MATLAB sessions...')
        evalc('matlabpool open');
        disp([char(8) 'Done '])
    end
catch  %#ok<*CTCH>
    disp('Parallel Processing toolbox not found');
end
try  %find coefficients of dynamical model
    if(tau~=0)
        sysc=tf(1,[tau 1]);
        sysd=c2d(sysc,1/Fs,'zoh');
        [b,a]=tfdata(sysd,'v');
    else
        b=1;
        a=[1 0];
    end
catch 
    disp('Warning:Control Systems Toolbox not found: Using unoptimized discretization for dynamics representation!');
    b=1/Fs/tau;
    a=[1 1/Fs/tau-1];
end
if(verbose)
    disp(['b: [' num2str(b) ']'])
    disp(['a: [' num2str(a) ']'])
end
[h,w]=freqz(b,a,1000);
noiseamp=sqrt(sum((abs(h)).^2)/length(h));
noise=max(0,(outputnoise-measnoise)/noiseamp); % this is thermal noise SD
est=y;
envelopeflag=0;
endbounding=0;
for pass=1:passes
    est=vitpass(y,est);
    steps=diff(est);
    if(verbose)
        disp( [char(8) ', ' num2str(sum(logical(steps))) ' steps']);
    end
end
if(nargout>1)
    stepfitconfidence(y,est);
end
%---------------------------------------------------------------------
function est=vitpass(y,estin)
global userargs outputnoise measnoise pass wt noise 
global verbose resolution lowsave highsave b a
% partition data into blocks and run viterbi on each block
% global trackcost
%% windowsize
N=length(y);
imp=zeros(1,1000);
imp(1)=1;
impf=filter(b,a,imp);
bw=find(impf>0.2*max(impf),1,'last');
window=max(10,min(length(y),bw));
%% noiseamp
if(pass>=1)%noiseamp
    [h,w]=freqz(b,a,10000);
    noiseamp=sqrt(sum((abs(h)).^2)/length(h));
    noise=(outputnoise-measnoise)/noiseamp;
end
%% ncores
try
    ncores=matlabpool('size');
catch
    ncores=1;
end%ncores
%% Blocksize

steps=diff(estin);
ind=find(steps);
if(pass==1)
    meandwell=500;
else
    meandwell=max(500,round(mean(diff(ind))));
end
blocksize=max(10*meandwell,ceil(N/ncores/ceil(N/10/meandwell)));
r=floor((N-1)/blocksize);
overlap=max(.1*blocksize,2*meandwell);
est=zeros(1,length(y));
%% noise reestimation
if(pass>1 && isempty(userargs{3}))%outputnoise 
    outputnoise=noise_std(y-filter(b,a,estin)); %re-estimate noise
end
%% envelope
if(pass==1)%envelope
    [low high maxstep minstep miny maxy]=findenvelope2(estin,window);
    lowsave=low;
    highsave=high;
    resolution=mean(high-low)/50; % refine resolution 
else
[low high maxstep minstep miny maxy resolution ]=findenvelope3(estin,window,resolution,lowsave,highsave);
lowsave=low;
highsave=high;
end
low=round((low-miny)/resolution)+1;
high=round((high-miny)/resolution)+1;
message=['Pass:' num2str(pass) ', Resolution=' sprintf('%4.3f',resolution) ];
if(verbose)
    fprintf(message);
end
%% compute weight 

if(pass<=1 )
    range=minstep-1:resolution:maxstep;
    [val in]=min(abs(range));
    range=range-range(in);
    wt.zeroindex=in;
    wt.pdf=abs(sign(range));
else
   wt=computeweight(estin,minstep,maxstep);
end%compute weight
start=zeros(1,r+1);
stop=zeros(1,r+1);
eststart=zeros(1,r+1);
eststop=zeros(1,r+1);

%% parallel operation
if(verbose)
    msg=[char(zeros(1,r+1)+'|')  ' '];
    fprintf(msg);
end
for j=1:r+1
    i=j-1;
    start(j)=max(blocksize*i+1-overlap,1);
    stop(j)=min(blocksize*(i+1)+overlap,N);
    eststart(j)=blocksize*i+1;
    eststop(j)=min(N,blocksize*(i+1));
    yme{j}=y(start(j):stop(j));
    lowtmp{j}=low(start(j):stop(j));
    hightmp{j}=high(start(j):stop(j));
end %%% slice data for parallel operation
param.noise=noise;
param.measnoise=measnoise;
param.outputnoise=outputnoise;
param.wt=wt;
param.resolution=resolution;
param.b=b;
param.a=a;
param.pass=pass;
parfor j=1:r+1  %*****************parfor**************
    [tmp costtmp]=viterbistepdetector(yme{j},miny,lowtmp{j},hightmp{j},param);
    myest{j}=tmp;
    if(verbose)
        disp([char(8) char(8)  ]);
    end
end
disp([char(8)]);
for j=1:r+1
    tempest=myest{j};
    est(eststart(j):eststop(j))=tempest(eststart(j)-start(j)+1:eststop(j)-start(j)+1);
end  %% merge sliced data 

%% display progress
if(verbose)
    figure(4);clf;
    %         nr=floor(sqrt(passes));
    %         nc=ceil(passes/nr);
    subplot(2,1,1);
    v=1:N;
    rv=N:-1:1;
    polyx=[v rv];
    polyenv=[(low-1)*resolution+miny (fliplr(high)-1)*resolution+miny];
    polyy=[y rv*0];
    py=fill(polyx,polyy,'k');hold on;set(py,'facealpha',0,'edgealpha',0.1);
    pf=fill(polyx,polyenv,'r');set(pf,'facealpha',0.2,'edgealpha',0);hold on;
    plot(est,'r','linewidth',2);
    xlabel('Samples');
    ylabel('Position');
%     legend('Data','Envelope','Fit','location','best');
    axis tight;
    
    subplot(2,1,2);
    cla
    x=resolution*((1:length(wt.pdf))-wt.zeroindex);
    smoothedhistogram=exp(-wt.pdf);
    smoothedhistogram=smoothedhistogram-min(smoothedhistogram);
    smoothedhistogram(wt.zeroindex)=0;        
    [bins nhist]=stephist(est,linspace(min(x),max(x),25),gcf,'normalized','bar');    
    if(max(nhist))
        bar(bins,nhist*max(smoothedhistogram)/max(nhist),'b');
    end        
    hold on;
    plot(x-resolution/2,smoothedhistogram,'r','linewidth',1);
    xlabel('Step-size');       
    ylabel('Probability');
    axis tight;
    drawnow();
    if(verbose==2)
        str=sprintf('StepFitIterationProgress %d .fig',pass);
        hgsave(str);
    end
end

function wt=computeweight(est,minstep,maxstep)
global resolution noise pass passes
redfactor=1-0.8*pass/passes;
smooth=max(3,round(redfactor*noise/resolution));
steps=diff(est);
s=minstep-2*smooth*resolution:resolution:maxstep+2*smooth*resolution;
[s0 zindex]=min(abs(s));
s=s-s(zindex);
N=histc(steps,s-resolution/2);

if(smooth)
t=N(zindex);
Ntmp=N;
Ntmp(zindex)=0;
Ntmp(zindex+1)=0 ;
Ntmp(zindex-1)=0 ;
d=smooth;
H=gaussfir(1/d,d,1);
Ntmpright=[Ntmp(zindex+1:end) zeros(1,2*d+1)];
Ntmpleft=[Ntmp(1:zindex-1) zeros(1,2*d+1)];
Ntmpright=filter(H,1,Ntmpright);
Ntmpleft=filter(H,1,Ntmpleft);
Ntmpright(1:end-d)=Ntmpright(d+1:end);
Ntmpleft(1:end-d)=Ntmpleft(d+1:end);
Ntmpleft=Ntmpleft(1:end-2*d-1);
Ntmpright=Ntmpright(1:end-2*d-1);
Ntmp=[Ntmpleft 0 Ntmpright];
Ntmp(zindex)=t*max(H);
N=Ntmp+1e-50;
end
N(zindex+1)=0 ;
N(zindex-1)=0 ;
N=N/sum(N);
N=-log(N);
wt.zeroindex=zindex;
wt.pdf=N;
%---------------------------------------------------------------------
function [est costtrack]=viterbistepdetector(y,miny,low,high,param)
noise=param.noise;
measnoise=param.measnoise;
wt=param.wt;
resolution=param.resolution;
b=param.b;
a=param.a;
pass=param.pass;
if(pass>1 ); %%sigmabar definition
    sigmabar=2*(sum(b.^2)*(noise.^2)+sum(a.^2)*(measnoise^2));
else
    sigmabar=9*(noise.^2+measnoise^2+resolution^2);
end
maxstep=max(high(2:end)-low(1:end-1));
minstep=min(low(2:end)-high(1:end-1));
maxband=max(high-low+1); % spread of envelope
S=(low(1):low(1)+maxband-1)'*ones(1,length(y)); % Initialize memory for surviving states.
cost=((S(:,1)-1)*resolution+miny-y(1)).^2; % Initialize cost associated with each end state.
levelsup=0:maxstep;
levelsdown=minstep:1:-1;
levels=max(1,[levelsdown levelsup]'+wt.zeroindex);
levels=min(levels,length(wt.pdf));
pdfcost1=wt.pdf(levels);
pdfcostoffset1=length(levelsdown)+1;
zhat=S*0+y(1); % Initialize memory for surviving output.
costtrack=0;
for t=max(length(a)):length(y)%% Viterbi
    lastrow=(low(t-1):high(t-1))';
    llastrow=length(lastrow);
    row=(low(t):high(t));
    lrow=length(row);
    vr=ones(llastrow,1);
    vc=ones(1,lrow);
    dx=-lastrow*vc+vr*row+pdfcostoffset1;
    z=b(1)*vr*((row-1)*resolution+miny); %************************************
    lastrowtmp=lastrow;
    for i=2:max(length(a),length(b))
        if(t-i>=2)
            idx=lastrowtmp-low(t-i+1)+1;
            if(i<=numel(b)); z=z+b(i)*((lastrowtmp-1)*resolution+miny)*vc; end %********************
            if(i<=numel(a))
                if(pass==1 ); 
                    z=z-a(i)*zhat(idx,t-i+1)*vc;
                else
                    z=z-a(i)*y(t-i+1)*vr*vc ;
                end
            end
            lastrowtmp=low(t-i)-1+S(idx,t-i+1);
        end
    end
    if(min(size(dx))==1)
        noisecost=reshape(pdfcost1(dx),size(dx))*sigmabar;
    else
        noisecost=pdfcost1(dx)*sigmabar;
    end
    measurementcost=(y(t)-z).^2;
    tempcost= cost(1:llastrow)*vc+measurementcost+noisecost;
    [mincost,ind]=min(tempcost);     %select candidate with least cost
    zhat(1:lrow,t)=z(ind+(0:length(ind)-1)*llastrow)';
    cost(1:lrow)=mincost';
    S(1:lrow,t)=ind;
end
[mincost,ind]=min(cost(1:high(end)-low(end)+1));
est=reconstruction(S,low,ind,resolution,miny);

function [low high maxstep minstep miny maxy]=findenvelope2(signal,window)

low=minfilter(signal,window);
high=maxfilter(signal,window);
mid=(low+high)/2;
low=2*low-mid-2;
high=2*high-mid+2;
maxstep=max(high-low);
minstep=min(low-high);
miny=min(low);
maxy=max(high);
%---------------------------------------------------------------------
function [low high maxstep minstep miny maxy resolution]=findenvelope3(signal,window,resolution,low,high)
global envelopeflag pass topthickness bottomthickness
if(pass<=2)
topthickness=max(5*resolution,mean(high-signal));
bottomthickness=max(5*resolution,mean(signal-low));
end
if(~(min(abs(high-signal))<=resolution || min(abs(low-signal))<=resolution))
    envelopeflag=1;
        topthickness=topthickness/2;
        bottomthickness=bottomthickness/2;
        low=signal-bottomthickness;
high=signal+topthickness;
else
    envelopeflag=0;
%     disp('Envelope too narrow');
    topper=topthickness*(logical(abs(high-signal)<=6*resolution));
    lower=bottomthickness*(logical(abs(signal-low)<=6*resolution));
    low=low-lower;
    high=high+topper;
end
low=min(low,minfilter(signal-2*resolution,window));
high=max(high,maxfilter(signal+2*resolution,window));
finelevel=max(100,1);
resolution=min(resolution,.5*(min(high-low)/10 +max(high-low)/finelevel));

high=max(low+2*resolution,high);
maxstep=max(high(2:end)-low(1:end-1));
minstep=min(low(2:end)-high(1:end-1));
miny=min(low);
maxy=max(high);
%---------------------------------------------------------------------
function y=maxfilter(x,window)
N=length(x);
y=x;
for i=1:N
    start=max(1,i-window);
    stop=min(i+window,N);
    y(i)=max(x(start:stop));
end
%---------------------------------------------------------------------
function y=minfilter(x,window)
N=length(x);
y=x;
for i=1:N
    start=max(1,i-window);
    stop=min(i+window,N);
    y(i)=min(x(start:stop));
end
%---------------------------------------------------------------------
function noiseD=noise_std(y)
y=y(1:min(length(y),10000));
N=floor(length(y)/10);
w=floor(logspace(log10(1),log10(N),50)');
w=unique(w);
m=zeros(length(w),1);
for j=1
    for i=1:length(w);
        v=movingvar(y,w(i)); %compute average noise variance for different window lengths
        m(i,j)=mean(v(w+1:end-w));
    end;
    st=sqrt(mean(m,2));
end;

dst=abs(diff(st)); %compute derivative of variance with respect to window size
[b,a]=butter(1,.6);
dst1=.5*filter(b,a,dst,dst(1))+.5*flipud(filter(b,a,flipud(dst),dst(end)));
dst1=.5*filter(b,a,dst1,dst1(1))+.5*flipud(filter(b,a,flipud(dst1),dst1(end)));
dst1=.5*filter(b,a,dst1,dst1(1))+.5*flipud(filter(b,a,flipud(dst1),dst1(end)));
% [val,in]=min(dst1);%find window size for which corresponds to minimum variation in noise
[val,in]=lmax(-dst1);
noiseD=min(st(in+1));
%---------------------------------------------------------------------
function [v,s]=movingvar(x,m)
offset=floor(m/2);
x=x(:);
s=movingmean(x,m);
N=length(x);
npts=floor(N/m);
locations=round(linspace(1,N,npts));
xt=spline(locations,s(locations),1:N);
xt=xt(:);
x=x-xt;
x=x.^2;
b=ones(1,m)/m;
a=1;
v=filter(b,a,x);
v(offset+1:end-m+offset)=v(m+1:end);
v(1:offset)=v(offset+1);
v(end-m+offset+1:end)=v(end-m+offset);
%---------------------------------------------------------------------
function xm=movingmean(x,m)
x=x(:);
b=ones(1,m)/m;
a=1;
xm=filter(b,a,x);
%---------------------------------------------------------------------
function [y i]=lmax(x)
der=diff(x);
der=sign(der);
i=find(diff(der)<0);
i=i+1;
y=x(i);
%---------------------------------------------------------------------
function c=stepfitconfidence(y,fit)
sm=fit;
sm1=fit;
steplocations=[0 find(diff(fit)) length(fit)];
ind=zeros(length(steplocations)-1,1);
for i=1:length(steplocations)-1
    ind(i)=round((steplocations(i)+1+steplocations(i+1))/2);
    sm(ind(i))=fit(ind(i));
end
sm1(ind(1):ind(end))=interp1(ind,sm(ind),ind(1):ind(end));
c=(-.5*(sum((y-fit).^2)-sum((sm1-y).^2))/var(y-fit))/length(y);

function est=reconstruction(S,low,lastindex,resolution,miny) 
N=size(S,2);
est(N)=(low(N)+lastindex-1-1)*resolution+miny;
ind=lastindex;
for k=N:-1:2
    ind=S(ind,k);
    est(k-1)=(ind+low(k-1)-1-1)*resolution+miny;    
end

function LLR=stepfitLLR(y,fit) % find likelihood ratio of stepping vs. smooth
y=y(:);
fit=fit(:);
sm=fit;
sm1=fit;
steplocations=[0 ;find(diff(fit)) ;length(fit)];
ind=zeros(length(steplocations)-1,1);
for i=1:length(steplocations)-1
    ind(i)=round((steplocations(i)+1+steplocations(i+1))/2);
    sm(ind(i))=fit(ind(i));
end
sm1(ind(1):ind(end))=interp1(ind,sm(ind),ind(1):ind(end));
LLR=(-.5*(sum((y-fit).^2)-sum((sm1-y).^2))/var(y-fit))/length(y);


function [bincenters Nhist]=stephist(stepsignal,bincenters,fignum,varargin)
s=size(stepsignal);
l=max(s);
m=min(s);
temp=reshape(stepsignal,l,m);
figure(fignum);
Nhist=zeros(length(bincenters),m);
for j=1:m
    stepsignal=temp(:,j);
    steps=diff(stepsignal);
    steps=steps(logical(steps));
    Nhist(:,j)=histc(steps,bincenters-(bincenters(2)-bincenters(1))/2);
    
    if(nargin>4)
        if(strcmp(varargin{1},'normalized') || strcmp(varargin{1},'Normalized'))
            Nhist(:,j)=Nhist(:,j)/sum(Nhist(:,j));
        elseif(strcmp(varargin{1},'smooth') || strcmp(varargin{1},'Smooth'))
            Nhist(:,j)=gaussfir_unbiased(Nhist(:,j),5);

        end
    end
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           