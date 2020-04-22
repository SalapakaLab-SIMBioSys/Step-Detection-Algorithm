%Program StepMaker
%this program creates an array with typical steps in it, drenched in noise

stochastic=1;
samples=500;
noise=5;
growthrate=20/25 %stepsize/window
Step1=20;
SigmaStep1=5;
nw1=Step1/growthrate;
%nw1=10;
Step2=0;
nw2=80;
Step3=0;
nw3=80;

data=zeros(samples,2);
data(:,2)=noise*randn(samples,1);        %gaussian noise;
data(:,1)=0.04*(1:1:samples)';
for i=1:samples
    
   
    if mod(i,nw1)==0 &  stochastic==0
        data(i:samples,2)=data(i:samples,2)+Step1;
    else
      if stochastic==1  
        teken=1;  %sign(randn(1,1));
            step=Step1*ceil(rand(1,1)-(1-1/nw1))+Step2*ceil(rand(1,1)-(1-1/nw2))+Step3*ceil(rand(1,1)-(1-1/nw3));
    		data(i:samples,2)=data(i:samples,2)+teken*step;
            %chance of 1/nw that it is 'step'; otherwise 0;
      end
    end
 
    if stochastic ==2  %make a distribution of steps, Gaussian size sigma=100, constant event chance
            stepstog=Step1+SigmaStep1*(randn(1,1));
            step=stepstog*ceil(rand(1,1)-(1-1/nw1));
    		data(i:samples,2)=data(i:samples,2)+step;
    end
end
plot(data(:,1),data(:,2));
data