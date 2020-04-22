% 03/31/07 by YJ SUN: use Kerssemaker's stepfinder algorithm to fit individual stepping traces in a concatenated data set that has multiple traces.
% Because the motor fluorescence intensity and even motor located actin's flexibility are different for different molecules,
% the algorithm does not work very well to assign appropriate number of steps to some of the noiser traces by sorting the quality factor throughout the concatenated data set.
% This function is used to call steps_find_ys2.m and steps_evaluate_ys.m to fit each individual trace in a multiple-trace data set.

% Run stepfindercaller.m to load a log.txt, an 2D array which holds the x,y coordinates of the fluorescence centroids and maybe other information such as standard deviation and intensity.
% The program will call the function xy2displacement(data) to convert x,y coordinates to displacement for each molecule.
% The program then will split the concatenated data into individual ones, and display two side by side plots. 
% The one on the left is the displacement plot. The one on the right is the S vs nst plot. Click in the S-vs-nst plot to choose NST for steps_evaluate.m.
% When the S-vs-nst plot shows a clear peak, click on the pick. Otherwise, refer to the displacement plot to choose a reasonable NST. 
% After NST is chosen by clicking, steps_evaluate.m will run and then update left plot with fitted results of the current molecule and the displacement plot of next molecule. 
% The right plot will show the S-vs-NST curve for next molecule, waiting for the click to choose NST. It repeats until the last molecule.
% Input of stepfindercaller.m: log.txt
% Output of stepfindercaller.m: 
% Log_Fit.csv: displacement (distance vs. time)
% Logpolyfitdata.csv: x-y by a 3rd order polynomial fit.
% Logsteps.csv: step dwell time and amplitude etc.
% Logstepsize.csv: displacement of all molecules.
% Logdwelltime.csv: dwell time of all molecules. The incomplete steps, i.e. the first and last steps, are not included.


function [stepfinderfit, polyfitdata, Steps_all]=stepfindercaller(datafromF2DG16_4,numbersteps)

% load data
if nargin==0
   [FileName,PathName] = uigetfile('*.txt','Select the directory and datafile');
   source=sprintf('%s%s',PathName,FileName);
   datafromF2DG16_4=load(source);
   file_name=(FileName(1:length(FileName)-4));
   outpath=PathName;
end

%call the function xy2displacement(data) to convert x,y coordinates to displacement for each molecule.
[data,polyfitdata]=xy2displacement(datafromF2DG16_4);

%find the indexes to split the data in the concatenated data set.
lijstall=Check_Multiple_Data(data);

%reserve space for output
stepfinderfit=zeros(size(data,1),3);
stepfinderfit(:,1:2)=data(:,1:2);
Steps_all=[];
stepsize=[]; dwelltime=[]; 

%figure(1) is for displaying displacement plot and figure(2) is for displaying the S vs NST plot.
figure(1), set(gcf,'position',[0 200 500 450]);
figure(2), set(gcf,'position',[550 200 500 450]);

% Display the S vs NST curve and allow users to click the approriate NST for the molecule
% When the S-vs-nst plot shows a clear peak, click on the pick. Otherwise, refer to the displacement plot to choose a reasonable NST. 
% After NST is chosen by clicking, steps_evaluate.m will run and then update left plot with fitted results of the current molecule and the displacement plot of next molecule. 
% The right plot will show the S-vs-NST curve for next molecule, waiting for the click to choose NST. It repeats until the last molecule.
for i=1:length(lijstall)-1
   singledata=data(lijstall(i)+1:lijstall(i+1),:);
   figure(1); plot(singledata(:,1),singledata(:,2),'bx'); hold on;
   [singledata,indexes,lijst,properties,initval]=Steps_Find_ys2(singledata);
   disp('click on the S vs NST plot to pick NST.')
   [nst,S,but] = ginput(1);
   [Fit,Steps,histossteps,All_Steppedness]=Steps_Evaluate(singledata,indexes,lijst,properties,initval,nst);     %the Steps_Evaluate function is changed from Dummy.
   figure(1); hold off;
   plot(singledata(:,1), Fit,singledata(:,1),singledata(:,2),'ro'); hold on;
   stepfinderfit(lijstall(i)+1:lijstall(i+1),3)=Fit;
   Steps_all=[Steps_all; Steps'];
   stepsizetmp=[]; dwelltimetmp=[];
   j=1;k=0;
   while j<lijstall(i+1)-lijstall(i)
       if Fit(j+1)==Fit(j)
           j=j+1
       else
           k=k+1
           stepsizetmp(k,1)=Fit(j+1)-Fit(j);
           dwelltimetmp(k,1)=j;
           j=j+1;
       end
   end
   stepsize=[stepsize; stepsizetmp];
   dwelltime=[dwelltime; diff(dwelltimetmp)];
end
  
fitname=strcat(outpath, file_name,'_Fit', '.csv');
csvwrite(fitname, stepfinderfit);                

fitname=strcat(outpath, file_name,'polyfitdata', '.csv');
csvwrite(fitname, polyfitdata);                

fitname=strcat(outpath, file_name,'steps', '.csv');
csvwrite(fitname, Steps_all);                

fitname=strcat(outpath, file_name,'stepsize', '.csv');
csvwrite(fitname, stepsize);                

fitname=strcat(outpath, file_name,'dwelltime', '.csv');
csvwrite(fitname, dwelltime);                


%Kerssemaker's subfunction to find the breaking points in glued data sets.
function lijst=Check_Multiple_Data(data)
  %This section determines if, and where, input data is divided in different
  %data segments (Meaning, whether it is data glued together). It returns the
  %indexes where segments stop, a new one starts on the next index. The last point is
  %followed by a zero
     teller=1;
	 ld=length(data);
	 lijst=zeros(ld,1); 
	 for i=2:ld
       if data(i,1)==0                          %datajump
          teller = teller+1;
          lijst(teller)=i-1;
       end
	 end
    lijst=lijst(1:teller+1);
    lijst(1)=0;
    lijst(teller+1)=ld;

