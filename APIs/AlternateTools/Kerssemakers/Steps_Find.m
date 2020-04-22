
%This program splits a time+data set in succesive stepfits
%Jacob Kerssemakers, 21-8-2005

%Run first from command line: [data, indexes,lijst,properties,initval]=Steps_Find;
%User actions: 
	%1) Click a most likely noise level in a presented graph of
	%noise level as function of pair distance. (The choice is not affecting the
	%stepfit itself)
    %2) Graphs of a steppedness-function are presented. Choose from this
    %a reasonable number of steps 'Nst' (for a pure step train, this is where the
    %function peaks).
%After that, run from the command line:
%dummy=Steps_Evaluate(data,indexes,lijst,properties,initval,Nst) with Nst
%the chosen number
    %The program determines the stepfit associated with the user choice,
    %collects step info and presents the final fits. 'Steps_Evaluate' can be
    %run many times (if you want to try different stepsizes) without
    %re-running 'Steps-Find'. Saves time.

%output:
%'indexes' contains values found during the splits:
%1=split-round, 2=index,3=stepsize,4=rank1=step*srqt(N)/noise,
%5=%rank2=step*sqrt(N/Var),6=Chi, 7=worstindex,8=worstChi

%'_Stepp' contains values found during ranked re-buildup of the curves:
%1=windowsize, 2=Steppedness, 3=LastStep, 4=Stepnumber,
%5=DerivativeSteppedness, 6=LastStep, 7=AvStep 8 normalized steppedness

function [data, indexes,lijst,properties,initval, Steppedness, selectie]=Steps_Find(data, verbosity);
	%Settings
	%------------------------------------------
	
    initval.showfits=1;          %shows the fitting, 
    initval.makemovie=0 ;       %Make a movie of the fits (also switch 'showfits on')
    initval.ra=20;                %for noise determination: max distance between data points
    initval.condense=0;             %condenses original data with this factor by skipping points
	initval.cutoff=200;             %max aantal steps fitted {to gain time}
	initval.deelfactor=1;           %pre-divisions ; per segment!
	initval.treshold = 3;           %number_of sigmas that a spike is allowed to rise
	initval.standardinput=0;    
    
    % Laptop
    %initval.inpath='C:\Users\siva_\Desktop\Workspace\SDA\Kerssemakers\'%'C:\MATLAB6p5\work\';        
    %initval.outpath='C:\Users\siva_\Desktop\Workspace\SDA\Kerssemakers\Out\'%'C:\MATLAB6p5\work\';
    
    % Desktop
    %initval.inpath = 'C:\Siva\Engineering\Workspace\SDA\Kerssemakers\';
    %initval.outpath = 'C:\Siva\Engineering\Workspace\SDA\Kerssemakers\Out\';
    
    % General
    initval.outpath = strcat(pwd, '\Out\')
    %---------------------------------------------
	
    %Main Program
    %-------------------------------------------
    % Sivaraman - edited this to accept array data
    %[data,initval.filename]=Select_data(initval);
    initval.filename = 'Default'
    %-------------------------------------------
    
    if initval.condense==1, data=Skip_Points(data,initval);,end;
        %Optional: Condense data to gain time
        
	lijst=Check_Multiple_Data(data);
        %This function looks for breaks in the time-axis, which indicate
        %breaks in combined data-sets. These are each treated separately
        %output: list of indexes where these breaks are.
    
    [noizes, properties,data]=Remove_Outliers_Get_Noise(data,initval,lijst);
        %This function removes outlying spikes and determines after that the Gaussian
        %background noise via a nearest-neighbour analysis. To check for
        %low-pass filtering, the analysis is also repeated as a function of
        %datapoint-distance. The result is plotted to allow the user to
        %click a most likely true noise level. If wished, the data is also condensed here by
	    %skipping points.
    [indexes,properties]=Iterative_Splitting(data,lijst,properties, initval);
        %Split Section. Output: a ranked list of potential step locations, ...
        %ordered by their expected or their measured error. A specific fit ...
        %to the whole dataset now just consists of the picking of these locations 
        %up to a certain number. This is user-tuned in Steps_Evaluate' 
        
        %'indexes' contains columns:
        %-----------------------------------
        %iteration round: the split round where this location was found
        %index: the location of the step, expressed as data point number
        %(last point of the left plateau)
        %stepsize: the fitted size of the step
        %step*srqt(N)/noise: A measure for the prominence of the step, equal to the inverted expected error. 
        %step*sqrt(N/Var):   The same value, but based on measured residual variance. 
                %...Also ranks how 'ideal' steps are but is more prone to statistical variation. 
        %minimum Chi: The minimum residual RMS Chi-square after the split.
                %Should for an ideal step approach the noise 
        %----------------------------------------------------------------
    
     if(verbosity == 1)
        Write_files(indexes,noizes,initval);
        %Storage of data
     end
    
    [dummy, Steppedness, selectie, Fit] = Steps_Evaluate(data,indexes,lijst,properties,initval,0, verbosity); % first time: full run
        %This function creates and analyzes a whole series of full step fits. Every fit
        %contains one step location more than the former, with the step
        %locations ranked by their 'prominence' in the 'indexes'-output
    
    
    function Write_files(indexes,noizes, initval);
    
    fitname=strcat(initval.outpath, initval.filename, '_Splitlist', '.csv');
    csvwrite(fitname, indexes);  
    fitname=strcat(initval.outpath, initval.filename, '_Noise', '.csv');
    csvwrite(fitname, noizes);  
    
    
    
    function [indexes,properties]=Iterative_Splitting(data,lijst,properties, initval);
         
        %Split Section. Output: a ranked list of potential step locations, ...
        %ordered by their expected or their measured error. A specific fit ...
        %to the whole dataset now just consists of the picking of these locations 
        %up to a certain number. This is user-tuned in Steps_Evaluate' 
        
        %'indexes' contains columns:
        %-----------------------------------
        %iteration round: the split round where this location was found
        %index: the location of the step, expressed as data point number
        %(last point of the left plateau)
        %stepsize: the fitted size of the step
        %step*srqt(N)/noise: A measure for the prominence of the step, equal to the inverted expected error. 
        %step*sqrt(N/Var):   The same value, but based on measured residual variance. 
                %...Also ranks how 'ideal' steps are but is more prone to statistical variation. 
        %minimum Chi: The minimum residual RMS Chi-square after the split.
                %Should for an ideal step approach the noise 
        %----------------------------------------------------------------
        
                
         lijst_teller=1;
        indexes=[0,0,0,0,0,0];
        while lijst_teller<length(lijst)   
            %1) Set some counters
            ronde=0;
            start=0;
            segment_aantal=0;
            found_last_round=initval.deelfactor-1; 
            
            %2) Pick a data segment (if the data consists of different data
            %sets glued together
            segment_data=data(lijst(lijst_teller)+1:lijst(lijst_teller+1),:);
            N1=length(segment_data);
            segment_indexes=zeros(N1,6);
            properties.growth_range=properties.growth_range+max(segment_data(:,2))-min(segment_data(:,2));

            %3) Make_Pre_locations;  
            %this first pre-split serves only to speed up the fits for very large
            %datasets; If iterations are done on segments with too many
            %steps, the initial splits might not locate true steps
            %(rememmber that the locations are not changed anymore after being found
            %------------------------------------------------------
            tres0=(max(segment_data(:,2))-min(segment_data(:,2)))*sqrt(N1);
            if initval.deelfactor > 1 
                segment_indexes([1:initval.deelfactor+1]',2)=ceil([0,[1:initval.deelfactor]*N1/initval.deelfactor]');
            else
                segment_indexes(1:2,2)=[0,N1]';
                segment_indexes(1:2,4)=[tres0,tres0]';
            end
            sorteer=sort(segment_indexes(1:initval.deelfactor+1,2));
 
            %4) Pre_Split; the used locations are not stored later on 
            %-------------------------------------------------------
            disp('splitting....')
            segment_indexes(1,:)=[0,0,0,tres0,tres0,tres0];
            segment_indexes(initval.deelfactor+2,:)=[0,N1,0,tres0,tres0,tres0];
            for c=1:initval.deelfactor 
                stop=sorteer(c+1);
                %----------------------------------------------------------
                %this section adresses a two-dim array 'segment_data'in a specific segmentand determines the best step-fit there
                spl=Split2(segment_data,start,stop);
  
                if spl~=[0,0,0,0,0]
                    segment_aantal=segment_aantal+1;
                    segment_indexes(c+1,:)=[1,spl(2:5) ,spl(1)];
                    %segment-indexes contains: 
                    %round, index, stepsize, step*srqt(N),%step*sqrt(N/Var), minimum Chi
                end
                start=stop;
            end %
            found_last_round=segment_aantal+2;       %includes begin and end
            start=0;   %resetting of a loop
            segment_aantal=segment_aantal+1;
            sorteer=sort(segment_indexes(1:segment_aantal+1,2));

            %5) Iterative_Split; 
            % Each new series of splitted steps is added to
            %the existing, earlier determined locations. The collection is regularly
            %sorted. The result is a ranked list of locations where steps can be
            %added; 
            %------------------------------------------------------------------
            while found_last_round ~=0   %nested2 continue until no more steps can be fitted
                ronde=ronde+1;
                oudsegment_aantal=segment_aantal;
                for c=1:oudsegment_aantal % nest3, run through the plateaus found in the former round
                    stop=sorteer(c+1);
                    %----------------------------------------------------------
                    %this section adresses a two-dim array 'segment_data'in a specific segment
                    %and determines the best step-fit there
                    spl=Split2(segment_data,start,stop);
                
                    %spl=minimum Chi, index, stepsize,  step*srqt(N), step*sqrt(N/Var)
                    if spl~=[0,0,0,0,0]    %add the results of the split to the collection
                        segment_aantal=segment_aantal+1;
                        segment_indexes(segment_aantal+1,:)=[ronde,spl(2:5),spl(1)]; %ronde, index, stepsize,  rank1, rank2, Chi, worstindex,worstChi
                    end
                    start=stop;
                end   % nest3......all plateaus were split, if possible
                found_last_round=segment_aantal-oudsegment_aantal;
                sorteer=sort(segment_indexes(1:segment_aantal+1,2));
                lijst_teller=lijst_teller;
                start=0;   %resetting of a loop; go back the beginning of the loop
            end  
            disp('all steps found')
            disp('sorting.....')    
            segment_indexes=sort_on_key(segment_indexes(1:segment_aantal+1,:),2); %cuts off initializing and unused parts and sorts on any column
            segment_indexes(:,2)=segment_indexes(:,2)+lijst(lijst_teller); %to change the relative index back to the absolute one
            indexes=[[indexes]',[segment_indexes(2:segment_aantal+1,:)]']' ;
            properties.aantal=properties.aantal+segment_aantal;
            lijst_teller=lijst_teller+1;    
     end 
	
    %6) some cleaning up
     
    indexes=indexes(2:properties.aantal+1,:);  %get the zero out of the system!  
	indexes(:,4)=indexes(:,4)/properties.noise; 
        %convert column 4 (=step*sqrt(N)) to a relative error that can be compared with other datasets with different noise levels 
    indexes=(fliplr((sort_on_key(indexes,1))'))'; 
    

function [data,file_name]=Select_data(initval);
	if initval.standardinput==1
        PathName=initval.inpath;
        %FileName='Simpledata.txt';
        FileName='Fakestepper1000.TXT';
        source=sprintf('%s%s',PathName,FileName);
        data=load(source);
        file_name=(FileName(1:length(FileName)-4));
    else
        [FileName,PathName] = uigetfile('*.txt','Select the directory and file');
	    source=sprintf('%s%s',PathName,FileName);
	    data=load(source);
	    file_name=(FileName(1:length(FileName)-4));
    end
    cd(initval.inpath);

   function lijst=Check_Multiple_Data(data);
    %This section determines if, and where, input data is divided in different
	%data segments (Meaning, whether it is data glued together). It returns the
	%indexes where segments stop, a new one starts on the next index. The last point is
	%followed by a zero
        teller=1;
		ld=length(data);
		lijst=ones(200,1); % maximaal 100 segmenten
        dt=median(data(2:ld,1)-data(1:ld-1,1));
		for i=1:ld-2
            d=abs(data(i+2,1)-data(i+1,1));
            if d>5*dt                          %datajump
                teller = teller+1;
                lijst(teller)=i+1;
            end
		end
		lijst=lijst(1:teller+1);
		lijst(1)=0;
		lijst(teller+1)=ld;

    function data=Skip_Points(data,initval);        
        %If wished, the data is condensed here by
	    %skipping points.
		N0=length(data);
        resampled_data=[];
		tel=0;
			for j=1:N0 
                if mod(j,initval.condense) ==0
                    tel=tel+1;
                    resampled_data(:,tel)=data(j,:)';
                end
            end
		
         data=resampled_data';
        
         
   function [far_neighbour, properties, data]=Remove_Outliers_Get_Noise(data,initval,lijst);
	%Section: properties.noise_Spike_Finder:
	%This program determines noise, then spots the outliers, rmoves them and
	%recalculates the noise. 
		N0=length(data);
        properties.aantal=0;
        properties.growth_range=0; %this parameter collects the total growth or shrinkage of all segments
        properties.N0=length(data);
		hold off;
   
    %1) Determine noise including outliers
        le=length(data);
        neighbour_dif=[];
        n=0;
        for i=1:le-1
            if sum(find(lijst==i-1))==0 & sum(find(lijst==i))==0 & sum(find(lijst==i+1))==0 
                n=n+1;
                neighbour_dif(:,n)=(data(i+1,2)-data(i,2))^2;
            end
        end
        properties.noise=(mean(neighbour_dif))^0.5/2^0.5 ;   
        
    %2) Discard spikes
        spikes=0;
        for i=2:le-1
            if sum(find(lijst==i-1))==0 & sum(find(lijst==i))==0 & sum(find(lijst==i+1))==0  %do not consider border points or next to those
                if data(i,2)-data(i-1,2)>initval.treshold*properties.noise & data(i,2)-data(i+1,2)>initval.treshold*properties.noise
                    data(i,2)=(data(i-1,2)+data(i+1,2))/2; %positive spike
                    spikes=spikes+1; 
                end
                if data(i,2)-data(i-1,2)<-initval.treshold*properties.noise & data(i,2)-data(i+1,2)<-initval.treshold*properties.noise
                    data(i,2)=(data(i-1,2)+data(i+1,2))/2; 
                    spikes=spikes+1;
                end %negative spike
            end
        end
        
        %3) Re-measure noise; also as function of further neighbours (up to
        %30 points away
        lel=length(lijst);
        far_neighbour=[];
        k=0;
        for delta_i=1:initval.ra                   
            neighbour_dif=[];
            n=0;
            for j=1:lel-1                   %for each segment separate
                for i=lijst(j)+1:lijst(j+1)-delta_i  %Pick_a_Segment; 
                    n=n+1 ;   
                    neighbour_dif(:,n)=data(i+delta_i,2)-data(i,2);
                end    
            end
            properties.noise=(std(neighbour_dif))/2^0.5;
            k=k+1;
            far_neighbour(:,k)=[k,properties.noise,0]';
        end    
        far_neighbour=far_neighbour';
        % Sivaraman commented
        %subplot(1,2,1);
        %hold;
        %title('Noise vs. Neighbour_distance','FontSize',12)
        %plot(far_neighbour(:,1),far_neighbour(:,2));
        %axis([0 initval.ra  0 ceil(max(far_neighbour(:,2)))]);
        %disp('Left mouse button picks points.')
        %[rolloff,n,but] = ginput(1);
        %properties.noise=n;
        %far_neighbour(:,3)=properties.noise;
    

function sorteer=sort_on_key(rij,key)
	%this function reads a (index,props) array ands sorts along the index with one of the
	%props as sort key
	size=length(rij(:,1));
	sorteer=0*rij;
	buf=rij;
	for i=1:size
        [g,h]=min(buf(:,key));
        sorteer(i,:)=buf(h,:);
        buf(h,key)=max(buf(:,key))+1;
	end


function spl=Split2(rij,i1,i2)
	%this function adresses a two-dim array 'rij'in a specific segment
	%and determines the best step-fit there
	window=i2-i1;
	if window>2
        Chisq=(1:window-1)*0;
        for t=2:window-2;
            left=rij(i1+1:i1+t,2);
            right=rij(i1+t+1:i2,2);
            left_t=rij(i1+1:i1+t,1);
            right_t=rij(i1+t+1:i2,1);
            dcleft=mean(left);
            dcright=mean(right);
            Chisq(t)=(sum((left-dcleft).^2)+sum((right-dcright).^2))/(window-1);
        end
        Chisq(1)=Chisq(2)+1;
        Chisq(window-1)=Chisq(window-2)+1;
        [g,h]=min(Chisq);
        r=rij(i1+h+1:i2,2);
        l=rij(i1+1:i1+h,2);   
        stp=mean(r)-mean(l);           
        rank1=abs(stp)*sqrt(window);        %expected rel.step accuracy relative to background noise
        rank2=abs(stp)*sqrt(window/g);       %measured rel. step accuracy  
        spl=[sqrt(g),h+i1,stp,rank1, rank2];      %minimum Chi, index, stepsize,  step*srqt(N), step*sqrt(N/Var)
	else
        spl=[0,0,0,0,0];
	end
    
    
    
