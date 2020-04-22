
%This program splits a time+data set in succesive stepfits
%Jacob Kerssemakers, 21-8-2005

%(1)Run first from command line: [data, indexes,lijst,properties,initval]=Steps_Find;
%User actions: 
	%1) Click a most likely noise level in a presented graph of
	%noise level as function of pair distance. (The choice is not affecting the
	%stepfit itself)
    %2) Graphs of a steppedness-function are presented. Choose from this
    %a reasonable number of steps 'Nst' (for a pure step train, this is where the
    %function peaks).

%(2) After that, run from the command line:
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

% Sivaraman - added extra output variables 'Steppedness', 'selectie', and
% 'Fit'
% 'Steppedness' and 'selectie' are used to determine the expected number of 
% steps in the original signal by maximizing the 'Steppedness' variable. 
% The correct number of steps is then 'selectie(argmax(Steppedness))'. 
% These two variables only have meaningful values 
% when doitforthisstepnumber == 0.
% 'Fit' contains the estimated signal.
% Verbosity = 0 avoid all the plots and saving data to the 'Out' folder


function [dummy, Steppedness, selectie, Fit]=Steps_Evaluate(data, indexes,lijst, properties,initval,doitforthisstepnumber, verbosity)
    Steppedness = -1;
    selectie = -1;
    %EVALUATION
	initval.ranking=1;              % 1: order steps by expected accuracy; 2: by measured accuracy 
    initval.showfits=0;
    initval.lo_hist=-120;            %nm
	initval.hi_hist=120;             %nm
	initval.bin_hist=6;             %nm
	initval.fast_steptest=4;        %number of points left and right of a step used for evaluating how instant a step is.
	initval.fast_steptest=4
    initval.Nst=doitforthisstepnumber;
    %-----------------------------------------------------------------------
	%end of splitting part; Up to now, all possible steps has been found.
	figcnt=0; %for counting frames for movie
    c=0;        %for collecting fits
    shft=0;
    CollectFits=[];
	lei=length(indexes);
	lel=length(lijst);
	%Now, rank them by their expected error of their measured one.
	disp('sorting...');
	if initval.ranking == 1, indexes(1:lei-lel,:)=sort_on_key(indexes(1:lei-lel,:),4); , end  %sort on expected error (proportional with [noise*stepsize*sqrt(windowsize)]^-1 )
	if initval.ranking == 2, indexes(1:lei-lel,:)=sort_on_key(indexes(1:lei-lel,:),5); , end  %sort on measured error (proportional with [Chisq^0.5*stepsize*sqrt(windowsize)]^-1 )

	
    %Initialize; %this section also uses the borderpoints to make a first fit
			%-----------------------------------------------------------------------
			Steppedness=zeros(properties.aantal,4);
			Fit=ones(properties.N0,1);
			le=length(lijst);
			for tel=1:le-1 
                Fit(lijst(tel)+1:lijst(tel+1))=mean(data(lijst(tel)+1:lijst(tel+1),2)); % Mean
			end
			Anti_Fit=Fit;
			Var=sum((data(:,2)-Fit).^2);
			Anti_Var=sum((data(:,2)-Anti_Fit).^2);
			Chi=Var/properties.N0;
			Anti_Chi=Anti_Var/properties.N0;
			teller=0;
			selected_aantal=1;
			
			if initval.Nst>0  
                stophere=initval.Nst+le;        %add number of segments
			else
                stophere=initval.cutoff;
			end
			if stophere > properties.aantal, stophere = properties.aantal; , end;
			%------------------------------------------------------------------------
	
    
    
    %Analyze
    while selected_aantal < stophere+1   %aantal is the total number after  full splitting
       teller=teller+1;
        if indexes(properties.aantal-selected_aantal+1,1)~=0  %skip 'border' points of segments
            new_index=indexes(properties.aantal+1-selected_aantal,2); 
            only_indexes=sort(indexes(properties.aantal+1-selected_aantal:properties.aantal,2));  %with new index(es) added from the main list
            only_indexes=[0,only_indexes']';
            
            %Pick_a_Segment(lijst, new_index);
            %--------------------------------
            %Find the segment where the new step-location is residing
            lel=length(lijst);
            pick=[0,0];
            for i=1:lel-1 
                if new_index > lijst(i) & new_index <lijst(i+1) 
                    pick=[lijst(i),lijst(i+1)];
                end 
            end
            
            %Change_Fits_Where_Needed;
            %the largest stretch affected by the new step is not that of the 'good' indexes on
            %the left and the right, but on the 'bad' indexes even more to the
            %left and the right of those. Meaning, three indexes determine
            %where the best fit should be changed, four indexes determine where
            %the 'worst' fit should be changed.
            %This speeds up the whole analysis especially in the
            %over-fit-range
            %------------------------------------------                      
            zoek=find(only_indexes==new_index);  %this is where a new step is in the indexrow
            gim=new_index;
            gil=only_indexes(zoek-1); % good index on left 
            gir=only_indexes(zoek+1); % ...... right.......
            %check for borders
            border=1;
            if gil== pick(1), border =2; , end %border on left side
            if gir == pick(2), border =3; , end %border on right side
            if gir == pick(2) & gil==pick(1), border =4; , end
            %determine left, right and middle bad indexes
            a=ceil((gil+gim)/2); % badindexmiddle left
            dummy=Split2(data,gil+1,gim);
            b=dummy(2);
            %bilm=max(a,b);
            if b==0, bilm=a;, else bilm=b;,end
            a=ceil((gim+gir)/2); % badindexmiddle right
            dummy=Split2(data,gim+1,gir);
            b=dummy(2);
            %birm=max(a,b);
            if b==0, birm=a;, else birm=b;,end
            %determine outer limits 
            switch border 
                case 1  %no borders
                    a=ceil((gil+only_indexes(zoek-2))/2);
                    dummy=Split2(data,only_indexes(zoek-2)+1,gil);
                    b=dummy(2);
                    %bil=max(a,b);
                    if b==0, bil=a;, else bil=b;,end
                    a=ceil((gir+only_indexes(zoek+2))/2);
                    dummy=Split2(data,gir+1,only_indexes(zoek+2));
                    b=dummy(2);
                    %bir=max(a,b);
                    if b==0, bir=a;, else bir=b;,end 
                case 2 %border on left
                    bil=pick(1);
                    a=ceil((gir+only_indexes(zoek+2))/2);
                    dummy=Split2(data,gir+1,only_indexes(zoek+2));
                    b=dummy(2);
                    %bir= max(a,b);
                    if b==0, bir=a;, else bir=b;,end
                case 3 %border on right
                    a=ceil((gil+only_indexes(zoek-2))/2);
                    dummy=Split2(data,only_indexes(zoek-2)+1,gil);
                    b=dummy(2);
                    %bil=max(a,b);
                    if b==0, bil=a;, else bil=b;,end
                    bir=pick(2);
                case 4 % both borders
                    bil=pick(1);
                    bir=pick(2);
                end
	
                %store values from parts that are going to be changed
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
	
                Fit(gil+1:gim)=mean(data(gil+1:gim,2)); 
                Fit(gim+1:gir)=mean(data(gim+1:gir,2)); 
	
                if gil ~= pick(1) %border not on left side
                    Anti_Fit(bil+1:bilm)=mean(data(bil+1:bilm,2));
                else
                    Anti_Fit(bil+1:bilm)=mean(data(bil+1:bil+3,2));
                end
                Anti_Fit(bilm+1:birm)=mean(data(bilm+1:birm,2));
                if gir ~= pick(2)%border not on right side
                    Anti_Fit(birm+1:bir)=mean(data(birm+1:bir,2));
                else
                    Anti_Fit(birm+1:bir)=mean(data(bir-2:bir,2));
                end
	            
                Chi=mean((data(:,2)-Fit).^2);
                Anti_Chi=mean((data(:,2)-Anti_Fit).^2);
                S=Anti_Chi/Chi;           %
                

                %Collect_Fit_Properties;
                %-----------------------------------------------
                %1) Determine stepsizes and distinguish between fast and
                %slow steps
                leo=length(only_indexes);
                lel=length(lijst);
                segfinder=[lijst(1:lel-1) lijst(2:lel)]';
                Steps_all=[];
                Steps_fast=[];
                FastFit=0*Fit;
                SlowFit=0*Fit;
                tel1=0;
                tel2=0;
                for j=2:leo-1
                    if any(lijst-only_indexes(j)==0)                  %skip border points
                    else
                        tel1=tel1+1; 
                        nw=only_indexes(j+1)-only_indexes(j-1);
                        step=Fit(only_indexes(j)+1)-Fit(only_indexes(j));
                        %keep track of segment
                        which_seg=find(sum(sign(segfinder-only_indexes(j)))==0);
                        %keep track of value (prior to step)
                        Steps_all(:,tel1)=[which_seg; only_indexes(j) ; only_indexes(j)-lijst(which_seg); Fit(only_indexes(j)); step; nw; 0];
                        %contains: segment totalindex segmentindex val step nw fast
                        %
                        if initval.Nst ~= 0
                            %Test how quick a step is
                            %-------------------------------

                            if      only_indexes(j)-only_indexes(j-1)>2*initval.fast_steptest -1& ...
                                    only_indexes(j+1)-only_indexes(j)>2*initval.fast_steptest-1  %left and right large enough
                                Faststep=mean(data(only_indexes(j)+1:only_indexes(j)+3,2))-mean(data(only_indexes(j)-2:only_indexes(j),2)); 
                                if step-Faststep < properties.noise*(2/initval.fast_steptest)^0.5  
                                    %the chance that an instant step yields
                                    %this much lower closely around the step location is only 15% (half of outside one-sigma confidence interval
                                    tel2=tel2+1;
                                    Steps_fast(:,tel2)=[step;nw];
                                    Steps_all(7,tel1)=1;
                                    FastFit(only_indexes(j-1)+1:only_indexes(j+1))=Fit(only_indexes(j-1)+1:only_indexes(j+1));
                                else
                                    SlowFit(only_indexes(j-1)+1:only_indexes(j+1))=Fit(only_indexes(j-1)+1:only_indexes(j+1));;
                                end
                            end
                            %_-----------------------
                        end  
                    end
                end
                                           
               
                average_step=properties.growth_range/selected_aantal;
                 
                Steppedness(teller,1)=properties.N0/selected_aantal;  %windowsize
                Steppedness(teller,2)=S;        %steppedness chibad/chigood    
                Steppedness(teller,3)=selected_aantal-length(lijst);      %stepnumber exclusief border points
                Steppedness(teller,4)=average_step;
    
                if initval.showfits==1
               %plot fits
                    if initval.Nst~=0 &selected_aantal<7
                        c=c+1;
                        shft=shft+1.5*max(Steps_all(5,:));
                        subplot(1,1,1)
                        plot(data(:,2),'LineWidth',1);
                        hold on;
                        plot( Fit-shft,'r','LineWidth',2);
                        figcnt=figcnt+1;
                        title(num2str(selected_aantal));
                        F(figcnt) = getframe;  
                       hold on;
                       CollectFits(:,c)=Fit-shft;
                   end
                   if initval.Nst~=0 &selected_aantal>5 
                       if mod(selected_aantal,5)==0
                        c=c+1;
                        shft=shft+1.5*max(Steps_all(5,:));
                        subplot(1,1,1)
                        plot(data(:,2),'LineWidth',1);
                        hold on;
                         title(num2str(selected_aantal));
                        plot( Fit-shft,'r','LineWidth',2);
                        figcnt=figcnt+1;
                        F(figcnt) = getframe;  
                       hold on;
                       CollectFits(:,c)=Fit-shft;
                        end
                   end
               end
	
            end
            selected_aantal=selected_aantal+1;
        end
        %----------------------------------------------------------------------
        if initval.showfits==1&initval.makemovie==1&initval.Nst~=0
            movie(F,1,3)
            m=strcat(initval.outpath, initval.filename,'_mov',int2str(initval.Nst), '.avi');
            movie2avi(F,m);
        end

        
        if initval.Nst ~= 0
            Bins=(initval.lo_hist:initval.bin_hist:initval.hi_hist);
            Step_hist_all=hist(Steps_all(5,:),Bins);
            Step_hist_fast=[];
            if sum(Steps_fast)>0, Step_hist_fast=hist(Steps_fast(1,:),Bins);, end;
        end
 
        Steppedness=Steppedness(1:teller,:);
        le=length(lijst);   %skip border points: 
        
        disp('Nw   S 0 Nst 0 0 AvStep 0 ')
        
        selectie=find(Steppedness(:,1)~=0);
        ls=length(selectie);
        
        subplot(1,1,1)
        hold ;
        if initval.Nst == 0
            selectie=find(Steppedness(:,1)~=0);
            subplot(3,1,1);
            hold on;
            title('Chibad/Chigood  vs. nst','FontSize',8)
            plot(Steppedness(selectie,3),Steppedness(selectie,2));  %
        
            subplot(3,1,2);
            hold on;
            title('(Chibad/Chigood  vs. nw ','FontSize',8);
            plot(Steppedness(selectie,1),Steppedness(selectie,2));  %4
       
            subplot(3,1,3);
            hold on;
            title('Chibad/Chigood-1 vs. st )','FontSize',8);
            plot(abs(Steppedness(selectie,4)),Steppedness(selectie,2)-1);  %.
            
            %figname=strcat(initval.filename,'_Stepp' , '.fig');
            %saveas(gcf,figname);
        end
        
        if initval.Nst == 0
            All_Steppedness=Steppedness; %keep for later 
            All_selectie=selectie;
            histossteps=[];
        else
            All_Steppedness=[];
            subplot(1,2,1);
            hold;
            title('Reconstructed steps','FontSize',8)
            %plot(data(:,1),data(:,2),data(:,1), FastFit,'r',data(:,1), SlowFit,'g');  
            plot(data(:,1),data(:,2),data(:,1), Fit,'r');  
            subplot(1,2,2);
            hold;
            title('Histogram of all steps: green; fast steps: red','FontSize',8)
            bar(Bins, Step_hist_all, 'g');
            if length(Step_hist_fast) > 0, bar(Bins, Step_hist_fast, 'r');,end;
            histossteps=[Bins', Step_hist_all'];  %,Step_hist_fast'
            %figname=strcat(initval.filename, '_Fit' , int2str(initval.Nst), '.fig');
            %saveas(gcf,figname);
            
        end
  
      if(verbosity == 1)
            dummy2=Write_files(initval,data,Fit,Steps_all,histossteps,All_Steppedness)  
      end
      
      disp('ready');
        

        dummy=1;


        
function dummy2=Write_files(initval,data,Fit,Steps_all,histossteps,All_Steppedness)
	if initval.Nst ~= 0
             fitname=strcat(initval.outpath,initval.filename, '_Hist' , int2str(initval.Nst), '.csv');
             csvwrite(fitname, histossteps);
             uit=[data(:,1),data(:,2), Fit];  
             fitname=strcat(initval.outpath, initval.filename,'_Fit' , int2str(initval.Nst), '.csv');
             csvwrite(fitname, uit);                
             fitname=strcat(initval.outpath, initval.filename , '_',int2str(initval.Nst), 'stps','.csv');
             csvwrite(fitname, Steps_all');  
	else
        fitname=strcat(initval.outpath, initval.filename,'_Stepp' , '.csv');
        csvwrite(fitname, All_Steppedness);
	end
	dummy2=1;

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
        rank1=abs(stp)*sqrt(window); %expected rel.step accuracy relative to background noise
        rank2=abs(stp)*sqrt(window/g);       %measured rel. step accuracy  
        spl=[g,h+i1,stp,rank1, rank2];   %minimum Chi, index, stepsize,  step*srqt(N), step*sqrt(N/Var)
	else
        spl=[0,0,0,0,0];
	end
    
        