function InverterpillarProcessing

%clear all; 
%uiopen
cd 'C:\Users\79829\Documents\Data';
cluster_class = cluster_class .';           %must have already run WaveClus spike sorting.

%%Raw data filtering 
%simple average
simpleavg = movmean(data(1,:),100);   %Smoothing filter design, normally 100

%High pass filter 
fc = 50; % Cut off frequency
fs = 10000; % Sampling rate


[b,a] = butter(6,fc/(fs/2),'high'); % Butterworth filter of order 6
HPB = filter(b,a,simpleavg); % Will be the filtered signal, high pass butter
%fvtool(HPB);      %Visualizes the high pass filter%

%% Trial segmenting
%
thresh = 1.508;                         %establish threshold
positrig = find(simpleavg>thresh);      %Find elements with value more than threshold
myDiscontinuityIndex = [zeros(1,length(positrig))];
myDiscontinuityValues = [zeros(1,length(positrig))];
myDiscontinuityValues = [positrig(1)];  %adding first point
for i=2:length(positrig)
    if positrig(i) - positrig(i-1) > 7000
        myDiscontinuityIndex = [myDiscontinuityIndex i];
        myDiscontinuityValues = [myDiscontinuityValues positrig(i-1) positrig(i)];
    end
end
myDiscontinuityValues(end+1) = [positrig(end)];    %adding last point

ii = 1;
for ii = 1:length(cluster_class(2,:))
    cluster_class(2,ii) = cluster_class(2,ii).*(10);
end    

currcyc = 1;
cycles = (numel(myDiscontinuityValues)/2);
%Trials = [];
Trials = cell(cycles,15);           %preallocate trials matrix in cell array
for i = 1:cycles
    tempstart = myDiscontinuityValues(currcyc)-50000;
    tempend = myDiscontinuityValues(currcyc+1)+50000;
    Trials{i,15} = i;                                        %inserting trial number so I can reference original data and not smoothed trace.  
    Trials{i,1} = simpleavg((tempstart):(tempend));         %work with the smoothed trace
    Trials{i,6} = HPB((tempstart):(tempend));               %add in High pass trace
    Trials{i,10} = diff(Trials{i,1});               %derivative of positon to get AC signal, run RMS on ac to get power
    clear movRMS;
    movRMS = dsp.MovingRMS(100);                 %so it is working now but i dont know what is correct. scalling shift with change of variable
    Trials{i,11} = movRMS(Trials{i,10});                         %RMS, to run on position i,1, run on HPB trace i,6, run on AC signal i,10
    Trials{i,12} = movmean(Trials{i,11},500);                    %smooth the RMS, normally set to around 500
    Trials{i,2} = data(2,(tempstart):(tempend));
    Trialcluster = find(cluster_class(2,:)>=tempstart & cluster_class(2,:)<=tempend);
    TCfirst(i) =  Trialcluster(1);
    TClast(i) = Trialcluster(end);
    Trials{i,3} = cluster_class(2,TCfirst(i):TClast(i));
    ii = 1;
    for ii=1:length(Trials{i,3})
        clustoffset(ii) = Trials{i,3}(ii)-tempstart;
    end
    Trials{i,3} = clustoffset(1,:);
    currcyc = currcyc+2;
    beep
end
clear thresh myDiscontinuityIndex myDiscontinuityValues;
clear cluster_class positrig currcyc TCfirst TClast;
clear Trialcluster tempstart tempend clustoffset;
clear simpleavg HPB a b


%% Binning
i = 1;
for i = 1:cycles
    bin1 = 1;
    bin2 = 1000;            %100 ms bin windows
    for ii = 1:((length(Trials{i,2}))./1000)
        Trials{i,4}(1,ii) = length(find(Trials{i,3}>=bin1 & Trials{i,3}<=bin2));
        bin1 = bin1+1000;
        bin2 = bin2+1000;
    end
end
clear bin1 bin2;

%% Fast Fourier Transform
i=1;
for i = 1:cycles               %FFT on position is i,1, If done on HPB then i,6
    Trialfft = fft(Trials{i,1});
    lengthdata=length(Trials{i,1});
    P2 = abs(Trialfft/lengthdata);
    P1 = P2(1:round(lengthdata/2)+1);
    P1 = P1(1:40000);
    P2 = P2(1:40000);
    P1(2:end-1)=2*P1(2:end-1);
    Trials{i,5} = P1;       %saves fft for each trial in 5th cell.
    %clear fft lengthdata P1 P2;
end
clear P1 P2 lengthdata fft;

%% Amplitude sorting
%first get slopes, then get heights
i=1;
for i = 1:cycles
    temptrial = Trials{i,1};
    amplitude{i} = movmean(temptrial,4000);
    ampmax{i} = max(amplitude{i});
    ampmax{i} = round(ampmax{i},2);
    ampmin{i} = min(amplitude{i});
    ampmin{i} = round(ampmin{i},2);
end

i = 1;
for i = 1:cycles
    if ampmax{i} == 3.39
        ampgroup1{i} = Trials(i,:);
    elseif ampmax{i} == 2.14
        ampgroup2{i} = Trials(i,:);
    end
end

ampgroup1 = ampgroup1(~cellfun('isempty',ampgroup1));
ampgroup2 = ampgroup2(~cellfun('isempty',ampgroup2));

clear temptrial amplitude ampmax ampmin;
    
%% Velocity grouping
%first get slopes, then get heights
i=1;
for i = 1:(length(ampgroup1))
    slopes{i} = ampgroup1{1,i};
    slopes{i} = movmean(slopes{i}{1,1},4500);
    slopesdiff{i} = diff(slopes{i});
    slopemax{i} = max(slopesdiff{i});
    slopemax{i} = round(slopemax{i},6);
    slopemin{i} = min(slopesdiff{i});
    slopemin{i} = round(slopemin{i},6);
end

i = 1;
for i = 1:(length(ampgroup1))
    if slopemax{i} == 3.4e-05
        group1(i) = ampgroup1(1,i);
    elseif slopemax{i} == 1.3e-04
        group2(i) = ampgroup1(1,i);
    elseif slopemax{i} == 6.6e-05
        group3(i) = ampgroup1(1,i);
    end
end

i=1;
for i = 1:(length(ampgroup2))
    slopes{i} = ampgroup2{1,i};
    slopes{i} = movmean(slopes{i}{1,1},4500);
    slopesdiff{i} = diff(slopes{i});
    slopemax{i} = max(slopesdiff{i});
    slopemax{i} = round(slopemax{i},6);
    slopemin{i} = min(slopesdiff{i});
    slopemin{i} = round(slopemin{i},6);
end

i = 1;
for i = 1:(length(ampgroup2))
    if slopemax{i} == 3.4e-05
        group4(i) = ampgroup2(1,i);
    elseif slopemax{i} == 7.1e-05
        group5(i) = ampgroup2(1,i);
    elseif slopemax{i} == 6.5e-05
        group6(i) = ampgroup2(1,i);
    end
end

group1 = group1(~cellfun('isempty',group1));
group2 = group2(~cellfun('isempty',group2));
group3 = group3(~cellfun('isempty',group3));
group4 = group4(~cellfun('isempty',group4));
group5 = group5(~cellfun('isempty',group5));
group6 = group6(~cellfun('isempty',group6));

clear slopes slopesdiff slopemax slopemin;


%crop data to ten trials for each group
i=1;
ii=11;          %cropping start trial
for i=ii:(length(group1))
    group1{1,i} = [];
end

i=1;
for i=ii:(length(group2))
    group2{1,i} = [];
end

i=1;
for i=ii:(length(group3))
    group3{1,i} = [];
end

i=1;
for i=ii:(length(group4))
    group4{1,i} = [];
end

i=1;
for i=ii:(length(group5))
    group5{1,i} = [];
end

i=1;
for i=ii:(length(group6))
    group6{1,i} = [];
end

group1 = group1(~cellfun('isempty',group1));
group2 = group2(~cellfun('isempty',group2));
group3 = group3(~cellfun('isempty',group3));
group4 = group4(~cellfun('isempty',group4));
group5 = group5(~cellfun('isempty',group5));
group6 = group6(~cellfun('isempty',group6));


%% Counting
%group1
i = 1;                                  %counting spike numbers during each phase of the trial
for i = 1:length(group1)
    group1{1,i}{1,9}(1,1) = (sum(group1{1,i}{1,4}(1,30:50))./2);      %pre
    group1{1,i}{1,9}(1,2) = (sum(group1{1,i}{1,4}(1,51:70))./1.9);     %rise1
    group1{1,i}{1,9}(1,3) = (sum(group1{1,i}{1,4}(1,71:117))./4.6);        %hold1
    group1{1,i}{1,9}(1,4) = (sum(group1{1,i}{1,4}(1,118:157))./3.9);       %rise2
    group1{1,i}{1,9}(1,5) = (sum(group1{1,i}{1,4}(1,158:205))./4.7);       %hold2
    group1{1,i}{1,9}(1,6) = (sum(group1{1,i}{1,4}(1,206:264))./5.8);       %fall
    group1{1,i}{1,9}(1,7) = (sum(group1{1,i}{1,4}(1,265:285))./2);       %post
end

%group2
i = 1;
for i = 1:length(group2)       
    group2{1,i}{1,9}(1,1) = (sum(group2{1,i}{1,4}(1,30:50))./(2));
    group2{1,i}{1,9}(1,2) = (sum(group2{1,i}{1,4}(1,51:55))./(0.4));
    group2{1,i}{1,9}(1,3) = (sum(group2{1,i}{1,4}(1,56:105))./(4.9));
    group2{1,i}{1,9}(1,4) = (sum(group2{1,i}{1,4}(1,106:115))./(0.9));
    group2{1,i}{1,9}(1,5) = (sum(group2{1,i}{1,4}(1,116:164))./(4.8));
    group2{1,i}{1,9}(1,6) = (sum(group2{1,i}{1,4}(1,165:180))./(1.5));
    group2{1,i}{1,9}(1,7) = (sum(group2{1,i}{1,4}(1,181:201))./(2));
end

%group3
i = 1;
for i = 1:length(group3)               
    group3{1,i}{1,9}(1,1) = (sum(group3{1,i}{1,4}(1,30:50))./(2));
    group3{1,i}{1,9}(1,2) = (sum(group3{1,i}{1,4}(1,51:60))./(0.9));
    group3{1,i}{1,9}(1,3) = (sum(group3{1,i}{1,4}(1,61:109))./(4.8));
    group3{1,i}{1,9}(1,4) = (sum(group3{1,i}{1,4}(1,110:129))./(1.9));
    group3{1,i}{1,9}(1,5) = (sum(group3{1,i}{1,4}(1,130:178))./(4.8));
    group3{1,i}{1,9}(1,6) = (sum(group3{1,i}{1,4}(1,179:208))./(2.9));
    group3{1,i}{1,9}(1,7) = (sum(group3{1,i}{1,4}(1,208:228))./(2));
end

%group4
 i = 1;
for i = 1:length(group4)                 
    group4{1,i}{1,9}(1,1) = (sum(group4{1,i}{1,4}(1,30:50))./(2));
    group4{1,i}{1,9}(1,2) = (sum(group4{1,i}{1,4}(1,51:60))./(0.9));
    group4{1,i}{1,9}(1,3) = (sum(group4{1,i}{1,4}(1,61:107))./(4.6));
    group4{1,i}{1,9}(1,4) = (sum(group4{1,i}{1,4}(1,108:118))./(1));
    group4{1,i}{1,9}(1,5) = (sum(group4{1,i}{1,4}(1,119:165))./(4.6));
    group4{1,i}{1,9}(1,6) = (sum(group4{1,i}{1,4}(1,166:185))./(1.9));
    group4{1,i}{1,9}(1,7) = (sum(group4{1,i}{1,4}(1,186:206))./(2));
end

%group5
i = 1;
for i = 1:length(group5)              
    group5{1,i}{1,9}(1,1) = (sum(group5{1,i}{1,4}(1,30:50))./(2));
    group5{1,i}{1,9}(1,2) = (sum(group5{1,i}{1,4}(1,51:53))./(0.2));
    group5{1,i}{1,9}(1,3) = (sum(group5{1,i}{1,4}(1,54:102))./(4.8));
    group5{1,i}{1,9}(1,4) = (sum(group5{1,i}{1,4}(1,103:105))./(0.2));
    group5{1,i}{1,9}(1,5) = (sum(group5{1,i}{1,4}(1,106:154))./(4.8));
    group5{1,i}{1,9}(1,6) = (sum(group5{1,i}{1,4}(1,155:159))./(0.4));
    group5{1,i}{1,9}(1,7) = (sum(group5{1,i}{1,4}(1,160:180))./(2));
end

%group6
i = 1;
for i = 1:length(group6)                   %change time points
    group6{1,i}{1,9}(1,1) = (sum(group6{1,i}{1,4}(1,30:50))./(2));
    group6{1,i}{1,9}(1,2) = (sum(group6{1,i}{1,4}(1,51:55))./(0.4));
    group6{1,i}{1,9}(1,3) = (sum(group6{1,i}{1,4}(1,56:103))./(4.7));
    group6{1,i}{1,9}(1,4) = (sum(group6{1,i}{1,4}(1,104:109))./(0.5));
    group6{1,i}{1,9}(1,5) = (sum(group6{1,i}{1,4}(1,110:157))./(4.7));
    group6{1,i}{1,9}(1,6) = (sum(group6{1,i}{1,4}(1,158:168))./(1));
    group6{1,i}{1,9}(1,7) = (sum(group6{1,i}{1,4}(1,169:189))./(2));
end

%% RMS averaging
%group1
i = 1;                                  %counting spike numbers during each phase of the trial
for i = 1:length(group1)
    group1{1,i}{1,7} = ((group1{1,i}{1,10}).^(2));              %run square of all data
    group1{1,i}{1,8}(1,1) = (sum(group1{1,i}{1,7}(1,30000:50000))./(20000));        %average data
    group1{1,i}{1,8}(1,1) = (sqrt(group1{1,i}{1,8}(1,1)));                 %RMS value for time sequence
    group1{1,i}{1,8}(1,2) = (sum(group1{1,i}{1,7}(1,51000:70000))./(19000));     %rise1
    group1{1,i}{1,8}(1,2) = (sqrt(group1{1,i}{1,8}(1,2)));
    group1{1,i}{1,8}(1,3) = (sum(group1{1,i}{1,7}(1,71000:117000))./(46000));        %hold1
    group1{1,i}{1,8}(1,3) = (sqrt(group1{1,i}{1,8}(1,3)));
    group1{1,i}{1,8}(1,4) = (sum(group1{1,i}{1,7}(1,118000:157000))./(39000));       %rise2
    group1{1,i}{1,8}(1,4) = (sqrt(group1{1,i}{1,8}(1,4)));
    group1{1,i}{1,8}(1,5) = (sum(group1{1,i}{1,7}(1,158000:205000))./(47000));       %hold2
    group1{1,i}{1,8}(1,5) = (sqrt(group1{1,i}{1,8}(1,5)));
    group1{1,i}{1,8}(1,6) = (sum(group1{1,i}{1,7}(1,206000:264000))./(58000));       %fall
    group1{1,i}{1,8}(1,6) = (sqrt(group1{1,i}{1,8}(1,6)));
    group1{1,i}{1,8}(1,7) = (sum(group1{1,i}{1,7}(1,265000:285000))./(20000));       %post
    group1{1,i}{1,8}(1,7) = (sqrt(group1{1,i}{1,8}(1,7)));
end

%group2
 i = 1;
for i = 1:length(group2)       
    group2{1,i}{1,7} = ((group2{1,i}{1,10}).^(2));         
    group2{1,i}{1,8}(1,1) = (sum(group2{1,i}{1,7}(1,30000:50000))./(20000));     
    group2{1,i}{1,8}(1,1) = (sqrt(group2{1,i}{1,8}(1,1)));                
    group2{1,i}{1,8}(1,2) = (sum(group2{1,i}{1,7}(1,51000:55000))./(4000));
    group2{1,i}{1,8}(1,2) = (sqrt(group2{1,i}{1,8}(1,2)));
    group2{1,i}{1,8}(1,3) = (sum(group2{1,i}{1,7}(1,56000:105000))./(49000));        
    group2{1,i}{1,8}(1,3) = (sqrt(group2{1,i}{1,8}(1,3)));
    group2{1,i}{1,8}(1,4) = (sum(group2{1,i}{1,7}(1,106000:115000))./(9000));      
    group2{1,i}{1,8}(1,4) = (sqrt(group2{1,i}{1,8}(1,4)));
    group2{1,i}{1,8}(1,5) = (sum(group2{1,i}{1,7}(1,116000:164000))./(48000));      
    group2{1,i}{1,8}(1,5) = (sqrt(group2{1,i}{1,8}(1,5)));
    group2{1,i}{1,8}(1,6) = (sum(group2{1,i}{1,7}(1,165000:180000))./(15000));    
    group2{1,i}{1,8}(1,6) = (sqrt(group2{1,i}{1,8}(1,6)));
    group2{1,i}{1,8}(1,7) = (sum(group2{1,i}{1,7}(1,181000:201000))./(20000));   
    group2{1,i}{1,8}(1,7) = (sqrt(group2{1,i}{1,8}(1,7)));
end


% %group3
i = 1;
for i = 1:length(group3)               
    group3{1,i}{1,7} = ((group3{1,i}{1,10}).^(2));        
    group3{1,i}{1,8}(1,1) = (sum(group3{1,i}{1,7}(1,30000:50000))./(20000));     
    group3{1,i}{1,8}(1,1) = (sqrt(group3{1,i}{1,8}(1,1)));                
    group3{1,i}{1,8}(1,2) = (sum(group3{1,i}{1,7}(1,51000:60000))./(9000));
    group3{1,i}{1,8}(1,2) = (sqrt(group3{1,i}{1,8}(1,2)));
    group3{1,i}{1,8}(1,3) = (sum(group3{1,i}{1,7}(1,61000:109000))./(48000));        
    group3{1,i}{1,8}(1,3) = (sqrt(group3{1,i}{1,8}(1,3)));
    group3{1,i}{1,8}(1,4) = (sum(group3{1,i}{1,7}(1,110000:129000))./(19000));      
    group3{1,i}{1,8}(1,4) = (sqrt(group3{1,i}{1,8}(1,4)));
    group3{1,i}{1,8}(1,5) = (sum(group3{1,i}{1,7}(1,130000:178000))./(48000));      
    group3{1,i}{1,8}(1,5) = (sqrt(group3{1,i}{1,8}(1,5)));
    group3{1,i}{1,8}(1,6) = (sum(group3{1,i}{1,7}(1,179000:208000))./(29000));    
    group3{1,i}{1,8}(1,6) = (sqrt(group3{1,i}{1,8}(1,6)));
    group3{1,i}{1,8}(1,7) = (sum(group3{1,i}{1,7}(1,208000:228000))./(20000));   
    group3{1,i}{1,8}(1,7) = (sqrt(group3{1,i}{1,8}(1,7)));
end

%group4
i = 1;
for i = 1:length(group4)                 
    group4{1,i}{1,7} = ((group4{1,i}{1,10}).^(2));         
    group4{1,i}{1,8}(1,1) = (sum(group4{1,i}{1,7}(1,30000:50000))./(20000));     
    group4{1,i}{1,8}(1,1) = (sqrt(group4{1,i}{1,8}(1,1)));                
    group4{1,i}{1,8}(1,2) = (sum(group4{1,i}{1,7}(1,51000:60000))./(9000));
    group4{1,i}{1,8}(1,2) = (sqrt(group4{1,i}{1,8}(1,2)));
    group4{1,i}{1,8}(1,3) = (sum(group4{1,i}{1,7}(1,61000:107000))./(46000));        
    group4{1,i}{1,8}(1,3) = (sqrt(group4{1,i}{1,8}(1,3)));
    group4{1,i}{1,8}(1,4) = (sum(group4{1,i}{1,7}(1,108000:118000))./(10000));      
    group4{1,i}{1,8}(1,4) = (sqrt(group4{1,i}{1,8}(1,4)));
    group4{1,i}{1,8}(1,5) = (sum(group4{1,i}{1,7}(1,119000:165000))./(46000));      
    group4{1,i}{1,8}(1,5) = (sqrt(group4{1,i}{1,8}(1,5)));
    group4{1,i}{1,8}(1,6) = (sum(group4{1,i}{1,7}(1,166000:185000))./(19000));    
    group4{1,i}{1,8}(1,6) = (sqrt(group4{1,i}{1,8}(1,6)));
    group4{1,i}{1,8}(1,7) = (sum(group4{1,i}{1,7}(1,181000:201000))./(20000));   
    group4{1,i}{1,8}(1,7) = (sqrt(group4{1,i}{1,8}(1,7)));
end

%group5
i = 1;
for i = 1:length(group5)              
    group5{1,i}{1,7} = ((group5{1,i}{1,10}).^(2));         
    group5{1,i}{1,8}(1,1) = (sum(group5{1,i}{1,7}(1,30000:50000))./(20000));     
    group5{1,i}{1,8}(1,1) = (sqrt(group5{1,i}{1,8}(1,1)));                
    group5{1,i}{1,8}(1,2) = (sum(group5{1,i}{1,7}(1,51000:53000))./(2000));
    group5{1,i}{1,8}(1,2) = (sqrt(group5{1,i}{1,8}(1,2)));
    group5{1,i}{1,8}(1,3) = (sum(group5{1,i}{1,7}(1,54000:102000))./(48000));        
    group5{1,i}{1,8}(1,3) = (sqrt(group5{1,i}{1,8}(1,3)));
    group5{1,i}{1,8}(1,4) = (sum(group5{1,i}{1,7}(1,103000:105000))./(2000));      
    group5{1,i}{1,8}(1,4) = (sqrt(group5{1,i}{1,8}(1,4)));
    group5{1,i}{1,8}(1,5) = (sum(group5{1,i}{1,7}(1,106000:154000))./(48000));      
    group5{1,i}{1,8}(1,5) = (sqrt(group5{1,i}{1,8}(1,5)));
    group5{1,i}{1,8}(1,6) = (sum(group5{1,i}{1,7}(1,155000:159000))./(4000));    
    group5{1,i}{1,8}(1,6) = (sqrt(group5{1,i}{1,8}(1,6)));
    group5{1,i}{1,8}(1,7) = (sum(group5{1,i}{1,7}(1,160000:180000))./(20000));   
    group5{1,i}{1,8}(1,7) = (sqrt(group5{1,i}{1,8}(1,7)));
end

%group6
i = 1;
for i = 1:length(group6)                  
    group6{1,i}{1,7} = ((group6{1,i}{1,10}).^(2));         
    group6{1,i}{1,8}(1,1) = (sum(group6{1,i}{1,7}(1,30000:50000))./(20000));     
    group6{1,i}{1,8}(1,1) = (sqrt(group6{1,i}{1,8}(1,1)));                
    group6{1,i}{1,8}(1,2) = (sum(group6{1,i}{1,7}(1,51000:55000))./(4000));
    group6{1,i}{1,8}(1,2) = (sqrt(group6{1,i}{1,8}(1,2)));
    group6{1,i}{1,8}(1,3) = (sum(group6{1,i}{1,7}(1,56000:103000))./(47000));        
    group6{1,i}{1,8}(1,3) = (sqrt(group6{1,i}{1,8}(1,3)));
    group6{1,i}{1,8}(1,4) = (sum(group6{1,i}{1,7}(1,104000:109000))./(5000));      
    group6{1,i}{1,8}(1,4) = (sqrt(group6{1,i}{1,8}(1,4)));
    group6{1,i}{1,8}(1,5) = (sum(group6{1,i}{1,7}(1,110000:157000))./(47000));      
    group6{1,i}{1,8}(1,5) = (sqrt(group6{1,i}{1,8}(1,5)));
    group6{1,i}{1,8}(1,6) = (sum(group6{1,i}{1,7}(1,158000:168000))./(10000));    
    group6{1,i}{1,8}(1,6) = (sqrt(group6{1,i}{1,8}(1,6)));
    group6{1,i}{1,8}(1,7) = (sum(group6{1,i}{1,7}(1,169000:189000))./(20000));   
    group6{1,i}{1,8}(1,7) = (sqrt(group6{1,i}{1,8}(1,7)));
end

%%Correlation
i = 1;
for i = 1:length(group1) 
    tempcorrcoef = corrcoef((group1{1,i}{1,8}),(group1{1,i}{1,9}));
    group1{1,i}{1,13} = tempcorrcoef(1,2);
end

i = 1;
for i = 1:length(group2) 
    tempcorrcoef = corrcoef((group2{1,i}{1,8}),(group2{1,i}{1,9}));
    group2{1,i}{1,13} = tempcorrcoef(1,2);
end

i = 1;
for i = 1:length(group3) 
    tempcorrcoef = corrcoef((group3{1,i}{1,8}),(group3{1,i}{1,9}));
    group3{1,i}{1,13} = tempcorrcoef(1,2);
end

i = 1;
for i = 1:length(group4) 
    tempcorrcoef = corrcoef((group4{1,i}{1,8}),(group4{1,i}{1,9}));
    group4{1,i}{1,13} = tempcorrcoef(1,2);
end

i = 1;
for i = 1:length(group5) 
    tempcorrcoef = corrcoef((group5{1,i}{1,8}),(group5{1,i}{1,9}));
    group5{1,i}{1,13} = tempcorrcoef(1,2);
end

i = 1;
for i = 1:length(group6) 
    tempcorrcoef = corrcoef((group6{1,i}{1,8}),(group6{1,i}{1,9}));
    group6{1,i}{1,13} = tempcorrcoef(1,2);
end
clear tempcorrcoef


i = 1;
figure;
hold on;
for i = 1:length(group4)
    clear f x y p;
    x = group4{1,i}{1,8};
    y = group4{1,i}{1,9};
    p = polyfit(x,y,1);
    f = polyval(p,x);
    plot (x,y,'o',x,f,'r');
    meanx(1,1) = mean(x(1,[1 3 5 7]));     
    meanx(1,2) = mean(x(1,[2 4 6]));
    g4meany(i,1) = mean(y(1,[1]));
    g4meany(i,2) = mean(y(1,[2 4  ]));
    %g4meany(i,3) = g4meany(i,2)./g4meany(i,1);
    g4meany(i,:) = (g4meany(i,:)).*(g4meany(1,1)./g4meany(i,1));
    %g4meany(i,:) = g4meany(i,1:2)-(g4meany(i,1)); 
    meanp = polyfit(meanx,g4meany(i,:),1);
    meanf = polyval(meanp,meanx);
    plot (meanx,g4meany(i,:),meanx,meanf,'b');
end

i = 1;
%figure;
hold on;
for i = 1:length(group5)
    clear f x y p;
    x = group5{1,i}{1,8};
    y = group5{1,i}{1,9};
    p = polyfit(x,y,1);
    f = polyval(p,x);
    plot (x,y,'o',x,f,'r');
    meanx(1,1) = mean(x(1,[1 3 5 7]));     
    meanx(1,2) = mean(x(1,[2 4 6]));
    g5meany(i,1) = mean(y(1,[1]));
    g5meany(i,2) = mean(y(1,[2 4  ]));
    %g5meany(i,3) = g5meany(i,2)./g5meany(i,1);
    g5meany(i,:) = (g5meany(i,:)).*(g5meany(1,1)./g5meany(i,1));
    %g5meany(i,:) = g5meany(i,1:2)-(g5meany(i,1)); 
    meanp = polyfit(meanx,g5meany(i,:),1);
    meanf = polyval(meanp,meanx);
    plot (meanx,g5meany(i,:),meanx,meanf,'r');
end

i = 1;
%figure;
hold on;
for i = 1:length(group6)
    clear f x y p;
    x = group6{1,i}{1,8};
    y = group6{1,i}{1,9};
    p = polyfit(x,y,1);
    f = polyval(p,x);
    plot (x,y,'o',x,f,'r');
    meanx(1,1) = mean(x(1,[1 3 5 7]));     
    meanx(1,2) = mean(x(1,[2 4 6]));
    g6meany(i,1) = mean(y(1,[1]));
    g6meany(i,2) = mean(y(1,[2 4  ]));
    %g6meany(i,3) = g6meany(i,2)./g6meany(i,1);
    g6meany(i,:) = (g6meany(i,:)).*(g6meany(1,1)./g6meany(i,1));
    %g6meany(i,:) = g6meany(i,1:2)-(g6meany(i,1)); 
    meanp = polyfit(meanx,g6meany(i,:),1);
    meanf = polyval(meanp,meanx);
    plot (meanx,g6meany(i,:),meanx,meanf,'g');
end

figure;
hold on;
i=1;
for i = 1:length(group1)
    %figure;
    %hold on;
    %group1
    meanx(1,1) = mean(group1{1,i}{1,8}(1,[1 3 5 7]));     
    meanx(1,2) = mean(group1{1,i}{1,8}(1,[2 4 6]));
    g1meany(i,1) = mean(group1{1,i}{1,9}(1,[1]));
    g1meany(i,2) = mean(group1{1,i}{1,9}(1,[2 4 6 ]));
    g1meany(i,:) = (g1meany(i,:)).*(g1meany(1,1)./g1meany(i,1));
    meanp = polyfit(meanx,g1meany(i,:),1);
    meanf = polyval(meanp,meanx);
    plot (meanx,g1meany(i,:),meanx,meanf);
    %group2
    meanx(1,1) = mean(group2{1,i}{1,8}(1,[1 3 5 7]));     
    meanx(1,2) = mean(group2{1,i}{1,8}(1,[2 4 6]));
    g2meany(i,1) = mean(group2{1,i}{1,9}(1,[1]));
    g2meany(i,2) = mean(group2{1,i}{1,9}(1,[2 4 6 ]));
    g2meany(i,:) = (g2meany(i,:)).*(g2meany(1,1)./g2meany(i,1));
    meanp = polyfit(meanx,g2meany(i,:),1);
    meanf = polyval(meanp,meanx);
    plot (meanx,g2meany(i,:),meanx,meanf);
    %group3
    meanx(1,1) = mean(group3{1,i}{1,8}(1,[1 3 5 7]));     
    meanx(1,2) = mean(group3{1,i}{1,8}(1,[2 4 6]));
    g3meany(i,1) = mean(group3{1,i}{1,9}(1,[1]));
    g3meany(i,2) = mean(group3{1,i}{1,9}(1,[2 4 6]));
    g3meany(i,:) = (g3meany(i,:)).*(g3meany(1,1)./g3meany(i,1));
    meanp = polyfit(meanx,g3meany(i,:),1);
    meanf = polyval(meanp,meanx);
    plot (meanx,g3meany(i,:),meanx,meanf);
    %pause(2);
    %hold off;
    %close;
end

%figure;
%hold on;
i=1;
for i = 1:length(group4)
    %figure;
    %hold on;
    %group4
    meanx(1,1) = mean(group4{1,i}{1,8}(1,[1 3 5 7]));     
    meanx(1,2) = mean(group4{1,i}{1,8}(1,[2 4 6]));
    g4meany(i,1) = mean(group4{1,i}{1,9}(1,[1]));
    g4meany(i,2) = mean(group4{1,i}{1,9}(1,[2 4 6 ]));
    g4meany(i,:) = (g4meany(i,:)).*(g4meany(1,1)./g4meany(i,1));
    meanp = polyfit(meanx,g4meany(i,:),1);
    meanf = polyval(meanp,meanx);
    plot (meanx,g4meany(i,:),'o',meanx,meanf);
    %group5
    meanx(1,1) = mean(group5{1,i}{1,8}(1,[1 3 5 7]));     
    meanx(1,2) = mean(group5{1,i}{1,8}(1,[2 4 6]));
    g5meany(i,1) = mean(group5{1,i}{1,9}(1,[1]));
    g5meany(i,2) = mean(group5{1,i}{1,9}(1,[2 4 6 ]));
    g5meany(i,:) = (g5meany(i,:)).*(g5meany(1,1)./g5meany(i,1));
    meanp = polyfit(meanx,g5meany(i,:),1);
    meanf = polyval(meanp,meanx);
    plot (meanx,g5meany(i,:),'o',meanx,meanf);
    %group6
    meanx(1,1) = mean(group6{1,i}{1,8}(1,[1 3 5 7]));     
    meanx(1,2) = mean(group6{1,i}{1,8}(1,[2 4 6]));
    g6meany(i,1) = mean(group6{1,i}{1,9}(1,[1]));
    g6meany(i,2) = mean(group6{1,i}{1,9}(1,[2 4 6]));
    g6meany(i,:) = (g6meany(i,:)).*(g6meany(1,1)./g6meany(i,1));
    meanp = polyfit(meanx,g6meany(i,:),1);
    meanf = polyval(meanp,meanx);
    plot (meanx,g6meany(i,:),'o',meanx,meanf);
    %pause(2);
    %hold off;
    %close;
end
clear f x y p meanx meanp meanf;

% %% displacement plotting
% i = 1;
% disp = [];
% figure;
% hold on;
% for i = 1:length(disp);
%     clear f x y p s;
%     x = [0 250 500];        %measured in microns, 0, 250, and 500 microns
%     y = disp(i,:);
%     s = scatter(x,y,10,'b');
%     %p = polyfit(x,y,1);
%     %f = polyval(p,x);
%     %plot (x,y,'o',x,f,'r');
% end
% 
% %%velocity plotting
% i = 1;
% figure;
% %velo = [];
% hold on;
% for i = 1:length(velo);
%     clear f x y p s;
%     x = [0.5 1 2];        %measured in mm/s, 0.5, 1, and 2 mm/s
%     y = velo(i,:);
%     s = scatter(x,y,10,'b');
%     h1 = lsline;
% end

%% PLOTTING
%
group1fft = plus(group1{1,1}{1,5},group1{1,2}{1,5});
for i = 3:length(group1)
    group1fft = plus(group1fft,group1{1,i}{1,5});
end
group1fft = tsmovavg(group1fft(1,:),'s',50);   %Smoothing filter design
%
group2fft = plus(group2{1,1}{1,5},group2{1,2}{1,5});
for i = 3:length(group2)
    group2fft = plus(group2fft,group2{1,i}{1,5});
end
group2fft = tsmovavg(group2fft(1,:),'s',50);   %Smoothing filter design
%
group3fft = plus(group3{1,1}{1,5},group3{1,2}{1,5});
for i = 3:length(group3)
    group3fft = plus(group3fft,group3{1,i}{1,5});
end
group3fft = tsmovavg(group3fft(1,:),'s',50);   %Smoothing filter design
%
group4fft = plus(group4{1,1}{1,5},group4{1,2}{1,5});
for i = 3:length(group4)
    group4fft = plus(group4fft,group4{1,i}{1,5});
end
group4fft = tsmovavg(group4fft(1,:),'s',50);   %Smoothing filter design
%
group5fft = plus(group5{1,1}{1,5},group5{1,2}{1,5});
for i = 3:length(group5)
    group5fft = plus(group5fft,group5{1,i}{1,5});
end
group5fft = tsmovavg(group5fft(1,:),'s',50);   %Smoothing filter design
%
group6fft = plus(group6{1,1}{1,5},group6{1,2}{1,5});
for i = 3:length(group6)
    group6fft = plus(group6fft,group6{1,i}{1,5});
end
group6fft = tsmovavg(group6fft(1,:),'s',50);   %Smoothing filter design

figure, 
subplot(3,1,1);
hold on;
plot(group1fft);
plot(group4fft,'r');
hold off;
subplot(3,1,2);
hold on;
plot(group2fft);
plot(group5fft,'r');
hold off;
subplot(3,1,3);
hold on;
plot(group3fft);
plot(group6fft,'r');
hold off;
close;

%%
ii = 1;         %FIGURE 1
figure;
%subplot(2,1,1);
for ii = 1:(length(group1))
    plot((group1{1,ii}{1,1}*0.4)+2.2); %group1{1,ii{1,1} for smoothed trace, group1{1,ii{1,6} for HPB trace
    hold on  
    plot((group1{1,ii}{1,6}*100)+2.3);    %plot high pass trace   
end


ii = 1;
L1 = 1.35;
L2 = 1.4;
for ii = 1:(length(group1))
    clear t1;
    hold on;
    iii = 1;
    t1 = group1{1,ii}{1,3};
     for iii = 1:length(t1)
         line([t1(iii) t1(iii)],[L1 L2]);
     end
    L1 = L1 - 0.05;
    L2 = L2 - 0.05;
    xlim([1 (length(group1{1,1}{1,2}))]);
end


group1avg = plus(group1{1,1}{1,4},group1{1,2}{1,4});        %averaged together firing freq.
for i = 3:length(group1)
    group1avg = plus(group1avg,group1{1,i}{1,4});
end
group1avg = resample(group1avg,1000,1);
plot((group1avg.*0.005)-0.4);
hline1 = refline(0,-0.4); %horizontal baseline
hline1.Color = 'k';
hline1.LineStyle = '-';

group1RMSavg = plus(group1{1,1}{1,12}(1:313000),group1{1,2}{1,12}(1:313000));
for i = 3:length(group1)       %since we already added the first two together we then add the rest starting at the 3rd trial.
    group1RMSavg = plus(group1RMSavg,group1{1,i}{1,12}(1:313000));
end
plot((group1RMSavg*200)+1.6);
hline2 = refline(0,1.6); %horizontal baseline
hline2.Color = 'k';
ylim([-1 4]);

%labeling panel
x1 = 0;      
y1 = 2.9;
txt1 = 'Position trace';
text(x1,y1,txt1)

x1 = 0;
y1 = 2.4;
txt2 = 'High pass trace (50 hz cutoff)';
text(x1,y1,txt2)

x1 = 0;
y1 = 1.7;
txt3 = 'RMS average';
text(x1,y1,txt3)

x1 = 0;
y1 = 1.5;
txt4 = 'Raster plot';
text(x1,y1,txt4)

x1 = 0;
y1 = -0.3;
txt5 = 'Averaged firing frequency (hz)';
text(x1,y1,txt5)


%%
ii = 1;
figure;
%subplot(2,1,1);
for ii = 1:(length(group2))
    plot((group2{1,ii}{1,1}*0.4)+2.2);
    hold on
    plot((group2{1,ii}{1,6})*100+2.3); %highpass trace plotting
end
ii = 1;
L1 = 1.35;
L2 = 1.4;
for ii = 1:(length(group2))
    clear t1;
    hold on;
    iii = 1;
    t1 = group2{1,ii}{1,3};
    for iii = 1:length(t1)
        line([t1(iii) t1(iii)],[L1 L2]);
    end
    L1 = L1 - 0.05;
    L2 = L2 - 0.05;
    xlim([1 (length(group2{1,1}{1,2}))]);
end

group2avg = plus(group2{1,1}{1,4},group2{1,2}{1,4});
for i = 3:length(group2)
    group2avg = plus(group2avg,group2{1,i}{1,4});
end
group2avg = resample(group2avg,1000,1);
plot((group2avg.*0.005)-0.4);
hline = refline(0,-0.4); %horizontal baseline
hline.Color = 'k';

group2RMSavg = plus(group2{1,1}{1,12}(1:228000),group2{1,2}{1,12}(1:228000));
for i = 3:length(group2)       %since we already added the first two together we then add the rest starting at the 3rd trial.
    group2RMSavg = plus(group2RMSavg,group2{1,i}{1,12}(1:228000));
end
plot((group2RMSavg*200)+1.6);
hline = refline(0,1.6); %horizontal baseline
hline.Color = 'k';
%labeling panel
x1 = 0;      
y1 = 2.9;
txt1 = 'Position trace';
text(x1,y1,txt1)

x1 = 0;
y1 = 2.4;
txt2 = 'High pass trace (50 hz cutoff)';
text(x1,y1,txt2)

x1 = 0;
y1 = 1.7;
txt3 = 'RMS average';
text(x1,y1,txt3)

x1 = 0;
y1 = 1.5;
txt4 = 'Raster plot';
text(x1,y1,txt4)

x1 = 0;
y1 = -0.3;
txt5 = 'Averaged firing frequency (hz)';
text(x1,y1,txt5)
ylim([-1 4]);

%%
figure;
ii = 1;
%subplot(2,1,2);
for ii = 1:(length(group3))
    plot((group3{1,ii}{1,1}*0.4)+2.2); %group1{1,ii{1,1} for smoothed trace, group1{1,ii{1,6} for HPB trace
    hold on
    plot((group3{1,ii}{1,6}*100)+2.3);
    %ylim([1.3 3.7]);
end
ii = 1;
L1 = 1.35;
L2 = 1.4;
for ii = 1:(length(group3))
    clear t1;
    hold on;
    iii = 1;
    t1 = group3{1,ii}{1,3};
    for iii = 1:length(t1)
        line([t1(iii) t1(iii)],[L1 L2]);
    end
    L1 = L1 - 0.05;
    L2 = L2 - 0.05;
    xlim([1 (length(group3{1,1}{1,2}))]);
end

group3avg = plus(group3{1,1}{1,4},group3{1,2}{1,4});
for i = 3:length(group3)
    group3avg = plus(group3avg,group3{1,i}{1,4});
end
group3avg = resample(group3avg,1000,1);
plot((group3avg.*0.005)-0.4);
hline = refline(0,-0.4); %horizontal baseline
hline.Color = 'k';

group3RMSavg = plus(group3{1,1}{1,12}(1:257000),group3{1,2}{1,12}(1:257000));
for i = 3:length(group3)       %since we already added the first two together we then add the rest starting at the 3rd trial.
    group3RMSavg = plus(group3RMSavg,group3{1,i}{1,12}(1:257000));
end
plot((group3RMSavg*200)+1.6);
hline = refline(0,1.6); %horizontal baseline
hline.Color = 'k';

ylim([-1 4]);

%labeling panel
x1 = 0;      
y1 = 2.9;
txt1 = 'Position trace';
text(x1,y1,txt1)

x1 = 0;
y1 = 2.4;
txt2 = 'High pass trace (50 hz cutoff)';
text(x1,y1,txt2)

x1 = 0;
y1 = 1.7;
txt3 = 'RMS average';
text(x1,y1,txt3)

x1 = 0;
y1 = 1.5;
txt4 = 'Raster plot';
text(x1,y1,txt4)

x1 = 0;
y1 = -0.3;
txt5 = 'Averaged firing frequency (hz)';
text(x1,y1,txt5)

%%
ii = 1;
figure;
%subplot(2,1,1);
for ii = 1:(length(group4))
    plot((group4{1,ii}{1,1}*0.4)+2.2);
    hold on
    plot((group4{1,ii}{1,6})*100+2.3);
end
ii = 1;
L1 = 1.35;
L2 = 1.4;
for ii = 1:(length(group4))
    clear t1;
    hold on;
    iii = 1;
    t1 = group4{1,ii}{1,3};
    for iii = 1:length(t1)
        line([t1(iii) t1(iii)],[L1 L2]);
    end
    L1 = L1 - 0.05;
    L2 = L2 - 0.05;
    xlim([1 (length(group4{1,1}{1,2}))]);
end

group4avg = plus(group4{1,1}{1,4}(1:234),group4{1,2}{1,4}(1:234));
for i = 3:length(group4)  %normally this is length but there is a problem with one not agreeing so its custom.
    group4avg = plus(group4avg,group4{1,i}{1,4}(1:234));
end
group4avg = resample(group4avg,1000,1);
plot((group4avg.*0.005)-0.4);
hline = refline(0,-0.4); %horizontal baseline
hline.Color = 'k';

group4RMSavg = plus(group4{1,1}{1,12}(1:234000),group4{1,2}{1,12}(1:234000));
for i = 3:length(group4)       %since we already added the first two together we then add the rest starting at the 3rd trial.
    group4RMSavg = plus(group4RMSavg,group4{1,i}{1,12}(1:234000));
end
plot((group4RMSavg*200)+1.6);
hline2 = refline(0,1.6); %horizontal baseline
hline2.Color = 'k';
ylim([-1 4]);

%labeling panel
x1 = 0;      
y1 = 2.9;
txt1 = 'Position trace';
text(x1,y1,txt1)

x1 = 0;
y1 = 2.4;
txt2 = 'High pass trace (50 hz cutoff)';
text(x1,y1,txt2)

x1 = 0;
y1 = 1.7;
txt3 = 'RMS average';
text(x1,y1,txt3)

x1 = 0;
y1 = 1.5;
txt4 = 'Raster plot';
text(x1,y1,txt4)

x1 = 0;
y1 = -0.3;
txt5 = 'Averaged firing frequency (hz)';
text(x1,y1,txt5)

%%
ii = 1;
%subplot(2,1,2);
figure
for ii = 1:(length(group5))
    plot((group5{1,ii}{1,1}*0.4)+2.2);
    hold on
    plot(((group5{1,ii}{1,6})*100)+2.3);
    %ylim([1.3 3.7]);
end
ii = 1;
L1 = 1.35;
L2 = 1.4;
for ii = 1:(length(group5))
    clear t1;
    hold on;
    iii = 1;
    t1 = group5{1,ii}{1,3};
    for iii = 1:length(t1)
        line([t1(iii) t1(iii)],[L1 L2]);
    end
    L1 = L1 - 0.05;
    L2 = L2 - 0.05;
    xlim([1 (length(group5{1,1}{1,2}))]);
end

group5avg = plus(group5{1,1}{1,4},group5{1,2}{1,4});
for i = 3:length(group5)
    group5avg = plus(group5avg,group5{1,i}{1,4});
end
group5avg = resample(group5avg,1000,1);
plot((group5avg.*0.005)-0.4);
hline = refline(0,-0.4); %horizontal baseline
hline.Color = 'k';

group5RMSavg = plus(group5{1,1}{1,12}(1:208000),group5{1,2}{1,12}(1:208000));
for i = 3:length(group5)       %since we already added the first two together we then add the rest starting at the 3rd trial.
    group5RMSavg = plus(group5RMSavg,group5{1,i}{1,12}(1:208000));
end
plot((group5RMSavg*200)+1.6);
hline2 = refline(0,1.6); %horizontal baseline
hline2.Color = 'k';
ylim([-1 4]);

%labeling panel
x1 = 0;      
y1 = 2.9;
txt1 = 'Position trace';
text(x1,y1,txt1)

x1 = 0;
y1 = 2.4;
txt2 = 'High pass trace (50 hz cutoff)';
text(x1,y1,txt2)

x1 = 0;
y1 = 1.7;
txt3 = 'RMS average';
text(x1,y1,txt3)

x1 = 0;
y1 = 1.5;
txt4 = 'Raster plot';
text(x1,y1,txt4)

x1 = 0;
y1 = -0.3;
txt5 = 'Averaged firing frequency (hz)';
text(x1,y1,txt5)
%pause(1)

%%
ii = 1;
%subplot(2,1,2);
figure
for ii = 1:(length(group6))
    plot((group6{1,ii}{1,1}*0.4)+2.2);
    hold on
    plot((group6{1,ii}{1,6}*100)+2.3);
end
ii = 1;
L1 = 1.35;
L2 = 1.4;
for ii = 1:(length(group6))
    clear t1;
    hold on;
    iii = 1;
    t1 = group6{1,ii}{1,3};
    for iii = 1:length(t1)
        line([t1(iii) t1(iii)],[L1 L2]);
    end
    L1 = L1 - 0.05;
    L2 = L2 - 0.05;
    xlim([1 (length(group6{1,1}{1,2}))]);
end

group6avg = plus(group6{1,1}{1,4},group6{1,2}{1,4});
for i = 3:length(group6)
    group6avg = plus(group6avg,group6{1,i}{1,4});
end
group6avg = resample(group6avg,1000,1);
plot((group6avg.*0.005)-0.4);
hline = refline(0,-0.4); %horizontal baseline
hline.Color = 'k';

group6RMSavg = plus(group6{1,1}{1,12}(1:217000),group6{1,2}{1,12}(1:217000));
for i = 3:length(group6)       %since we already added the first two together we then add the rest starting at the 3rd trial.
    group6RMSavg = plus(group6RMSavg,group6{1,i}{1,12}(1:217000));
end
plot((group6RMSavg*200)+1.6);
hline2 = refline(0,1.6); %horizontal baseline
hline2.Color = 'k';
ylim([-1 4]);

%labeling panel
x1 = 0;      
y1 = 2.9;
txt1 = 'Position trace';
text(x1,y1,txt1)

x1 = 0;
y1 = 2.4;
txt2 = 'High pass trace (50 hz cutoff)';
text(x1,y1,txt2)

x1 = 0;
y1 = 1.7;
txt3 = 'RMS average';
text(x1,y1,txt3)

x1 = 0;
y1 = 1.5;
txt4 = 'Raster plot';
text(x1,y1,txt4)

x1 = 0;
y1 = -0.3;
txt5 = 'Averaged firing frequency (hz)';text(x1,y1,txt5)

%% Ending pause
pause(1)

