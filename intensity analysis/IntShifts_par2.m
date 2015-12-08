% Description 
% Paper
%
% First execute matlabpool to activate parallel computing
% Code written for Parallel Computing: n CPUs will decrease execution time
% by a factor n.
% 
% Run together with:
% - IntShifts_Algorithm.m
% - IntShifts_binints.m
% - IntShifts_CR.m
% - IntShifts_P4vars.m
% - IntShifts_TestSM.m
% - IntShifts_TrimTrace.m
% - DrawIntLevels.m
% 
% Tjaart Krüger
% Vrije Universiteit Amsterdam
% Aug 2011
%
% Improvements: 
% - Don't specify all k's separately but calculate all for 1
% user-defined sigma deviation
% - fprintf: 1...n possible to display?
% - Don't need all temp parameters with variable sizes: calculations cost time.
% Better to define large vectors/matrices and remove all the zeros in the end!

%% Initialisation
tic;
clear all;
startfile=1; endfile=300; % files must be numbered consecutively!
skipfiles=[13 44 49 51 57 59 60 63 69 81 83 130 133 141 168 175 176 177 185 189 193 197 201 203 209 220 223 225 226 228 246 259 269 273 278 280 284 287 291 293 294 130 180 186 218];
usefiles=[1:100];
readdir='D:\LHC2\TROLOX\030912 GOC RT';
writedir='D:\LHC2\TROLOX\030912 GOC RT\ana1';

ki = [3.5 1.96 1.37 1]; % sigma-deviation for 1..4 consecutive trace points, resp.
Nf = 1.1; % noise multiplication factor (1.0 - 1.1 for good SNR)

thrliveI = 5000;  % threshold for final intensity level (in cps, background included)
thrlivet = 1;  % minimal survival time (in seconds)
thr2c = 12000;  % intensity threshold when definitely >1 complex (in cps, background included)
intbin = 1;   % number of intensity values to be binned together
maxbg = 1400; % max estimated background; used to test for SM (in cps)
testSM = true; % make sure this is set to "false" when your complexes exhibit very steady fluorescence (i.e. hardly switch to background level)
testphotonburst = true;
photonburstfactor = 10;   % factor by which dwelling time in Q state exceeds dwelling time in unQ states to define photonburst (typically 10)
intermedfactor = 1; % parameter used to test if there are too many unnatural fluctuations (signifying unstable complex)

binsize = 200;  % binning intensity levels for plots (in cps)
Qextent = 0.2;  % extent of quenching (for panel 4)

figtitle = 'LHCII, pH 5.5';

drawlevels = true;  % draw intensity levels of the first indexed complex
drawti = 0;    % start time (in seconds)
drawtf = 60;    % finish time (in seconds)

maxsize = 30/0.01/2;  % max size of any of the data sets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           User may change parameter values until here                 %%%                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate matrix sizes
nfiles=endfile-startfile+1;
sortedtraces=zeros(nfiles,5);  % [1..5] = [allusedtraces alldead alldble allphotonburst allusedtracebg]
allnc=zeros(nfiles,1); % number of intensity levels for each complex
alltauon=zeros(maxsize,nfiles);
alltauoff=zeros(maxsize,nfiles);
allinton=zeros(maxsize,nfiles);
allintoff=ones(maxsize,maxsize)*10*thr2c; % since int can be 0 after bg subtraction
allint=ones(maxsize,nfiles)*10*thr2c;
alldwelltime=zeros(maxsize,nfiles);
allstarttime=zeros(maxsize,nfiles);

alltlive=0;
%readfile = fullfile(readdir,['trace' int2str(startfile)]);  % why - only for timeres in end?????
%firsttrace = dlmread(readfile);                             % in that case you can use last file that was opened

%% Core
for i=1:nfiles

    
       j=1;

    goodfile=true;
    while (j<=length(skipfiles))&&(skipfiles(j)<=i)&&(goodfile)
        if i==skipfiles(j)
            goodfile=false;
        end
        j=j+1;
    end
%         while (j<=length(usefiles))&&(usefiles(j)<=i)
%         if i==usefiles(j)
%             goodfile=true;
%         end
%         j=j+1;
%     end
    if goodfile
    tr = startfile+i-1;
    tempindex = zeros(1,5);
    if exist(fullfile(readdir,['trace' int2str(tr)]));
        
    readfile = fullfile(readdir,['trace' int2str(tr)]);
    trace=dlmread(readfile);
   % trace=dlmread(strcat(readdir,'\trace (',int2str(tr),')'));
    if intbin>1     % bin intensity values
        trace=IntShifts_binints(trace,intbin);
    end
    temptimeres=trace(1,1);
    trace(:,2)=trace(:,2)*temptimeres;
    
    tlive=IntShifts_TrimTrace(trace,thrliveI); % trim trace into surviving time
    if tlive<thrlivet/temptimeres
        tempindex(2)=tr; 
    else
       trace=IntShifts_CR(trace,thr2c*temptimeres,tlive); % remove excessively large intensities
       dble=false; photonburst=false;
       if testSM
           SM=IntShifts_testSM(trace,maxbg,tlive,thrliveI);
           if ~SM
              dble=true;
              tempindex(3)=tr;
           end
       end
       if (testphotonburst)&&((length(find(trace(1:tlive,2)<thrliveI*temptimeres))>photonburstfactor*length(find(trace(1:tlive,2)>thrliveI*temptimeres)))...
                            ||(sum(trace(1:tlive,2)<intermedfactor*thrliveI*temptimeres & trace(1:tlive,2)>maxbg*temptimeres)... % intermediate intensities
                            > sum(trace(1:tlive,2)>intermedfactor*thrliveI*temptimeres | trace(1:tlive,2)<maxbg*temptimeres)))       % are more than Q and unQ
            photonburst=true;
            tempindex(4)=tr; 
       end
        
        if (~photonburst)&&(~dble)
                if i==14
    
                end
            [intlevels,inttimes,intstart,SM] = IntShifts_Algorithm(trace,ki,Nf,tlive,thr2c*temptimeres);
            if ~isempty(intlevels)
            if ~SM
                tempindex(3)=tr; 
                dble=true;
            end

            % capture also very last int level??
            
            % estimate background
%             bg=1;
%             bgi=1;

            if i==49
                
            end

            if length(intlevels)==1
                bg=sum(trace(tlive+1:tlive+20,2))/21;
                bgi=1;
            else
            [bg,bgi] = min(intlevels);
            end
            
            if inttimes(bgi)<=2    % short dwell times can represent bad estimations
                levels2 = intlevels;
                times2 = inttimes;
                while (bgi>1)&&(bgi<length(levels2))&&(length(levels2)>3)&&(times2(bgi)<=2)
                    levels2 = [levels2(1:bgi-1) levels2(bgi+1:length(levels2))];
                    times2 = [times2(1:bgi-1) times2(bgi+1:length(times2))];
                    [bg,bgi] = min(levels2);
                end
            end
            
                   
            tempindex(5) = bg;
            intlevels = intlevels-bg;
           
                            
            inttimes = inttimes.*temptimeres;
            intstart = intstart.*temptimeres;
            alltlive = alltlive+tlive;
            
            % calculations for 4th panel of figure
            if ~dble
                tempalltauon = zeros(maxsize,1);
                tempalltauoff = zeros(maxsize,1);
                tempallinton = zeros(maxsize,1);
                tempallintoff = zeros(maxsize,1);
                tempallint = ones(maxsize,1)*10*thr2c;
                tempalldwelltime = zeros(maxsize,1);
                tempallstarttime = zeros(maxsize,1);
                
                P4vars = IntShifts_P4vars(intlevels,inttimes,thrliveI*temptimeres-bg,Qextent);
                l = size(P4vars,1);
                if l>0
                tempalltauon(1:l,1) = P4vars(1:end,1);     
                alltauon(:,i) = tempalltauon;
                tempalltauoff(1:l,1) = P4vars(1:end,2);     
                alltauoff(:,i) = tempalltauoff;
                tempallinton(1:l,1) = P4vars(1:end,3);     
                allinton(:,i) = tempallinton;
                tempallintoff(1:l,1) = P4vars(1:end,4);     
                allintoff(:,i) = tempallintoff;
                end
                
                l = length(intlevels);
                tempallint(1:l,1) = intlevels;
                allint(:,i) = tempallint;
                tempalldwelltime(1:l,1) = inttimes;
                alldwelltime(:,i) = tempalldwelltime;
                tempallstarttime(1:l,1) = intstart(1:l);
                allstarttime(:,i) = tempallstarttime;
                
                tempindex(1) = tr; 
                allnc(i) = length(intlevels);
            end
            end
        end
    end
    sortedtraces(i,:) = tempindex;
    fprintf(1,'%6.2f\n',tr);
    else
    end
    end
end
% remove zeros from matrices
allusedtraces=sortedtraces(sortedtraces(:,1)>0,1);
alldead=sortedtraces(sortedtraces(:,2)>0,2);
alldble=sortedtraces(sortedtraces(:,3)>0,3);
allphotonburst=sortedtraces(sortedtraces(:,4)>0,4);
allusedtracebg=sortedtraces(sortedtraces(:,5)>0,5);
allint=allint(allint<(10*thr2c));
alldwelltime=alldwelltime(alldwelltime>0);
allstarttime=allstarttime(allstarttime>0);
allnc=allnc(allnc>0);
alltauon=alltauon(alltauon>0);
alltauoff=alltauoff(alltauoff>0);
allinton=allinton(allinton>0);
allintoff=allintoff(allintoff<10*thr2c);

%% Plot output
if ~isempty(allint)
    tracenr=6;
    firstfile=allusedtraces(tracenr);
    readfile = fullfile(readdir,['trace' int2str(firstfile)]);
    firsttrace = dlmread(readfile);                           
    timeres=firsttrace(1,1);
    
    IntLim1=sum(allnc(1:tracenr-1));
    IntLim2=sum(allnc(1:tracenr-1))+allnc(tracenr);

    % Draw intensity levels
    if drawlevels
        if intbin>1     % bin intensity values
            firsttrace=IntShifts_binints(firsttrace,intbin);
        end
        firsttrace(:,2)=firsttrace(:,2)*timeres-sortedtraces(firstfile,5);
        DrawIntLevels(drawti,drawtf,firsttrace,allint(IntLim1:IntLim2),alldwelltime(IntLim1:IntLim2),allstarttime(IntLim1:IntLim2));
    end

    % Global overview of switches
    figure; subplot(2,2,1); semilogx(alldwelltime,allint,'.');
    xlabel('Dwell time (s)'); ylabel(['Intensity (c/',int2str(timeres*1000),' ms)']);
    title(figtitle);

    binsize=binsize*timeres; % Take note!
    maxbin=ceil(max(allint)/binsize)*binsize;
    %maxbin=ceil(thr2c*timeres/binsize)*binsize;
    minbin=ceil(min(allint)/binsize)*binsize;

    markers=minbin:binsize:maxbin;
    m=length(markers);
    bin.times=zeros(1,m);
    bin.n2=zeros(1,m);
    for i=1:length(allint) 
        p=ceil((allint(i)-minbin)./binsize);
        if (p<=0)||(isnan(p)), p=1; end
%        if p>maxbin, p=maxbin; end
        bin.times(p)=bin.times(p)+alldwelltime(i);
        bin.n2(p)=bin.n2(p)+1;
    end
    bin.times=bin.times./(alltlive*timeres)*100;
    bin.n2=bin.n2./(alltlive*timeres)*60;
    subplot(2,2,2); bar(markers,bin.n2); 
    xlabel(['Intensity (c/',int2str(timeres*1000),' ms)']); ylabel('Access rate (min^{-1})');
    subplot(2,2,3); bar(markers,bin.times);
    xlabel(['Intensity (c/',int2str(timeres*1000),' ms)']); ylabel('Total dwell time (%)');

    edges=timeres*(exp(0:log(2):9));
    tausortoff=histc(alltauoff,edges)./(alltlive*timeres)*60;
    if ~isempty(tausortoff)
        subplot(2,2,4); semilogx(edges,tausortoff,'.-'); ylabel('Switching frequency (min^{-1})'); xlabel('Dwell time (s)');
    end
    tausorton=histc(alltauon,edges)./(alltlive*timeres)*60;
    if ~isempty(tausorton)
        hold on;
        subplot(2,2,4); semilogx(edges,tausorton,'o-');
    end
end

%% File output
%dlmwrite(strcat(writedir,'\','parameters'),['bg  ',int2str(bg)],'');    
if ~isempty(allusedtraces)
    dlmwrite(strcat(writedir,'\','parameters'),['k1 k2 k3 k4 Nf ',int2str(ki(1)),' ',int2str(ki(2)),' ',int2str(ki(3)),' ',int2str(ki(4)),' ',int2str(Nf)],'');
    dlmwrite(strcat(writedir,'\','parameters'),['thrliveI thrlivet  ',int2str(thrliveI),' ',int2str(thrlivet)],'-append','delimiter','');
    dlmwrite(strcat(writedir,'\','parameters'),['thr2c  ',int2str(thr2c)],'-append','delimiter','');
    dlmwrite(strcat(writedir,'\','parameters'),['binsize Qextent  ',int2str(binsize),' ',int2str(Qextent)],'-append','delimiter','');

    dlmwrite(strcat(writedir,'\','panel1'),[alldwelltime;allint]',' ');  
    dlmwrite(strcat(writedir,'\','panel23'),[markers;bin.n2;bin.times]',' ');
    dlmwrite(strcat(writedir,'\','panel4'),[edges;tausortoff';tausorton']',' ');

    dlmwrite(strcat(writedir,'\','complexes used'),allusedtraces',''); 
    dlmwrite(strcat(writedir,'\','complexes dead'),alldead',''); 
    dlmwrite(strcat(writedir,'\','complexes double'),alldble',''); 

    save(strcat(writedir,'\analysis.mat'));
end
toc