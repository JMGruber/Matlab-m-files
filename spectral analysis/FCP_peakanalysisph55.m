% USAGE: Just run, but don't forget to change the parameters first
%        Use together with:
%        - skewgaussfit4.m
%        - dblgaussfit.m
%        - skewgaussian1.m, skewgaussian2.m, skewgaussian3.m
%        - gaussfit.m, mingauss.m, gaussian.m
%        - bgsubtr.m
%        - subtrCR.m
%        - adjavg.m
%        - dblgaussian2.m - ???
%        - 
%
% This script resolves the spectra from a spectral sequence, currently tuned for Palustris
%
% OUTPUT 
% All intensities, peak positions, fwhm, skewness of double- or single-band spectra
% as well as averages of spectral peaks and spectral jumps
% 
% ALGORITHM
% The main problem is how to accurately identify and resolve a double-band spectrum
% 1. We first fit the average spectrum of a particular complex. If the width > fwhmdbl or skewness > skewdbl, go to 3
% 2. Then we run through all spectra of a particular complex. If the width > fwhmdbl or skewness > skewdbl, go to 3
%    If not, there's probably only one band, which makes life much easier
% 3. So we assume we're dealing with a double-band spectrum. 
%    There could be ealier double-band spectra which we missed. So we start from the beginning.
%    A double Gaussian is fitted as a first estimate.
% 4. The larger peak is resolved first by cutting off the other peak
% 4a. If the smaller peak of the double Gaussian fit is pushed against its
%     boundary values, this is usually a sign that it's too small to be
%     resolved well. In such a case the tail (with the small peak) is cut off.
% 4b. If 4a does not apply, the smaller peak (as resolved from the double
%     Gaussian fit) is subtracted from the data and the larger peak is resolved.
% 5. The smaller peak is subsequently resolved. If the intensity of the small peak 
%    is below a certain threshold we presume the complex jumped to a single-band state.
% 6. The misfits are refitted by using the data of good fits in the same spectral sequence.

% Output of fit results:
% -2: single peak
% -1: misfit
% 0: intensity too small

% Feb. 2011, Tjaart Kr�ger


clear all;

startfile=2; endfile=220;
skipfiles=[37 38 39 40 42 43 44 49 50 53 54 55 94 95 96 125 126 127 148 149 150 153 154 155 156 157 158 159 160 161 183 184 185 192 193 194 198 199 200 203 204 205 206 207 208 209 210 217 218 219 220 221 222 223 224 225]; %37 38 39 40 42 43 44 49 50 53 54 55 94 95 96 125 126 127 148 149 150 153 154 155 156 157 158 159 160 161 183 184 185 192 193 194 198 199 200 203 204 205 206 207 208 209 210 217 218 219 220 221 222 223 224 225numbers don't need to be in any particular order; matrix can also be empty
readdir='H:\Windows\Desktop\Maximes SMS DATA TREATED\fcpa pH5.5\FCPa pH55';
writedir='H:\Windows\Desktop\Maximes SMS DATA TREATED\fcpa pH5.5\FCPa pH55\pH5.5_1_Treated';

dblpeaks=true; %checks for double bands
fittype=3; %type of skewed Gaussian used for fitting
% fittype=3: FS69 function (Frazer and Suzuki 1969)
% fittype=2: Mikas's function (more flexible)
% fittype=1: Composite function

% first estimations for fitting
fwhm=18; 
skewness=0.2;
peak1=678; peak2=700; %typical peak positions for double-band spectra;

% boundary conditions for fitting
xmin=630; xmax=800; % maximum wavelength window considered in all calculations. Set 0 if full range is requireed.
Qthr=800; % if intensity is below this threshold, chances are good for a misfit. 
          %Also used in double-Gaussian fits to identify a single-peak spectrum
deadthr=800; %if the intensity of all spectra is below this threshold, we could consider it dead
minwavint=600; maxwavint=800; %wavelengths considered together with deadthr. Should be within (xmin,xmax)

fwhmdbl=30; %max fwhm of a single peak before the script considers double bands
fwhmblue=20; %used together with skewdbl for blue spectra to differentiate between a red tail and an extra red peak
skewdbl=0.5; %max skewness of a single peak before the script considers double bands
intdbl=2000; %intensity above which the script considers a double band
amplsingle=3;    %minimum ampl of one of the peaks of a double band, otherwise it's regarded as a single-band spectrum

peakmin=650; peakmax=800;
fwhmmin=15; fwhmmax=25;
skewmin=-.5; skewmax=2;

%bins
wavbin=0;   % wavelength binning: number of adjacent datapoints
wavbinsmall=3; % binning for small peaks
timebin=0;  % number of consecutive spectra to be binned  (not working for my version of Matlab)***
binFLP=1; %bin size of FL peak position (in nm)
binJumps=1; %bin size of jump sizes (in nm)
jumpthr=1; %jump size threshold (in nm) - two consecutive spectral peaks that differ by more than this number are considered a jump

% ***user can change values until here***

lbs = [0 fwhmmin peakmin skewmin];
ubs = [Inf fwhmmax peakmax skewmax]; 
ubs1 = [Inf fwhmdbl peakmax skewmax]; % only to choose between a single and double band 
lbd = [0 fwhmmin peakmin 0 fwhmmin (peak1+peak2)/2];
ubd = [Inf fwhmmax (peak1+peak2)/2 Inf fwhmmax peakmax];

deads=[];
doubles=[]; 
avg.int=[]; avg.fwhm=[]; avg.skew=[]; avg.peaks=[];
stdev.int=[]; stdev.fwhm=[]; stdev.skew=[]; stdev.peaks=[];
all.jumps=[]; all.jumps2=[];
markers=(xmin:binFLP:xmax);
histpeaks=zeros(size(markers));

for i=startfile:1:endfile
    all.int=[]; all.fwhm=[]; all.skew=[]; all.peaks=[]; pp=0;pp2=0;
    all.fwhm2=[]; all.skew2=[]; all.peaks2=[]; all.int2=[]; all.spec=[];
    j=1;
    goodfile=true;                                  %sort first......
    while (j<length(skipfiles))&&(goodfile)
        if i==skipfiles(j), goodfile=false; end
        j=j+1;
    end
    if goodfile
        mat=dlmread(strcat(readdir,'\spec',int2str(i)));
        nmin = find(mat(:,1)>xmin, 1);
        nmax = find(mat(:,1)>xmax, 1);
        datamat=[];
        for j=2:size(mat,2)
           datamat=[datamat mat(nmin:nmax,j)];
           mat(:,j)=subtrCR(mat(:,j),150);  %subtract intensity peaks due to cosmic rays
        end
        wav=mat(nmin:nmax,1);
        mat=[wav datamat];
        mat=bgsubtr(mat,30);    %subtract background
        nminwavint = find(mat(:,1)>minwavint, 1);
        nmaxwavint = find(mat(:,1)>maxwavint, 1);
        dead=false;
        if max(sum(mat(nminwavint:nmaxwavint,2:size(mat,2))))<deadthr
            dead=true;
            deads=[deads i];
        end
        if ~dead
            kk=size(mat,2);
            while (kk>2)&&(sum(mat(:,kk))<deadthr)
                kk=kk-1;
            end
            lastgood=kk;
            single=true; singlebroad=false;
            if lastgood>2
                avgdata=mean(mat(:,2:lastgood)')';
            else
                avgdata=mat(:,2);
            end
            single=true; singlebroad=false;
            avgdata=mean(mat(:,2:lastgood-timebin)')';
            matbinned=wav;
            avgfit=skewgaussfit4([wav,avgdata],[],fittype,lbs,ubs);
            if (dblpeaks)&&(sum(avgdata(nminwavint:nmaxwavint))>intdbl)
                if (avgfit(4)<-skewdbl)||((avgfit(2)>fwhmblue)&&(avgfit(4)>skewdbl))
                  single=false; 
                  avgdblref=dblgaussfit([wav avgdata],fwhm,peak1,fwhm,peak2,lbd,ubd);
                elseif avgfit(2)>fwhmdbl-.1 %differentiate between double band and broad single band
                    ub = [Inf 40 peakmax skewmax];
                    avgfit=skewgaussfit4([wav,avgdata],[],fittype,lbs,ub);
                    if (avgfit(4)<-skewdbl)||((avgfit(2)>fwhmblue)&&(avgfit(4)>skewdbl))
                      single=false; 
                      avgdblref=dblgaussfit([wav avgdata],fwhm,peak1,fwhm,peak2,lbd,ubd);
                    else
                        singlebroad=true;
                    end
                end
            end
            spec=1;
            singlespecdata=[];
            while (single)&&(spec<lastgood-timebin)
                data=mat(:,spec+1);
                int=sum(data);
                if int<Qthr*2
                    data=adjavg(data,wavbinsmall);
                else
                    data=adjavg(data,wavbin);
                end
                matbinned=[matbinned data];
                all.int=[all.int int];
                all.int2=[all.int2 -2];
                if (int<Qthr)
                    specdata=[0 0 0 0 0];
                    singlespecdata=[singlespecdata specdata'];
                else
                    if singlebroad
                        ub=[Inf 100 peakmax skewmax];
                        specdata=skewgaussfit4([wav,data],[],fittype,lbs,ub);
                    else
                        specdata=skewgaussfit4([wav,data],[],fittype,lbs,ubs1);
                    end
                    if (specdata(3)<peakmin+1)||(specdata(3)>peakmax-1)   %misfit
                        specdata=skewgaussfit4([wav,data],[fwhm,peak1,skewness],fittype,lbs,ubs1);
                        if (specdata(3)<peakmin+1)||(specdata(3)>peakmax-1)||(specdata(2)<fwhmmin+.1)||(specdata(2)>fwhmmax-.1)   %misfit
                            specdata=[-1 -1 -1 -1 -1]; 
                            singlespecdata=[singlespecdata specdata'];
                        end
                    else
                        if (dblpeaks)&&(~singlebroad)&&((specdata(2)>fwhmdbl-.1)||(specdata(4)<-skewdbl))
                            if sum(data(nminwavint:nmaxwavint))>intdbl
                                single=false;
                                avgdblref=dblgaussfit([wav data],fwhm,peak1,fwhm,peak2,lbd,ubd);
                            else
                                ub=[Inf 100 peakmax skewmax];    %low-ampl single peaks can be pretty broad
                                specdata=skewgaussfit4([wav,data],[],fittype,lbs,ub);
                                if length(specdata)==4, specdata(5)=0; end                                
                                singlespecdata=[singlespecdata specdata'];
                            end
                        else
                            if length(specdata)==4, specdata(5)=0; end
                            singlespecdata=[singlespecdata specdata'];
                        end
                    end
                end
                spec=spec+1;
            end
            if ~single
                all.int=[]; all.int2=[];
                for spec=1:lastgood-1-timebin
                    data=mat(:,spec+1);
                    int=sum(data);
                    if int<Qthr*2
                        data=adjavg(data,wavbinsmall);
                    else
                        data=adjavg(data,wavbin);
                    end
                    matbinned=[matbinned data];

                    if (int<Qthr)
                        specdata1=[0 0 0 0 0];
                        specdata2=[0 0 0 0 0];
                        int1=0; int2=0;                        
                    else
                        dblspecdata=dblgaussfit([wav data],avgdblref(2),avgdblref(3),avgdblref(5),avgdblref(6),lbd,ubd);
                        dbl=true;
                        if dblspecdata(1)<amplsingle    
                            specdata=skewgaussfit4([wav,data],[],fittype,lbs,ubs1); 
                            if (specdata(2)>fwhmdbl-.1)   % single-band spectrum
                                 specdata2=[-2 -2 -2 -2 -2];     
                                 specdata1=skewgaussfit4([wav,data],[fwhm,peak1,skewness],fittype,lbs,ubs);
                                 if (specdata1(2)>fwhmmax-.1)
                                    ub=[Inf fwhmmax+10 peakmax skewmax];
                                    specdata1=skewgaussfit4([wav,data],[fwhm,peak1,skewness],fittype,lbs,ub);
                                 end
                                 if pp>0, all.jumps=[all.jumps specdata1(3)-pp]; end
                                 if pp2>0, all.jumps2=[all.jumps2 specdata1(3)-pp2]; end
                                 pp=specdata1(3);pp2=0;
                                 dbl=false;
                            end
                        end
                        if dbl
                            if dblspecdata(1)/dblspecdata(4)>3
                                midwavx=(dblspecdata(3)+2*dblspecdata(6))/3;
                            elseif dblspecdata(4)/dblspecdata(1)>3
                                midwavx=(2*dblspecdata(3)+dblspecdata(6))/3;
                            else    
                                midwavx=(dblspecdata(3)+dblspecdata(6))/2;
                            end
                            midwav=find(wav>midwavx,1);
                            A1=max(data(1:midwav));
                            A2=max(data(midwav:length(data)));
                            if A1>=A2
                                ub=[Inf fwhmmax midwavx skewmax];
                                if (dblspecdata(5)<fwhmmin+.1)||(dblspecdata(5)>fwhmmax-.1)
                                    specdata1=skewgaussfit4([wav(1:midwav),data(1:midwav)],[dblspecdata(2),dblspecdata(3),skewness],fittype,lbs,ub); 
                                    singlepeak1=gaussian(wav,specdata1(1),specdata1(2),specdata1(3));
                                else
                                    singlepeak1=data-gaussian(wav,dblspecdata(4),dblspecdata(5),dblspecdata(6));
                                    specdata1=skewgaussfit4([wav,singlepeak1],[dblspecdata(2),dblspecdata(3),skewness],fittype,lbs,ub);
                                end
                                int1=sum(singlepeak1);
                                int2=int-int1;
                                if int2<Qthr
                                    specdata1=[-1 -1 -1 -1 -1];
                                else
                                    switch fittype
                                        case 1
                                            singlepeak2=data-skewgaussian1([specdata1(1),specdata1(2),specdata1(3),specdata1(4)],wav);
                                        case 2
                                            singlepeak2=data-skewgaussian2([specdata1(1),specdata1(2),specdata1(3),specdata1(4)],wav);
                                        case 3
                                            singlepeak2=data-skewgaussian3([specdata1(1),specdata1(2),specdata1(3),specdata1(4)],wav);
                                    end
                                    if int2<Qthr*2
                                        singlepeak2=adjavg(singlepeak2,wavbinsmall);
                                    end
                                    lb=[0 fwhmmin midwavx skewmin];
                                    specdata2=skewgaussfit4([wav,singlepeak2],[dblspecdata(5),dblspecdata(6),skewness],fittype,lb,ubs);
                                    border=midwav;
                                    borderx=wav(border);
                                    while (border>midwav/2)&&((specdata2(1)<amplsingle)||(specdata2(2)<fwhmmin+.1)||(specdata2(2)>fwhmmax-.1)||(specdata2(3)>peakmax-1)||(specdata2(3)<borderx+1))
                                        lb=[0 fwhmmin borderx skewmin];
                                        specdata2=skewgaussfit4([wav(border:length(wav)),singlepeak2(border:length(wav))],[dblspecdata(2),dblspecdata(3),skewness],fittype,lb,ubs);
                                        border=border-10;
                                        borderx=wav(border);
                                    end

                                    if (border<midwav/2)||(abs(specdata2(3)-dblspecdata(6))>10)     % misfit
                                        specdata2=[-1 -1 -1 -1 -1];
                                        ub=[Inf fwhmmax midwavx skewmax];
                                        specdata1=skewgaussfit4([wav,data],[fwhm,peak1,skewness],fittype,lbs,ub);
                                        if pp>0, all.jumps=[all.jumps specdata1(3)-pp]; end
                                        if pp2>0, all.jumps2=[all.jumps2 specdata1(3)-pp2]; end
                                        pp=specdata1(3);pp2=0;
                                    else  %double-band spectrum
                                        if pp>0, all.jumps=[all.jumps specdata1(3)-pp]; end
                                        if pp2>0, all.jumps2=[all.jumps2 specdata2(3)-pp2]; end
                                        pp=specdata1(3);pp2=specdata2(3);
                                    end
                                end
                            else
                                lb=[0 fwhmmin midwavx skewmin];
                                if (dblspecdata(2)<fwhmmin+.1)||(dblspecdata(2)>fwhmmax-.1)
                                    specdata2=skewgaussfit4([wav(midwav:length(wav)),data(midwav:length(data))],[dblspecdata(5),dblspecdata(6),skewness],fittype,lb,ubs); 
                                    singlepeak2=gaussian(wav,specdata2(1),specdata2(2),specdata2(3));
                                else
                                    singlepeak2=data-gaussian(wav,dblspecdata(1),dblspecdata(2),dblspecdata(3));
                                    specdata2=skewgaussfit4([wav,singlepeak2],[dblspecdata(5),dblspecdata(6),skewness],fittype,lb,ubs);
                                end
                                int2=sum(singlepeak2);
                                int1=int-int2;
                                if int1<Qthr
                                    specdata1=[-1 -1 -1 -1 -1];
                                else
                                    switch fittype
                                        case 1
                                            singlepeak1=data-skewgaussian1([specdata2(1),specdata2(2),specdata2(3),specdata2(4)],wav);
                                        case 2
                                            singlepeak1=data-skewgaussian2([specdata2(1),specdata2(2),specdata2(3),specdata2(4)],wav);
                                        case 3
                                            singlepeak1=data-skewgaussian3([specdata2(1),specdata2(2),specdata2(3),specdata2(4)],wav);
                                    end
                                    if int1<Qthr*2
                                        singlepeak1=adjavg(singlepeak1,wavbinsmall);
                                    end
                                    ub=[Inf fwhmmax midwavx skewmax];
                                    specdata1=skewgaussfit4([wav,singlepeak1],[dblspecdata(2),dblspecdata(3),skewness],fittype,lbs,ub);
                                    border=midwav;
                                    borderx=wav(border);
                                    while (border<(length(wav)+midwav)/2)&&((specdata1(1)<amplsingle)||(specdata1(2)<fwhmmin+.1)||(specdata1(2)>fwhmmax-.1)||(specdata1(3)>borderx-1)||(specdata1(3)<peakmin+1))
                                        ub=[Inf fwhmmax borderx skewmax];
                                        specdata1=skewgaussfit4([wav(1:border),singlepeak1(1:border)],[dblspecdata(2),dblspecdata(3),skewness],fittype,lbs,ub);
                                        border=border+10;
                                        borderx=wav(border);
                                    end
                                    if (border>(length(wav)+midwav)/2)||(abs(specdata1(3)-dblspecdata(3))>10)     % misfit
                                        specdata1=[-1 -1 -1 -1 -1];
                                        ub=[Inf fwhmmax midwavx skewmax];
                                        specdata2=skewgaussfit4([wav,data],[],fittype,lbs,ub);
                                        if pp>0, all.jumps=[all.jumps specdata2(3)-pp]; end
                                        if pp2>0, all.jumps2=[all.jumps2 specdata2(3)-pp2]; end
                                        pp=0;pp2=specdata2(3);
                                    else
                                        if pp>0, all.jumps=[all.jumps specdata1(3)-pp]; end
                                        if pp2>0, all.jumps2=[all.jumps2 specdata2(3)-pp2]; end
                                        pp=specdata1(3);pp2=specdata2(3);
                                    end 
                                end
                            end
                        end
                    end
                    all.fwhm=[all.fwhm specdata1(2)]; 
                    all.fwhm2=[all.fwhm2 specdata2(2)];
                    all.peaks=[all.peaks specdata1(3)];
                    all.peaks2=[all.peaks2 specdata2(3)];
                    all.skew=[all.skew specdata1(4)];
                    all.skew2=[all.skew2 specdata2(4)];
                    if exist('int1'), all.int=[all.int int1]; else all.int=[all.int int]; end
                    all.int2=[all.int2 int2];
                    all.spec=[all.spec spec];
                end
                % try to improve the misfits
                if sum(all.peaks>0 & all.peaks2>0)>0
                    misfits=all.spec(all.peaks==-1 | all.peaks2==-1);    
                    goodpeaksavg=mean(all.peaks(all.peaks>0 & all.peaks2>0));
                    goodpeaksstd=std(all.peaks(all.peaks>0 & all.peaks2>0));
                    goodpeaks2avg=mean(all.peaks2(all.peaks>0 & all.peaks2>0));
                    goodpeaks2std=std(all.peaks2(all.peaks>0 & all.peaks2>0));
                    goodfwhmavg=mean(all.fwhm(all.peaks>0 & all.peaks2>0));
                    goodfwhmstd=std(all.fwhm(all.peaks>0 & all.peaks2>0));
                    goodfwhm2avg=mean(all.fwhm2(all.peaks>0 & all.peaks2>0));
                    goodfwhm2std=std(all.fwhm2(all.peaks>0 & all.peaks2>0));
                    if goodpeaksstd>0
                        lb=[0 goodfwhmavg-goodfwhmstd goodpeaksavg-goodpeaksstd 0 goodfwhm2avg-goodfwhm2std goodpeaks2avg-goodpeaks2std];
                        ub=[Inf goodfwhmavg+goodfwhmstd goodpeaksavg+goodpeaksstd Inf goodfwhm2avg+goodfwhm2std goodpeaks2avg+goodpeaks2std]; 
                    else
                        lb=[0 goodfwhmavg-1 goodpeaksavg-1 0 goodfwhm2avg-1 goodpeaks2avg-1];
                        ub=[Inf goodfwhmavg+1 goodpeaksavg+1 Inf goodfwhm2avg+1 goodpeaks2avg+1]; 
                    end
                    for k=misfits(1:length(misfits))
                        data=mat(:,k+1);
                        int=sum(data);
                        if int<Qthr*2
                            data=adjavg(data,wavbinsmall);
                        else
                            data=adjavg(data,wavbin);
                        end
                        dblspecdata=dblgaussfit([wav data],goodfwhmavg,goodpeaksavg,goodfwhm2avg,goodpeaks2avg,lb,ub);
                        all.fwhm(k)=dblspecdata(2);
                        all.fwhm2(k)=dblspecdata(5);
                        all.peaks(k)=dblspecdata(3);
                        all.peaks2(k)=dblspecdata(6);
                        all.skew(k)=0;
                        all.skew2(k)=0;
                    end
                end
                                    
            else % if single-band spectrum
                singlespecdata=singlespecdata';
                for k=1:size(singlespecdata,1)
                    all.fwhm=[all.fwhm singlespecdata(k,2)]; 
                    all.peaks=[all.peaks singlespecdata(k,3)];
                    all.skew=[all.skew singlespecdata(k,4)];
                    all.fwhm2=[all.fwhm2 -2];
                    all.peaks2=[all.peaks2 -2];
                    all.skew2=[all.skew2 -2];
                end
                for k=2:size(all.peaks)
                    all.jumps=[all.jumps all.peaks(k)-allpeaks(k-1)];
                end
            end
                    
        end
    end
    if ~isempty(all.int)
        %only values of fitted peaks; ****** If preferred, we could use fit values of the average of all specs
        avg.int=[avg.int mean(all.int(max(all.peaks,all.peaks2)>0)) mean(all.int(max(all.peaks,all.peaks2)>0))]; %write it twice for the sake of matrix dimensions
        avg.fwhm=[avg.fwhm mean(all.fwhm(all.peaks>0)) mean(all.fwhm2(all.peaks2>0))];
        avg.peaks=[avg.peaks mean(all.peaks(all.peaks>0)) mean(all.peaks2(all.peaks2>0))];
        avg.skew=[avg.skew mean(all.skew(all.peaks>0)) mean(all.skew2(all.peaks2>0))];
        stdev.int=[stdev.int std(all.int(all.peaks>0)) std(all.int(all.peaks>0))];
        stdev.fwhm=[stdev.fwhm std(all.fwhm(all.peaks>0)) std(all.fwhm2(all.peaks2>0))];
        stdev.peaks=[stdev.peaks std(all.peaks(all.peaks>0)) std(all.peaks2(all.peaks2>0))];
        stdev.skew=[stdev.skew std(all.skew(all.peaks>0)) std(all.skew(all.peaks>0))];

        histpeaks=histpeaks+histc([all.peaks(all.peaks>0) all.peaks2(all.peaks2>0)],markers);
        %write data to matrix
        dlmwrite(strcat(writedir,'\','dataspec',int2str(i)),[all.int' all.peaks' all.fwhm' all.skew' all.int2' all.peaks2' all.fwhm2' all.skew2'],'\t');
        dlmwrite(strcat(writedir,'\','allints'),[i all.int],'-append','delimiter','\t');    
        dlmwrite(strcat(writedir,'\','allpeaks'),[i all.peaks],'-append','delimiter','\t');    
        dlmwrite(strcat(writedir,'\','allfwhm'),[i all.fwhm],'-append','delimiter','\t');    
        dlmwrite(strcat(writedir,'\','allskew'),[i all.skew],'-append','delimiter','\t');   
        dlmwrite(strcat(writedir,'\','allints2'),[i all.int2],'-append','delimiter','\t');            
        dlmwrite(strcat(writedir,'\','allpeaks2'),[i all.peaks2],'-append','delimiter','\t');    
        dlmwrite(strcat(writedir,'\','allfwhm2'),[i all.fwhm2],'-append','delimiter','\t');    
        dlmwrite(strcat(writedir,'\','allskew2'),[i all.skew2],'-append','delimiter','\t');        
        dlmwrite(strcat(writedir,'\','specminbg',int2str(i)),matbinned,' ');       

        %write data to single column
        dlmwrite(strcat(writedir,'\','convert_allints'),all.int','-append');    
        dlmwrite(strcat(writedir,'\','convert_allpeaks'),all.peaks','-append');    
        dlmwrite(strcat(writedir,'\','convert_allfwhm'),all.fwhm','-append');    
        dlmwrite(strcat(writedir,'\','convert_allskew'),all.skew','-append');   
        dlmwrite(strcat(writedir,'\','convert_allints2'),all.int2','-append');            
        dlmwrite(strcat(writedir,'\','convert_allpeaks2'),all.peaks2','-append');    
        dlmwrite(strcat(writedir,'\','convert_allfwhm2'),all.fwhm2','-append');    
        dlmwrite(strcat(writedir,'\','convert_allskew2'),all.skew2','-append');
    end
    disp(i);
end
if ~isempty(all.int)
    a=dlmread(strcat(writedir,'\','allints'))'; dlmwrite(strcat(writedir,'\','allints'),a,'\t');
    a=dlmread(strcat(writedir,'\','allpeaks'))'; dlmwrite(strcat(writedir,'\','allpeaks'),a,'\t');
    a=dlmread(strcat(writedir,'\','allfwhm'))'; dlmwrite(strcat(writedir,'\','allfwhm'),a,'\t');
    a=dlmread(strcat(writedir,'\','allskew'))'; dlmwrite(strcat(writedir,'\','allskew'),a,'\t');
    a=dlmread(strcat(writedir,'\','allints2'))'; dlmwrite(strcat(writedir,'\','allints2'),a,'\t');
    a=dlmread(strcat(writedir,'\','allpeaks2'))'; dlmwrite(strcat(writedir,'\','allpeaks2'),a,'\t');
    a=dlmread(strcat(writedir,'\','allfwhm2'))'; dlmwrite(strcat(writedir,'\','allfwhm2'),a,'\t');
    a=dlmread(strcat(writedir,'\','allskew2'))'; dlmwrite(strcat(writedir,'\','allskew2'),a,'\t');

    dlmwrite(strcat(writedir,'\','avgspec'),[avg.int' avg.peaks' avg.fwhm' avg.skew'],'\t');
    dlmwrite(strcat(writedir,'\','stdspec'),[stdev.int' stdev.peaks' stdev.fwhm' stdev.skew'],'\t');

    peakmin=floor(min(avg.peaks)/binFLP)*binFLP;
    peakmax=ceil(max(avg.peaks)/binFLP)*binFLP;
    markavgpeak=peakmin:binFLP:peakmax;
    figure; subplot(2,2,1); hist(avg.peaks,markavgpeak);
    xlabel('avg FLP (nm)'); ylabel('# spectra');
    dlmwrite(strcat(writedir,'\','histavgpeaks'),[markavgpeak' histc(avg.peaks,markavgpeak)'],'\t');

    markstdpeak=0:2:max(stdev.peaks);
    subplot(2,2,2); hist(stdev.peaks,10);
    xlabel('avg SDFP (nm)'); ylabel('# spectra');
    dlmwrite(strcat(writedir,'\','histstdpeaks'),[markstdpeak' histc(stdev.peaks,markstdpeak)'],'\t');

    subplot(2,2,3); bar(markers,histpeaks); axis tight;
    xlabel('all FLP (nm)'); ylabel('# spectra');
    dlmwrite(strcat(writedir,'\','histallpeaks'),[markers' histpeaks'],'\t');

    all.jumps=[all.jumps all.jumps2];

    if ~isempty(all.jumps)
        minjump=floor(min(all.jumps)/binJumps)*binJumps-jumpthr/2; %last term for symmetry around 0
        maxjump=ceil(max(all.jumps)/binJumps)*binJumps+jumpthr/2;
        markjumps=minjump:binJumps:maxjump;
        subplot(2,2,4); hist(all.jumps(abs(all.jumps)>jumpthr),markjumps);
        xlabel('Jump size (nm)'); ylabel('Occurrence');
        dlmwrite(strcat(writedir,'\','histalljumps'),[markjumps' histc(all.jumps,markjumps)'],'\t');
    end

    allints=dlmread(strcat(writedir,'\','convert_allints'))';
    allpeaks=dlmread(strcat(writedir,'\','convert_allpeaks'))';
    allfwhm=dlmread(strcat(writedir,'\','convert_allfwhm'))';
    allskew=dlmread(strcat(writedir,'\','convert_allskew'))';
    allints2=dlmread(strcat(writedir,'\','convert_allints2'))';
    allpeaks2=dlmread(strcat(writedir,'\','convert_allpeaks2'))';
    allfwhm2=dlmread(strcat(writedir,'\','convert_allfwhm2'))';
    allskew2=dlmread(strcat(writedir,'\','convert_allskew2'))';

% colour code:
% red: redder peak of double-band spectra
% blue: bluer peak of double-band spectra
% black: single-band spectra
    figure; 
    subplot(2,2,1); plot(allfwhm(allpeaks>0),allskew(allpeaks>0),'b.'); 
    hold on; plot(allfwhm2(allpeaks2>0),allskew2(allpeaks2>0),'r.');
    plot(allfwhm(allpeaks2==-2&allpeaks>0),allskew(allpeaks2==-2&allpeaks>0),'k.');
    plot(allfwhm2(allpeaks==-2&allpeaks2>0),allskew2(allpeaks==-2&allpeaks2>0),'k.');
    xlabel('fwhm (nm)'); ylabel('skewness');
    subplot(2,2,2); plot(allpeaks(allpeaks>0),allskew(allpeaks>0),'b.'); 
    hold on; plot(allpeaks2(allpeaks2>0),allskew2(allpeaks2>0),'r.');
    plot(allpeaks(allpeaks2==-2&allpeaks>0),allskew(allpeaks2==-2&allpeaks>0),'k.');
    plot(allpeaks2(allpeaks==-2&allpeaks2>0),allskew2(allpeaks==-2&allpeaks2>0),'k.');
    xlabel('FLP (nm)'); ylabel('skewness');
    subplot(2,2,3); plot(allpeaks(allpeaks>0),allfwhm(allpeaks>0),'b.'); 
    hold on; plot(allpeaks2(allpeaks2>0),allfwhm2(allpeaks2>0),'r.');
    plot(allpeaks(allpeaks2==-2&allpeaks>0),allfwhm(allpeaks2==-2&allpeaks>0),'k.');
    plot(allpeaks2(allpeaks==-2&allpeaks2>0),allfwhm2(allpeaks==-2&allpeaks2>0),'k.');
    xlabel('FLP (nm)'); ylabel('fwhm (nm)');
    subplot(2,2,4); plot(allpeaks(allpeaks>0),allints(allpeaks>0),'b.'); 
    hold on; plot(allpeaks2(allpeaks2>0),allints2(allpeaks2>0),'r.');
    plot(allpeaks(allpeaks2==-2&allpeaks>0),allints(allpeaks2==-2&allpeaks>0),'k.');
    plot(allpeaks2(allpeaks==-2&allpeaks2>0),allints2(allpeaks==-2&allpeaks2>0),'k.');
    xlabel('FLP (nm)'); ylabel('intensity (cps)');
end
save(strcat(writedir,'\','allvariables.mat'));