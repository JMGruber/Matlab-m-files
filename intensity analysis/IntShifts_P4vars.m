function P4vars = IntShifts_P4vars(intlevels,inttimes,thrliveI,Qextent)
% Usage: Calculate intensity fluctuations that exceed the factor "Qextent".
% Called by IntShifts.m and plotted as Panel 4
% Output: P4vars = [tauon' tauoff' inton' intoff'];
%                = ["on" dwelltimes, "off" dwelltimes, corresponding intensities]

avgallint=sum(intlevels)/length(intlevels);
tauon=[]; tauoff=[]; inton=[]; intoff=[];
tau_Qint=find(intlevels<min(avgallint,thrliveI)); %starting point
j=1; prev=0;
while j<length(tau_Qint)
    Qs=intlevels(tau_Qint(j)); % Q levels
    Qsno=tau_Qint(j);           % indices of Q levels
    while (length(tau_Qint)>j)&&(tau_Qint(j+1)-tau_Qint(j)==1)
        Qs=[Qs intlevels(tau_Qint(j+1))]; % for avg Q levels
        Qsno=[Qsno tau_Qint(j+1)];
        j=j+1;
    end
    if Qsno(1)>1
        Qsmin=min(Qs);      % lowest intensity level
        unQs=intlevels(prev+1:Qsno(1)-1);
        Qsmax=max(unQs);    % highest intensity level

        Qtime=0;
        Qintall=[]; tauoffall=[];
        for i=1:length(Qs)
            if Qs(i)/Qsmax<1-Qextent
                Qtime=Qtime+inttimes(Qsno(i));
                tauoffall=[tauoffall inttimes(Qsno(i))];
                Qintall=[Qintall Qs(i)];
            end
        end
        if Qtime>0
            Qint=sum(tauoffall.*Qintall)./Qtime;   % weighted intensity
            tauoff=[tauoff Qtime];
            intoff=[intoff Qint];
        end

        unQtime=0;
        unQintall=[]; tauonall=[];
        for i=1:length(unQs)
            if Qsmin/unQs(i)<1-Qextent
                unQtime=unQtime+inttimes(prev+1);
                tauonall=[tauonall inttimes(prev+1)];
                unQintall=[unQintall unQs(i)];
            end
        end
        if unQtime>0
            unQint=sum(tauonall.*unQintall)./unQtime;   % weighted intensity
            tauon=[tauon unQtime];
            inton=[inton unQint];
        end
    end
    prev=tau_Qint(j);
    j=j+1;

end
P4vars=[tauon' tauoff' inton' intoff'];