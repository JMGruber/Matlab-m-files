function [intlevels,inttimes,intstart,SM] = IntShifts_Algorithm(trace,ki,Nf,tlive,thr2c)

SM=true;
intstart=1;
intlevels=[];
inttimes=[];
for t=2:tlive
    for i=t:-1:max(intstart(length(intstart))+1,t-8)  % 2nd part to enable fewer iterations
        k=ki(4);
        I1v=trace(intstart(length(intstart)):i-1,2);
        I1=sum(I1v)/length(I1v); % faster calculation of the mean
        sigma=sqrt(I1)*Nf; % noise
        I2v=trace(i:t,2); % values in I2
        I2=sum(I2v)/length(I2v);
        if min(length(I1),length(I2))<4
            if min(length(I1),length(I2))<3
                if min(length(I1),length(I2))<2
                    k=ki(1);
                else k=ki(2);
                end
            else k=ki(3);
            end
        end
        if abs(I2-I1)>k*sigma      % intensity jump identified
            if I2>thr2c     % check if intensity represents 2 complexes
                SM=false;
                break
            else
                % check transitory point;
                transit=false;  %0: in favour of I1; 1: in favour of I2
                if length(I1v)>1
                    if (k==ki(3))||(k==ki(4))
                        if ((I2<I1)&&(I1v(length(I1v))<0.5*(I1-I2)))||((I2>I1)&&(I1v(length(I1v))>0.5*(I2-I1)))
                            transit=true;
                        end
                    else
                        sigma2=sqrt(I2);
                        if ((I2<I1)&&(I1v(length(I1v))-I2<ki(2)*sigma2))||((I2>I1)&&(I2-I1v(length(I1v))<ki(3)*sigma2))
                            transit=true;
                        end
                    end
                end
                if transit
                    I2v=[I1v(length(I1v)) I2v'];
                    I2=sum(I2v)/length(I2v);
                    I1v=I1v(1:length(I1v)-1);
                    I1=sum(I1v)/length(I1v);
                end
                % save I1
                intlevels=[intlevels I1];
                inttimes=[inttimes i-intstart(length(intstart))];
                intstart=[intstart i];
                I1=I2;
                break
            end
        end
    end
end