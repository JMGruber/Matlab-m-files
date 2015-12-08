function DrawIntLevels(ti,tf,trace,intlevels,inttimes,intstart)
% Usage: Draw intensity levels of "trace" within the time range ti:tf (in
% seconds) as calculated by IntShifts.m
% intlevels = all resolved intensity levels
% inttimes = dwell time corersponding to each resolved level
% intstart = starting indices in "trace" of resolved levels

timeres=trace(1,1);
if ti < timeres
    ti = timeres;
end
if tf > intstart(end)
    tf = intstart(end);
end
ri = find(intstart<=ti,1,'last');        
rf = find(intstart<tf,1,'last');
range = round(intstart(ri)/timeres):round((intstart(rf)+inttimes(rf))/timeres);
figure; plot(trace(range,1),trace(range,2),'g');
hold on;

for j = ri:rf
  x = [intstart(j) intstart(j)+inttimes(j)]; 
  y = [intlevels(j) intlevels(j)];
  plot(x,y,'LineWidth',4,'Color','k');
end
hold off;