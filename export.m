function [] = export(writedir,parameter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dlmwrite(strcat(writedir,'\',parameter),parameter,'');

end

