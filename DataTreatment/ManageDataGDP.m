% =================================================================== %
%  
% =================================================================== %
clear
close all 
clc

mfilename = pwd;
mfilename = [mfilename '\Population'];
addpath(genpath(mfilename));
mfilename = pwd;
mfilename = [mfilename '\RealGDP'];
addpath(genpath(mfilename));


GDP_tab = readtable('GDPC1.csv'); 
CNP16OV_tab = readtable('CNP16OV.csv'); % population
CNP16OV_tab(1:47, :) = []; CNP16OV_tab(end, :) = [];
CNP16OV = table2array(CNP16OV_tab(:,2));
GPDP = table2array(GDP_tab(:,2));
GDP = diff(log(GPDP./CNP16OV));
dates_GDP = [1960:0.25:2021.50]';
GDP_tab = table(dates_GDP, GDP);


% SAVE
filename = [pwd '\OUTPUT' '\GDP'];
save ( filename,...
'GDP_tab'...
);