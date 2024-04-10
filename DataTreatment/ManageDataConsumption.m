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
mfilename = [mfilename '\Consumption'];
addpath(genpath(mfilename));

CNP16OV_tab = readtable('CNP16OV.csv'); CNP16OV_tab(end, :) = [];
PCECC96_tab = readtable('PCECC96_new_new.csv');
PCECC96_tab(1:4,:) = []; % remove the first year (consumption not available for that)
CNP16OV = table2array(CNP16OV_tab(:,2));
PCECC96 = table2array(PCECC96_tab(:,2));
cc = diff(log(PCECC96./CNP16OV));
dates_cc = [1948.25:0.25:2021.50]';
cc_tab = table(dates_cc, cc);

% SAVE
filename = [pwd '\OUTPUT' '\cc'];
save ( filename,...
'cc_tab'...
);