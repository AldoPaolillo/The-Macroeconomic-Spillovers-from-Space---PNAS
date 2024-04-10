% =================================================================== %
% 
% 
% =================================================================== %
clear
close all
clc

% Generate the series for Aerospace industrial production

mfilename = pwd;
mfilename = [mfilename '\AerospaceIndustrialProduction'];
addpath(genpath(mfilename));
IPG3364S_tab = readtable('IPG3364S.csv');
IPG3364S = str2double(table2array(IPG3364S_tab(:,2)))/100;
load('ASindprod_imp.mat')
dates_ASindprod = [1960:0.25:2022]';

ASindprod = NaN(size(dates_ASindprod));

index_beginning = find(dates_ASindprod == 1972.25 );
index_end = find(dates_ASindprod == 2022 );

ASindprod(1:index_beginning-1) = table2array(ASindprod_imp_tab(1:index_beginning-1, 2));
ASindprod(index_beginning:index_end) = IPG3364S;

ASindprod_tab = table(dates_ASindprod, ASindprod);
ASindprod_tab(end, :) = [];


figure()
plot(dates_ASindprod, ASindprod, '-o', 'linewidth',1.4)



% SAVE
filename = [pwd '\OUTPUT' '\ASindprod'];



save ( filename,...
'ASindprod_tab'...
);