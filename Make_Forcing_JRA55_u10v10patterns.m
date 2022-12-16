% This script makes forcing files for MOM-AnENSO runs.
%
% This version makes forcing files from JRA55 for JRA55 v1.3 idealized Nino 3.4 IAF run.

clear all;
close all;

addpath(genpath('/g/data/e14/rmh561/software/matlab-utilities/'));
startup;

% Source data:
source_type = 1; % 0 = ERA, 1 = JRA
if (source_type)
    base = 'JRAdata/';
    label = 'JRA55-do';
else
    base = 'ERAdata/';
    label = 'ERA Interim';
end

% Output data:
out_type = 1; % 0 = ERA, 1 = JRA

%%% Get SST indices:

% $$$ % Nino 3.4:
% $$$ % OISST:
% $$$ DATA = load('index_data/NinoIndices.txt');
% $$$ n34 = DATA(:,10);
% $$$ n34yr = DATA(:,1);
% $$$ n34mn = DATA(:,2);
% $$$ n34yrd = n34yr + n34mn/12;

% $$$ % ERSST:
% $$$ DATA = load('index_data/nino34.data_ERSSTraw');
% HadISST:
DATA = load('index_data/nino34.long.data_HadISSTraw');
n34nyrs = length(DATA(:,1));
n34 = reshape(DATA(:,2:end)',[n34nyrs*12 1]);
n34(abs(n34)>50) = NaN;
n34yr = reshape(repmat(DATA(:,1),[1 12])',[n34nyrs*12 1]);
n34mn = repmat([1:12]',[n34nyrs 1]);
n34yrd = n34yr + n34mn/12;
n34cli = zeros(12,1);
for mi = 1:12
    n34cli(mi) = nanmean(n34(n34mn == mi & (n34yr >= 1981 & n34yr <=2010)));
    n34(n34mn == mi) = n34(n34mn == mi) - n34cli(mi);
end

%%% Get wind data:

Uname = [base 'U10_anom.nc'];
Vname = [base 'V10_anom.nc'];
if (source_type)
    U10 = ncread(Uname,'uas_10m');
    V10 = ncread(Vname,'vas_10m');
    lat = ncread(Uname,'latitude');
    lon = ncread(Uname,'longitude');
    time = ncread(Uname,'time');
    dnum = datenum([1900 1 1 0 0 0])+time;
else
    U10 = ncread(Uname,'10U_GDS0_SFC');
    V10 = ncread(Vname,'10V_GDS0_SFC');
    lat = ncread(Uname,'g0_lat_1');
    lon = ncread(Uname,'g0_lon_2');
    time = ncread(Uname,'initial_time0_hours');
    dnum = datenum([1800 1 1 0 0 0])+time/24;
end
[xL,yL,tL] = size(U10);

[X,Y] = ndgrid(lon,lat);
dvec = datevec(dnum);
yr = dvec(:,1);
mn = dvec(:,2);

%%% Calculate N34 regression on source data grid:
minyr = max([min(yr) min(n34yr)]);
maxyr = min([max(yr) max(n34yr)]);
minyr = 1982;

n34i = find(n34yr==minyr,1,'first');
n34f = find(n34yr==maxyr,1,'last');
ERAi = find(yr==minyr,1,'first');
ERAf = find(yr==maxyr,1,'last');

tL = ERAf-ERAi+1;

% Full year:
U10reg = reshape(reshape(U10(:,:,ERAi:ERAf),[xL*yL tL])*(n34(n34i: ...
                                                  n34f)-mean(n34(n34i:n34f)))/tL/std(n34(n34i:n34f)),[xL yL]);
V10reg = reshape(reshape(V10(:,:,ERAi:ERAf),[xL*yL tL])*(n34(n34i: ...
                                                  n34f)-mean(n34(n34i:n34f)))/tL/std(n34(n34i:n34f)),[xL yL]);
