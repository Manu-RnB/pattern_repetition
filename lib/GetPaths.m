function GetPaths(Computer_username)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% [outputArg1,outputArg2]
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;


clear all; clc

computer_username='emmanuelcoulon'; %Mac 13"
% computer_username='manu'; %Mac 15"

cd (['/Users/',computer_username,'/Documents/MATLAB/PitchChange'])
 
letswave_path = ['/Users/',computer_username,'/Documents/MATLAB/Toolboxes/letswave7-master'];
if ~exist(letswave_path); warning('provide valid path to letswave folder'); end
addpath(genpath(letswave_path))

data_path = ['/Users/',computer_username,'/Documents/MATLAB/PitchChange']; % Where the data is stored (sub-001, sub002, etc folders)
if ~exist(data_path); warning('provide valid path to data folder'); end
addpath(genpath(data_path))

figure_path = ['/Users/',computer_username,'/Documents/MATLAB/PitchChange/figures']; % Where I want the figures to be saved
 
if ~exist('PitchChange_Data.mat'); Data=[]; 
else; load('PitchChange_Data.mat')
end

end

