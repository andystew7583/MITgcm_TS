%%% 
%%% backupAverages.m
%%%
%%% Creates a .mat backup file of all of the time-averaged model fields,
%%% for convenient storage.
%%%

loadexp;
sdir = './backups';
sfname = fullfile(sdir,[expname,'_backup.mat']);
load(sfname);