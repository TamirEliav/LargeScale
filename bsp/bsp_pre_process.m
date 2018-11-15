function bsp_pre_process(main_dir)


%%
clear all
clc

%% read files
main_dir = 'L:\DATA\9892_Rocky\20160306\Bsp';
load( fullfile(main_dir,'bsp_pos.mat') );


%% pre-process
% (1) clear artifacts.
% (2) smooth
% (3) clac velocity
% (4) project to 1D


%% save mat files


end