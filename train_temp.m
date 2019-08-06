cd RFECS_software
addpath(pwd)
cd ..
%matlabpool
load set_p300_tssbg
load vals
forest_training_part1_parallel(set_p300_tssbg,vals,15,65,20,'forest_p300_tssbg')
%matlabpool close
exit
