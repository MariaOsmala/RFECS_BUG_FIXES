cd RFECS_software
addpath(pwd)
cd ..
%set(0,'RecursionLimit', 5000);
%matlabpool close force local
%matlabpool

%matlabpool open 12 %THIS WAS UNCOMMENTED
extn ='_enhancers';
threshold=0.5 %change this
len_var=20;
chr_set={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX'};
load forest_p300_tssbg

%forest_predict(forest_p300_tssbg_all,		...
%		nvar,			...
%		ntrees,			...
%		[1:nvar],		...
%		extn,			...
%		prob_bg_p300_dist,			...
%		20,		...
%		'test_set/all_mods_',		...
%		'test_set/',		...
%		chr_set)		
%
%forest_predict(forest_all,		forest_p300_tssbg_all
%		nvar,			nvar
%		ntrees,			ntrees
%		mod_vec,		[1:nvar]
%		extn,			extn
%		prob_bg_p300_dist,	?????
%		peak_filt_dist,		20
%		input_path,
%		output_path,
%		chr_set)


forest_predict_par_threshold(forest_p300_tssbg_all, ... 
                   nvar,                  ... 
                   ntrees,                ... 
                   len_var,               ... 
                   [1:nvar],              ... 
                   extn,                  ... 
                   20,                    ... 
                   'test_set/all_mods_',  ... 
                   'test_set/',           ... 
                   chr_set, threshold)			



%forest_predict_par(forest_all,
%			nvar,
%			ntrees,
%			len_var,
%			mod_vec,		[1:nvar]
%			extn,
%			peak_filt_dist,
%			input_path,
%			output_path,
%			chr_set)
%


%matlabpool close %THIS WAS UNCOMMENTED
exit
