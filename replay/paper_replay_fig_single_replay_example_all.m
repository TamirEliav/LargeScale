%%
clear
clc

%%
examples_list = [
    % sleep
    struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',36),
	struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',109),
    struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',138),
	struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',139),
    struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',143),
    struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',152),
	struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',178),
    struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',226),
    struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',279),
    struct('exp_ID','b0184_d191201', 'epoch_type','sleep', 'event_num',75),
	struct('exp_ID','b0184_d191201', 'epoch_type','sleep', 'event_num',100),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',10),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',26),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',59),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',60),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',79), % merged seqs
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',124),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',137),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',138),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',154),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',159),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',184),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',187),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',204),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',206),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',209),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',217),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',231),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',234),
	struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',241),
	struct('exp_ID','b0184_d191203', 'epoch_type','sleep', 'event_num',31),
	struct('exp_ID','b0184_d191203', 'epoch_type','sleep', 'event_num',34),
	struct('exp_ID','b0184_d191203', 'epoch_type','sleep', 'event_num',44),
	struct('exp_ID','b0184_d191203', 'epoch_type','sleep', 'event_num',60),
	struct('exp_ID','b0184_d191203', 'epoch_type','sleep', 'event_num',66),
	struct('exp_ID','b0184_d191203', 'epoch_type','sleep', 'event_num',89),
	struct('exp_ID','b0184_d191203', 'epoch_type','sleep', 'event_num',90),
	struct('exp_ID','b0184_d191203', 'epoch_type','sleep', 'event_num',93),
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',14),
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',17),
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',32),
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',35),
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',45),
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',84), % merged seqs
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',96),
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',98),
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',104),
	struct('exp_ID','b0184_d191204', 'epoch_type','sleep', 'event_num',120),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',4),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',63),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',64),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',73),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',78),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',83),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',98),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',125), % merged seqs
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',174),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',195),
	struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',227),
	struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',36),
	struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',42),
	struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',47),
    struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',51),
    struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',64),
	struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',70),
    struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',73),
	struct('exp_ID','b0184_d191209', 'epoch_type','sleep', 'event_num',24),
	struct('exp_ID','b0184_d191209', 'epoch_type','sleep', 'event_num',63),
	struct('exp_ID','b0184_d191209', 'epoch_type','sleep', 'event_num',66),
	struct('exp_ID','b0184_d191209', 'epoch_type','sleep', 'event_num',78), % merged seqs
	struct('exp_ID','b0184_d191210', 'epoch_type','sleep', 'event_num',42),
	struct('exp_ID','b0184_d191210', 'epoch_type','sleep', 'event_num',70),
	struct('exp_ID','b0184_d191210', 'epoch_type','sleep', 'event_num',76),
	struct('exp_ID','b0184_d191210', 'epoch_type','sleep', 'event_num',80),
	struct('exp_ID','b0184_d191211', 'epoch_type','sleep', 'event_num',36),
	struct('exp_ID','b0184_d191211', 'epoch_type','sleep', 'event_num',92),
	struct('exp_ID','b0184_d191212', 'epoch_type','sleep', 'event_num',22),
	struct('exp_ID','b0184_d191212', 'epoch_type','sleep', 'event_num',27),
	struct('exp_ID','b0184_d191212', 'epoch_type','sleep', 'event_num',56),
	struct('exp_ID','b0184_d191216', 'epoch_type','sleep', 'event_num',1),
	struct('exp_ID','b0184_d191220', 'epoch_type','sleep', 'event_num',13),
	struct('exp_ID','b0184_d191220', 'epoch_type','sleep', 'event_num',14),
	struct('exp_ID','b0184_d191220', 'epoch_type','sleep', 'event_num',18),
	
	struct('exp_ID','b2382_d190623', 'epoch_type','sleep', 'event_num',36),
	struct('exp_ID','b2382_d190623', 'epoch_type','sleep', 'event_num',42),
	struct('exp_ID','b2382_d190623', 'epoch_type','sleep', 'event_num',77),
	struct('exp_ID','b2382_d190627', 'epoch_type','sleep', 'event_num',43),
	struct('exp_ID','b2382_d190627', 'epoch_type','sleep', 'event_num',50),
	struct('exp_ID','b2382_d190712', 'epoch_type','sleep', 'event_num',24),
	struct('exp_ID','b2382_d190712', 'epoch_type','sleep', 'event_num',57),
	struct('exp_ID','b2382_d190716', 'epoch_type','sleep', 'event_num',126),
	struct('exp_ID','b2382_d190716', 'epoch_type','sleep', 'event_num',138),
	struct('exp_ID','b2382_d190716', 'epoch_type','sleep', 'event_num',139),
	struct('exp_ID','b2382_d190718', 'epoch_type','sleep', 'event_num',69),
	struct('exp_ID','b2382_d190718', 'epoch_type','sleep', 'event_num',85),
	struct('exp_ID','b2382_d190718', 'epoch_type','sleep', 'event_num',87),
	struct('exp_ID','b2382_d190728', 'epoch_type','sleep', 'event_num',472),
	struct('exp_ID','b2382_d190729', 'epoch_type','sleep', 'event_num',67),
	struct('exp_ID','b2382_d190729', 'epoch_type','sleep', 'event_num',76),
	struct('exp_ID','b2382_d190729', 'epoch_type','sleep', 'event_num',83),
	struct('exp_ID','b2382_d190729', 'epoch_type','sleep', 'event_num',84),
	struct('exp_ID','b2382_d190730', 'epoch_type','sleep', 'event_num',211),
	struct('exp_ID','b2382_d190730', 'epoch_type','sleep', 'event_num',213),
	struct('exp_ID','b2382_d190730', 'epoch_type','sleep', 'event_num',216),
	struct('exp_ID','b2382_d190801', 'epoch_type','sleep', 'event_num',71),
	struct('exp_ID','b2382_d190804', 'epoch_type','sleep', 'event_num',27),
	struct('exp_ID','b2382_d190804', 'epoch_type','sleep', 'event_num',69),
	struct('exp_ID','b2382_d190807', 'epoch_type','sleep', 'event_num',28),
	struct('exp_ID','b2382_d190808', 'epoch_type','sleep', 'event_num',61),
	struct('exp_ID','b2382_d190812', 'epoch_type','sleep', 'event_num',60),
	struct('exp_ID','b2382_d190813', 'epoch_type','sleep', 'event_num',33),
	struct('exp_ID','b2382_d190813', 'epoch_type','sleep', 'event_num',51),
	struct('exp_ID','b2382_d190813', 'epoch_type','sleep', 'event_num',119),
	struct('exp_ID','b2382_d190814', 'epoch_type','sleep', 'event_num',44),
	struct('exp_ID','b2382_d190814', 'epoch_type','sleep', 'event_num',52),
	
	
	
	% need to go over 9861 34 148 194 2289
    struct('exp_ID','b9861_d180526', 'epoch_type','sleep', 'event_num',49),
    struct('exp_ID','b0194_d180614', 'epoch_type','sleep', 'event_num',2),
    struct('exp_ID','b0194_d180614', 'epoch_type','sleep', 'event_num',3),
    struct('exp_ID','b0194_d180606', 'epoch_type','sleep', 'event_num',1),
    struct('exp_ID','b0194_d180520', 'epoch_type','sleep', 'event_num',8),
	struct('exp_ID','b0034_d180313', 'epoch_type','sleep', 'event_num',21),
    struct('exp_ID','b0034_d180313', 'epoch_type','sleep', 'event_num',39),
    struct('exp_ID','b0148_d170613', 'epoch_type','sleep', 'event_num',8),
    struct('exp_ID','b0148_d170614', 'epoch_type','sleep', 'event_num',5),
    struct('exp_ID','b0148_d170625', 'epoch_type','sleep', 'event_num',11),
    struct('exp_ID','b0148_d170625', 'epoch_type','sleep', 'event_num',18),
    struct('exp_ID','b0148_d170625', 'epoch_type','sleep', 'event_num',25),
    struct('exp_ID','b0148_d170628', 'epoch_type','sleep', 'event_num',11),
    struct('exp_ID','b0148_d170703', 'epoch_type','sleep', 'event_num',3),
    struct('exp_ID','b0148_d170718', 'epoch_type','sleep', 'event_num',2),
    
    

];

%% params
win_s = 0.5;
params_opt = 11;

%%
for ii_example = 1:length(examples_list)
    example = examples_list(ii_example);
    paper_replay_fig_single_replay_example(example.exp_ID,example.epoch_type,params_opt,example.event_num,win_s);
    close all
end
