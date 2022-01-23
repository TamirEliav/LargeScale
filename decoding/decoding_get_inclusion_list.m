function [exp_list, T] = decoding_get_inclusion_list()

%%
inc_list_filename = 'F:\sequences\inclusion\exp_inc_list.xlsx';
T = readtable(inc_list_filename, 'ReadRowNames',1,'ReadVariableNames',1);
exp_list = T.exp_ID(T.included_auto);

end