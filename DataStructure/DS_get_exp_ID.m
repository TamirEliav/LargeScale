function exp_ID = DS_get_exp_ID(str)

exp_ID = regexp(str, '(b[\d]+_d[\d]+)', 'tokens');
exp_ID = exp_ID{1};
exp_ID = exp_ID{1};

%%
end