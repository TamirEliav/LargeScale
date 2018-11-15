function exp_ID = DS_get_exp_ID_from_cell_ID(cell_ID)

exp_ID = regexprep(cell_ID, '_TT([\d]+)_SS_([\d]+)', '');

end