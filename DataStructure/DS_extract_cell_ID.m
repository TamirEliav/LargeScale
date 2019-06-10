function cell_ID = DS_extract_cell_ID(str)


%%
% str = 'b0034_d180312_TT4_SS01_FR_map'; % example
cell_ID = regexp(str, '(b[\d]+_d[\d]+_TT[\d]+_SS[\d]+)*','tokens');
cell_ID = cell_ID{1}{1};

end