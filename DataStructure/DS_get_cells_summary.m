function cells_t = DS_get_cells_summary()

    %%
    cells_data_file = 'L:\Analysis\Code\inclusion_lists\SfN_2018_Shir_sorting';
%     opt = 'ReadFromExcel';
    opt = 'ReadFromMat';
    
    %%
    switch opt
        case 'ReadFromExcel'
            cells_t = readtable([cells_data_file '.xlsx'], 'ReadRowNames',1);
            save(cells_data_file, 'cells_t');
        case 'ReadFromMat'
            load(cells_data_file, 'cells_t');
    end
    
end


