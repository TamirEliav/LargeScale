function exp_create_details(exp_ID)

%% load cells summary
exp_t = DS_get_exp_summary();
details = table2struct( exp_t(exp_ID,:) );

%% add some details...
details.exp_ID = regexprep(exp_ID, '_TT\d+_SS\d+','');
details.session_names = eval(details.session_names);
details.session_ts = eval(details.session_ts);
details.activeChannels = eval(details.activeChannels);
details.refCh = eval(details.refCh);
details.TT_to_use = eval(details.TT_to_use);

%% save exp details
filename_exp_details = ['L:\Analysis\Results\exp\details\' exp_ID '_exp_details' ];
save(filename_exp_details, 'details');

end