function exp_sync_bsp2nlg(exp_ID)

exp=exp_load_data(exp_ID,'details','path');
PRE_sync_bsp_to_nlg(exp.path.bsp, exp.path.nlx, exp.path.sync, exp.details.sync_jump_ts);

end