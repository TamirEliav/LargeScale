function paramset = PARAMS_GetParamset()

%% global param set definition
global paramset_global_var;
if isempty(paramset_global_var)
    paramset_global_var = 0;
end
paramset = paramset_global_var;