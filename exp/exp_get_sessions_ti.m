function sessions_ti = exp_get_sessions_ti(exp_ID, varargin)

exp = exp_load_data(exp_ID,'details');
for ii_arg = 1:length(varargin)
    IX = find(strcmp(exp.details.session_names, varargin{ii_arg}));
    if isempty(IX)
        sessions_ti(ii_arg,:) = [nan nan];
    else
        sessions_ti(ii_arg,:) = exp.details.session_ts(IX,:);
    end
    
end

end