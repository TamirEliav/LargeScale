function decode = decoding_load_data(exp_ID, epoch_type, params_opt)
arguments
    exp_ID = 'b9861_d180526'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
    params_opt = 11
end

dir_IN = 'F:\sequences\decoded';
decode_filename = fullfile(dir_IN, epoch_type, exp_ID, sprintf('%s_%s_opt_%d.nc',exp_ID,epoch_type,params_opt));
decode = decoding_read_decoded_file(decode_filename);
decode.exp_ID = exp_ID;
decode.epoch_type = epoch_type;
decode.params_opt = params_opt;

end