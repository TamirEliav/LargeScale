function nlx_ntt_write(filename_out, Timestamps, Samples, CellNumbers, fs)
    
header_file = 'Nlx_header_NTT.txt';
header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');

if isempty(CellNumbers)
    CellNumbers = zeros(size(Timestamps));
end
%     csc header
ADMaxValue = 32767;
InputRange = max(abs(Samples(:)));
ADC = InputRange / ADMaxValue / 1e6;
Samples = Samples ./ ADC ./ 1e6;
ADC_str = sprintf('%.24f',ADC);
InputRange_str = sprintf('%g',InputRange);
ADC_str_IX = contains(header, 'ADBitVolts');
InputRange_str_IX = contains(header, 'InputRange');
header{ADC_str_IX} = sprintf('-ADBitVolts %s %s %s %s', ADC_str, ADC_str, ADC_str, ADC_str);
header{InputRange_str_IX} = sprintf('-InputRange %s %s %s %s', InputRange_str, InputRange_str, InputRange_str, InputRange_str);
header{contains(header, '-SamplingFrequency')} = sprintf('-SamplingFrequency %g',fs);

Mat2NlxSpike(filename_out, 0, 1, [], [1 0 1 0 1 1], ...
    Timestamps, CellNumbers, Samples, header);

end
