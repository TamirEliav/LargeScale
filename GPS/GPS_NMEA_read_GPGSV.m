function data = GPS_NMEA_read_GPGSV(filename)


%%
clear 
clc
% Filename = 'C:\Users\Administrator\Downloads\nmea-GeoTag (1).txt';
Filename = 'L:\GPS\WILD_BATS_YOSSI_DATA\nmea-GeoTag_1.nmea';
% Filename = 'L:\GPS\WILD_BATS_YOSSI_DATA\nmea-GeoTag_test.nmea';

%%
fid = fopen(Filename, 'rt');

%%
sat_numbers = [];
SNRs = [];
ii_entry = 1;
while ~feof(fid)
    
    %% read line
    nline = fgetl(fid);    
    %%
    if isempty(nline)
        continue;
    end
    
    %% Find which string we're dealing with
    fields = textscan(nline,'%s','delimiter',',');
    % pull the checksum out of the last field and make a new one for it
    fields{1}{end+1} = fields{1}{end}(end-1:end);
    % cut off the old last field at the chksum delimiter
    fields{1}(end-1) = strtok(fields{1}(end-1), '*');
    nmea_sentence = fields{1}{1};
%     fields = char(fields{1});
%     fields([1 end]) = [];
    fields = fields{1};
    fields([1 end]) = [];
    switch nmea_sentence 
        case '$GPGSV'
            fields([1:3]) = [];
            nsat = length(fields) / 4;
            for ii_sat = 1:nsat
                sat_numbers(ii_entry) = str2num(fields{4*(ii_sat-1)+1});
                SNRs(ii_entry) =     str2num(fields{4*(ii_sat-1)+4});
                ii_entry = ii_entry + 1;
            end
        otherwise
    end
end
warning on
fclose(fid);


%% change to a table
sat_numbers = sat_numbers';
SNRs = SNRs';
data = table(sat_numbers, SNRs);

%%
all_sat = unique(data.sat_numbers);
nsat = length(all_sat);
figure
hold on
for ii_sat = 1:nsat
    sat_ID = all_sat(ii_sat);
    IX = find(data.sat_numbers==sat_ID);
%     plot(IX, data.SNRs(IX), '.')
%     histogram( data.SNRs(IX) , linspace(30,50,20));
end


%%
% h=histogram( data.SNRs , linspace(30,50,20));
h=histogram( data.SNRs );
h.Normalization = 'probability';
xlabel('SNR (dB)')
ylabel('Prob.')

end



    










