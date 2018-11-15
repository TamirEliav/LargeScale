function nlg2nlx_audio(main_dir)
% convert audio logger file to nlx format

%% arrange files/folders
disp('1. constract/verify file structure...');
header_file = 'Nlg2Nlx_header.txt';
% if ~exist(main_dir,'dir')
%     mkdir(main_dir)
% end

%% parameters setting
DATA_file_prefix = 'NEUR0';
% DATA_file_prefix = 'BACK';
zero_DC_level_bit = 2048;
is_invert_data = true;
is_remove_DC = true;
use_clock_diff_correction = true;
num_channels = 16;

%% read EVENT file
event_file_name_xlsx = fullfile(main_dir, 'EVENTLOG.csv');
[NUM,TXT,RAW]=xlsread(event_file_name_xlsx);

% extract recording details from event file header
file_header_lines = TXT(1:3,1);
[splitstr] = regexp(file_header_lines{2}, '[;]+', 'split'); % 2nd header row
firmware_ver = regexp(splitstr{1}, '\d*[.]\d*','match');
serial_number = regexp(splitstr{2}, '\d*','match');
time = regexp(splitstr{3}, '\d*:\d*:\d*','match');
date = regexp(splitstr{4}, '\d*/\d*/\d*','match');
% ADC_period_usec = regexp(splitstr{5}, '\d*[.]\d*','match');
ADC_period_usec = regexp(splitstr{5}, '\d*[.]?\d*','match');
% ADC_period_usec = regexp(splitstr{5}, '\d*','match');
[splitstr] = regexp(file_header_lines{3}, '[;]+', 'split'); % 3rd header row
ADC_resolution_uVolt = regexp(splitstr{1}, '\d*[.]\d*','match');

ADC_SAMPLE_PERIOD = str2num(cell2mat(ADC_period_usec))/num_channels*1e-6;
fs =  1/ADC_SAMPLE_PERIOD;
% uVolt_per_bit = str2num(cell2mat(ADC_resolution_uVolt));
block_period_time_usec = (512/fs) * 1e6;
digital_data_block_period_time_usec = (512/fs) * 1e6;
block_period_time_usec_2 = 512 * (ADC_SAMPLE_PERIOD * num_channels) * 1e6;
file_len_time_usec = block_period_time_usec*1024; % TODO: calc
fs_acc = 1e6/(file_len_time_usec / 2048);
ts_offset_acc_neur_usec = file_len_time_usec / 2048;

%% extract events details
events_IX = NUM(1:end,1);
events_TS = NUM(1:end,3).*1e3;
events_TS_source = TXT(5:end,4);
events_type = TXT(5:end,5);
events_details = TXT(5:end,6);

%% join '...Continued' events
continued_event_lines_IX = find(isnan(events_IX))';
for ii_continued_event_line = continued_event_lines_IX
    % take the last line with a valid number in the event index column as
    % the event index
    last_valid_line = find( ~isnan(events_IX(1:ii_continued_event_line)), 1, 'last');
    events_details{last_valid_line} = [events_details{last_valid_line} events_details{ii_continued_event_line}];
end
events_IX(continued_event_lines_IX) = [];
events_TS(continued_event_lines_IX) = [];
events_TS_source(continued_event_lines_IX) = [];
events_type(continued_event_lines_IX) = [];
events_details(continued_event_lines_IX) = [];

%% join event type with event details to a single string
EventStrings = {};
for ii_event = 1:length(events_type)
    event_str = [events_type{ii_event} ' - ' events_details{ii_event}];
    EventStrings{ii_event} = event_str;
end

%% extract recording parameteres from the event file

% % % % % disp('2. reading+saving parameters from event file...');
% % % % % 
% % % % % for i=1:size(events_type,1)
% % % % %     if strcmp(events_type{i},'Recording parameters'),break;end
% % % % % end
% % % % % 
% % % % % paramst = strsplit(events_details{i},';');
% % % % % for i=1:length(paramst)
% % % % %     if length(paramst{i})>1
% % % % %         t = strsplit(paramst{i},' = ');
% % % % %         t1 = t{1};
% % % % %         t1= t1(find(~isspace(t1)));
% % % % %         if strfind(t1,'-'), t1(find(t1=='-'))='';end
% % % % %         if strfind(t1,'/'), t1(find(t1=='/'))='';end
% % % % %         t{1} = t1;
% % % % %         param.(t1) = t{2};
% % % % %     end
% % % % % end
% % % % % 
% % % % % % param.ChannelOverwrittenbyMotionSensor = str2double(char(regexp(param.ChannelOverwrittenbyMotionSensor,'(\d*)','match')));
% % % % % % param.GyroscopeRangeIndex = str2double(char(regexp(param.GyroscopeRangeIndex,'(\d*)','match')));
% % % % % % param.AccelerometerRangeIndex = str2double(char(regexp(param.AccelerometerRangeIndex,'(\d*)','match')));
% % % % % param.LowCutoffFrequency = str2double(char(regexp(param.LowCutoffFrequency,'(\d*)','match')));
% % % % % 
% % % % % save( fullfile(Nlx_OutDir,'param.mat'), 'param' );


%% Apply clock difference correction (Transceiver vs. logger time)
% TODO: plot figure showing the clocks drift

% CD = clock difference
disp('3. correcting for clock diff...');

if use_clock_diff_correction
    
    PC_gen_events_IX = find(strcmp('PC-generated comment', events_type));
    CD_str = 'CD=';
    CD_values = [];
    CD_timestamps = [];
    for ii_event = 1:length(PC_gen_events_IX)
        curr_event_details = events_details{PC_gen_events_IX(ii_event)};
        CD_str_pos = strfind(curr_event_details, CD_str);
        if isempty(CD_str_pos)
            continue;
        end
        str = curr_event_details(CD_str_pos+3:end);
        CD = str2num(str(1:min(strfind(str,' '))));
        CD_values(end+1) = CD;
        CD_timestamps(end+1) = events_TS(PC_gen_events_IX(ii_event));
    end
    CD_values = CD_values.*1e6; % change from sec to usec
    
    % TODO: we need to check what is the meaning of CD (clock difference). Logger-Tx
    % or Tx-Logger?
    % Is this old comment relevant? Didi
    % TODO: best to do is to check this again after Jacob bring us the new
    % version that fix the problem from version 1.69
    
    CD_event_ts__logger_time = CD_timestamps-CD_timestamps(1);        % removed 1st timestamp to balance the fit data closer to zero (otheriwse it gives a warning message)
    CD_event_ts__Tx_time = CD_event_ts__logger_time - CD_values;
    transceiver_2_logger_time_fit = polyfit(CD_event_ts__Tx_time, CD_event_ts__logger_time , 1);
    ts_source_transceiver_events_IX = find(strcmp('Transceiver', events_TS_source)); % TODO: add Transceiver (fine) - validate string !!
    events_TS(ts_source_transceiver_events_IX) = polyval(transceiver_2_logger_time_fit , events_TS(ts_source_transceiver_events_IX) - CD_timestamps(1)) + CD_timestamps(1);
    for ii_event = 1:length(ts_source_transceiver_events_IX)
        events_TS_source{ts_source_transceiver_events_IX(ii_event)} = 'Logger';
    end
    
    % TODO: plot some figure, and make a check that the clock dcorrection
    % was OK
end

%% write event file in Nlx format

disp('4. writing event file...');
NlxEventFile = fullfile(main_dir, 'EVENTS.nev');
Mat2NlxEV(NlxEventFile, 1, 1, [], [1 0 0 0 1 0] , events_TS', EventStrings');

%% create sperate event file for each event category
event_type_list = unique(events_type);
for ii_event_type = 1:length( event_type_list )
    event_type_string = event_type_list{ii_event_type};
    event_type_events_IX = find(strcmp( event_type_string , events_type));
    Nlx_event_type_file = fullfile(main_dir, ['EVENTS__' event_type_string '.nev']);
    Mat2NlxEV(Nlx_event_type_file , 1, 1, [], [1 0 0 0 1 0] , events_TS(event_type_events_IX)', EventStrings(event_type_events_IX)');
end

%% create battery discharge plot
disp('5. calculating battery discharge...');

PC_gen_events_IX = find(strcmp('PC-generated comment', events_type));
mode_changed_events_IX = find(strcmp('Mode change', events_type));
BV_str = 'BV=';
BV_values = [];
BV_timestamps = [];
for ii_event = 1:length(PC_gen_events_IX)
    curr_event_details = events_details{PC_gen_events_IX(ii_event)};
    BV_str_pos = strfind(curr_event_details, BV_str);
    if isempty(BV_str_pos)
        continue;
    end
    
    BV_values_str_interval = (BV_str_pos:BV_str_pos+4) + length(BV_str);
    BV_values(end+1) = str2num( curr_event_details(BV_values_str_interval) );
    BV_timestamps(end+1) = events_TS(PC_gen_events_IX(ii_event));
end
usec_2_min = 1/60*1e-6;
% BV_timestamps = BV_timestamps(4:end);
% BV_values = BV_values(4:end)
figure
plot( usec_2_min.*(BV_timestamps - BV_timestamps(1)), BV_values , '-o')
title('Battery discharge')
xlabel('time (minutes)')
ylabel('Voltage (V)')
grid on

hold on
mode_changed_events_ts = events_TS(mode_changed_events_IX) - BV_timestamps(1);
ylimits = get(gca, 'ylim');
plot( repmat(usec_2_min .* mode_changed_events_ts,1,2)' , ylimits', '--m' )

legend({'battey voltage','record mode change event'})

fig_file = fullfile(main_dir, 'Battery_discharge');
saveas(gcf, fig_file , 'jpg')
saveas(gcf, fig_file , 'fig')


%% create Nlx data files - continuous sampling files called (.ncs)

% identify  'File started' events and take timestamps
FileStarted_IX = strcmp('File started', events_type);
FileStarted_TS = events_TS(FileStarted_IX);
FileStarted_details = events_details(FileStarted_IX);

% read header template
header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
% change fields that are recording specific
header_audio = strrep(header,'-SamplingFrequency', ['-SamplingFrequency ' char(vpa(fs))]);
disp('6. Nlg -> Nlx...');

for ii_file_start_entry = 1:length(FileStarted_TS)
    
    %%
    file_str = FileStarted_details{ii_file_start_entry};
%     file_name = [DATA_file_prefix '0' file_str(end-2:end) '.DAT'];
%     file_name = [DATA_file_prefix '_' file_str(end-2:end) '.DAT'];
    file_num_str = regexp(file_str, '([\d]+)','tokens');
    file_num_str = file_num_str{1}; file_num_str = file_num_str{1};
    file_name = [DATA_file_prefix file_num_str '.DAT'];
    fid = fopen(fullfile(main_dir, file_name));
    filedata = fread(fid, 'int16', 0, 'l');
    fclose(fid);
    filedata = double(filedata);
    data = reshape(filedata, 1, length(filedata)/1);
    file_TS_usec = FileStarted_TS(ii_file_start_entry);
    
    %% write data chunks to Nlx file format
    %     file_TS_usec = FileStarted_TS(ii_file_start_entry);
    disp(['writing file ' file_str])
    cnl_data = data;
    cnl_data_blocks = vec2mat(cnl_data,512)';
    num_blocks = size(cnl_data_blocks, 2);
    blocks_timestamps_usec = file_TS_usec + (0:block_period_time_usec:(block_period_time_usec*(num_blocks-1)));
    file_name = ['audio.ncs'];
    out_file = fullfile(main_dir, file_name);

    %         we do this because Nlx functions have bugs with writing header
    %         when using append... (this way works...)
    if exist(out_file, 'file')
        append_flag = 1;
    else
        append_flag = 0;
    end
    %append_flag

    Mat2NlxCSC(out_file, append_flag, 1, 0, [1 0 0 0 1 1], blocks_timestamps_usec, cnl_data_blocks, header_audio );
end

%%
disp('Finished converting all files from Nlg format to Nlx format!')









%%