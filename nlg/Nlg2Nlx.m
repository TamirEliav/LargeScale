function Nlg2Nlx(main_dir,forcecalc)

% function Nlg2Nlx(main_dir)
%
%   main_dir: data dir
%
%
% Tamir Eliav
% revised - Didi Omer, May  2015
% revised - Didi Omer, June 2015
%

%% defaults
if nargin==1
    forcecalc = 0;
end

%% arrange files/folders
disp('1. constract/verify file structure...');
Nlg_InDir = fullfile(main_dir, 'nlg');
Nlx_OutDir = fullfile(main_dir, 'nlx');
if ~exist(Nlg_InDir,'dir')
    error('NLG input folder does not exist');
end
if exist(Nlx_OutDir,'dir')
    if forcecalc
%         % rename existing output dir
%         Nlx_OutDir_OLD_rename = sprintf([Nlx_OutDir '_OLD_%s'],datetime(datetime,'Format','yyyymmdd_HHmmss'));
%         movefile(Nlx_OutDir, Nlx_OutDir_OLD_rename);
        % delete existing output dir
        warning('NLX output dir already existing and you chose to override it, deleting old NLX dir!')
        rmdir(Nlx_OutDir,'s')
    else
        error('NLX output folder already exist, use forcecalc to override it!');
    end
end
% at this point we should not have the nlx output dir, so let's create it!
mkdir(Nlx_OutDir)

%% parameters setting
num_channels = 16;
data_cnl_ind = [0:15]; % note that this numbering system is of the neurologger which means channels 0-15
DATA_file_prefix = 'NEUR0';
% DATA_file_prefix = 'BACK';
zero_DC_level_bit = 2048;
is_invert_data = true;
is_remove_DC = true;
is_remove_flash_write_artifact = false;
use_clock_diff_correction = true;
use_post_rec_ref_channel = false; % use this if you recorded with GND as ref channel
post_rec_ref_channel = 1; % (1-16) If the above is true - choose the channel you want to substract from all the other channels
motion_channel = 1;
header_file = 'Nlg2Nlx_header.txt';

%% read EVENT file
event_file_name_xlsx = fullfile(Nlg_InDir, 'EVENTLOG.csv');
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
fs =  1/(ADC_SAMPLE_PERIOD * num_channels);
uVolt_per_bit = str2num(cell2mat(ADC_resolution_uVolt));
block_period_time_usec = (512/fs) * 1e6;
digital_data_block_period_time_usec = (512/fs) * 1e6;
block_period_time_usec_2 = 512 * (ADC_SAMPLE_PERIOD * num_channels) * 1e6;
file_len_time_usec = block_period_time_usec*1024;
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

disp('2. reading+saving parameters from event file...');

is_recording = false;
for i=1:size(events_type,1)
    if strcmp(events_type{i},'Recording parameters')
        is_recording = true;
        paramst = strsplit(events_details{i},';');
        for i=1:length(paramst)
            if length(paramst{i})>1
                t = strsplit(paramst{i},' = ');
                t1 = t{1};
                t1= t1(find(~isspace(t1)));
                if strfind(t1,'-'), t1(find(t1=='-'))='';end
                if strfind(t1,'/'), t1(find(t1=='/'))='';end
                t{1} = t1;
                param.(t1) = t{2};
            end
        end

        % param.ChannelOverwrittenbyMotionSensor = str2double(char(regexp(param.ChannelOverwrittenbyMotionSensor,'(\d*)','match')));
        % param.GyroscopeRangeIndex = str2double(char(regexp(param.GyroscopeRangeIndex,'(\d*)','match')));
        % param.AccelerometerRangeIndex = str2double(char(regexp(param.AccelerometerRangeIndex,'(\d*)','match')));
        param.LowCutoffFrequency = str2double(char(regexp(param.LowCutoffFrequency,'(\d*)','match')));
        save( fullfile(Nlx_OutDir,'param.mat'), 'param' );
        break
    end
end


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
    
    CD_event_ts__logger_time = CD_timestamps;
    CD_event_ts__Tx_time = CD_event_ts__logger_time - CD_values;
    tx2logger_time_lm = fitlm(CD_event_ts__Tx_time, CD_event_ts__logger_time , 'linear');
    ts_src_tx_events_IX = find(contains(events_TS_source, 'Transceiver'));
    events_TS(ts_src_tx_events_IX) = predict(tx2logger_time_lm, events_TS(ts_src_tx_events_IX));
    for ii_event = 1:length(ts_src_tx_events_IX)
        events_TS_source{ts_src_tx_events_IX(ii_event)} = 'Logger';
    end
    
    % plot figure to make a check that the clock correction was OK
    figure
    subplot(1,2,1)
    plot(CD_event_ts__logger_time, CD_values.*1e-3, 'o-');
    xlabel('Logger time (usec)')
    ylabel('Clock Diff (ms)')
    title('logger vs. tx clock drift')
    subplot(1,2,2)
    plot(tx2logger_time_lm)
    text(0.1,0.8, sprintf('clock gain (logger vs. tx):\n%.24f',tx2logger_time_lm.Coefficients.Estimate(2)), 'Units','normalized','FontSize',12)
    axis equal
    xlabel('Logger time (usec)')
    ylabel('Transciever time (usec)')
    suptitle(main_dir)
    fig_file = fullfile(Nlx_OutDir, 'Clock_diff_correction');
    saveas(gcf, fig_file , 'tif')
end

%% write event file in Nlx format

disp('4. writing event file...');
NlxEventFile = fullfile(Nlx_OutDir, 'EVENTS.nev');
Mat2NlxEV(NlxEventFile, 1, 1, [], [1 0 0 0 1 0] , events_TS', EventStrings');

%% create sperate event file for each event category
event_type_list = unique(events_type);
for ii_event_type = 1:length( event_type_list )
    event_type_string = event_type_list{ii_event_type};
    event_type_events_IX = find(strcmp( event_type_string , events_type));
    Nlx_event_type_file = fullfile(Nlx_OutDir, ['EVENTS__' event_type_string '.nev']);
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
if is_recording
    mode_changed_events_ts = events_TS(mode_changed_events_IX) - BV_timestamps(1);
    ylimits = get(gca, 'ylim');
    plot( repmat(usec_2_min .* mode_changed_events_ts,1,2)' , ylimits', '--m' )
end

legend({'battey voltage','record mode change event'})
suptitle(main_dir)

fig_file = fullfile(Nlx_OutDir, 'Battery_discharge');
saveas(gcf, fig_file , 'tif')
% saveas(gcf, fig_file , 'fig')


%% create Nlx data files - continuous sampling files called (.ncs)

% identify  'File started' events and take timestamps
FileStarted_IX = strcmp('File started', events_type);
FileStarted_TS = events_TS(FileStarted_IX);
FileStarted_details = events_details(FileStarted_IX);

% read header template
header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
header_neural = header;
header_motion = header;
% change fields that are recording specific
ADMaxValue = 32767;
InputRange = 7000;
ADBitVolts = InputRange / ADMaxValue / 1e6;
if is_invert_data
    InputInvertedStr = 'True';
else
    InputInvertedStr = 'False';
end
header_neural{contains(header_neural,'SamplingFrequency')}  = sprintf('-SamplingFrequency %g', fs);
header_neural{contains(header_neural,'ADMaxValue')}         = sprintf('-ADMaxValue %d', ADMaxValue);
header_neural{contains(header_neural,'InputRange')}         = sprintf('-InputRange %g', InputRange);
header_neural{contains(header_neural,'ADBitVolts')}         = sprintf('-ADBitVolts %.24f', ADBitVolts);
header_neural{contains(header_neural,'InputInverted')}      = sprintf('-InputInverted %s', InputInvertedStr);
header_neural{contains(header_neural,'SamplingFrequency')}  = sprintf('-SamplingFrequency %g', fs);
header_neural{contains(header_neural,'SamplingFrequency')}  = sprintf('-SamplingFrequency %g', fs);
header_motion{contains(header_motion,'SamplingFrequency')}  = sprintf('-SamplingFrequency %g', fs_acc);

disp('6. Nlg -> Nlx...');
for ii_file_start_entry = 1:length(FileStarted_TS)
    
    %%
    file_str = FileStarted_details{ii_file_start_entry};
%     file_name = [DATA_file_prefix '0' file_str(end-2:end) '.DAT'];
%     file_name = [DATA_file_prefix '_' file_str(end-2:end) '.DAT'];
    file_num_str = regexp(file_str, '([\d]+)','tokens');
    file_num_str = file_num_str{1}; file_num_str = file_num_str{1};
    file_name = [DATA_file_prefix file_num_str '.DAT'];
    fid = fopen(fullfile(Nlg_InDir, file_name));
    filedata = fread(fid, 'uint16', 0, 'l');
    fclose(fid);
    filedata = double(filedata);
    data = reshape(filedata, num_channels, length(filedata)/num_channels);
    file_TS_usec = FileStarted_TS(ii_file_start_entry);
    
    
    if isfield(param, 'MotionSensorLogging');
    if strcmp(param.MotionSensorLogging,'Enabled')
        
        
        %% parse accelerometer data (digital)
        block_data = data(str2num(param.ChannelOverwrittenbyMotionSensor)+1,:);
        block_data = reshape(block_data,256,2048);
        acc_data = [];
        gyro_data = [];
        magnet_data = [];
        digital_data_ts = [];
        for ii_rec = 1:2048
            %% check barker
            rec_data = block_data(:,ii_rec);
            barker1 = rec_data(1);
            barker2 = rec_data(2);
            if (barker1~=13579) || (barker2~=24680)
                acc_data    = [acc_data,    [0 0 0]'];
                gyro_data   = [gyro_data,   [0 0 0]'];
                magnet_data = [magnet_data, [0 0 0]'];
                warning('wrong BARKER')
                continue;
            end
            acc_offset = rec_data(3);
            gyro_offset = rec_data(4);
            magnet_offset = rec_data(5);
            acc_len = rec_data(7);
            gyro_len = rec_data(8);
            magnet_len = rec_data(9);
            
            % Didi - for each of the 2048 sequences take the mean for each sensor.
            acc_data_seg = rec_data( (acc_offset+1) : (acc_offset+acc_len));
            gyro_data_seg = rec_data( (gyro_offset+1) : (gyro_offset+gyro_len));
            magnet_data_seg = rec_data( (magnet_offset+1) : (magnet_offset+magnet_len));
            if ~isempty(acc_data_seg)
                acc_data_seg = acc_data_seg(1:floor(length(acc_data_seg)/3)*3);
                gyro_data_seg = gyro_data_seg(1:floor(length(gyro_data_seg)/3)*3);
                magnet_data_seg = magnet_data_seg(1:floor(length(magnet_data_seg)/3)*3);
                
                acc_data_seg = median(reshape(acc_data_seg,[3,length(acc_data_seg)/3]),2);
                gyro_data_seg = median(reshape(gyro_data_seg,[3,length(gyro_data_seg)/3]),2);
                magnet_data_seg = median(reshape(magnet_data_seg,[3,length(magnet_data_seg)/3]),2);
                
                acc_data    = [acc_data,    acc_data_seg];
                gyro_data   = [gyro_data,   gyro_data_seg];
                magnet_data = [magnet_data, magnet_data_seg];
                
            else
                % Didi - case the motion data is empty (first sample, first
                % trial
                acc_data    = [acc_data,    [0;0;0]];
                gyro_data   = [gyro_data,   [0;0;0]];
                magnet_data = [magnet_data, [0;0;0]];
                
            end
            
        end
        %     digital_data_ts = digital_data_ts';
        %     digital_data_ts_usec = digital_data_ts .* 1e3;
        acc_data_X = acc_data(1,:);
        acc_data_Y = acc_data(2,:);
        acc_data_Z = acc_data(3,:);
        gyro_data_X = gyro_data(1,:);
        gyro_data_Y = gyro_data(2,:);
        gyro_data_Z = gyro_data(3,:);
        magnet_data_X = magnet_data(1,:);
        magnet_data_Y = magnet_data(2,:);
        magnet_data_Z = magnet_data(3,:);
        
        % convert to signed values
        uint16_to_int16 = @(a)(a.*(a<=(2^15-1)) + (a-2^16).*(a>(2^15-1)));
        acc_data_X = uint16_to_int16(acc_data_X);
        acc_data_Y = uint16_to_int16(acc_data_Y);
        acc_data_Z = uint16_to_int16(acc_data_Z);
        gyro_data_X = uint16_to_int16(gyro_data_X);
        gyro_data_Y = uint16_to_int16(gyro_data_Y);
        gyro_data_Z = uint16_to_int16(gyro_data_Z);
        magnet_data_X = uint16_to_int16(magnet_data_X);
        magnet_data_Y = uint16_to_int16(magnet_data_Y);
        magnet_data_Z = uint16_to_int16(magnet_data_Z);
        
        acc_data_X_blocks = reshape(acc_data_X,[],4)';
        acc_data_Y_blocks = reshape(acc_data_Y,[],4)';
        acc_data_Z_blocks = reshape(acc_data_Z,[],4)';
        gyro_data_X_blocks = reshape(gyro_data_X,[],4)';
        gyro_data_Y_blocks = reshape(gyro_data_Y,[],4)';
        gyro_data_Z_blocks = reshape(gyro_data_Z,[],4)';
        magnet_data_X_blocks = reshape(magnet_data_X,[],4)';
        magnet_data_Y_blocks = reshape(magnet_data_Y,[],4)';
        magnet_data_Z_blocks = reshape(magnet_data_Z,[],4)';
        
        temp = linspace(0,file_len_time_usec,2049);
        temp = temp(1:end-1);
        digital_data_ts_usec = file_TS_usec + temp;
        digital_data_ts_usec = digital_data_ts_usec - ts_offset_acc_neur_usec;
        digital_data_block_ts_usec = digital_data_ts_usec(1:512:end);
        digital_data_block_ts_usec =digital_data_block_ts_usec+(8.7381333/1000);
        digital_data_block_ts_usec = round(digital_data_block_ts_usec);          % usec is good enough... (and we can have double/integer problems...)
        
        acc_parsed_data = cat(3,...
            acc_data_X_blocks,...
            acc_data_Y_blocks,...
            acc_data_Z_blocks,...
            gyro_data_X_blocks,...
            gyro_data_Y_blocks,...
            gyro_data_Z_blocks,...
            magnet_data_X_blocks,...
            magnet_data_Y_blocks,...
            magnet_data_Z_blocks...
            );
        
        acc_ch_names = {...
            'acc_X',...
            'acc_Y',...
            'acc_Z',...
            'gyro_X',...
            'gyro_Y',...
            'gyro_Z',...
            'magnet_X',...
            'magnet_Y',...
            'magnet_Z',...
            };
        
        % write motion data to nlx file format
        for ii_acc_chnl = 1:9
            file_name = ['Motion_' acc_ch_names{ii_acc_chnl} '.ncs'];
            out_file = fullfile(Nlx_OutDir, file_name);
            
            % we do this because Nlx functions have bugs with writing header
            % when using append... (this way works...)
            if exist(out_file, 'file')
                append_flag = 1;
            else
                append_flag = 0;
            end
            
            blocks_to_write = squeeze(acc_parsed_data(:,:,ii_acc_chnl))';
            Mat2NlxCSC(out_file, append_flag, 1, 0, [1 0 1 0 1 1], digital_data_block_ts_usec, repmat(fs_acc,1,4), blocks_to_write , header_motion );
            
%             for ii_block = 2:2
%                 data_chunk = squeeze(acc_parsed_data(ii_block,:,ii_acc_chnl))';
%                 Mat2NlxCSC(out_file, append_flag, 1, 0, [1 0 0 0 1 1], digital_data_block_ts_usec(ii_block), data_chunk, header_motion );
%             end
        end
        
    end
    end
    
    %% remove DC
    if is_remove_DC
        for ii_channel = 1:size(data,1)
            %         data(ii_channel,:) = data(ii_channel,:) - mean(data(ii_channel,:));  % TAMIR: we prefer not to use this since data could be biased around real zero voltage (artifacts for example...)
            data(ii_channel,:) = data(ii_channel,:) - zero_DC_level_bit;         % There should be a theoretical constant representing the real zero voltage, but because the INTAN have a DC dependency on freq it is not perfect
        end
    end
    
    %% if we recorded with GND as ref. channel we want to have a
    % "post-recording ref. channel"
    if use_post_rec_ref_channel
        gain_relative_to_ref_chnl = 1.*ones(1,16);
        for ii_channel = 1:size(data,1)
            if ii_channel == post_rec_ref_channel
                continue;
            end
            data(ii_channel,:) = data(ii_channel,:) - gain_relative_to_ref_chnl(ii_channel).*data(post_rec_ref_channel,:);
        end
    end
    
    %% change to uVolt units and invert
    data = data.*uVolt_per_bit;
    if is_invert_data
        % invert data
        data = -data;
    end
    
    
    %% write data chunks to Nlx file format
    %     file_TS_usec = FileStarted_TS(ii_file_start_entry);
    disp(['writing file ' file_str])
    for cnl = data_cnl_ind
        cnl_data = data(cnl+1,:);
        cnl_data_blocks = vec2mat(cnl_data,512)';
        num_blocks = size(cnl_data_blocks, 2);
        blocks_timestamps_usec = file_TS_usec + (0:block_period_time_usec:(block_period_time_usec*(num_blocks-1)));
        file_name = ['CSC' num2str(cnl) '.ncs'];
        out_file = fullfile(Nlx_OutDir, file_name);
        
        %         we do this because Nlx functions have bugs with writing header
        %         when using append... (this way works...)
        if exist(out_file, 'file')
            append_flag = 1;
        else
            append_flag = 0;
        end
        %append_flag
        
        Mat2NlxCSC(out_file, append_flag, 1, 0, [1 0 0 0 1 1], blocks_timestamps_usec, cnl_data_blocks, header_neural );
    end
end

%%
disp('Finished converting all files from Nlg format to Nlx format!')









%%