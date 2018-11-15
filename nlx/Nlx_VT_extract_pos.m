function [VT_pos, VT_ts] = Nlx_VT_extract_pos(dir_in)


%% load cameras data

% Extract the video files (x,y,z) from both cameras:
% --------------------------------------------------

Extract_Fields = [1 1 1 0 1 0] ; % For each frame extract all variables except ExtractedAngle & Points
Extract_Header = 1 ; % Extract the Header as well
Extraction_Mode = 1 ; % Extract all video frames
ExtractionModeArray = []; % Timestamp Range to Extract (it is defined above)

VT1_file = fullfile(dir_in,'VT1.nvt');
[Timestamps_v1, ExtractedX_VT1, ExtractedY_VT1, Targets_VT1, NlxHeader] = ...
    Nlx2MatVT( VT1_file, Extract_Fields, Extract_Header,Extraction_Mode,ExtractionModeArray) ; % Extract data - camera 1

VT2_file = fullfile(dir_in,'VT2.nvt');
[Timestamps_v2, ExtractedX_VT2, ExtractedY_VT2, Targets_VT2, NlxHeader] = ...
    Nlx2MatVT( VT2_file, Extract_Fields, Extract_Header, Extraction_Mode,ExtractionModeArray ) ; % Extract data - camera 2

%% 3D re-construction
% -------------------------
% Re-order the time frames:
% We do this to avoid cases where frames between two camera's will be
% shifted relative to one another due to events such as a miss of a frame
% of temporal jitter in the sampling between the two cameras:
% The algorithm will run relative to the First camera (VT1):

% We allow the frames from the two video camera's be apart by no more then
% half the sampling rate:
frames_per_second = 25;
min_time_diff = 1000/(frames_per_second/2);
% Initilalize:
Timestamps_v2_corrected = [];
Targets_VT2_corrected = [];
ExtractedX_VT2_corrected = [];
ExtractedY_VT2_corrected = [];
frame_counter = 0;

disp(['Adjusting video frames'])
for ii_frame = 1:length(Timestamps_v1)
    current_frame_timestamp = Timestamps_v1(ii_frame);
    % find the nearest (temporally) frame in the second camera:
    [val_nearest_VT2,nearest_VT2_IX] = min(abs(Timestamps_v2 - current_frame_timestamp));
    if (val_nearest_VT2/10^3 < min_time_diff)
        frame_counter = frame_counter + 1;
        Timestamps_v2_corrected(frame_counter) = Timestamps_v2(nearest_VT2_IX);
        Targets_VT2_corrected(:,frame_counter) = Targets_VT2(:,nearest_VT2_IX);
        Targets_VT2_corrected(:,frame_counter) = Targets_VT2(:,nearest_VT2_IX);
        
        ExtractedX_VT2_corrected(frame_counter) = ExtractedX_VT2(nearest_VT2_IX);
        ExtractedY_VT2_corrected(frame_counter) = ExtractedY_VT2(nearest_VT2_IX);
    else end
end
disp(['Finished adjusting video frames'])

% [x_det_v1,y_det_v1 location_LED_interp_v1 idx_notzero_v1 idx_IS_zero_v1] = Get_Luminance_Values_Michael(Targets_VT1,Timestamps_v1);
% [x_det_v2,y_det_v2 location_LED_interp_v2 idx_notzero_v2 idx_IS_zero_v2] = Get_Luminance_Values_Michael(Targets_VT2_corrected,Timestamps_v2_corrected);

% Get 3D locations - Koby's DLT code
% --------------------------------------------------
% Synchronize the two Video timestamps -
% VT marks which timestamp is selected to "lead" -
% VT=0 - camera 1, VT=1 - camera 2
diff_start =  Timestamps_v1(1) - Timestamps_v2_corrected(1);
if diff_start>0 % Meaning if VT1 started working AFTER VT2
    Timestamps_v = Timestamps_v1;
    start_idx = min(find(Timestamps_v2_corrected > Timestamps_v1(1))); % The first timestamp of VT2 after VT1 started working
    VT = 0;
elseif diff_start<0
    Timestamps_v = Timestamps_v2_corrected;
    start_idx = min(find(Timestamps_v1 > Timestamps_v2_corrected(1)));
    VT = 1;
elseif length(Timestamps_v1) < length(Timestamps_v2_corrected)
    Timestamps_v = Timestamps_v1;
    start_idx = 1;
    VT = 0;
else
    Timestamps_v = Timestamps_v1;
    start_idx = 1;
    VT = 1;
end

disp(['Run 3D recunstruction (find_xyz func)'])
tic
% [ x1, y1, z1] = NLX_VT_find_xyz(location_LED_interp_v1(1,:), location_LED_interp_v1(2,:), location_LED_interp_v2(1,:), location_LED_interp_v2(2,:));
% [ x1, y1, z1] = NLX_VT_find_xyz(x_det_v1, y_det_v1, x_det_v1, y_det_v1);
[ x1, y1, z1] = NLX_VT_find_xyz(ExtractedX_VT1, ExtractedY_VT1, ExtractedX_VT2_corrected, ExtractedY_VT2_corrected);
toc
disp(['Finished 3D recunstruction (find_xyz func)'])

VT_pos = [ x1; y1; z1];
VT_ts = Timestamps_v;

[C,IA,IC] = unique(VT_ts);
VT_ts = VT_ts(IA);
VT_pos = VT_pos(:,IA);

VT_pos = VT_pos.*1e-3; % change units to meters (from mm)




