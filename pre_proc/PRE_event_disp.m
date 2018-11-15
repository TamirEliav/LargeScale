function PRE_event_disp(exp_ID)



%%
[exp_path exp_info] = DS_get_path(exp_ID);

%%
T = table();

%% free text
EV_filename = fullfile(exp_path.nlx,'EVENTS__Free text.nev');
FieldSelectionFlags = [1 0 0 0 1];
[Timestamps,EventStrings] = Nlx2MatEV(EV_filename,FieldSelectionFlags,0,1,[]);
% EventStrings = regexprep(EventStrings, 'Free text - ', '');
T1 = table();
T1.ts = Timestamps';
T1.str = EventStrings;

%% mode change
EV_filename = fullfile(exp_path.nlx,'EVENTS__Mode change.nev');
FieldSelectionFlags = [1 0 0 0 1];
[Timestamps,EventStrings] = Nlx2MatEV(EV_filename,FieldSelectionFlags,0,1,[]);
T2 = table();
T2.ts = Timestamps';
T2.str = EventStrings;

%% merge 
T = [T1;T2];
T = sortrows(T, 'ts');

%% disp all
disp(['exp ID: ' exp_ID])
sdf=get(0,'Format');
format long
T
eval(['format ' sdf]);



end