% making order in Shir's data
%%
clear
clc
close all

%% 
% bat = 184;
% bat = 2382;
bat = 2311;
% bat = 194;
% bat = 9845;
switch bat
    %% bat 0184
    case 184
        exp_ID = 'b0184_d200102';
        dir_IN = 'Z:\students_personal_backup_space\Shir\Shir\Project\bat0184_exp_structs';
        load(fullfile(dir_IN,'pos',[exp_ID '_exp_pos.mat']))
        load(fullfile(dir_IN,'flight',[exp_ID '_exp_flight.mat']))
        FE = flight.FE;

        % ts = [42996716000 43398344000; 43448478283 46758852405; 46833804000 47374483000]; % b0184_d191128
        % ts = [48990278000 49360176767; 49385702555 52410745642; 52453715000 52924663000]; % b0184_d191129
        % ts = [33421026000 33743453000; 33814579705 36571517679; 36638738000 37377320000]; % b0184_d191130
        % ts = [47651774000 48084419000; 48145514378 50853468195; 50916426000 51462587000]; % b0184_d191201
        % ts = [43297999000 43676890000; 43703190036 46198524689; 46248112000 47441443000]; % b0184_d191202
        % ts = [42008508000 42509030916; 42549382153 45847232000; 45891038000 46349301000]; % b0184_d191203
        % ts = [53072884000 53525554248; 53561649000 56039684731; 56073933000 56700724000]; % b0184_d191204
        % ts = [45451593000 45959709000; 46014682197 48249405865; 48306648000 49505503000]; % b0184_d191205
        % ts = [58751273000 59087557000; 59125261260 61614240000; 61553064299 62173276000]; % b0184_d191208
        % ts = [42250304000 42603202524; 42679783380 44921609000; 44987515000 45444983000]; % b0184_d191209
        % ts = [53267796000 53844871000; 53899448613 56206936743; 56261380000 57003580000]; % b0184_d191210
        % ts = [48806414000 49409970000; 49445590711 51784103472; 51846449000 52225157000]; % b0184_d191212
        % ts = [44592759000 45045998000; 45075851632 46992051329; 47063238000 47623872000]; % b0184_d191215
        % ts = [49925120000 50551000000; 50620905164 52988270184; 53054221000 54394434000]; % b0184_d191216
        % ts = [35936083000 36452649000; 36518390854 38668695143; 38731961000 39212575000]; % b0184_d191217
        % ts = [42533979000 42859936000; 45294018000 45297903000; 45529647000 45529647000]; % b0184_d191218 - no behavior... ignore
        % ts = [41969900000 42287386000; 42322454185 44507265110; 44582143000 45129013000]; % b0184_d191220
        % ts = [57783002000 58687385000; 58728865921 60602453000; 60627515000 60978019000]; % b0184_d191222
        % ts = [40889960000 41276030000; 41308462435 43418558820; 43474109000 43928114000]; % b0184_d191224
        % ts = [47714400000 48126668590; 48189894153 49640148066; 49682755000 50216473000]; % b0184_d191225
        % ts = [39927708000 40275994244; 40367561459 41924687825; 41976757000 42307975000]; % b0184_d191226
        % ts = [55929622000 56553213000; 55929622000 56553213000; 58927717000 59452227000]; % b0184_d191229 - bs sync problem... ignore
        % ts = [40328462000 40669666471; 40740890384 42966271431; 43055435000 43631985000]; % b0184_d200101
        % ts = [39954098000 40362586512; 40486117434 42625785856; 42710397336 43214325000]; % b0184_d200102
        labels = {'Sleep1 start','Sleep1 end';'Behave start','Behave end';'Sleep2 start','Sleep2 end'};
        
    %% bat 2382
    case 2382
        ii_exp = 34;
        exp_list = {
            'b2382_d190623';
            'b2382_d190624';
            'b2382_d190625';
            'b2382_d190626';
            'b2382_d190627';
            'b2382_d190628';
            'b2382_d190630';
            'b2382_d190701';
            'b2382_d190703';
            'b2382_d190707';
            'b2382_d190708';
            'b2382_d190709';
            'b2382_d190712';
            'b2382_d190714';
            'b2382_d190715';
            'b2382_d190716';
            'b2382_d190718';
            'b2382_d190721';
            'b2382_d190722';
            'b2382_d190724';
            'b2382_d190725';
            'b2382_d190728';
            'b2382_d190729';
            'b2382_d190730';
            'b2382_d190731';
            'b2382_d190801';
            'b2382_d190804';
            'b2382_d190805';
            'b2382_d190807';
            'b2382_d190808';
            'b2382_d190811';
            'b2382_d190812';
            'b2382_d190813';
            'b2382_d190814';
        };
        exp_ID = exp_list{ii_exp};
        exp_create_details(exp_ID);
        exp=exp_load_data(exp_ID,'details','path','pos','flight');
        ts = exp.details.session_ts;
        labels = exp.details.session_names;
        labels = [labels+" start" labels+" end"]';
        labels = reshape({labels{:}},2,[])';
        pos = exp.pos;
        FE = exp.flight.FE;
        
    %% bat 2311
    case 2311
        ii_exp = 1;
        exp_list = {
            'b2311_d191218';
            'b2311_d191219';
            'b2311_d191220';
            'b2311_d191222';
            'b2311_d191223';
            'b2311_d191224';
            'b2311_d191225';
            'b2311_d191226';
            'b2311_d191229';
            'b2311_d191230';
            'b2311_d191231';
            'b2311_d200101';
            'b2311_d200102';
        };
        exp_ID = exp_list{ii_exp};
        exp_create_details(exp_ID);
        exp=exp_load_data(exp_ID,'details','path','pos','flight');
        ts = exp.details.session_ts;
        labels = exp.details.session_names;
        labels = [labels+" start" labels+" end"]';
        labels = reshape({labels{:}},2,[])';
        pos = exp.pos;
        FE = exp.flight.FE;
     
    %% bat 194
    case 194
        ii_exp = 1
        exp_list = {
            'b0194_d180429'; % TODO
            'b0194_d180501';
            'b0194_d180502';
            'b0194_d180503';
            'b0194_d180505';
            'b0194_d180507';
            'b0194_d180508';
            'b0194_d180509';
            'b0194_d180510';
            'b0194_d180513';
            'b0194_d180514';
            'b0194_d180515';
            'b0194_d180516';
            'b0194_d180520';
            'b0194_d180521';
            'b0194_d180522';
            'b0194_d180528';
            'b0194_d180531';
            'b0194_d180604';
            'b0194_d180605';
            'b0194_d180606';
            'b0194_d180607';
            'b0194_d180610';
            'b0194_d180611';
            'b0194_d180612';
            'b0194_d180614';
                };
        exp_ID = exp_list{ii_exp};
        exp_create_details(exp_ID);
        exp=exp_load_data(exp_ID,'details','path','pos','flight');
        ts = exp.details.session_ts;
        labels = exp.details.session_names;
        labels = [labels+" start" labels+" end"]';
        labels = reshape({labels{:}},2,[])';
        pos = exp.pos;
        FE = exp.flight.FE;
    %% bat 9845
    case 9845
        ii_exp = 22
        exp_list = {
            'b9845_d170212'; % TODO: pre-proc pos+flight
            'b9845_d170213'; % TODO: pre-proc pos+flight
            'b9845_d170215'; % TODO: pre-proc pos+flight
            'b9845_d170216'; % TODO: problem with bespoon. I think sync is bad (timestemps are wrong)
            'b9845_d170217';
            'b9845_d170219';
            'b9845_d170220';
            'b9845_d170221';
            'b9845_d170222';
            'b9845_d170223';
            'b9845_d170516';
            'b9845_d170517';
            'b9845_d170518';
            'b9845_d170525';
            'b9845_d170527'; % TODO: the tag (bsp_pos is moving) during sleep1, but I think this is shir walking because the seed is rather slow...
            'b9845_d170528';
            'b9845_d170529';
            'b9845_d170603';
            'b9845_d170605';
            'b9845_d170606'; % TODO: pre-proc pos+flight
            'b9845_d170608'; % TODO: pre-proc pos+flight
            'b9845_d170612';
            'b9845_d170614';
            'b9845_d170615'; % TODO: pre-proc pos+flight
            'b9845_d170622';
            'b9845_d170628'; % TODO: pre-proc pos+flight
            'b9845_d170709'; % TODO: pre-proc pos+flight
                };
        exp_ID = exp_list{ii_exp}
        exp_create_details(exp_ID);
        exp=exp_load_data(exp_ID,'details','path','pos','flight');
        ts = exp.details.session_ts;
        labels = exp.details.session_names;
        labels = [labels+" start" labels+" end"]';
        labels = reshape({labels{:}},2,[])';
        pos = exp.pos;
        FE = exp.flight.FE;
end
%%
exp_ID
% lp = {'r-','r--';'g-','g--';'b-','b--';'c-','c--'};
lp = {'r-','r--';'g-','g--';'b-','b--';'c-','c--';'m-','m--'};
lp = [lp;lp];
lp = lp(1:size(labels,1),:);
fig=figure;
fig.WindowState = 'maximized';
hold on
plot(pos.proc_1D.ts, pos.proc_1D.pos, '-k')
plot([FE.ts],[FE.pos],'.r')
scale_bar_min = 10;
scale_bar_x = ts(1)+[0 scale_bar_min*60*1e6];
scale_bar_y = [0 0];
plot(scale_bar_x, scale_bar_y, 'k','Clipping','off','LineWidth',2);
text(scale_bar_x(1),3,scale_bar_min+" min",'FontSize',14);
cellfun(@(x,p,lbl)(xline(x,p,lbl)), num2cell(ts), lp, labels)
% for ii=1:size(labels,1)
%     xline(ts(ii,1),'-')
% end
title(exp_ID)
zoom on






%%
