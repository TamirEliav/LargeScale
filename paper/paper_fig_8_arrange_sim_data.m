%% ========================================================================
% Code loads input and output binary maps
% Naming convention:
% variables ending in S : Single field CA3
% variables ending in M : Multi field CA3
% variables ending in P : Periodic MEC input (1D slices of 2D hexagonal grid)
% variables ending in Xw: Input from both CA3 and MEC, varying relative synaptic weight
% variables ending in Xs: Input from both CA3 and MEC, varying relative sparsity

% variable dimensions: 
% L        = 1000 spatial bins
% repangle = 30   1D slice angles
% nrCA3    = 20   values of the relative MEC/CA3 input
% Nmax     = 2000 input neurons 
% rep      = 1000 output neurons for each parameter combination

% for example, binary map variable names and dimensions
% input:
% fS: Single field CA3, 1000 x 2000      (spatial bins x examples)
% fM: Multi  field CA3, 1000 x 2000      (spatial bins x examples)
% fP: Periodic        , 1000 x 2000 x 30 (spatial bins x examples x angles)

% output:
% fResS : Single field CA3,            1000 x 1000      (spatial bins x examples)
% fResM : Multi  field CA3,            1000 x 1000      (spatial bins x examples)
% fResXs: Periodic x Single field CA3, 1000 x 1000 x 20 (spatial bins x examples x input ratios) [results for each angle saved in a different file] 
% fResXw: Periodic x Single field CA3, 1000 x 1000 x 20 (spatial bins x examples x input ratios) [results for each angle saved in a different file] 


% Based on these maps, the code computes 
% field size/number of fields/max min ratio distributions, 
% spectra 
% and produces a draft of Figure 8 (FF part) and corresponding supplementary
% note that the plots should include distributions of quantities from data.
% to ensure that these match the distributions shown elsewhere in the paper
% I did not include plotting these data here in case I missed some
% exclusion criterion.

%% data folder
sim_data_dir = 'L:\Yonatan_theory\20200706__FF_model\Data';

%% params and definitions
HS = @(x) (sign(x)+1)/2 ;
L = 20000 ; % size of environment in cm
repangle = 30 ; % number of 1D slice angles
ds      = 20 ;  % downsample factor
L       = L/ds ;
xFF     = linspace(0,200,L) ; % linear coordinate 

Nmax    = 2000      ; % number of input neurons of each type

aMEC   = 25*[1 sqrt(2) sqrt(5) sqrt(13) sqrt(23)] ;                                     % module sizes 
A2 = [1 0 ; 1/2 sqrt(3)/2 ; -1/2 sqrt(3)/2 ; -1 0 ; -1/2 -sqrt(3)/2 ; 1/2 -sqrt(3)/2] ; % hexagonal coordinates to create 2D grid 
f2MEC = @(x,a) HS(sum(abs(mod(x*A2'/a,1)-0.5).^(1/2),2)-3.45) ;                         % binary hexagonal map function 
kPow = (0:(L-1))*ds/L/4 ;                                                               % spatial frequency vector for Fourier plots
kPowmin = 1 ;                                                                           % first valid position for max (dPower/dk)
kMEC = 1./(aMEC/2.5) ;                                                                  % spatial frequencies of modules for 0 degree slice
naMEC  = length(aMEC) ;                                                                 % number of modules

r2  = (1:L) ;
xi  = zeros(L,2) ;

thMEC2 = pi/180*(0:29) ;                                                                % slice angles

[x1,x2] = meshgrid(r2,r2) ;
xplot = zeros(L^2,2) ;
xplot(:,1) = x1(:) ;
xplot(:,2) = x2(:) ;

rCA3 = 0:0.05:0.95 ;  % relative values of CA3/MEC input. 0 means only MEC, 1 only CA3
nrCA3 = length(rCA3) ;

rep    = 1000 ;
%%
nFieldS  = zeros(rep,1) ;
nFieldM  = zeros(rep,1) ;
nFieldXw = zeros(rep,nrCA3,repangle) ;
nFieldXs = zeros(rep,nrCA3,repangle) ;

FieldSizeS  = zeros(rep,200) ;
FieldSizeM  = zeros(rep,200) ;
FieldSizeXw = cell(nrCA3,repangle)  ;
FieldSizeXs = cell(nrCA3,repangle)  ;

GapSizeS  = zeros(rep,200) ;
GapSizeM  = zeros(rep,200) ;
GapSizeXw = cell(nrCA3,repangle) ;
GapSizeXs = cell(nrCA3,repangle) ;
    
for iang = 1:repangle
    for i3 = 1:nrCA3
        FieldSizeXw{i3,iang} = zeros(rep,200) ;
        FieldSizeXs{i3,iang} = zeros(rep,200) ;
        GapSizeXw{i3,iang}   = zeros(rep,200) ;
        GapSizeXs{i3,iang}   = zeros(rep,200) ;
    end
end

MaxMinRatioS = zeros(rep,1) ;
MaxMinRatioM = zeros(rep,1) ;
MaxMinRatioXw = zeros(rep,nrCA3,repangle) ;
MaxMinRatioXs = zeros(rep,nrCA3,repangle) ;

fiS  = zeros(1,L+2) ;
fiM  = zeros(1,L+2) ;
fiXw = zeros(1,L+2) ;
fiXs = zeros(1,L+2) ;

%%
load(fullfile(sim_data_dir,'CA3MEC_FFModel_Angles_1_30_2.mat'));
load(fullfile(sim_data_dir,'CA3MEC_FFModel_fP_2.mat'));

%%
for ir = 1:rep
    fiS(2:L+1) = fResS(:,ir) ;
    dfiS = diff(fiS) ;
    js = find(dfiS==+1) ;
    je = find(dfiS==-1) ;
    nFieldS(ir) = length(js) ;
    FieldSizeS(ir,1:nFieldS(ir)) = je-js ;
    GapSizeS(ir,1:nFieldS(ir)-1) = js(2:end)-je(1:end-1) ;
    MaxMinRatioS(ir) = max(FieldSizeS(ir,1:nFieldS(ir)))/min(FieldSizeS(ir,1:nFieldS(ir))) ;

    fiM(2:L+1) = fResM(:,ir) ;
    dfiM = diff(fiM) ;
    js = find(dfiM==+1) ;
    je = find(dfiM==-1) ;
    nFieldM(ir) = length(js) ;
    FieldSizeM(ir,1:nFieldM(ir)) = je-js ;
    GapSizeM(ir,1:nFieldM(ir)-1) = js(2:end)-je(1:end-1) ;
    MaxMinRatioM(ir) = max(FieldSizeM(ir,1:nFieldM(ir)))/min(FieldSizeM(ir,1:nFieldM(ir))) ;
end

for iang = 1:repangle
    load(fullfile(sim_data_dir,['CA3MEC_FFModel_fResXw_iangle_' num2str(iang) '.mat']));
    load(fullfile(sim_data_dir,['CA3MEC_FFModel_fResXs_iangle_' num2str(iang) '.mat']));
    for ir = 1:rep
        for i3 = 1:nrCA3
            fiXw(2:L+1) = fResXw(:,ir,i3) ;
            dfiXw = diff(fiXw) ;
            js = find(dfiXw==+1) ;
            je = find(dfiXw==-1) ;
            nFieldXw(ir,i3,iang) = length(js) ;
            FieldSizeXw{i3,iang}(ir,1:nFieldXw(ir,i3,iang)) = je-js ;
            GapSizeXw{i3,iang}(ir,1:nFieldXw(ir,i3,iang)-1) = js(2:end)-je(1:end-1) ;
            MaxMinRatioXw(ir,i3,iang) = max(FieldSizeXw{i3,iang}(ir,1:nFieldXw(ir,i3,iang)))/min(FieldSizeXw{i3,iang}(ir,1:nFieldXw(ir,i3,iang))) ;

            fiXs(2:L+1) = fResXs(:,ir,i3) ;
            dfiXs = diff(fiXs) ;
            js = find(dfiXs==+1) ;
            je = find(dfiXs==-1) ;
            nFieldXs(ir,i3,iang) = length(js) ;
            FieldSizeXs{i3,iang}(ir,1:nFieldXs(ir,i3,iang)) = je-js ;
            GapSizeXs{i3,iang}(ir,1:nFieldXs(ir,i3,iang)-1) = js(2:end)-je(1:end-1) ;
            MaxMinRatioXs(ir,i3,iang) = max(FieldSizeXs{i3,iang}(ir,1:nFieldXs(ir,i3,iang)))/min(FieldSizeXs{i3,iang}(ir,1:nFieldXs(ir,i3,iang))) ;

        end
    end
end
FieldSizeS(FieldSizeS(:)==0) = [] ;
FieldSizeM(FieldSizeM(:)==0) = [] ;
GapSizeS(GapSizeS(:)==0) = [] ;
GapSizeM(GapSizeM(:)==0) = [] ;
for i3 = 1:nrCA3
    for iang = 1:repangle
        FieldSizeXw{i3,iang}(FieldSizeXw{i3,iang}(:)==0) = [] ;
        FieldSizeXs{i3,iang}(FieldSizeXs{i3,iang}(:)==0) = [] ;
        GapSizeXw{i3,iang}(GapSizeXw{i3,iang}(:)==0) = [] ;
        GapSizeXs{i3,iang}(GapSizeXs{i3,iang}(:)==0) = [] ;
    end
end
%%
% clrM = [1 0 0 ; 0 0 1 ; 0.7 0.7 0 ; 0.5 0 0.5] ;
clrM = [0 0 1 ; 
        0.0742    0.6758    0.0195 ;
        1 0 0 ;
        0.5 0 0.5] ;
lw   = 4 ;

bField = 0:1:30 ;
xField = bField(1:end-1)+(bField(2)-bField(1))/2 ;
bGap = 0:1:100 ;
xGap = bGap(1:end-1)+(bGap(2)-bGap(1))/2 ;
bnField = 0:1:50 ;
xnField = bnField(1:end-1)+(bnField(2)-bnField(1))/2 ;
bMaxMinRatio = 1:4:100 ;
xMaxMinRatio = bMaxMinRatio(1:end-1)+(bMaxMinRatio(2)-bMaxMinRatio(1))/2 ;

powInS  = mean(abs(fft(fS,[],1)).^2,2)/L ;
powInM  = mean(abs(fft(fM,[],1)).^2,2)/L ;
powInP = zeros(L,repangle) ;
for iang = 1:repangle
    powInP(:,iang) = mean(abs(fft(fP(:,:,iang),[],1)).^2,2)/L ; 
end
%%
powOutS  = mean(abs(fft(fResS,[],1)).^2,2)/L ;
powOutM  = mean(abs(fft(fResM,[],1)).^2,2)/L ; 
powOutXw = zeros(L,nrCA3,repangle) ;
powOutXs = zeros(L,nrCA3,repangle) ;
for iang = 1:repangle
    load(fullfile(sim_data_dir,['CA3MEC_FFModel_fResXw_iangle_' num2str(iang) '.mat']));
    load(fullfile(sim_data_dir,['CA3MEC_FFModel_fResXs_iangle_' num2str(iang) '.mat']));
    for i3 = 1:nrCA3
        powOutXw(:,i3,iang) = mean(abs(fft(fResXw(:,:,i3),[],1)).^2,2)/L ; 
        powOutXs(:,i3,iang) = mean(abs(fft(fResXs(:,:,i3),[],1)).^2,2)/L ; 
    end
end

MaxdpowOutS = max(diff(powOutS(kPowmin:L/2)))/(kPow(2)-kPow(1)) ;
MaxdpowOutM = max(diff(powOutM(kPowmin:L/2)))/(kPow(2)-kPow(1)) ;
MaxdpowOutXw = squeeze(max(diff(powOutXw(kPowmin:L/2,:,:),[],1),[],1)/(kPow(2)-kPow(1))) ;
MaxdpowOutXs = squeeze(max(diff(powOutXs(kPowmin:L/2,:,:),[],1),[],1)/(kPow(2)-kPow(1))) ;


rMaxdpowOutXw = zeros(nrCA3+1,repangle) ;
rMaxdpowOutXs = zeros(nrCA3+1,repangle) ;
for iang = 1:repangle
    rMaxdpowOutXw(:,iang) = [MaxdpowOutXw(:,iang)' MaxdpowOutS] ;
    rMaxdpowOutXw(:,iang) = rMaxdpowOutXw(:,iang)/max(rMaxdpowOutXw(:,iang)) ;
    rMaxdpowOutXs(:,iang) = [MaxdpowOutXs(:,iang)' MaxdpowOutS] ;
    rMaxdpowOutXs(:,iang) = rMaxdpowOutXs(:,iang)/max(rMaxdpowOutXs(:,iang)) ;
end
%%
pFieldSizeS  = histcounts(FieldSizeS*ds/100,bField,'Normalization','pdf') ;
pFieldSizeM  = histcounts(FieldSizeM*ds/100,bField,'Normalization','pdf') ;
pFieldSizeXw = zeros(nrCA3,repangle,length(xField)) ;
pFieldSizeXs = zeros(nrCA3,repangle,length(xField)) ;
for iang = 1:repangle
    for i3 = 1:nrCA3
        pFieldSizeXw(i3,iang,:)  = histcounts(FieldSizeXw{i3,iang}*ds/100,bField,'Normalization','pdf') ;
        pFieldSizeXs(i3,iang,:)  = histcounts(FieldSizeXs{i3,iang}*ds/100,bField,'Normalization','pdf') ;
    end
end

pGapSizeS  = histcounts(GapSizeS*ds/100,bGap,'Normalization','pdf') ;
pGapSizeM  = histcounts(GapSizeM*ds/100,bGap,'Normalization','pdf') ;
pGapSizeXw = zeros(nrCA3,repangle,length(xGap)) ;
pGapSizeXs = zeros(nrCA3,repangle,length(xGap)) ;
for iang = 1:repangle
    for i3 = 1:nrCA3
        pGapSizeXw(i3,iang,:)  = histcounts(GapSizeXw{i3,iang}*ds/100,bGap,'Normalization','pdf') ;
        pGapSizeXs(i3,iang,:)  = histcounts(GapSizeXs{i3,iang}*ds/100,bGap,'Normalization','pdf') ;
    end
end

pnFieldS  = histcounts(nFieldS,bnField,'Normalization','pdf') ;
pnFieldM  = histcounts(nFieldM,bnField,'Normalization','pdf') ;
pnFieldXw = zeros(nrCA3,repangle,length(xnField)) ;
pnFieldXs = zeros(nrCA3,repangle,length(xnField)) ;
for iang = 1:repangle
    for i3 = 1:nrCA3
        pnFieldXw(i3,iang,:)  = histcounts(nFieldXw(:,i3,iang),bnField,'Normalization','pdf') ;
        pnFieldXs(i3,iang,:)  = histcounts(nFieldXs(:,i3,iang),bnField,'Normalization','pdf') ;
    end
end

pMaxMinRatioS  = histcounts(MaxMinRatioS,bMaxMinRatio,'Normalization','pdf') ;
pMaxMinRatioM  = histcounts(MaxMinRatioM,bMaxMinRatio,'Normalization','pdf') ;
pMaxMinRatioXw = zeros(nrCA3,repangle,length(xMaxMinRatio)) ;
pMaxMinRatioXs = zeros(nrCA3,repangle,length(xMaxMinRatio)) ;
for iang = 1:repangle
    for i3 = 1:nrCA3
        pMaxMinRatioXw(i3,iang,:)  = histcounts(MaxMinRatioXw(:,i3,iang),bMaxMinRatio,'Normalization','pdf') ;
        pMaxMinRatioXs(i3,iang,:)  = histcounts(MaxMinRatioXs(:,i3,iang),bMaxMinRatio,'Normalization','pdf') ;
    end
end
%%
clrrCA3 = zeros(nrCA3+2,3) ;
clrrCA3(:,1) = linspace(1,0,nrCA3+2) ;
clrrCA3(:,2) = linspace(0,0,nrCA3+2) ;
clrrCA3(:,3) = linspace(0,1,nrCA3+2) ;

nplot1 = 10 ;
nplot2 = 20 ;
iang = 1 ;
i3   = 1 ;
colormap('gray') ;