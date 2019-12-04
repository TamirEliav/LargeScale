%% arrange data for Fig. 4 - Theoretical analysis

%%
% note that the conversion of model names from files and code to figures is
% A --> A
% B --> B
% C --> C
% D --> D
% E --> E
% F --> D'
% G --> E'
% H --> C'
% I --> F (not currently included)
% J --> F' (not currently included)

%% variable names:
% me .. mean error
% pe .. percentile error
% ple .. probability of large error
% Nerr .. min N for absolute error
% Nperr .. min N for relative error


%% load data
switch dt
    case 0.5
        load('L:\Theory\Yonatan_code_data\Multiscale_PV_ML_decoder_2_Summary_1.mat');
        load('L:\Theory\Yonatan_code_data\Multiscale_PV_ML_decoder_2IJ_Summary_1.mat')
    case 0.2
        load('L:\Theory\Yonatan_code_data\Multiscale_PV_ML_decoder_4_Summary_1.mat');
        load('L:\Theory\Yonatan_code_data\Multiscale_PV_ML_decoder_4IJ_Summary_1.mat')
end
nL = length(L);
ds = 20 ;
Nint = 10:1:250 ;
 
%% arrange data
clr  = [1.0 0.0 0.5   ; ... % 1
        0.0 1.0 0.5   ; ... % 2
        1.0 0.5 0.0   ; ... % 3
        0.25 0.75 1.0 ; ... % 4
        0.5 0.0 1.0   ; ... % 5
        0.5 0.9 1.0   ; ... % 6
        0.85 0.4 1.0  ; ... % 7
        1.0 0.75 0.25 ; ... % 8
        0   0    0    ; ... % 9
        0.7 0.7 0.7 ] ;     % 10
    
errs  = [1 2 5] ;
perrs = [0.01 0.02 0.05] ;
nerrs = length(errs) ;

jLfit = 1:38 ;

NerrPVA = zeros(nL,nerrs) ;
NerrPVB = zeros(nL,nerrs) ;
NerrPVC = zeros(nL,nerrs) ;
NerrPVD = zeros(nL,nerrs) ;
NerrPVE = zeros(nL,nerrs) ;
NerrPVF = zeros(nL,nerrs) ;
NerrPVG = zeros(nL,nerrs) ;
NerrPVH = zeros(nL,nerrs) ;
NerrPVI = zeros(nL,nerrs) ;
NerrPVJ = zeros(nL,nerrs) ;

NerrMLA = zeros(nL,nerrs) ;
NerrMLB = zeros(nL,nerrs) ;
NerrMLC = zeros(nL,nerrs) ;
NerrMLD = zeros(nL,nerrs) ;
NerrMLE = zeros(nL,nerrs) ;
NerrMLF = zeros(nL,nerrs) ;
NerrMLG = zeros(nL,nerrs) ;
NerrMLH = zeros(nL,nerrs) ;
NerrMLI = zeros(nL,nerrs) ;
NerrMLJ = zeros(nL,nerrs) ;

NperrPVA = zeros(nL,nerrs) ;
NperrPVB = zeros(nL,nerrs) ;
NperrPVC = zeros(nL,nerrs) ;
NperrPVD = zeros(nL,nerrs) ;
NperrPVE = zeros(nL,nerrs) ;
NperrPVF = zeros(nL,nerrs) ;
NperrPVG = zeros(nL,nerrs) ;
NperrPVH = zeros(nL,nerrs) ;
NperrPVI = zeros(nL,nerrs) ;
NperrPVJ = zeros(nL,nerrs) ;

NperrMLA = zeros(nL,nerrs) ;
NperrMLB = zeros(nL,nerrs) ;
NperrMLC = zeros(nL,nerrs) ;
NperrMLD = zeros(nL,nerrs) ;
NperrMLE = zeros(nL,nerrs) ;
NperrMLF = zeros(nL,nerrs) ;
NperrMLG = zeros(nL,nerrs) ;
NperrMLH = zeros(nL,nerrs) ;
NperrMLI = zeros(nL,nerrs) ;
NperrMLJ = zeros(nL,nerrs) ;

NerrPVBfit  = cell(1,nerrs) ;
NerrPVDfit  = cell(1,nerrs) ;
NerrPVEfit  = cell(1,nerrs) ;

NerrMLBfit  = cell(1,nerrs) ;
NerrMLDfit  = cell(1,nerrs) ;
NerrMLEfit  = cell(1,nerrs) ;

for jerr = 1:nerrs
    for jL = 1:nL
        jjNPVA = find(interp1(N,mePVA(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVB = find(interp1(N,mePVB(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVC = find(interp1(N,mePVC(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVD = find(interp1(N,mePVD(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVE = find(interp1(N,mePVE(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVF = find(interp1(N,mePVF(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVG = find(interp1(N,mePVG(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVH = find(interp1(N,mePVH(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVI = find(interp1(N,mePVI(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVJ = find(interp1(N,mePVJ(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        
        jjNMLA = find(interp1(N,meMLA(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLB = find(interp1(N,meMLB(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLC = find(interp1(N,meMLC(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLD = find(interp1(N,meMLD(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLE = find(interp1(N,meMLE(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLF = find(interp1(N,meMLF(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLG = find(interp1(N,meMLG(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLH = find(interp1(N,meMLH(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLI = find(interp1(N,meMLI(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLJ = find(interp1(N,meMLJ(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        
        if isempty(jjNPVA) && jL>1 NerrPVA(jL,jerr) = NerrPVA(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrPVA(jL-1,jerr)-NerrPVA(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrPVA(jL,jerr) = Nint(jjNPVA) ; end ;
        if isempty(jjNPVB) && jL>1 NerrPVB(jL,jerr) = NerrPVB(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrPVB(jL-1,jerr)-NerrPVB(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrPVB(jL,jerr) = Nint(jjNPVB) ; end ;
        if isempty(jjNPVD) && jL>1 NerrPVD(jL,jerr) = NerrPVD(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrPVD(jL-1,jerr)-NerrPVD(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrPVD(jL,jerr) = Nint(jjNPVD) ; end ;
        if isempty(jjNPVE) && jL>1 NerrPVE(jL,jerr) = NerrPVE(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrPVE(jL-1,jerr)-NerrPVE(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrPVE(jL,jerr) = Nint(jjNPVE) ; end ;
        if isempty(jjNPVF) && jL>1 NerrPVF(jL,jerr) = NerrPVF(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrPVF(jL-1,jerr)-NerrPVF(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrPVF(jL,jerr) = Nint(jjNPVF) ; end ;
        if isempty(jjNPVG) && jL>1 NerrPVG(jL,jerr) = NerrPVG(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrPVG(jL-1,jerr)-NerrPVG(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrPVG(jL,jerr) = Nint(jjNPVG) ; end ;
        if isempty(jjNPVI) && jL>1 NerrPVI(jL,jerr) = NerrPVI(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrPVI(jL-1,jerr)-NerrPVI(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrPVI(jL,jerr) = Nint(jjNPVI) ; end ;
        
        if isempty(jjNPVC) NerrPVC(jL,jerr) = N(end)+1 ; else NerrPVC(jL,jerr) = Nint(jjNPVC) ; end ;
        if isempty(jjNPVH) NerrPVH(jL,jerr) = N(end)+1 ; else NerrPVH(jL,jerr) = Nint(jjNPVH) ; end ;
        if isempty(jjNPVJ) NerrPVJ(jL,jerr) = N(end)+1 ; else NerrPVJ(jL,jerr) = Nint(jjNPVJ) ; end ;
        
        if isempty(jjNMLA) && jL>1 NerrMLA(jL,jerr) = NerrMLA(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrMLA(jL-1,jerr)-NerrMLA(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrMLA(jL,jerr) = Nint(jjNMLA) ; end ;
        if isempty(jjNMLB) && jL>1 NerrMLB(jL,jerr) = NerrMLB(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrMLB(jL-1,jerr)-NerrMLB(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrMLB(jL,jerr) = Nint(jjNMLB) ; end ;
        if isempty(jjNMLD) && jL>1 NerrMLD(jL,jerr) = NerrMLD(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrMLD(jL-1,jerr)-NerrMLD(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrMLD(jL,jerr) = Nint(jjNMLD) ; end ;
        if isempty(jjNMLE) && jL>1 NerrMLE(jL,jerr) = NerrMLE(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrMLE(jL-1,jerr)-NerrMLE(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrMLE(jL,jerr) = Nint(jjNMLE) ; end ;
        if isempty(jjNMLF) && jL>1 NerrMLF(jL,jerr) = NerrMLF(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrMLF(jL-1,jerr)-NerrMLF(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrMLF(jL,jerr) = Nint(jjNMLF) ; end ;
        if isempty(jjNMLG) && jL>1 NerrMLG(jL,jerr) = NerrMLG(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrMLG(jL-1,jerr)-NerrMLG(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrMLG(jL,jerr) = Nint(jjNMLG) ; end ;
        if isempty(jjNMLI) && jL>1 NerrMLI(jL,jerr) = NerrMLI(jL-1,jerr) + (L(jL)-L(jL-1))*(NerrMLI(jL-1,jerr)-NerrMLI(jL-2,jerr))/(L(jL-1)-L(jL-2)) ; else NerrMLI(jL,jerr) = Nint(jjNMLI) ; end ;
    
        if isempty(jjNMLC) NerrMLC(jL,jerr) = N(end)+1 ; else NerrMLC(jL,jerr) = Nint(jjNMLC) ; end ;
        if isempty(jjNMLH) NerrMLH(jL,jerr) = N(end)+1 ; else NerrMLH(jL,jerr) = Nint(jjNMLH) ; end ;
        if isempty(jjNMLJ) NerrMLJ(jL,jerr) = N(end)+1 ; else NerrMLJ(jL,jerr) = Nint(jjNMLJ) ; end ;
    end
    jLfiti = intersect(jLfit,find(NerrPVB(:,jerr)<=N(end))) ; NerrPVBfit{jerr} = fit(L(jLfiti)'/100,NerrPVB(jLfiti,jerr),'poly1') ;
    jLfiti = intersect(jLfit,find(NerrPVD(:,jerr)<=N(end))) ; NerrPVDfit{jerr} = fit(L(jLfiti)'/100,NerrPVD(jLfiti,jerr),'poly1') ;
    jLfiti = intersect(jLfit,find(NerrPVE(:,jerr)<=N(end))) ; NerrPVEfit{jerr} = fit(L(jLfiti)'/100,NerrPVE(jLfiti,jerr),'poly1') ;

    jLfiti = intersect(jLfit,find(NerrMLB(:,jerr)<=N(end))) ; NerrMLBfit{jerr} = fit(L(jLfiti)'/100,NerrMLB(jLfiti,jerr),'poly1') ;
    jLfiti = intersect(jLfit,find(NerrMLD(:,jerr)<=N(end))) ; NerrMLDfit{jerr} = fit(L(jLfiti)'/100,NerrMLD(jLfiti,jerr),'poly1') ;
    jLfiti = intersect(jLfit,find(NerrMLE(:,jerr)<=N(end))) ; NerrMLEfit{jerr} = fit(L(jLfiti)'/100,NerrMLE(jLfiti,jerr),'poly1') ;
end



for jerr = 1:nerrs
    for jL = 1:nL
        jjNPVA = find(interp1(N,mePVA(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNPVB = find(interp1(N,mePVB(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNPVC = find(interp1(N,mePVC(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNPVD = find(interp1(N,mePVD(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNPVE = find(interp1(N,mePVE(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNPVF = find(interp1(N,mePVF(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNPVG = find(interp1(N,mePVG(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNPVH = find(interp1(N,mePVH(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNPVI = find(interp1(N,mePVI(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNPVJ = find(interp1(N,mePVJ(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        
        jjNMLA = find(interp1(N,meMLA(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLB = find(interp1(N,meMLB(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLC = find(interp1(N,meMLC(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLD = find(interp1(N,meMLD(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLE = find(interp1(N,meMLE(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLF = find(interp1(N,meMLF(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLG = find(interp1(N,meMLG(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLH = find(interp1(N,meMLH(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLI = find(interp1(N,meMLI(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLJ = find(interp1(N,meMLJ(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        
        if isempty(jjNPVA) NperrPVA(jL,jerr) = N(end)+1 ; else NperrPVA(jL,jerr) = Nint(jjNPVA) ; end ;
        if isempty(jjNPVB) NperrPVB(jL,jerr) = N(end)+1 ; else NperrPVB(jL,jerr) = Nint(jjNPVB) ; end ;
        if isempty(jjNPVC) NperrPVC(jL,jerr) = N(end)+1 ; else NperrPVC(jL,jerr) = Nint(jjNPVC) ; end ;
        if isempty(jjNPVD) NperrPVD(jL,jerr) = N(end)+1 ; else NperrPVD(jL,jerr) = Nint(jjNPVD) ; end ;
        if isempty(jjNPVE) NperrPVE(jL,jerr) = N(end)+1 ; else NperrPVE(jL,jerr) = Nint(jjNPVE) ; end ;
        if isempty(jjNPVF) NperrPVF(jL,jerr) = N(end)+1 ; else NperrPVF(jL,jerr) = Nint(jjNPVF) ; end ;
        if isempty(jjNPVG) NperrPVG(jL,jerr) = N(end)+1 ; else NperrPVG(jL,jerr) = Nint(jjNPVG) ; end ;
        if isempty(jjNPVH) NperrPVH(jL,jerr) = N(end)+1 ; else NperrPVH(jL,jerr) = Nint(jjNPVH) ; end ;
        if isempty(jjNPVI) NperrPVI(jL,jerr) = N(end)+1 ; else NperrPVI(jL,jerr) = Nint(jjNPVI) ; end ;
        if isempty(jjNPVJ) NperrPVJ(jL,jerr) = N(end)+1 ; else NperrPVJ(jL,jerr) = Nint(jjNPVJ) ; end ;
        
        if isempty(jjNMLA) NperrMLA(jL,jerr) = N(end)+1 ; else NperrMLA(jL,jerr) = Nint(jjNMLA) ; end ;
        if isempty(jjNMLB) NperrMLB(jL,jerr) = N(end)+1 ; else NperrMLB(jL,jerr) = Nint(jjNMLB) ; end ;
        if isempty(jjNMLC) NperrMLC(jL,jerr) = N(end)+1 ; else NperrMLC(jL,jerr) = Nint(jjNMLC) ; end ;
        if isempty(jjNMLD) NperrMLD(jL,jerr) = N(end)+1 ; else NperrMLD(jL,jerr) = Nint(jjNMLD) ; end ;
        if isempty(jjNMLE) NperrMLE(jL,jerr) = N(end)+1 ; else NperrMLE(jL,jerr) = Nint(jjNMLE) ; end ;
        if isempty(jjNMLF) NperrMLF(jL,jerr) = N(end)+1 ; else NperrMLF(jL,jerr) = Nint(jjNMLF) ; end ;
        if isempty(jjNMLG) NperrMLG(jL,jerr) = N(end)+1 ; else NperrMLG(jL,jerr) = Nint(jjNMLG) ; end ;
        if isempty(jjNMLH) NperrMLH(jL,jerr) = N(end)+1 ; else NperrMLH(jL,jerr) = Nint(jjNMLH) ; end ;
        if isempty(jjNMLI) NperrMLI(jL,jerr) = N(end)+1 ; else NperrMLI(jL,jerr) = Nint(jjNMLI) ; end ;
        if isempty(jjNMLJ) NperrMLJ(jL,jerr) = N(end)+1 ; else NperrMLJ(jL,jerr) = Nint(jjNMLJ) ; end ;
    end
end

%% calc slopes of minN required for error vs. L
lmNerrMLA = {};
lmNerrMLB = {};
lmNerrMLC = {};
lmNerrMLD = {};
lmNerrMLE = {};
lmNerrMLF = {};
lmNerrMLG = {};
lmNerrMLH = {};
lmNerrMLI = {};
lmNerrPVA = {};
lmNerrPVB = {};
lmNerrPVC = {};
lmNerrPVD = {};
lmNerrPVE = {};
lmNerrPVF = {};
lmNerrPVG = {};
lmNerrPVH = {};
lmNerrPVI = {};
L_thr = 50;
L_IX = L./100 >= L_thr;
L_IX = true(size(L));
for jerr = 1:nerrs
	IX = L_IX' & NerrMLA(:,jerr)<max(N); lmNerrMLA{jerr} = fitlm(L(IX)./100, NerrMLA(IX,jerr));
	IX = L_IX' & NerrMLB(:,jerr)<max(N); lmNerrMLB{jerr} = fitlm(L(IX)./100, NerrMLB(IX,jerr));
% 	IX = L_IX' & NerrMLC(:,jerr)<max(N); lmNerrMLC{jerr} = fitlm(L(IX)./100, NerrMLC(IX,jerr));
% 	IX = L_IX' & NerrMLD(:,jerr)<max(N); lmNerrMLD{jerr} = fitlm(L(IX)./100, NerrMLD(IX,jerr));
% 	IX = L_IX' & NerrMLE(:,jerr)<max(N); lmNerrMLE{jerr} = fitlm(L(IX)./100, NerrMLE(IX,jerr));
	IX = L_IX' & NerrMLF(:,jerr)<max(N); lmNerrMLF{jerr} = fitlm(L(IX)./100, NerrMLF(IX,jerr));
	IX = L_IX' & NerrMLG(:,jerr)<max(N); lmNerrMLG{jerr} = fitlm(L(IX)./100, NerrMLG(IX,jerr));
% 	IX = L_IX' & NerrMLH(:,jerr)<max(N); lmNerrMLH{jerr} = fitlm(L(IX)./100, NerrMLH(IX,jerr));
	IX = L_IX' & NerrMLI(:,jerr)<max(N); lmNerrMLI{jerr} = fitlm(L(IX)./100, NerrMLI(IX,jerr));

	IX = L_IX' & NerrPVA(:,jerr)<max(N); lmNerrPVA{jerr} = fitlm(L(IX)./100, NerrPVA(IX,jerr));
	IX = L_IX' & NerrPVB(:,jerr)<max(N); lmNerrPVB{jerr} = fitlm(L(IX)./100, NerrPVB(IX,jerr));
% 	IX = L_IX' & NerrPVC(:,jerr)<max(N); lmNerrPVC{jerr} = fitlm(L(IX)./100, NerrPVC(IX,jerr));
% 	IX = L_IX' & NerrPVD(:,jerr)<max(N); lmNerrPVD{jerr} = fitlm(L(IX)./100, NerrPVD(IX,jerr));
% 	IX = L_IX' & NerrPVE(:,jerr)<max(N); lmNerrPVE{jerr} = fitlm(L(IX)./100, NerrPVE(IX,jerr));
	IX = L_IX' & NerrPVF(:,jerr)<max(N); lmNerrPVF{jerr} = fitlm(L(IX)./100, NerrPVF(IX,jerr));
	IX = L_IX' & NerrPVG(:,jerr)<max(N); lmNerrPVG{jerr} = fitlm(L(IX)./100, NerrPVG(IX,jerr));
% 	IX = L_IX' & NerrPVH(:,jerr)<max(N); lmNerrPVH{jerr} = fitlm(L(IX)./100, NerrPVH(IX,jerr));
	IX = L_IX' & NerrPVI(:,jerr)<max(N); lmNerrPVI{jerr} = fitlm(L(IX)./100, NerrPVI(IX,jerr));
end











