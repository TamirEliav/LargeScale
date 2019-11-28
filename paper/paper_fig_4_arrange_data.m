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
load('L:\Theory\Yonatan_code_data\Multiscale_PV_ML_decoder_2_Summary_1.mat')
nL = length(L);
ds = 20 ;
Nint = 10:1:250 ;
 
%% arrange data
clr  = [1.0 0.0 0.5   ; ...
        0.0 1.0 0.5   ; ...
        1.0 0.5 0.0   ; ...
        0.25 0.75 1.0 ; ...
        0.5 0.0 1.0   ; ...
        0.5 0.9 1.0   ; ...
        0.85 0.4 1.0  ; ...
        1.0 0.75 0.25] ;
    
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

NerrMLA = zeros(nL,nerrs) ;
NerrMLB = zeros(nL,nerrs) ;
NerrMLC = zeros(nL,nerrs) ;
NerrMLD = zeros(nL,nerrs) ;
NerrMLE = zeros(nL,nerrs) ;
NerrMLF = zeros(nL,nerrs) ;
NerrMLG = zeros(nL,nerrs) ;
NerrMLH = zeros(nL,nerrs) ;

NperrPVA = zeros(nL,nerrs) ;
NperrPVB = zeros(nL,nerrs) ;
NperrPVC = zeros(nL,nerrs) ;
NperrPVD = zeros(nL,nerrs) ;
NperrPVE = zeros(nL,nerrs) ;
NperrPVF = zeros(nL,nerrs) ;
NperrPVG = zeros(nL,nerrs) ;
NperrPVH = zeros(nL,nerrs) ;

NperrMLA = zeros(nL,nerrs) ;
NperrMLB = zeros(nL,nerrs) ;
NperrMLC = zeros(nL,nerrs) ;
NperrMLD = zeros(nL,nerrs) ;
NperrMLE = zeros(nL,nerrs) ;
NperrMLF = zeros(nL,nerrs) ;
NperrMLG = zeros(nL,nerrs) ;
NperrMLH = zeros(nL,nerrs) ;

NerrPVBfit  = cell(1,nerrs) ;
NerrPVDfit  = cell(1,nerrs) ;
NerrPVEfit  = cell(1,nerrs) ;

NerrMLBfit  = cell(1,nerrs) ;
NerrMLDfit  = cell(1,nerrs) ;
NerrMLEfit  = cell(1,nerrs) ;

for jerr = 1:nerrs
    for jL = 1:nL
%         jjNPVA = find(mePVA(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNPVB = find(mePVB(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNPVC = find(mePVC(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNPVD = find(mePVD(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNPVE = find(mePVE(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNPVF = find(mePVF(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNPVG = find(mePVG(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         
%         jjNMLA = find(meMLA(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNMLB = find(meMLB(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNMLC = find(meMLC(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNMLD = find(meMLD(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNMLE = find(meMLE(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNMLF = find(meMLF(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         jjNMLG = find(meMLG(:,jL)*ds/100<errs(jerr),1,'first') ; 
%         
%         if isempty(jjNPVA) NerrPVA(jL,jerr) = N(end)+1 ; else NerrPVA(jL,jerr) = N(jjNPVA) ; end ;
%         if isempty(jjNPVB) NerrPVB(jL,jerr) = N(end)+1 ; else NerrPVB(jL,jerr) = N(jjNPVB) ; end ;
%         if isempty(jjNPVC) NerrPVC(jL,jerr) = N(end)+1 ; else NerrPVC(jL,jerr) = N(jjNPVC) ; end ;
%         if isempty(jjNPVD) NerrPVD(jL,jerr) = N(end)+1 ; else NerrPVD(jL,jerr) = N(jjNPVD) ; end ;
%         if isempty(jjNPVE) NerrPVE(jL,jerr) = N(end)+1 ; else NerrPVE(jL,jerr) = N(jjNPVE) ; end ;
%         if isempty(jjNPVF) NerrPVF(jL,jerr) = N(end)+1 ; else NerrPVF(jL,jerr) = N(jjNPVF) ; end ;
%         if isempty(jjNPVG) NerrPVG(jL,jerr) = N(end)+1 ; else NerrPVG(jL,jerr) = N(jjNPVG) ; end ;
%         
%         if isempty(jjNMLA) NerrMLA(jL,jerr) = N(end)+1 ; else NerrMLA(jL,jerr) = N(jjNMLA) ; end ;
%         if isempty(jjNMLB) NerrMLB(jL,jerr) = N(end)+1 ; else NerrMLB(jL,jerr) = N(jjNMLB) ; end ;
%         if isempty(jjNMLC) NerrMLC(jL,jerr) = N(end)+1 ; else NerrMLC(jL,jerr) = N(jjNMLC) ; end ;
%         if isempty(jjNMLD) NerrMLD(jL,jerr) = N(end)+1 ; else NerrMLD(jL,jerr) = N(jjNMLD) ; end ;
%         if isempty(jjNMLE) NerrMLE(jL,jerr) = N(end)+1 ; else NerrMLE(jL,jerr) = N(jjNMLE) ; end ;
%         if isempty(jjNMLF) NerrMLF(jL,jerr) = N(end)+1 ; else NerrMLF(jL,jerr) = N(jjNMLF) ; end ;
%         if isempty(jjNMLG) NerrMLG(jL,jerr) = N(end)+1 ; else NerrMLG(jL,jerr) = N(jjNMLG) ; end ;

        jjNPVA = find(interp1(N,mePVA(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVB = find(interp1(N,mePVB(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVC = find(interp1(N,mePVC(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVD = find(interp1(N,mePVD(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVE = find(interp1(N,mePVE(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVF = find(interp1(N,mePVF(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVG = find(interp1(N,mePVG(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNPVH = find(interp1(N,mePVH(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        
        jjNMLA = find(interp1(N,meMLA(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLB = find(interp1(N,meMLB(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLC = find(interp1(N,meMLC(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLD = find(interp1(N,meMLD(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLE = find(interp1(N,meMLE(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLF = find(interp1(N,meMLF(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLG = find(interp1(N,meMLG(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        jjNMLH = find(interp1(N,meMLH(:,jL)*ds/100,Nint,'spline')<errs(jerr),1,'first') ; 
        
        if isempty(jjNPVA) NerrPVA(jL,jerr) = N(end)+1 ; else NerrPVA(jL,jerr) = Nint(jjNPVA) ; end ;
        if isempty(jjNPVB) NerrPVB(jL,jerr) = N(end)+1 ; else NerrPVB(jL,jerr) = Nint(jjNPVB) ; end ;
        if isempty(jjNPVC) NerrPVC(jL,jerr) = N(end)+1 ; else NerrPVC(jL,jerr) = Nint(jjNPVC) ; end ;
        if isempty(jjNPVD) NerrPVD(jL,jerr) = N(end)+1 ; else NerrPVD(jL,jerr) = Nint(jjNPVD) ; end ;
        if isempty(jjNPVE) NerrPVE(jL,jerr) = N(end)+1 ; else NerrPVE(jL,jerr) = Nint(jjNPVE) ; end ;
        if isempty(jjNPVF) NerrPVF(jL,jerr) = N(end)+1 ; else NerrPVF(jL,jerr) = Nint(jjNPVF) ; end ;
        if isempty(jjNPVG) NerrPVG(jL,jerr) = N(end)+1 ; else NerrPVG(jL,jerr) = Nint(jjNPVG) ; end ;
        if isempty(jjNPVH) NerrPVH(jL,jerr) = N(end)+1 ; else NerrPVH(jL,jerr) = Nint(jjNPVH) ; end ;
        
        if isempty(jjNMLA) NerrMLA(jL,jerr) = N(end)+1 ; else NerrMLA(jL,jerr) = Nint(jjNMLA) ; end ;
        if isempty(jjNMLB) NerrMLB(jL,jerr) = N(end)+1 ; else NerrMLB(jL,jerr) = Nint(jjNMLB) ; end ;
        if isempty(jjNMLC) NerrMLC(jL,jerr) = N(end)+1 ; else NerrMLC(jL,jerr) = Nint(jjNMLC) ; end ;
        if isempty(jjNMLD) NerrMLD(jL,jerr) = N(end)+1 ; else NerrMLD(jL,jerr) = Nint(jjNMLD) ; end ;
        if isempty(jjNMLE) NerrMLE(jL,jerr) = N(end)+1 ; else NerrMLE(jL,jerr) = Nint(jjNMLE) ; end ;
        if isempty(jjNMLF) NerrMLF(jL,jerr) = N(end)+1 ; else NerrMLF(jL,jerr) = Nint(jjNMLF) ; end ;
        if isempty(jjNMLG) NerrMLG(jL,jerr) = N(end)+1 ; else NerrMLG(jL,jerr) = Nint(jjNMLG) ; end ;
        if isempty(jjNMLH) NerrMLH(jL,jerr) = N(end)+1 ; else NerrMLH(jL,jerr) = Nint(jjNMLH) ; end ;
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
        
        jjNMLA = find(interp1(N,meMLA(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLB = find(interp1(N,meMLB(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLC = find(interp1(N,meMLC(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLD = find(interp1(N,meMLD(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLE = find(interp1(N,meMLE(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLF = find(interp1(N,meMLF(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLG = find(interp1(N,meMLG(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        jjNMLH = find(interp1(N,meMLH(:,jL)*ds/100,Nint,'spline')<L(jL)/100*perrs(jerr),1,'first') ; 
        
        if isempty(jjNPVA) NperrPVA(jL,jerr) = N(end)+1 ; else NperrPVA(jL,jerr) = Nint(jjNPVA) ; end ;
        if isempty(jjNPVB) NperrPVB(jL,jerr) = N(end)+1 ; else NperrPVB(jL,jerr) = Nint(jjNPVB) ; end ;
        if isempty(jjNPVC) NperrPVC(jL,jerr) = N(end)+1 ; else NperrPVC(jL,jerr) = Nint(jjNPVC) ; end ;
        if isempty(jjNPVD) NperrPVD(jL,jerr) = N(end)+1 ; else NperrPVD(jL,jerr) = Nint(jjNPVD) ; end ;
        if isempty(jjNPVE) NperrPVE(jL,jerr) = N(end)+1 ; else NperrPVE(jL,jerr) = Nint(jjNPVE) ; end ;
        if isempty(jjNPVF) NperrPVF(jL,jerr) = N(end)+1 ; else NperrPVF(jL,jerr) = Nint(jjNPVF) ; end ;
        if isempty(jjNPVG) NperrPVG(jL,jerr) = N(end)+1 ; else NperrPVG(jL,jerr) = Nint(jjNPVG) ; end ;
        if isempty(jjNPVH) NperrPVH(jL,jerr) = N(end)+1 ; else NperrPVH(jL,jerr) = Nint(jjNPVH) ; end ;
        
        if isempty(jjNMLA) NperrMLA(jL,jerr) = N(end)+1 ; else NperrMLA(jL,jerr) = Nint(jjNMLA) ; end ;
        if isempty(jjNMLB) NperrMLB(jL,jerr) = N(end)+1 ; else NperrMLB(jL,jerr) = Nint(jjNMLB) ; end ;
        if isempty(jjNMLC) NperrMLC(jL,jerr) = N(end)+1 ; else NperrMLC(jL,jerr) = Nint(jjNMLC) ; end ;
        if isempty(jjNMLD) NperrMLD(jL,jerr) = N(end)+1 ; else NperrMLD(jL,jerr) = Nint(jjNMLD) ; end ;
        if isempty(jjNMLE) NperrMLE(jL,jerr) = N(end)+1 ; else NperrMLE(jL,jerr) = Nint(jjNMLE) ; end ;
        if isempty(jjNMLF) NperrMLF(jL,jerr) = N(end)+1 ; else NperrMLF(jL,jerr) = Nint(jjNMLF) ; end ;
        if isempty(jjNMLG) NperrMLG(jL,jerr) = N(end)+1 ; else NperrMLG(jL,jerr) = Nint(jjNMLG) ; end ;
        if isempty(jjNMLH) NperrMLH(jL,jerr) = N(end)+1 ; else NperrMLH(jL,jerr) = Nint(jjNMLH) ; end ;
    end
end
