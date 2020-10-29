%% arrange data for Fig. 7 - Theoretical analysis (decoder)

%% variable names:
% merr  .. mean error
% perr  .. percentile error
% plerr .. probability of large error
% Nerr  .. min N for absolute error
% Nperr .. min N for relative error

%% load data
% load('C:\Tamir\work\Projects\LargeScale\Yonatan_theory\20200630__new_simulations_data+script\SummaryDecoderResults_VariableEnvironmentSize.mat');
% load('C:\Tamir\work\Projects\LargeScale\Yonatan_theory\20200630__new_simulations_data+script\SummaryDecoderResults_VariableEnvironmentSize_delta0.3.mat');
load('L:\Yonatan_theory\20200710__simulations_updated\SummaryDecoderResults_VariableEnvironmentSize_delta0.3_2.mat');

%%
% L = Lv*ds/100 ; % environment size variable in meters
nL = length(L) ;
ndt = length(dt) ; % all decoder results were computed for integration times dt = [0.2 0.5]
jLfit = 1:nL;
% Nint = 10:1:250;
% errs  = [2];
% nerrs = length(errs);
% perrs = [0.05];

%%
% % % % % variables holding minimal number of neruons required for error for each
% % % % % scheme
% % % % errs = 2 ; % error in meters for which the "minimum number of required neurons is computed.
% % % % Nint = 10:1:250 ; % interpolation of N values, for computation of minimal number of neurons required for specific (2m) error
% % % % 
% % % % Nerr_ML_S1  = zeros(nL,ndt) ;
% % % % Nerr_ML_S2  = zeros(nL,ndt) ;
% % % % Nerr_ML_S3  = zeros(nL,ndt) ;
% % % % Nerr_ML_S4  = zeros(nL,ndt) ;
% % % % Nerr_ML_S5  = zeros(nL,ndt) ;
% % % % Nerr_ML_S6  = zeros(nL,ndt) ;
% % % % Nerr_ML_S5v  = zeros(nL,ndt) ;
% % % % Nerr_ML_S6v  = zeros(nL,ndt) ;
% % % % 
% % % % Nerr_PV_S1  = zeros(nL,ndt) ;
% % % % Nerr_PV_S2  = zeros(nL,ndt) ;
% % % % Nerr_PV_S3  = zeros(nL,ndt) ;
% % % % Nerr_PV_S4  = zeros(nL,ndt) ;
% % % % Nerr_PV_S5  = zeros(nL,ndt) ;
% % % % Nerr_PV_S6  = zeros(nL,ndt) ;
% % % % Nerr_PV_S5v  = zeros(nL,ndt) ;
% % % % Nerr_PV_S6v  = zeros(nL,ndt) ;
% % % % 
% % % % for jdt = 1:ndt
% % % %     for jL = 1:nL
% % % %         jjNPV_S1 = find(interp1(N,merr_PV_S1(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNPV_S2 = find(interp1(N,merr_PV_S2(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNPV_S3 = find(interp1(N,merr_PV_S3(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNPV_S4 = find(interp1(N,merr_PV_S4(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNPV_S5 = find(interp1(N,merr_PV_S5(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNPV_S6 = find(interp1(N,merr_PV_S6(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNPV_S5v = find(interp1(N,merr_PV_S5v(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNPV_S6v = find(interp1(N,merr_PV_S6v(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         
% % % %         jjNML_S1 = find(interp1(N,merr_ML_S1(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNML_S2 = find(interp1(N,merr_ML_S2(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNML_S3 = find(interp1(N,merr_ML_S3(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNML_S4 = find(interp1(N,merr_ML_S4(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNML_S5 = find(interp1(N,merr_ML_S5(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNML_S6 = find(interp1(N,merr_ML_S6(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNML_S5v = find(interp1(N,merr_ML_S5v(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         jjNML_S6v = find(interp1(N,merr_ML_S6v(:,jL,jdt)*ds/100,Nint,'spline')<errs,1,'first') ;
% % % %         
% % % %         if isempty(jjNPV_S1) && jL>1 Nerr_PV_S1(jL,jdt) = Nerr_PV_S1(jL-1,jdt) + (L(jL)-L(jL-1))*(Nerr_PV_S1(jL-1,jdt)-Nerr_PV_S1(jL-2,jdt))/(L(jL-1)-L(jL-2)) ; else Nerr_PV_S1(jL,jdt) = Nint(jjNPV_S1) ; end ;
% % % %         if isempty(jjNPV_S2) && jL>1 Nerr_PV_S2(jL,jdt) = Nerr_PV_S2(jL-1,jdt) + (L(jL)-L(jL-1))*(Nerr_PV_S2(jL-1,jdt)-Nerr_PV_S2(jL-2,jdt))/(L(jL-1)-L(jL-2)) ; else Nerr_PV_S2(jL,jdt) = Nint(jjNPV_S2) ; end ;
% % % %         if isempty(jjNPV_S3) && jL>1 Nerr_PV_S3(jL,jdt) = Nerr_PV_S3(jL-1,jdt) + (L(jL)-L(jL-1))*(Nerr_PV_S3(jL-1,jdt)-Nerr_PV_S3(jL-2,jdt))/(L(jL-1)-L(jL-2)) ; else Nerr_PV_S3(jL,jdt) = Nint(jjNPV_S3) ; end ;
% % % %         if isempty(jjNPV_S4) && jL>1 Nerr_PV_S4(jL,jdt) = Nerr_PV_S4(jL-1,jdt) + (L(jL)-L(jL-1))*(Nerr_PV_S4(jL-1,jdt)-Nerr_PV_S4(jL-2,jdt))/(L(jL-1)-L(jL-2)) ; else Nerr_PV_S4(jL,jdt) = Nint(jjNPV_S4) ; end ;
% % % %         if isempty(jjNPV_S5) && jL>1 Nerr_PV_S5(jL,jdt) = Nerr_PV_S5(jL-1,jdt) + (L(jL)-L(jL-1))*(Nerr_PV_S5(jL-1,jdt)-Nerr_PV_S5(jL-2,jdt))/(L(jL-1)-L(jL-2)) ; else Nerr_PV_S5(jL,jdt) = Nint(jjNPV_S5) ; end ;
% % % %         if isempty(jjNPV_S6) && jL>1 Nerr_PV_S6(jL,jdt) = Nerr_PV_S6(jL-1,jdt) + (L(jL)-L(jL-1))*(Nerr_PV_S6(jL-1,jdt)-Nerr_PV_S6(jL-2,jdt))/(L(jL-1)-L(jL-2)) ; else Nerr_PV_S6(jL,jdt) = Nint(jjNPV_S6) ; end ;
% % % %         if isempty(jjNPV_S5v) && jL>1 Nerr_PV_S5v(jL,jdt) = Nerr_PV_S5v(jL-1,jdt) + (L(jL)-L(jL-1))*(Nerr_PV_S5v(jL-1,jdt)-Nerr_PV_S5v(jL-2,jdt))/(L(jL-1)-L(jL-2)) ; else Nerr_PV_S5v(jL,jdt) = Nint(jjNPV_S5v) ; end ;
% % % %         if isempty(jjNPV_S6v) && jL>1 Nerr_PV_S6v(jL,jdt) = Nerr_PV_S6v(jL-1,jdt) + (L(jL)-L(jL-1))*(Nerr_PV_S6v(jL-1,jdt)-Nerr_PV_S6v(jL-2,jdt))/(L(jL-1)-L(jL-2)) ; else Nerr_PV_S6v(jL,jdt) = Nint(jjNPV_S6v) ; end ;
% % % %         
% % % %         if isempty(jjNML_S1) && jL>1 Nerr_ML_S1(jL,jdt) = Nerr_ML_S1(jL-1,1) + (L(jL)-L(jL-1))*(Nerr_ML_S1(jL-1,1)-Nerr_ML_S1(jL-2,1))/(L(jL-1)-L(jL-2)) ; else Nerr_ML_S1(jL,jdt) = Nint(jjNML_S1) ; end ;
% % % %         if isempty(jjNML_S2) && jL>1 Nerr_ML_S2(jL,jdt) = Nerr_ML_S2(jL-1,1) + (L(jL)-L(jL-1))*(Nerr_ML_S2(jL-1,1)-Nerr_ML_S2(jL-2,1))/(L(jL-1)-L(jL-2)) ; else Nerr_ML_S2(jL,jdt) = Nint(jjNML_S2) ; end ;
% % % %         if isempty(jjNML_S3) && jL>1 Nerr_ML_S3(jL,jdt) = Nerr_ML_S3(jL-1,1) + (L(jL)-L(jL-1))*(Nerr_ML_S3(jL-1,1)-Nerr_ML_S3(jL-2,1))/(L(jL-1)-L(jL-2)) ; else Nerr_ML_S3(jL,jdt) = Nint(jjNML_S3) ; end ;
% % % %         if isempty(jjNML_S4) && jL>1 Nerr_ML_S4(jL,jdt) = Nerr_ML_S4(jL-1,1) + (L(jL)-L(jL-1))*(Nerr_ML_S4(jL-1,1)-Nerr_ML_S4(jL-2,1))/(L(jL-1)-L(jL-2)) ; else Nerr_ML_S4(jL,jdt) = Nint(jjNML_S4) ; end ;
% % % %         if isempty(jjNML_S5) && jL>1 Nerr_ML_S5(jL,jdt) = Nerr_ML_S5(jL-1,1) + (L(jL)-L(jL-1))*(Nerr_ML_S5(jL-1,1)-Nerr_ML_S5(jL-2,1))/(L(jL-1)-L(jL-2)) ; else Nerr_ML_S5(jL,jdt) = Nint(jjNML_S5) ; end ;
% % % %         if isempty(jjNML_S6) && jL>1 Nerr_ML_S6(jL,jdt) = Nerr_ML_S6(jL-1,1) + (L(jL)-L(jL-1))*(Nerr_ML_S6(jL-1,1)-Nerr_ML_S6(jL-2,1))/(L(jL-1)-L(jL-2)) ; else Nerr_ML_S6(jL,jdt) = Nint(jjNML_S6) ; end ;
% % % %         if isempty(jjNML_S5v) && jL>1 Nerr_ML_S5v(jL,jdt) = Nerr_ML_S5v(jL-1,1) + (L(jL)-L(jL-1))*(Nerr_ML_S5v(jL-1,1)-Nerr_ML_S5v(jL-2,1))/(L(jL-1)-L(jL-2)) ; else Nerr_ML_S5v(jL,jdt) = Nint(jjNML_S5v) ; end ;
% % % %         if isempty(jjNML_S6v) && jL>1 Nerr_ML_S6v(jL,jdt) = Nerr_ML_S6v(jL-1,1) + (L(jL)-L(jL-1))*(Nerr_ML_S6v(jL-1,1)-Nerr_ML_S6v(jL-2,1))/(L(jL-1)-L(jL-2)) ; else Nerr_ML_S6v(jL,jdt) = Nint(jjNML_S6v) ; end ;
% % % %         
% % % %         if isempty(jjNPV_S1) Nerr_PV_S1(jL,jdt) = N(end)+1 ; else Nerr_PV_S1(jL,jdt) = Nint(jjNPV_S1) ; end ;
% % % %         if isempty(jjNPV_S2) Nerr_PV_S2(jL,jdt) = N(end)+1 ; else Nerr_PV_S2(jL,jdt) = Nint(jjNPV_S2) ; end ;
% % % %         if isempty(jjNPV_S3) Nerr_PV_S3(jL,jdt) = N(end)+1 ; else Nerr_PV_S3(jL,jdt) = Nint(jjNPV_S3) ; end ;
% % % %         if isempty(jjNPV_S4) Nerr_PV_S4(jL,jdt) = N(end)+1 ; else Nerr_PV_S4(jL,jdt) = Nint(jjNPV_S4) ; end ;
% % % %         
% % % %         if isempty(jjNML_S1) Nerr_ML_S1(jL,jdt) = N(end)+1 ; else Nerr_ML_S1(jL,jdt) = Nint(jjNML_S1) ; end ;
% % % %         if isempty(jjNML_S2) Nerr_ML_S2(jL,jdt) = N(end)+1 ; else Nerr_ML_S2(jL,jdt) = Nint(jjNML_S2) ; end ;
% % % %         if isempty(jjNML_S3) Nerr_ML_S3(jL,jdt) = N(end)+1 ; else Nerr_ML_S3(jL,jdt) = Nint(jjNML_S3) ; end ;
% % % %         if isempty(jjNML_S4) Nerr_ML_S4(jL,jdt) = N(end)+1 ; else Nerr_ML_S4(jL,jdt) = Nint(jjNML_S4) ; end ;
% % % %     end
% % % % end

%% calc slopes of minN required for 2m error vs. L
lmNerr_ML_S1 = {};
lmNerr_ML_S2 = {};
lmNerr_ML_S3 = {};
lmNerr_ML_S4 = {};
lmNerr_ML_S5 = {};
lmNerr_ML_S6 = {};
lmNerr_ML_S5v = {};
lmNerr_ML_S6v = {};

lmNerr_PV_S1 = {};
lmNerr_PV_S2 = {};
lmNerr_PV_S3 = {};
lmNerr_PV_S4 = {};
lmNerr_PV_S5 = {};
lmNerr_PV_S6 = {};
lmNerr_PV_S5v = {};
lmNerr_PV_S6v = {};

% L_thr = 50;
% L_IX = L >= L_thr;
L_IX = true(size(L));
for jdt = 1:ndt
	IX = L_IX' & Nerr_ML_S1(:,jdt)<max(N); lmNerr_ML_S1{jdt} = runFit(L(IX), Nerr_ML_S1(IX,jdt));
    IX = L_IX' & Nerr_ML_S2(:,jdt)<max(N); lmNerr_ML_S2{jdt} = runFit(L(IX), Nerr_ML_S2(IX,jdt));
    IX = L_IX' & Nerr_ML_S3(:,jdt)<max(N); lmNerr_ML_S3{jdt} = runFit(L(IX), Nerr_ML_S3(IX,jdt));
    IX = L_IX' & Nerr_ML_S4(:,jdt)<max(N); lmNerr_ML_S4{jdt} = runFit(L(IX), Nerr_ML_S4(IX,jdt));
    IX = L_IX' & Nerr_ML_S5(:,jdt)<max(N); lmNerr_ML_S5{jdt} = runFit(L(IX), Nerr_ML_S5(IX,jdt));
    IX = L_IX' & Nerr_ML_S6(:,jdt)<max(N); lmNerr_ML_S6{jdt} = runFit(L(IX), Nerr_ML_S6(IX,jdt));
    IX = L_IX' & Nerr_ML_S5v(:,jdt)<max(N); lmNerr_ML_S5v{jdt} = runFit(L(IX), Nerr_ML_S5v(IX,jdt));
    IX = L_IX' & Nerr_ML_S6v(:,jdt)<max(N); lmNerr_ML_S6v{jdt} = runFit(L(IX), Nerr_ML_S6v(IX,jdt));
    
	IX = L_IX' & Nerr_PV_S1(:,jdt)<max(N); lmNerr_PV_S1{jdt} = runFit(L(IX), Nerr_PV_S1(IX,jdt));
    IX = L_IX' & Nerr_PV_S2(:,jdt)<max(N); lmNerr_PV_S2{jdt} = runFit(L(IX), Nerr_PV_S2(IX,jdt));
    IX = L_IX' & Nerr_PV_S3(:,jdt)<max(N); lmNerr_PV_S3{jdt} = runFit(L(IX), Nerr_PV_S3(IX,jdt));
    IX = L_IX' & Nerr_PV_S4(:,jdt)<max(N); lmNerr_PV_S4{jdt} = runFit(L(IX), Nerr_PV_S4(IX,jdt));
    IX = L_IX' & Nerr_PV_S5(:,jdt)<max(N); lmNerr_PV_S5{jdt} = runFit(L(IX), Nerr_PV_S5(IX,jdt));
    IX = L_IX' & Nerr_PV_S6(:,jdt)<max(N); lmNerr_PV_S6{jdt} = runFit(L(IX), Nerr_PV_S6(IX,jdt));
    IX = L_IX' & Nerr_PV_S5v(:,jdt)<max(N); lmNerr_PV_S5v{jdt} = runFit(L(IX), Nerr_PV_S5v(IX,jdt));
    IX = L_IX' & Nerr_PV_S6v(:,jdt)<max(N); lmNerr_PV_S6v{jdt} = runFit(L(IX), Nerr_PV_S6v(IX,jdt));
end

%% calc slopes of minN required for 5% error vs. L
lmNrerr_ML_S1 = {};
lmNrerr_ML_S2 = {};
lmNrerr_ML_S3 = {};
lmNrerr_ML_S4 = {};
lmNrerr_ML_S5 = {};
lmNrerr_ML_S6 = {};
lmNrerr_ML_S5v = {};
lmNrerr_ML_S6v = {};

lmNrerr_PV_S1 = {};
lmNrerr_PV_S2 = {};
lmNrerr_PV_S3 = {};
lmNrerr_PV_S4 = {};
lmNrerr_PV_S5 = {};
lmNrerr_PV_S6 = {};
lmNrerr_PV_S5v = {};
lmNrerr_PV_S6v = {};

% L_thr = 50;
% L_IX = L >= L_thr;
L_IX = true(size(L));
for jdt = 1:ndt
	IX = L_IX' & Nrerr_ML_S1(:,jdt)<max(N); lmNrerr_ML_S1{jdt} = runFit(L(IX), Nrerr_ML_S1(IX,jdt));
    IX = L_IX' & Nrerr_ML_S2(:,jdt)<max(N); lmNrerr_ML_S2{jdt} = runFit(L(IX), Nrerr_ML_S2(IX,jdt));
    IX = L_IX' & Nrerr_ML_S3(:,jdt)<max(N); lmNrerr_ML_S3{jdt} = runFit(L(IX), Nrerr_ML_S3(IX,jdt));
    IX = L_IX' & Nrerr_ML_S4(:,jdt)<max(N); lmNrerr_ML_S4{jdt} = runFit(L(IX), Nrerr_ML_S4(IX,jdt));
    IX = L_IX' & Nrerr_ML_S5(:,jdt)<max(N); lmNrerr_ML_S5{jdt} = runFit(L(IX), Nrerr_ML_S5(IX,jdt));
    IX = L_IX' & Nrerr_ML_S6(:,jdt)<max(N); lmNrerr_ML_S6{jdt} = runFit(L(IX), Nrerr_ML_S6(IX,jdt));
    IX = L_IX' & Nrerr_ML_S5v(:,jdt)<max(N); lmNrerr_ML_S5v{jdt} = runFit(L(IX), Nrerr_ML_S5v(IX,jdt));
    IX = L_IX' & Nrerr_ML_S6v(:,jdt)<max(N); lmNrerr_ML_S6v{jdt} = runFit(L(IX), Nrerr_ML_S6v(IX,jdt));
    
	IX = L_IX' & Nrerr_PV_S1(:,jdt)<max(N); lmNrerr_PV_S1{jdt} = runFit(L(IX), Nrerr_PV_S1(IX,jdt));
    IX = L_IX' & Nrerr_PV_S2(:,jdt)<max(N); lmNrerr_PV_S2{jdt} = runFit(L(IX), Nrerr_PV_S2(IX,jdt));
    IX = L_IX' & Nrerr_PV_S3(:,jdt)<max(N); lmNrerr_PV_S3{jdt} = runFit(L(IX), Nrerr_PV_S3(IX,jdt));
    IX = L_IX' & Nrerr_PV_S4(:,jdt)<max(N); lmNrerr_PV_S4{jdt} = runFit(L(IX), Nrerr_PV_S4(IX,jdt));
    IX = L_IX' & Nrerr_PV_S5(:,jdt)<max(N); lmNrerr_PV_S5{jdt} = runFit(L(IX), Nrerr_PV_S5(IX,jdt));
    IX = L_IX' & Nrerr_PV_S6(:,jdt)<max(N); lmNrerr_PV_S6{jdt} = runFit(L(IX), Nrerr_PV_S6(IX,jdt));
    IX = L_IX' & Nrerr_PV_S5v(:,jdt)<max(N); lmNrerr_PV_S5v{jdt} = runFit(L(IX), Nrerr_PV_S5v(IX,jdt));
    IX = L_IX' & Nrerr_PV_S6v(:,jdt)<max(N); lmNrerr_PV_S6v{jdt} = runFit(L(IX), Nrerr_PV_S6v(IX,jdt));
end


%% extrap (using the linear fit) minimal neurons for 2m error
for jdt = 1:ndt
    IX = Nerr_ML_S1(:,jdt)>=max(N); Nerr_ML_S1(IX,jdt) = predict(lmNerr_ML_S1{jdt}, L(IX)');
    IX = Nerr_ML_S2(:,jdt)>=max(N); Nerr_ML_S2(IX,jdt) = predict(lmNerr_ML_S2{jdt}, L(IX)');
    IX = Nerr_ML_S3(:,jdt)>=max(N); Nerr_ML_S3(IX,jdt) = predict(lmNerr_ML_S3{jdt}, L(IX)');
    IX = Nerr_ML_S4(:,jdt)>=max(N); Nerr_ML_S4(IX,jdt) = predict(lmNerr_ML_S4{jdt}, L(IX)');
    IX = Nerr_ML_S5(:,jdt)>=max(N); Nerr_ML_S5(IX,jdt) = predict(lmNerr_ML_S5{jdt}, L(IX)');
    IX = Nerr_ML_S6(:,jdt)>=max(N); Nerr_ML_S6(IX,jdt) = predict(lmNerr_ML_S6{jdt}, L(IX)');
    IX = Nerr_ML_S5v(:,jdt)>=max(N); Nerr_ML_S5v(IX,jdt) = predict(lmNerr_ML_S5v{jdt}, L(IX)');
    IX = Nerr_ML_S6v(:,jdt)>=max(N); Nerr_ML_S6v(IX,jdt) = predict(lmNerr_ML_S6v{jdt}, L(IX)');
    
    IX = Nerr_PV_S1(:,jdt)>=max(N); Nerr_PV_S1(IX,jdt) = predict(lmNerr_PV_S1{jdt}, L(IX)');
    IX = Nerr_PV_S2(:,jdt)>=max(N); Nerr_PV_S2(IX,jdt) = predict(lmNerr_PV_S2{jdt}, L(IX)');
    IX = Nerr_PV_S3(:,jdt)>=max(N); Nerr_PV_S3(IX,jdt) = predict(lmNerr_PV_S3{jdt}, L(IX)');
    IX = Nerr_PV_S4(:,jdt)>=max(N); Nerr_PV_S4(IX,jdt) = predict(lmNerr_PV_S4{jdt}, L(IX)');
    IX = Nerr_PV_S5(:,jdt)>=max(N); Nerr_PV_S5(IX,jdt) = predict(lmNerr_PV_S5{jdt}, L(IX)');
    IX = Nerr_PV_S6(:,jdt)>=max(N); Nerr_PV_S6(IX,jdt) = predict(lmNerr_PV_S6{jdt}, L(IX)');
    IX = Nerr_PV_S5v(:,jdt)>=max(N); Nerr_PV_S5v(IX,jdt) = predict(lmNerr_PV_S5v{jdt}, L(IX)');
    IX = Nerr_PV_S6v(:,jdt)>=max(N); Nerr_PV_S6v(IX,jdt) = predict(lmNerr_PV_S6v{jdt}, L(IX)');
end


%% extrap (using the linear fit) minimal neurons for 5% error
for jdt = 1:ndt
    IX = Nrerr_ML_S1(:,jdt)>=max(N); Nrerr_ML_S1(IX,jdt) = predict(lmNrerr_ML_S1{jdt}, L(IX)');
    IX = Nrerr_ML_S2(:,jdt)>=max(N); Nrerr_ML_S2(IX,jdt) = predict(lmNrerr_ML_S2{jdt}, L(IX)');
    IX = Nrerr_ML_S3(:,jdt)>=max(N); Nrerr_ML_S3(IX,jdt) = predict(lmNrerr_ML_S3{jdt}, L(IX)');
    IX = Nrerr_ML_S4(:,jdt)>=max(N); Nrerr_ML_S4(IX,jdt) = predict(lmNrerr_ML_S4{jdt}, L(IX)');
    IX = Nrerr_ML_S5(:,jdt)>=max(N); Nrerr_ML_S5(IX,jdt) = predict(lmNrerr_ML_S5{jdt}, L(IX)');
    IX = Nrerr_ML_S6(:,jdt)>=max(N); Nrerr_ML_S6(IX,jdt) = predict(lmNrerr_ML_S6{jdt}, L(IX)');
    IX = Nrerr_ML_S5v(:,jdt)>=max(N); Nrerr_ML_S5v(IX,jdt) = predict(lmNrerr_ML_S5v{jdt}, L(IX)');
    IX = Nrerr_ML_S6v(:,jdt)>=max(N); Nrerr_ML_S6v(IX,jdt) = predict(lmNrerr_ML_S6v{jdt}, L(IX)');
    
    IX = Nrerr_PV_S1(:,jdt)>=max(N); Nrerr_PV_S1(IX,jdt) = predict(lmNrerr_PV_S1{jdt}, L(IX)');
    IX = Nrerr_PV_S2(:,jdt)>=max(N); Nrerr_PV_S2(IX,jdt) = predict(lmNrerr_PV_S2{jdt}, L(IX)');
    IX = Nrerr_PV_S3(:,jdt)>=max(N); Nrerr_PV_S3(IX,jdt) = predict(lmNrerr_PV_S3{jdt}, L(IX)');
    IX = Nrerr_PV_S4(:,jdt)>=max(N); Nrerr_PV_S4(IX,jdt) = predict(lmNrerr_PV_S4{jdt}, L(IX)');
    IX = Nrerr_PV_S5(:,jdt)>=max(N); Nrerr_PV_S5(IX,jdt) = predict(lmNrerr_PV_S5{jdt}, L(IX)');
    IX = Nrerr_PV_S6(:,jdt)>=max(N); Nrerr_PV_S6(IX,jdt) = predict(lmNrerr_PV_S6{jdt}, L(IX)');
    IX = Nrerr_PV_S5v(:,jdt)>=max(N); Nrerr_PV_S5v(IX,jdt) = predict(lmNrerr_PV_S5v{jdt}, L(IX)');
    IX = Nrerr_PV_S6v(:,jdt)>=max(N); Nrerr_PV_S6v(IX,jdt) = predict(lmNrerr_PV_S6v{jdt}, L(IX)');
end


%% 
function mdl = runFit(x,y)
switch 1
    case 1
        mdl = fitlm(x,y);
    case 2
        mdl = fitlm(x,y,'Intercept',false);
    case 3
        % fit only using the last 5 points
        x = x(end-4:end);
        y = y(end-4:end);
        mdl = fitlm(x,y);
end
end


%%




%%

