function [f] = MultiscalePlace_GenerateTuning_AllModels_Rev1(L)
alp     = 0.3 ; % alpha, scaling exponent of environment size

ds      = 20 ; % downsampling factor
L       = L/ds ;
L0_S6   = 20000/ds  ; % anchoring environment size for scheme 6
L0_S4   = 5000/ds ;   % anchoring environment size for scheme 4
Nmax    = 1000      ;
Ninvgma = 100  ;
gma_k   = 3.1647    ; % gamma shape parameter
gma_th0 = 179.23/ds ; % gamma scale parameter
phi0    = 0.15 ;
phimax   = 0.5 ;
phi      = min(phi0*(L0_S6./L)^alp,phimax) ;

gmaAL_th0 = 100/7.75/ds ;
gmaAL_k   = 0.57 ;

l0_S1    = 100/ds ;
l0_S2    = round(phi0*L0_S6*(L/L0_S6)^alp) ;
li_S3    = sort(repmat(round(linspace(l0_S1,l0_S2,Ninvgma)),[1 Nmax/Ninvgma])) ;
l0_S4    = 100/ds ;
li_S5    = sort(repmat(round(gaminv(((Ninvgma-0.5):-1:0.5)/Ninvgma,gma_k,gma_th0)*(L/L0_S6)^alp),[1 Nmax/Ninvgma]),'Descend') ; li_S5(li_S5==0) = 1 ;

f_S1  = zeros(L,Nmax) ;
f_S2  = zeros(L,Nmax) ;
f_S3  = zeros(L,Nmax) ;
f_S4  = zeros(L,Nmax) ;
f_S5  = zeros(L,Nmax) ;
f_S6  = zeros(L,Nmax) ;
f_S5v  = zeros(L,Nmax) ;
f_S6v  = zeros(L,Nmax) ;

for i = 1:Nmax
    f_S1(randi(L-l0_S1)+(1:l0_S1),i) = 1 ;
end

for i = 1:Nmax
    i0 = min(randi(L+l0_S2,1)-l0_S2,L) ;
    f_S2(max(i0,1):1:min(i0+l0_S2-1,L),i) = 1 ;
end

for i = 1:Nmax
    i0 = min(randi(L+li_S3(i),1)-li_S3(i),L) ;
    f_S3(max(i0,1):1:min(i0+li_S3(i)-1,L),i) = 1 ;
end

for i = 1:Nmax
    scf = 0 ;
    while ~scf
        ri  = random('Gamma',gmaAL_k,gmaAL_th0*(L0_S4/L)^alp,1)  ;
        f_S4i = zeros(1,L) ;
        i0  = random('Exp',1/ri*(100/ds)^2,1,1000) ;
        i0  = ceil(cumsum(i0)) ;
        i0(i0>L) = [] ;
        f_S4i(i0)  = 1 ;
        if sum(f_S4i)>0
            f_S4i = conv(f_S4i,ones(1,l0_S4),'same') ;
            f_S4i(f_S4i>1) = 1 ;
            f_S4(:,i) = circshift(f_S4i,randi(L)) ;
            scf = 1 ;
        end
    end
end

i = 1 ;
while i <= Nmax
    Ki = max(1,round(phi*L/li_S5(i))) ;
    i0  = sort(randi(L-min(L-1,li_S5(i)),Ki,1),'Ascend') ;
    if all(i0<=(L-min(L-1,li_S5(i))))
        for k = 1:Ki
            f_S5(i0(k)+(1:min(L-1,li_S5(i))),i) = 1 ;
        end
        f_S5(:,i) = circshift(f_S5(:,i),randi(L)) ;
        i = i+1 ;
    end
    
end

i = 1 ;
while i <= Nmax
    lii = li_S5(randperm(Nmax)) ;
    [~,Ki] = min(abs(cumsum(lii)/L-phi)) ;
    lii = lii(1:Ki) ;
    i0 = sort(randi(L-min(L-1,lii(Ki)),1,Ki)) ;
    di0 = diff(i0) ;
    if all(lii(1:(Ki-1))<di0)
        for k = 1:Ki
            f_S6(i0(k)+(0:1:min(L-i0(k),lii(k)-1)),i) = 1 ;
        end
        f_S6(:,i) = circshift(f_S6(:,i),randi(L)) ;
        i = i+1 ;
    end
end

i = 1 ;
phii = min(random('Exp',phi,1,Nmax),phimax) ;
while i <= Nmax
    Ki = max(1,round(phii(i)*L/li_S5(i))) ;
    i0  = sort(randi(L-min(L-1,li_S5(i)),Ki,1),'Ascend') ;
    if all(i0<=(L-min(L-1,li_S5(i))))
        for k = 1:Ki
            f_S5v(i0(k)+(1:min(L-1,li_S5(i))),i) = 1 ;
        end
        f_S5v(:,i) = circshift(f_S5v(:,i),randi(L)) ;
        i = i+1 ;
    end
end

i = 1 ;
ii = 0 ;
while i <= Nmax
    ii = ii + 1 ;
    if ii > 10000
        phii(i) = min(random('Exp',phi),phimax) ;
    end
    lii = li_S5(randperm(Nmax)) ;
    [~,Ki] = min(abs(cumsum(lii)/L-phii(i))) ;
    lii = lii(1:Ki) ;
    i0 = sort(randi(L-min(L-1,lii(Ki)),1,Ki)) ;
    di0 = diff(i0) ;
    if all(lii(1:(Ki-1))<di0)
        for k = 1:Ki
            f_S6v(i0(k)+(0:1:min(L-i0(k),lii(k)-1)),i) = 1 ;
        end
        f_S6v(:,i) = circshift(f_S6v(:,i),randi(L)) ;
        i = i+1 ;
        ii = 0 ;
    end
end

%% arrange output
f(1).maps = f_S1;
f(2).maps = f_S2;
f(3).maps = f_S3;
f(4).maps = f_S4;
f(5).maps = f_S5;
f(6).maps = f_S6;
f(7).maps = f_S5v;
f(8).maps = f_S6v;

f(1).name_short = '1';
f(2).name_short = '2';
f(3).name_short = '3';
f(4).name_short = '4';
f(5).name_short = '5';
f(6).name_short = '6';
f(7).name_short = '5v';
f(8).name_short = '6v';

f(1).name_long = '1: Single small field';
f(2).name_long = '2: Single large field';
f(3).name_long = '3: Single field D-V';
f(4).name_long = '4: Multiple small fields (Rich et al. 2014)';
f(5).name_long = '5: Multi-scale (population)';
f(6).name_long = '6: Multi-scale (single-cell)';
f(7).name_long = '5v: Multi-scale (population) - variable coverage';
f(8).name_long = '6v: Multi-scale (single-cell) - variable coverage';



