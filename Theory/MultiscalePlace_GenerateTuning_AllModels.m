function [f] = MultiscalePlace_GenerateTuning_AllModels(L,coverage)

ds = 20 ;
L = L/ds ;
L0      = 20000/ds  ;
L0H     = 5000/ds   ;
L0I     = 5000/ds   ;
Nmax    = 1000      ;
Ninvgma = 100  ;
N       = [10:10:150 170:20:250] ;
nN      = length(N) ;
gma_th0 = 316.47/ds    ;
gma_k   = 1.7923    ;

gmaAL_th0 = 100/7.75/ds ;
gmaAL_k   = 0.57 ;

l0A    = 100/ds ;
l0B    = round(coverage*L0*sqrt(L/L0)) ;
l0C    = 100/ds ;
l0H    = 100/ds ;
l0I = 100/ds ;
liD    = sort(repmat(round(gaminv(((Ninvgma-0.5):-1:0.5)/Ninvgma,gma_k,gma_th0)),[1 Nmax/Ninvgma]),'Descend') ;
liF    = sort(repmat(round(gaminv(((Ninvgma-0.5):-1:0.5)/Ninvgma,gma_k,gma_th0)*sqrt(L/L0)),[1 Nmax/Ninvgma]),'Descend') ; liF(liF==0) = 1 ;
phi    = coverage*sqrt(L0./L)  ;
KC     = round(phi./(l0C./L)) ;
liH    = l0H * sqrt(L/L0H) ;
KH     = round(phi./(liH/L)) ;

fA  = zeros(L,Nmax) ;
fB  = zeros(L,Nmax) ;
fC  = zeros(L,Nmax) ;
fD  = zeros(L,Nmax) ;
fE  = zeros(L,Nmax) ;
fF  = zeros(L,Nmax) ;
fG  = zeros(L,Nmax) ;
fH  = zeros(L,Nmax) ;
fI  = zeros(L,Nmax) ;
fJ  = zeros(L,Nmax) ;

for i = 1:Nmax
    fA(randi(L-l0A)+(1:l0A),i) = 1 ;
end

for i = 1:Nmax
    i0 = min(randi(L+l0B,1)-l0B,L) ;
    fB(max(i0,1):1:min(i0+l0B-1,L),i) = 1 ;
end

i = 1 ;
while i <= Nmax
    di0 = round(random('Exp',max(1,(1-phi)*L/KC-l0C),1,KC+1)+(l0C+1)) ;

    di0 = round(((L-l0C)*di0- l0C)/(sum(di0)) ) ;
    i0  = cumsum(di0(1:KC)) ;
    if all(i0<=(L-l0C))
        for k = 1:KC
            fC(i0(k)+(1:l0C),i) = 1 ;
        end
        fC(:,i) = circshift(fC(:,i),randi(L)) ;
        i = i+1 ;
    end
end

i = 1 ;
while i <= Nmax
    Ki = max(1,round(phi*L/liD(i))) ;
    if Ki < 4000 % 250
        i0  = sort(randi(L-min(L-1,liD(i)),Ki,1),'Ascend') ;
    else
        di0 = round(random('Exp',(1-phi)*L/Ki-liD(i),1,Ki+1)+(liD(i)+1)) ;
        di0 = round((L-liD(i))*di0/(sum(di0))) ;
        i0  = cumsum(di0(1:Ki)) ;
    end
    if all(i0<=(L-min(L-1,liD(i))))
        for k = 1:Ki
            fD(i0(k)+(1:min(L-1,liD(i))),i) = 1 ;
        end
        i = i+1 ;
    end
end

i = 1 ;
while i <= Nmax
    lii = liD(randperm(Nmax)) ;
    [~,Ki] = min(abs(cumsum(lii)/L-phi)) ;
    lii = lii(1:Ki) ;
    i0 = sort(randi(L-min(L-1,lii(Ki)),1,Ki)) ;
    di0 = diff(i0) ;
    if all(lii(1:(Ki-1))<di0)
        for k = 1:Ki
            fE(i0(k)+(0:1:min(L-i0(k),lii(k)-1)),i) = 1 ;
        end
        fE(:,i) = circshift(fE(:,i),randi(L)) ;
        i = i+1 ;
    end
end

i = 1 ;
while i <= Nmax
    Ki = max(1,round(phi*L/liF(i))) ;
    if Ki < 250
        i0  = sort(randi(L-min(L-1,liF(i)),Ki,1),'Ascend') ;
    else
        di0 = round(random('Exp',max(1,(1-phi)*L/Ki-liF(i)),1,Ki+1)+(liF(i)+1)) ;
        di0 = round((L-liF(i))*di0/(sum(di0))) ;
        i0  = cumsum(di0(1:Ki)) ;
    end
    if all(i0<=(L-min(L-1,liF(i))))
        for k = 1:Ki
            fF(i0(k)+(1:min(L-1,liF(i))),i) = 1 ;
        end
        i = i+1 ;
    end
end

i = 1 ;
while i <= Nmax
    lii = liF(randperm(Nmax)) ;
    [~,Ki] = min(abs(cumsum(lii)/L-phi)) ;
    lii = lii(1:Ki) ;
    i0 = sort(randi(L-min(L-1,lii(Ki)),1,Ki)) ;
    di0 = diff(i0) ;
    if all(lii(1:(Ki-1))<di0)
        for k = 1:Ki
            fG(i0(k)+(0:1:min(L-i0(k),lii(k)-1)),i) = 1 ;
        end
        fG(:,i) = circshift(fG(:,i),randi(L)) ;
        i = i+1 ;
    end
end

i = 1 ;
while i <= Nmax
    di0 = round(random('Exp',max(1,(1-phi)*L/KH-liH),1,KH+1)+(liH+1)) ;

    di0 = round(((L-liH)*di0- liH)/(sum(di0)) ) ;
    i0  = cumsum(di0(1:KH)) ;
    if all(i0<=(L-liH))
        for k = 1:KH
            fH(i0(k)+(1:liH),i) = 1 ;
        end
        fH(:,i) = circshift(fH(:,i),randi(L)) ;
        i = i+1 ;
    end
end

for i = 1:Nmax
    scf = 0 ;
    while ~scf
        ri  = random('Gamma',gmaAL_k,gmaAL_th0*sqrt(L0I/L),1)  ;
        fIi = zeros(1,L) ;
        i0  = random('Exp',1/ri*(100/ds)^2,1,1000) ;
        i0  = ceil(cumsum(i0)) ;
        i0(i0>L) = [] ;
        fIi(i0)  = 1 ;
        if sum(fIi)>0
            fIi = conv(fIi,ones(1,l0I),'same') ;
            fIi(fIi>1) = 1 ;
            fI(:,i) = circshift(fIi,randi(L)) ;
            scf = 1 ;
        end
    end
end

fJ = zeros(L,Nmax) ;
for i = 1:Nmax
    scf = 0 ;
    while ~scf
        ri  = random('Gamma',gmaAL_k,gmaAL_th0*(L0I/L),1)  ;
        fJi = zeros(1,L) ;
        i0  = random('Exp',1/ri*(100/ds)^2,1,1000) ;
        i0  = ceil(cumsum(i0)) ;
        i0(i0>L) = [] ;
        fJi(i0)  = 1 ;
        if sum(fJi)>0
            fJi = conv(fJi,ones(1,round(l0I*sqrt(L/L0I))),'same') ;
            fJi(fJi>1) = 1 ;
            fJ(:,i) = circshift(fJi,randi(L)) ;
            scf = 1 ;
        end
    end
end

f.fA = fA ;
f.fB = fB ;
f.fC = fC ;
f.fD = fD ;
f.fE = fE ;
f.fF = fF ;
f.fG = fG ;
f.fH = fH ;
f.fI = fI ;
f.fJ = fJ ;
