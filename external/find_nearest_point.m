function ib = find_nearest_point(a,b)
% this function finds for each element in a the closest element in b
% ib is a vector holding indeces of the closest elements in b (has same length as a)

a = shiftdim(a)';
b = shiftdim(b)';

m = length(a);
n = length(b);
[c,p] = sort([a,b]);
q = 1:m+n;
q(p) = q;
t = cumsum(p>m);
r = 1:n;
r(t(q(m+1:m+n))) = r;
s = t(q(1:m));
id = r(max(s,1));
iu = r(min(s+1,n));
[d,it] = min([abs(a-b(id));abs(b(iu)-a)]);
ib = id+(it-1).*(iu-id);

end