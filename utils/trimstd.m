function s = trimstd(x,p)

p = [0 100] +[1 -1].*p/2;
range = prctile(x,p);
x(x<range(1)) = [];
x(x>range(2)) = [];
s = std(x);

end