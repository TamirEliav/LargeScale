function sem = nansem(x)

sem = nanstd(x) ./ sqrt(sum(~isnan(x)));

end