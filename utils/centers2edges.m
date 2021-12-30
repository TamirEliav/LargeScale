function [e,w] = centers2edges(c)
c = c(:)'; % make row vector
e = (c(1:end-1)+c(2:end))/2;
w = median(diff(c));
e = [e(1)-w e e(end)+w];
end
