function a = soa2aos(s)
% convert struct of arrays to an array of structs

%%
fn=fieldnames(s);
n=numel(s.(fn{1}));

%%
a=repelem(struct(),n);
for ii = 1:length(fn)
    name = fn{ii};
    [a(:).(name)] = disperse(s.(name));
end

end