function util_fix_ncs(dir_IN)

%% create output folder
dir_OUT = [dir_IN '_NCS_fixed_header'];
mkdir(dir_OUT)

%% get files
files = dir(dir_IN);
files([files.isdir])=[];

[~,~,EXT] = arrayfun(@(name)(fileparts(name)), ({files.name}),'UniformOutput',0);
is_ncs = contains(EXT,'.ncs');
ncs_files = files(is_ncs);
other_files = files(~is_ncs);

%% copy other files
for ii_file = 1:length(other_files)
    file = other_files(ii_file);
    src_file = fullfile(file.folder, file.name);
    dst_file = fullfile(dir_OUT, file.name);
    copyfile( src_file, dst_file )
end

%% fix ncs files
for ii_file = 1:length(ncs_files)
    file = ncs_files(ii_file);
    file_IN = fullfile(file.folder, file.name);
    file_OUT = fullfile(dir_OUT, file.name);
    nlx_change_header(file_IN, file_OUT);
end



end