function util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list, ext)

%% parse inout
if nargin<4
    ext = '';
end

%% get input files
files = dir(dir_IN);
files(~contains({files.name},tmpl_list)) = [];
if ~isempty(ext)
    [~,~,extensions] = cellfun(@(x)(fileparts(x)),{files.name},'UniformOutput',false);
    files(~contains(extensions,ext)) = [];
end
disp('files to copy:')
{files.name}'

%% copy files
mkdir(dir_OUT);
for ii_file = 1:length(files)
    file_IN =  fullfile(files(ii_file).folder, files(ii_file).name);
    file_OUT = fullfile(dir_OUT, files(ii_file).name);
    [status,msg,msgID] = copyfile(file_IN,file_OUT);
    if ~status
        disp('problem with copying file:');
        disp(file_IN)
        disp(msgID)
        disp(msg)
    end
end






end