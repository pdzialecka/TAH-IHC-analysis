%% Find files to be renamed
% @author: pdzialecka

base_folder = 'D:\TAH';
% f_str = '12F4';
f_str = 'GFAP';
c_str = lower(f_str);

files = dir(fullfile(base_folder,'Data','**',strcat('*',f_str,'*')));
keep_idxs = find(contains({files.name}',caseSensitivePattern(f_str)));

files_to_rename = files(keep_idxs);

%% Rename files
for idx = 1:length(files_to_rename)
    file = files_to_rename(idx);
    fname = fullfile(file.folder,file.name);
    
    new_file_name = file.name;
    str_idxs = strfind(file.name,f_str):(strfind(file.name,f_str)+length(f_str)-1);
    new_file_name(str_idxs) = c_str;
    
    new_fname = fullfile(file.folder,new_file_name);
    tmp_new_fname = strcat(new_fname(1:end-4),'_tmp',new_fname(end-3:end));

    movefile(fname,tmp_new_fname);
    movefile(tmp_new_fname,new_fname);
end
