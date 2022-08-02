function savetofile(data,fullfilename,variblename,type)
% Save data inside a parfor loop
eval([variblename '=data;']);
if ~type
    save(fullfilename,variblename);
else
    save(fullfilename,variblename,'-append');
end
end