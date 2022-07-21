function savetofile(data,fullfilename,variblename,type)
eval([variblename '=data;']);
if ~type
    save(fullfilename,variblename);
else
    save(fullfilename,variblename,'-append');
end
end