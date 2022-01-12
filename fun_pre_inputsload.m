function [inputs] = fun_pre_inputsload(modelname)

fid = fopen(strcat(modelname,'.txt'));

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end % if end the file, break the while loop
    
    inputs.modelname = modelname;
    
    if tline(1:5) == 'penal'
        inputs.penal = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:5) == 'rmin '
        inputs.rmin = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:6) == 'initdv'
        inputs.initdv = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:2) == 'VT'
        inputs.VT = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:7) == 'volfrac'
        inputs.volfrac = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:5) == 'mprop'
        inputs.mprop = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:5) == 'force'
        inputs.force = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:4) == 'conv'
        inputs.conv = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:7) == 'bt_init'
        inputs.bt = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:5) == 'bt_ic'
        inputs.bt_ic = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:5) == 'bt_ns'
        inputs.bt_ns = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:5) == 'bt_fn'
        inputs.bt_fn = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:5) == 'MMA_c'
        inputs.MMA_c = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:5) == 'resol'
        inputs.resol = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:4) == 'xmir'
        inputs.xmir = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:4) == 'ymir'
        inputs.ymir = str2num(tline(strfind(tline,':')+1:end));
    elseif tline(1:4) == 'zmir'
        inputs.zmir = str2num(tline(strfind(tline,':')+1:end));        
    end
    
end

if size( struct2table(inputs),  2) ~= 18
    error('check inputs text file')
end

end