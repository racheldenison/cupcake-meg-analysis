function tag = meg_getTag(str, fileBase)
%
% function tag = meg_getTag(str, delimiter)
%
% the tag is the run string that comes after filebase and before the
% dot in the file extension

if nargin < 2
    delimiter = '_';
end

fileBaseIdx = strfind(str, fileBase) + length(fileBase);
dotIdx = strfind(str, '.');
tag = str(fileBaseIdx(end)+1:dotIdx(end)-1);