function myfs_addpath
[mypath,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(mypath,'external')))
addpath(fullfile(mypath,'utilities'))
end
