% REMRELPATH  Remove a relative path (which gets expanded into a absolute path)
%             from MATLABPATH.
%
%     remrelpath(path_string)
%
%     INPUT
%            path_string        : string containing relative path
%
%
%     Example : remrelpath('../foo');
%
%     This will first expand '../foo' to an absolute path and then remove it from
%     the current MATLABPATH
%
%   See also
%   ADDRELPATH
%
%   Copyright (c) Oregon Health & Science University (2006)
%
%   This file is part of the ReBEL Toolkit. The ReBEL Toolkit is available free for
%   academic use only (see included license file) and can be obtained from
%   http://choosh.csee.ogi.edu/rebel/.  Businesses wishing to obtain a copy of the
%   software should contact rebel@csee.ogi.edu for commercial licensing information.
%
%   See LICENSE (which should be part of the main toolkit distribution) for more
%   detail.

%=============================================================================================
function remrelpath(path_string)

if strncmp(computer,'PC',2),
   path_string = strrep(path_string,'/','\'); % Convert '/' to '\'
else
   path_string = strrep(path_string,'\','/'); % Convert '\' to '/'
end

current_dir = pwd;
cd(path_string);
absolute_path = pwd;
cd(current_dir);
rmpath(absolute_path);

