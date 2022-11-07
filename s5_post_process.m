% This programme aims to export the toolpath file, for the continuing post
% processing procedure in the IMSpost software to get the nc file.
isAPP = true;
if isAPP
else
    close all;
    addpath(genpath('funcs'));
end




%%
% rmpath(genpath('funcs'));