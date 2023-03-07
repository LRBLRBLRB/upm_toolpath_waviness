function axesPos = importNCFile(filename, dataLines)
% IMPORTFILE 从文本文件中导入数据
%  LENSARRAY210721RT0 = IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。
%  返回数值数据。
%
%  LENSARRAY210721RT0 = IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  Lensarray210721RT0 = importfile("D:\Code\2021-11_ToolWaviness\upm_toolpath_waviness\tool_path\Lens array_210721_RT0.2_tilt0_rot0_SF_s6a12d0.66.nc", [1, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2022-06-21 14:57:39 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [1, Inf];
end

%% 设置导入选项并导入数据
opts = delimitedTextImportOptions("NumVariables", 5);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = " ";

% 指定列名称和类型
opts.VariableNames = ["B", "C", "X", "Y", "Z"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% 指定变量属性
opts = setvaropts(opts, opts.VariableNames, "TrimNonNumeric", true);
opts = setvaropts(opts, opts.VariableNames, "ThousandsSeparator", ",");

opts.MissingRule = "fill";
% opts = setvaropts(opts, opts.VariableNames, "TreatAsMissing", "nan", "FillValue", z0);

% 导入数据
axesPos = readtable(filename, opts);

%% 转换为输出类型
axesPos = table2array(axesPos);
end