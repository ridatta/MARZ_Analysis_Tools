function LDLSshapeoct2015 = importLDLS(filename, dataLines)
%IMPORTFILE7 Import data from a text file
%  LDLSSHAPEOCT2015 = IMPORTFILE7(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  LDLSSHAPEOCT2015 = IMPORTFILE7(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  LDLSshapeoct2015 = importfile7("C:\Users\rdatta\Dropbox (MIT)\PUFFIN\Codes\MARZ\SVS\LDLS_shape_oct_2015.dat", [7, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 08-Sep-2022 14:50:50

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [7, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["ColumnDefinitions", "VarName2", "Var3", "Var4"];
opts.SelectedVariableNames = ["ColumnDefinitions", "VarName2"];
opts.VariableTypes = ["double", "double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var3", "Var4"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var3", "Var4"], "EmptyFieldRule", "auto");

% Import the data
LDLSshapeoct2015 = readtable(filename, opts);

end