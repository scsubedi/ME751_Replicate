%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads the input file from input.txt and saves them as variables in
% workspace

clear;clc;
close all;


% Reading the input file by avoiding white spaces and delimiters
syms t
filename  ='input';
str = fileread([filename '.txt']);
val = regexp(str,'^%\s*(.+?)\n(.+?)$','tokens','lineanchors');
val = strtrim(vertcat(val{:}));
num = cellfun(@(s)sscanf(s,'%f:'),val(:,2),'UniformOutput',false);
len = cellfun('length',num);
vec = cellfun(@num2cell,num(len>1),'UniformOutput',false);
num(len>1) = cellfun(@(v)colon(v{:}),vec,'UniformOutput',false);

i = str2num(val{3,2})';
j = str2num(val{4,2})';
L = str2num(val{5,2})';
s1_p = L*str2num(val{6,2})';
s2_q = L*str2num(val{7,2})';

ts = str2num(val{8,2});


joint = val{10,2};
saveResults = val{11,2};
type = val{12,2};

k = str2num(val{13,2});
c0 = str2num(val{14,2});
l0 = str2num(val{15,2});

dynamicAnalysis = val{16,2};

rho = str2num(val{17,2});
width = str2num(val{18,2});
Vol = 2*L*width*width;
len = 2*L;
M = rho*Vol;
J = 1/12*M*(width^2 + width^2);
Ix = 1/12*M * (len^2 + width^2);
Iy = 1/12*M * (len^2 + width^2);
Iz = 1/12*M * (2*width^2);
torque = zeros(length(ts),1);
