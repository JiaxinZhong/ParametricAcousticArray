%% add the subfolders to the path
addpath(sprintf('%s/absorp',pwd))
addpath(sprintf('%s/absorp/fig',pwd))
addpath(sprintf('%s/direct',pwd))
addpath(sprintf('%s/direct/fig',pwd))
addpath(sprintf('%s/kzk',pwd))
addpath(sprintf('%s/kzk/fig',pwd))
addpath(sprintf('%s/misc',pwd))
addpath(sprintf('%s/efficiency',pwd));
addpath(sprintf('%s/nah',pwd))

%% figure settings
fig_type = 'normal'; 
switch fig_type
    case 'normal' % normal type
        FONT_SIZE = 16;
        MARKER_SIZE = 10;
        LINE_WIDTH = 2;
    case 'subfig' % subfig. which means the fontsize should be larger
        FONT_SIZE = 24;
        MARKER_SIZE = 16;
        LINE_WIDTH = 3;
end
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');
set(groot, 'defaultLineLineWidth', LINE_WIDTH);
set(groot, 'defaultImplicitFunctionLineLineWidth', LINE_WIDTH);
% set(groot, 'DefaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultLegendLocation', 'best');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', FONT_SIZE);
set(groot, 'defaultPolarAxesFontSize', FONT_SIZE);
set(groot, 'defaultFigureColormap', jet)
% set(groot, 'defaultColorbarTextInterpreter', 'latex')
set(groot, 'defaultLineMarkerSize', MARKER_SIZE)
    
%% Color scheme
car = [.85, .33, .09];      ca1 = car;
cag = [.11, .83, 0];        ca2 = cag;
cab = [0, .45, .74];        ca3 = cab;
cap = [.5, 0, .5];          ca4 = cap;
cay = [.93, .70, .15];      ca5 = cay;
cao = [255,127,0]/255;      ca6 = cao;
cagr = [153,153,153]/255;   ca7 = cagr;
capink = [247,129,191]/255; ca8 = capink;
cat = [0, .5, .5];          ca9 = cat;
cacell = cell(9,1);
cacell{1} = ca1; cacell{2} = ca2; cacell{3} = ca3; cacell{4} = ca4;
cacell{5} = ca5; cacell{6} = ca6; cacell{7} = ca7; cacell{8} = ca8;
cacell{9} = ca9;


%% display settings
format('long')
