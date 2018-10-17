function figAddMarker(numMarker,runExample)
% =========================================================================
% FUNCTION 
%	- Add the markers to all the lines in gcf (i.e. current figure)
%	- the maximal number of the lines is restricted to EIGHT
% -------------------------------------------------------------------------
% TEST VERSIONS
%   Sucessful in:
%       - MATLAB R2018a macOS
% =========================================================================


%% Example
if ~exist('runExample','var')
    runExample = 0;
end
if runExample
    x = 1:300;
    y1 = sin(x*2*pi/200);
    y2 = cos(x*2*pi/200);
    y3 = x/300;
    figure
    plot(x,y1,'-o');
    hold on
    plot(x,y2,'-x');
    plot(x,y3,'-h');
    xlabel('$x$');
    ylabel('y');
    legend({'Line 1', 'Line 2', 'Line 3'})
end

%% Color scheme recommended by the author
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

%% Read data from the current figure
hLine = findobj(gca, 'type', 'line');
hLine = flipud(hLine); % the hangle of lines are reversely read
numLine = length(hLine); % number of the lines
x = cell(numLine,1);
y = cell(numLine,1);
for i = 1:numLine
	x{i} = get(hLine(i), 'xdata');
	y{i} = get(hLine(i), 'ydata');
end

%% Get the lengend from the current figure
hLegend = findobj(gcf, 'type', 'legend');
strLegend = get(hLegend, 'string');
locLegend = get(hLegend, 'location');

%% Get the labels from the current figure
xLabel = get(get(gca, 'xlabel'),'string');
yLabel = get(get(gca, 'ylabel'),'string');

%% Get the limits of the axe
xLim = get(gca, 'xlim');
yLim = get(gca, 'ylim');

%% close gcf
close(gcf);

%% Make the new figure
figure
% the number of the markers in a line
if ~exist('numMarker','var')
    numMarker = 8;
end
% the stepsize of the marker
mkrStep = zeros(numLine, 1);
for i = 1:numLine
	mkrStep(i) = ceil(length(y{i})/numMarker);
end
% the start point of the marker
mkrStart = ceil((1:numLine)'/numLine.*mkrStep);
mkrIndex = cell(numLine,1);
for i = 1:numLine
    mkrIndex{i} = mkrStart(i):mkrStep(i):length(y{i});
end

mkrType0 = {'x','*','o','^','d','o','^','d'};
mkrType = strcat('-',mkrType0);
mkrFacecolor = {'none', 'none', 'none', 'none', cacell{5:numLine}};
for i = 1:numLine
    plot(x{i}(mkrStart(i)), y{i}(mkrStart(i)), mkrType{i}, ...
        'color', cacell{i}, 'MarkerFaceColor', mkrFacecolor{i});
    hold on
end
for i = 1:numLine
    plot(x{i}, y{i}, '-', 'color', cacell{i});
    plot(x{i}(mkrIndex{i}), y{i}(mkrIndex{i}), mkrType0{i}, ...
        'color', cacell{i}, 'MarkerFaceColor', mkrFacecolor{i});
end

% Reserve labels
xlabel(xLabel)
ylabel(yLabel)

% reserve limits
xlim(xLim);
ylim(yLim);

if ~isempty(strLegend)
    legend(strLegend,'location',locLegend);
end
