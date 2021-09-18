%-----------------------------------------------------------------------------%
% Plot setup
imgdir = './fig/';

fscale = 0.8;
width = fscale*5.50;                 % width in inches
height = fscale*4.50;                % height in inches
fsz = 13;                            % font size
lw = 1.5;                            % line width
msz = 7;                             % marker size

% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;

% Lines and markers
set(0,'defaultLineLineWidth',lw);    % set the default line width to lw
set(0,'defaultLineMarkerSize',msz);  % set the default line marker size to msz

% Font sizes
set(0,'defaultAxesFontSize',fsz);
set(0,'defaultColorbarFontSize',fsz);
set(0,'defaultGeoaxesFontSize',fsz);
set(0,'defaultLegendFontSize',fsz);
set(0,'defaultPolaraxesFontSize',fsz);
set(0,'defaultTextboxshapeFontSize',fsz);
set(0,'defaultTextFontSize',fsz);

% Text interpreters
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultColorbarTickLabelInterpreter','latex');
set(0,'defaultGraphplotInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultPolaraxesTickLabelInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultTextarrowshapeInterpreter','latex');
set(0,'defaultTextboxshapeInterpreter','latex');

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on');  % This is the default anyway
set(0,'defaultFigurePaperUnits','inches');  % This is the default anyway
defsize = [8.5, 11];
left = (defsize(1) - width)/2;
bottom = (defsize(2) - height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
%-----------------------------------------------------------------------------%
