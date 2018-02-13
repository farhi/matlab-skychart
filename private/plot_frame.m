function [h,x,y,new] = plot_frame(Date, sc)
  % plot_frame: plot the chart frame (horizon), return the handle and x/y limits

  h = findobj('Tag','SkyChart');
  if isempty(h)
    h = figure('Tag','SkyChart', ...
      'Name', [ 'SkyChart: ' datestr(Date) ' (UTC)' ], ...
      'MenuBar','none', 'ToolBar','figure');
      
    %--- Horizon ---
    Theta = (0:5:360)';
    X     = cosd(Theta);
    Y     = sind(Theta);
    plot(X,Y,'LineWidth',3); title([ datestr(Date) ' (UTC)' ]); axis tight

    %--- Plot N/E/S/W ---
    Offset = 0.02;
    Letter = 0.05;
    text(0,                1+Offset,       'N', 'Color','b');
    text(-1-Offset-Letter, 0,              'E', 'Color','b');
    text(0,               -1-Offset-Letter,'S', 'Color','b');
    text(1+Offset,         0,              'W', 'Color','b');
    set(gca,'XTick',[],'YTick',[], 'Color','k');
    new = true;
    % now create some menu entries
    set(h, 'UserData', sc); % store the skychart handle in the figure
    m = uimenu(h, 'Label', 'SkyChart');
    uimenu(m, 'Label', 'Close',        ...
      'Callback', 'filemenufcn(gcbf,''FileClose'')','Accelerator','w');
    uimenu(m, 'Label', 'Save',        ...
      'Callback', 'filemenufcn(gcbf,''FileSave'')','Accelerator','s');
    uimenu(m, 'Label', 'Save As...',        ...
      'Callback', 'filemenufcn(gcbf,''FileSaveAs'')');
    uimenu(m, 'Label', 'Print',        ...
      'Callback', 'printdlg(gcbf)');  
    uimenu(m, 'Label', 'Update To Current Time',  'Separator','on', ...
      'Callback', @MenuCallback, 'Accelerator','u');
    uimenu(m, 'Label', 'Refresh Plot', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Reset Plot', ...
      'Callback', @MenuCallback);
    % bound listeners for gca:xlim/ylim and figure:resize actions
    propListener = addlistener(gca,'XLim','PostSet',@axesLimitsCallback);
    propListener = addlistener(gca,'YLim','PostSet',@axesLimitsCallback);
  else
    set(0,'CurrentFigure', h); % select but not raise
    set(h, 'Name', [ 'SkyChart: ' datestr(Date) ' (UTC)' ]);
    new = false;
  end
  x = xlim(gca);
  y = ylim(gca);
  set(gca, 'Tag', 'SkyChart_Axes');

end % plot_frame
