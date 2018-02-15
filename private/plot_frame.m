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
    % figure menus
    m = uimenu(h, 'Label', 'File');
    
    uimenu(m, 'Label', 'Save',        ...
      'Callback', 'filemenufcn(gcbf,''FileSave'')','Accelerator','s');
    uimenu(m, 'Label', 'Save As...',        ...
      'Callback', 'filemenufcn(gcbf,''FileSaveAs'')');
    uimenu(m, 'Label', 'Print',        ...
      'Callback', 'printdlg(gcbf)');
    uimenu(m, 'Label', 'Close',        ...
      'Callback', @MenuCallback, ...
      'Accelerator','w', 'Separator','on');
      
    m = uimenu(h, 'Label', 'SkyChart');
    uimenu(m, 'Label', 'Compute For Given Time', ...
      'Callback', @MenuCallback, 'Accelerator','t');
    uimenu(m, 'Label', 'Update To Current Time', ...
      'Callback', @MenuCallback, 'Accelerator','u');
    uimenu(m, 'Label', 'Refresh Plot', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Reset Plot', ...
      'Callback', @MenuCallback);
      
    m = uimenu(h, 'Label', 'Scope');
    uimenu(m, 'Label', 'Connect to Scope', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'GOTO Selected Object', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Add Selected Object to List', ...
      'Callback', @MenuCallback, 'Separator','on');
    uimenu(m, 'Label', 'Show List', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Set List Period', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Clear List', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Start/Stop GOTO List', ...
      'Callback', @MenuCallback);
      
    % bound listeners for gca:xlim/ylim and figure:resize actions
    propListener = addlistener(gca,'XLim','PostSet',@axesLimitsCallback);
    propListener = addlistener(gca,'YLim','PostSet',@axesLimitsCallback);
  else
    set(0,'CurrentFigure', h); % select but not raise
    set(h, 'Name', [ 'SkyChart: ' datestr(Date) ' (UTC)' ]);
    title([ datestr(Date) ' (UTC)' ]); 
    new = false;
  end
  x = xlim(gca);
  y = ylim(gca);
  set(gca, 'Tag', 'SkyChart_Axes', 'UserData', sc);

end % plot_frame

function axesLimitsCallback(src, evnt)
  % axesLimitsCallback: trigered when a zoom was used
  ax = evnt.AffectedObject;
  self=get(ax, 'UserData');
  plot(self, 1);
end

function MenuCallback(src, evnt)
  % MenuCallback: execute callback from menu.
  %   the action depends on the src Label (uimenu)
  sc = get(gcbf,'UserData');
  
  switch lower(get(src, 'Label'))
  case {'compute for given time'}
    % request Date/Time
    prompt = {'Enter Date Time (e.g. 14-Feb-2018 11:58:15)'};
    name = 'SkyChart: Set Date-Time';
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer=inputdlg(prompt,name, 1, {sc.utc}, options);
    if ~isempty(answer), date(sc, answer{1}); end
    compute(sc);
    plot(sc, 1);
  case {'update', 'update to current time'}
    compute(sc,'now');
    plot(sc, 1);
  case {'replot', 'refresh plot'}
    plot(sc, 1);
  case {'reset', 'reset plot'}
    figure(sc.figure);
    set(gca, 'XLim', [-1 1], 'YLim', [-1 1]);
  case 'find'
    % find an object from its name
  case 'connect to scope'
    % instantiate a StarBook object
    connect(sc);
  case 'send scope to selected object'
    goto(sc);
  case 'close'
    % close figure stop timer, etc
    filemenufcn(gcbf,'FileClose');
    stop(sc.timer);
    if ~isempty(sc.telescope) && isvalid(sc.telescope)
      if ~isempty(sc.telescope.timer) && isvalid(sc.telescope.timer)
        stop(sc.telescope.timer);
      end
      close(sc.telescope.figure);
      delete(sc.telescope);
    end
  case 'add selected object to list'
    listAdd(sc);
  case 'show list'
    listShow(sc);
  case 'clear list'
    listClear(sc);
  case 'start/stop goto list'
    if sc.list_start
      % already running: we stop execution
      sc.list_start = 0;
    else
      listRun(sc);
    end
  case 'set list period'
    listPeriod(sc);
  end
end % MenuCallback
