function [sc, new] = plot_frame(sc)
  % plot_frame: plot the chart frame (horizon), return the handle and x/y limits

  if ishandle(sc.figure)
    h = sc.figure;
  else 
    h = findall(0, 'Tag','SkyChart'); 
  end
  if numel(h) > 1, delete(h(2:end)); h = h(1); end
  if isempty(h) || ~ishandle(h)
    h = figure('Tag','SkyChart', ...
      'MenuBar','none', 'ToolBar','figure', ...
      'WindowScrollWheelFcn', @ScrollWheelCallback, ...
      'CloseRequestFcn',@MenuCallback, 'UserData', sc);     
    
    % now create some menu entries
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
    uimenu(m, 'Label', 'Compute For Given Location', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Update To Current Time', ...
      'Callback', @MenuCallback, 'Accelerator','u');
    uimenu(m, 'Label', 'Refresh Plot', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Reset Plot', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Help', 'Callback', @MenuCallback);
    uimenu(m, 'Label', 'About SkyChart', 'Callback', @MenuCallback, ...
      'Separator','on');
      
    m = uimenu(h, 'Label', 'Scope');
    uimenu(m, 'Label', 'Connect to Scope', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'GOTO Selected Object', ...
      'Callback', @MenuCallback);
      
    m = uimenu(h, 'Label', 'Planning');
    uimenu(m, 'Label', 'Find object...', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Add Selected Object', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Add Grid around Selected Object', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Edit/Show...', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Set Period...', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Clear', ...
      'Callback', @MenuCallback);
    uimenu(m, 'Label', 'Start/Stop Execution', ...
      'Callback', @MenuCallback, 'Separator','on');
      
    % plot horizon NSEW
    %--- Horizon ---
    Theta = (0:5:360)';
    X     = cosd(Theta);
    Y     = sind(Theta);
    plot(X,Y,'LineWidth',3); 
    axis tight

    %--- Plot N/E/S/W ---
    Offset = 0.02;
    Letter = 0.05;
    text(0,                1+Offset,       'N', 'Color','b');
    text(-1-Offset-Letter, 0,              'E', 'Color','b');
    text(0,               -1-Offset-Letter,'S', 'Color','b');
    text(1+Offset,         0,              'W', 'Color','b');
    set(gca,'XTick',[],'YTick',[], 'Color','k');
    new = true;
      
    % bound listeners for gca:xlim/ylim and figure:resize actions
    propListener = addlistener(gca,'XLim','PostSet',@axesLimitsCallback);
    propListener = addlistener(gca,'YLim','PostSet',@axesLimitsCallback);
    
    sc.figure = h;
    sc.axes   = gca; % where to send plots
    set(sc.axes, 'Tag','SkyChart_Axes');
    
    % start the update timer (when replotting after creation)
    if ~isempty(sc.timer) && isvalid(sc.timer) && strcmp(sc.timer.Running, 'off')
      start(sc.timer);
    end
  else
    new       = false;
  end
  if isempty(sc.axes) || ~ishandle(sc.axes)
    sc.axes = findall(0, 'Tag','SkyChart_Axes');
  end
  % select figure and axes, but not raise
  set(0, 'CurrentFigure', sc.figure);

end % plot_frame

function axesLimitsCallback(src, evnt)
  % axesLimitsCallback: trigered when a zoom/pan was used
  h    = findall(0, 'Tag','SkyChart');
  if numel(h) > 1, delete(h(2:end)); h=h(1); end
  if isempty(h) || ~ishandle(h)
    return
  end
  self = get(h, 'UserData');
  if ~isempty(self) && isobject(self)
    plot(self, 1);
  end
end

function MenuCallback(src, evnt)
  % MenuCallback: execute callback from menu.
  %   the action depends on the src Label (uimenu)
  sc = get(gcbf,'UserData');
  try
    lab = get(src, 'Label');
  catch
    % is this a Figure ?
    lab = get(src, 'Name');
    lab = 'close';
  end
  
  switch lower(lab)
  case {'compute for given time'}
    % request self.utc/Time
    prompt = {'{\color{blue}Enter Date/Time} (e.g. 14-Feb-2018 11:58:15)'};
    name = 'SkyChart: Set Date/Time';
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer=inputdlg(prompt,name, 1, {sc.utc}, options);
    if ~isempty(answer), date(sc, answer{1}); end
    compute(sc);
    plot(sc, 1);
  case {'compute for given location'}
    prompt = {'{\color{blue}Enter GPS location [deg]} (e.g. 5.5 45.2)'};
    name = 'SkyChart: Set GPS location';
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer=inputdlg(prompt,name, 1, {num2str(sc.place)}, options);
    if ~isempty(answer)
      sc.place = str2num(answer{1});
      compute(sc, 'now');
      plot(sc, 1);
    end
  case {'update', 'update to current time'}
    compute(sc,'now');
    plot(sc, 1);
  case {'replot', 'refresh plot'}
    plot(sc, 1);
  case {'reset', 'reset plot'}
    figure(sc.figure);
    set(sc.axes, 'XLim', [-1 1], 'YLim', [-1 1]);
  case {'find','find object...'}
    % find an object from its name
    prompt = {'{\color{blue}Enter a Star/Object Name} e.g. Betelgeuse, M 42, NGC 224. Use spaces between Catalog Name and ID. Known Catalogs include: StarID, HD, HR, M, NGC, IC, ...'};
    name = 'SkyChart: Find Object';
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer=inputdlg(prompt,name, 1, {'M 42'}, options);
    if ~isempty(answer)
      found = findobj(sc, answer{1});
      
      if ~isempty(found) && ishandle(sc.figure)
        % center plot on object and replot
        if isfield(found, 'X') && found.X^2+found.Y^2 < 1
          figure(sc.figure);
          set(sc.axes, 'XLim', [found.X-.1 found.X+.1], ...
                       'YLim', [found.Y-.1 found.Y+.1]);
        end
        plot(sc, 1); 
      end
    end
  case 'connect to scope'
    % instantiate a StarBook object
    connect(sc);
  case {'send scope to selected object','goto selected object'}
    goto(sc);
  case 'close'
    close(sc);
  case 'add selected object'
    listAdd(sc);
  case 'edit/show...'
    listShow(sc);
  case 'clear'
    listClear(sc);
  case 'start/stop execution'
    if sc.list_start
      % already running: we stop execution
      sc.list_start = 0;
    else
      listRun(sc);
    end
  case {'add grid around selected object','add grid around'}
    % request grid size / angular step
    prompt = {'{\color{blue}Enter Grid size} [n x n] or [n x m], e.g. "3 4" for DEC and RA', ...
              '{\color{blue}Enter Angular step} ([deg]). FOV=CameraSensor_{mm}/FocalLength_{mm}*57.3 e.g. "0.75" or " 0.75 1.2" for DEC and RA'};
    name = 'SkyChart: Create Grid';
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer=inputdlg(prompt,name, 1, {'3','0.75'}, options);
    if ~isempty(answer)
      n = str2num(answer{1});
      da= str2num(answer{2});
      listGrid(sc, sc.selected, n, da);
    end
  case 'about skychart'
    try
      im = imread(fullfile(fileparts(which(mfilename)),'..','doc','SkyChart.png'));
    catch
      im = '';
    end
    msg = { [ 'SkyChart <https://github.com/farhi/matlab-skychart>' ], ...
              'A Matlab interface to plot the sky', ...
              '(c) E. Farhi GPL2 2018' };
    if ~isempty(im)
      msgbox(msg,  'About SkyChart', 'custom', im);
    else
      helpdlg(msg, 'About SkyChart');
    end
  case 'help'
    help(sc);
  case 'set period...'
    listPeriod(sc);
  end
end % MenuCallback

function ScrollWheelCallback(src, evnt)
  % ScrollWheelCallback: callback to change speed/zoom with mouse wheel
  
  sc = get(gcbf,'UserData');
  
  if evnt.VerticalScrollCount > 0
    zoom(gcbf, .5);
  else
    zoom(gcbf, 2);
  end

end % ScrollWheelCallback
