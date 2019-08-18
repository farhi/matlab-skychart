classdef skychart < handle
  % SKYCHART: a class to plot a sky chart with stars/objects
  %
  % This class computes and plots the sky seen at given location and time. About
  % 43000 stars and 13000 deep sky objects are considered, as well as the Sun, the 
  % Moon and 7 planets. The actual number of rendered objects depends on the zoom 
  % level in the sky chart.
  %
  % You may zoom the plot using the Zoom tool (in the Toolbar). You may as well 
  % use the drag tool to move the visible area. Right-click shows a contextual
  % menu with the under-lying object properties (coordinates, type, ...).
  %
  % To use this code, type
  %
  % >> sc = skychart;
  %
  % displays the view at the current UTC and location. 
  % You can set the location (in deg) with:
  %
  % >> sc.place=[ 10 -40 ]; compute(sc,'force'); plot(sc, 1);
  %
  % Main Methods
  % ============
  % - skychart:   create the view.
  % - date:       set/get the date (UTC).
  % - getplace:   get the current GPS location from the network.
  % - plot:       plot/replot the view.
  % - connect:    connect to a scope controler.
  % - goto:       send connected scope to selected location.
  % - findobj:    search for a named object and select it.
  %
  % You may force a re-computation and replot of the sky view with:
  %
  % >> compute(sc, 'force')
  % >> plot(sc, 1)
  %
  % Connecting to a Scope
  % =====================
  % You may connect to a telescope mount using e.g.
  %
  % >> connect(sc, scope)
  %
  % where 'scope' should be an object with methods:
  % - goto(sc,RA,DEC):  sends the mount to location (RA,DEC)
  % - get_ra( sc):      returns the RA coordinate as [HH MM SS]
  % - get_dec(sc):      returns the DEC coordinate as [DD MM SS]
  % - get_ra( sc,'target'):      returns the target RA coordinate as [HH MM SS] (optional)
  % - get_dec(sc,'target'):      returns the target DEC coordinate as [DD MM SS] (optional)
  % - get_state(sc):    returns the mount state, e.g. MOVING, TRACKING (optional)
  %
  % Note: when behind a firewall, in order to get the initial GPS location, you may need to set
  %   ProxyHost='proxy.ill.fr'; % Proxy address if you are behind a proxy [e.g. myproxy.mycompany.com or empty]
  %   ProxyPort=8888;           % Proxy port if you are behind a proxy [8888 or 0 or empty]
  %   java.lang.System.setProperty('http.proxyHost', ProxyHost); 
  %   com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(true);
  %   com.mathworks.mlwidgets.html.HTMLPrefs.setProxyHost(ProxyHost);
  %   java.lang.System.setProperty('http.proxyPort', num2str(ProxyPort));
  %   com.mathworks.mlwidgets.html.HTMLPrefs.setProxyPort(num2str(ProxyPort));
  % otherwise, the default location will be used, but can be changed afterwards.
  %
  % Credits
  % =======
  % E. Ofek     MAAT            GPL3 2004
  %   http://weizmann.ac.il/home/eofek/matlab/ 
  % F. Glineur  parse_json      BSD  2009
  %   https://fr.mathworks.com/matlabcentral/fileexchange/23393--another--json-parser
  % P. Tenenbaum Local Time to UTC, BSD 2008
  %   https://fr.mathworks.com/matlabcentral/fileexchange/22295-local-time-to-utc
  % Stars (~46000) data base  from http://astrosci.scimuze.com/stellar_data.htm
  % Deep sky objects (~13000) from http://klima-luft.de/steinicke/ngcic/ngcic_e.htm
  % StarBook (Vixen) https://fr.mathworks.com/matlabcentral/fileexchange/65944-vixen-starbook-control

  
  properties

    catalogs  = [];       % Contains data bases as a cell of struct.
    utc       = [];       % The UTC offset (time-zone, daylight saving) in hours.
    place     = [];       % The current observation [ longitude latitude ] in deg.
    julianday = 0;        % The date/time in Julian calendar.
    update_time = 0;      % The local time of last computation.
    update_period = 120;  % The chart re-computation interval in seconds.
    figure    = [];       % The figure showing the Chart.
    figure_insert = false;% True when inserting the Chart in an existing figure.
    axes_insert   = false;% True when inserting the Chart in an existing axis.
    axes      = [];       % The axis showing the Chart.
    telescope = [];       % The handle to the connected Scope.
    xlim      = [0 0];    % The current X axis limits.
    ylim      = [0 0];    % The current Y axis limits.
    timer     = [];       % The update timer.
    selected  = [];       % The last selected/searched object.
    list      = [];       % The list of selected objects.
    list_start= 0;        % The start time when running the list.
    list_period = 1800;   % The time interval between list/planning GOTO actions in seconds.
    plotting  = false;    % True when plotting
    UserData  = [];       % Open for any further storage from User
    update_counter = 0;   % A counter to automatically replot all after e.g. 10 min.
    
    % catalogs is a struct array of single catalog entries.
    % Each named catalog entry has fields:
    %  RA:    Right Ascension in deg
    %  DEC:   Declinaison in deg
    %  DIST:  Distance pc for objects and stars, au for planets.
    %  MAG:   Visual magnitude
    %  TYPE:  star (O B A F G K M), planet, galaxy, open cluster, ...
    %  NAME:  Usual name; Can be a cellstr or string
    %  SIZE:  Visual size, X or [X Y], 0 for stars (in minutes)
    %  Description: Name and source of catalog
    %
    % Then, at compute, the catalogs are enriched with:
    %  Alt:   Altitude
    %  Az:    Azimuth
    %  X
    %  Y
    %  julianday: the Julian Day used for last computation
    %  visible = (X in xlim, Y in ylim, Alt > 0 and Alt,Az inpolygon(user))
    %
    % We use: catalogs = struct(stars, deepskyobjects, planets)
  end % properties
  
  methods
    function sc = skychart(varargin)
      % SKYCHART Create a SkyChart View.
      %   SKYCHART alone creates a view.
      %
      %   SKYCHART('Property', value, ...) specifies how to create the View.
      %   SKYCHART('figure',f) inserts the SkyChart view into the given figure.
      %   SKYCHART('axes',a)   inserts the SkyChart view into the given axes.
      %   SKYCHART('location',[ longitude latitude ]) specifies the GPS 
      %   location in [deg].
      
      % handle input name/value argument pairs
      flag_set_location = false;
      if mod(nargin,2) == 0 % name/value pairs
        for index=1:2:numel(varargin)
          switch varargin{index}
          case {'axis','axes'}
            sc.axes = varargin{index+1};
            sc.axes_insert   = true;
          case 'figure'
            sc.figure = varargin{index+1};
            sc.figure_insert = true;
          case 'catalogs'
            sc.catalogs = varargin{index+1};
          case {'location','site','place'}
            sc.place = varargin{index+1};
          end
        end
      end
      
      if isempty(sc.catalogs)
        sc.catalogs = getcatalogs;
      end
      
      % populate with starting stuff
      if isempty(sc.place), getplace(sc); end
      
      % update all and compute Alt/Az, X/Y coordinates
      compute(sc, 'now');
      
      % plot the sky view
      plot(sc);
      
      % attach a timer for regular updates (5 sec)
      sc.timer = timer('TimerFcn', @TimerCallback, ...
                'Period', 5.0, 'ExecutionMode', 'fixedDelay', 'Name', mfilename);
      set(sc.timer, 'UserData', sc);
      start(sc.timer);
    end % skychart
   
    function getplace(self, sb)
      % GETPLACE Get the current location.
      %   GETPLACE(sc) requests location from the Network.
      %
      %   GETPLACE(sc, [lon lat]) sets location given in degrees.

      if nargin > 1 && isnumeric(sb) && numel(sb) == 2
        self.place = sb;
      else % use ip-api to get the location from the IP
        try
          % could also use: https://api.ipdata.co/
          ip = urlread('http://ip-api.com/json');
          ip = parse_json(ip);  % into struct (private)
          self.place = [ ip.lon ip.lat ];
          disp([ mfilename ': You are located in ' ip.city ' ' ip.country ' [long lat]=' num2str(self.place) ]);
        catch
          self.place = [ 45.26 5.45 ];
        end
      end
    end % getplace
   
    function d=date(self, utc)
      % DATE Get the current UTC date.
      %   DATE(sc) sets the current time from the clock.
      %
      %   DATE(sc, UTC) sets the date as given UTC (string, number, vector).
      if nargin > 1 && ~isempty(utc)
        if strcmp(utc, 'now')
          self.utc = local_time_to_utc(now);
        else
          self.utc = utc;
        end
      else
        self.utc = local_time_to_utc(now); % in private, uses Java
        % could also use http://www.convert-unix-time.com/api?timestamp=now
      end
      % compute the Julian Day
      d = datenum(self.utc);
      Date = str2num(datestr(d,'dd mm yyyy HH MM SS')); % vector
      Frac = convertdms(Date(:,4:6),'H','f');
      self.julianday = julday([Date(:,1:3), Frac]);
     
      d = self.utc;
    end % date
    
    function ret = compute(self, utc)
      % COMPUTE Compute position of objects for all catalogs.
      %   COMPUTE(sc) updates only when time has changed significantly (2 min).
      %
      %   COMPUTE(sc, utc)     the same, but for a given UTC.
      %
      %   COMPUTE(sc, 'now')   update UTC to now, and force computation.
      %
      %   COMPUTE(sc, 'force') force to compute for previously set date/time.
      
      if nargin < 2, utc = []; end
    
      if strcmp(utc,'force')   force = true; utc = [];
      elseif strcmp(utc,'now') force = true;
      else force = false; end
      
      if ~isempty(utc) % we set a specific date before computation, else use date as stored
        self.date(utc); % set UTC
      end

      % catalogs -> store new stereographic polar coordinates and update time.
      for f=fieldnames(self.catalogs)'
        % update all valid catalogs
        catalog = self.catalogs.(f{1});
        if isempty(catalog) || ~isfield(catalog, 'RA')
          continue;
        end
        
        % do we need to compute/update ?
        if isfield(catalog,'julianday') && ~force && ...
          abs(self.julianday - catalog.julianday) < self.update_period/3600/24
          continue; % no need to update
        end

        % planets
        if strcmp(f{1}, 'planets') % compute new RA/DEC for Sun and Moon
          catalog = compute_planets(catalog, self.julianday, self.place);
        end
        
        % compute the horizontal coordinates (Alt-Az)
        [catalog.Az, catalog.Alt] = radec2altaz(catalog.RA, catalog.DEC, ...
          self.julianday, self.place);
        % compute the stereographic polar coordinates
        delta_az = 0;
        [catalog.X, catalog.Y]    = pr_stereographic_polar( ...
          catalog.Az+90-delta_az, catalog.Alt);
          
        % special case for constellations: compute lines for patterns
        if strcmp(f{1}, 'constellations')
          catalog = compute_constellations(catalog, ...
            self.julianday, self.place, delta_az);
        end
          
        % update catalog
        catalog.julianday    = self.julianday;
        self.catalogs.(f{1}) = catalog;
        
      end
      
      self.update_time = self.julianday;
      ret = self.update_time;
    end % compute
    
    function h = plot(self, force)
      % PLOT Plot the sky chart.
      %   PLOT(s) plot the skyview when axes limits have changed. Any attached 
      %   Scope coordinates are also shown (red cross/circle).
      %
      %   PLOT(s,'force') force to re-plot the skyview. Any defined list of objects
      %   is also shown (white X).
      
      h = [];
      if self.plotting, return; end
      if nargin < 2, force = false; end
      if isempty(self.figure) || ~ishandle(self.figure) || self.figure_insert, force=true; end
      
      self.plotting = true;
      % create or get current figure. Sets fig focus
      [self, new] = plot_frame(self);
      set(self.figure, 'HandleVisibility','on', 'NextPlot','add');
      hold on
      title([ datestr(self.utc) ' (UTC)' ]);
      set(self.figure, 'Name', [ 'SkyChart: ' datestr(self.utc) ' (UTC)' ]); 

      % when a scope is connected, replot its location
      plot_telescope(self);
      h = self.figure;

      % only plot if the figure was closed, or zoom/visible area has changed
      if ~isempty(self.figure) && (isempty(force) ...
        || (ischar(force) && ~strcmp(force,'force')) ...
        || (isscalar(force) && ~force)) ...
        && all(xlim(self.axes) == self.xlim) && all(ylim(self.axes) == self.ylim)
        hold off
        set(self.figure, 'HandleVisibility','off', 'NextPlot','new');
        self.plotting = false;
        return
      end
      
      % plot selection
      plot_list(self);

      % replot constellations and store current xlim/ylim
      plot_constellations(self);
      self.xlim = xlim(self.axes);
      self.ylim = ylim(self.axes);

      % replot catalogs, restricting to magnitude and xlim/ylim
      self = plot_catalogs(self);
      if new, plot_legend(self); end
      
      hold off
      set(self.figure, 'HandleVisibility','off', 'NextPlot','new');
      self.plotting = false;
      
    end % plot
    
    function close(self)
      % CLOSE Close the SkyChart.
      if isa(self.timer, 'timer') && isvalid(self.timer)
        stop(self.timer); delete(self.timer);
      end
      if ~isempty(self.figure) && ishandle(self.figure) && ~self.figure_insert, delete(self.figure); end
      self.figure = []; self.axes = []; self.plotting = false;
    end
    
    function delete(self)
      % DELETE Delete the SkyChart object.
      close(self);
    end % delete
    
    function found = findobj(self, name)
      % FINDOBJ Search for an object (star, DSO, planet) in catalogs.
      %   FINDOBJ(sc, name) search for a given object in catalogs. Select it.
      %   The found object is returned as a structure.
      %
      %   FINDOBJ(sc) opens a Dialogue to enter the nameof an object to search.
      catalogs = fieldnames(self.catalogs);
      found = [];
      
      if nargin < 2
        prompt = { '{\color{blue}Enter a Star/Object Name} e.g. Betelgeuse, M 42, NGC 224, Venus.', ...
          'Use spaces between Catalog Name and ID.', ...
          'Known Catalogs include: Planets, StarID, HD, HR, Messier, NGC, IC, ...' };
        name = 'SkyChart: Find Object';
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        answer=inputdlg(prompt,name, 1, {'M 42'}, options);
        if ~isempty(answer), name = answer{1}; 
        else return; end
      end
      
      % check first for name without separator
      if ~any(name == ' ')
        [n1,n2]  = strtok(name, '0123456789');
        found = findobj(self, [ n1 ' ' n2 ]);
        if ~isempty(found) return; end
      end
      
      namel= strtrim(lower(name));
      for f=catalogs(:)'
        catalog = self.catalogs.(f{1});
        if ~isfield(catalog, 'X') || ~isfield(catalog, 'MAG'), continue; end
        NAME = lower(catalog.NAME);
        NAME = regexprep(NAME, '\s*',' ');
        % search for name
        index = find(~cellfun(@isempty, strfind(NAME, [ ';' namel ';' ])));
        if isempty(index)
        index = find(~cellfun(@isempty, strfind(NAME, [ namel ';' ])));
        end
        if isempty(index)
        index = find(~cellfun(@isempty, strfind(NAME, [ ';' namel ])));
        end
        if isempty(index)
        index = find(~cellfun(@isempty, strfind(NAME, [ namel ])));
        end
        if ~isempty(index)
          found.index   = index(1);
          found.catalog = f{1};
          found.RA      = catalog.RA(found.index);
          found.DEC     = catalog.DEC(found.index);
          found.Alt     = catalog.Alt(found.index);
          found.Az      = catalog.Az(found.index);
          found.X       = catalog.X(found.index);
          found.Y       = catalog.Y(found.index);
          found.MAG     = catalog.MAG(found.index);
          found.TYPE    = catalog.TYPE{found.index};
          found.NAME    = catalog.NAME{found.index};
          found.DIST    = catalog.DIST(found.index);
          break;
        end
      end

      if ~isempty(found)
        disp([ mfilename ': Selecting object ' name ' as: ' found.NAME ])
        switch found.catalog
        case {'deep_sky_objects','stars'}, un='ly'; found.DIST = found.DIST*3.26;
        case 'planets', un='a.u.';
        otherwise, un='??';
        end
        if found.DIST > 0
          disp(sprintf('  %s: Magnitude: %.1f ; Type: %s ; Dist: %.3g [%s]', ...
            found.catalog, found.MAG, found.TYPE, found.DIST, un ));
        else
          disp(sprintf('  %s: Magnitude: %.1f ; Type: %s', ...
            found.catalog, found.MAG, found.TYPE ));
        end
        if found.X^2+found.Y^2 < 1
          self.selected = found;
          if ishandle(self.figure)
            % center plot on object and replot
            if isfield(found, 'X') && found.X^2+found.Y^2 < 1
              figure(self.figure);
              set(self.axes, 'XLim', [found.X-.1 found.X+.1], ...
                           'YLim', [found.Y-.1 found.Y+.1]);
            end
            plot(self, 1); 
          end
        else
          disp([ mfilename ': object ' name ' is not visible.' ])
        end
      else
        disp([ mfilename ': object ' name ' was not found.' ])
      end
      
    end
    
    function url=help(self)
      % HELP Open the SkyChart Help page.
      url = fullfile('file:///',fileparts(which(mfilename)),'doc','SkyChart.html');
      open_system_browser(url);
    end
    
    % scope methods ------------------------------------------------------------
    
    function connect(self, sb)
      % CONNECT Connect the SkyChart to a scope controller.
      if ~isempty(self.telescope) && isvalid(self.telescope)
        disp([ mfilename ': Scope is already connected.' ])
        plot(self.telescope);
        return
      end
      if nargin > 1
        disp([ mfilename ': Connecting Scope to SkyChart' ]);
        self.telescope = sb;
      end
    end % connect
    
    function goto(self, varargin)
      % GOTO Send the scope to given object/coordinates.
      %   GOTO(sc, ra, dec) send scope to given location
      %     ra  is given in hh:mm:ss or DEG
      %     dec is given in deg:min or DEG
      %
      %   GOTO(sc) send the scope to the last Selected object, e.g. after a search.
      if nargin ==1
        if isfield(self.selected, 'RA') && isfield(self.selected, 'DEC')
          RA = self.selected.RA;
          DEC= self.selected.DEC;
          varargin = { self.selected };
        else return; end
      end
      if isobject(self.telescope) && isvalid(self.telescope) && ismethod(self.telescope, 'goto')
        % send scope
        self.telescope.goto(varargin{:});
      else
        disp([ mfilename ': No Scope is Connected yet. Use "connect(' inputname(1) ', scope)" first' ]);
      end
    end % goto
    
    % list/planning methods ----------------------------------------------------
    
    function l = listAdd(self, RA, DEC, name)
      % listAdd Add an object to the List of planned observations.
      %   listAdd(sc) add the last selected object to the List.
      %
      %   listAdd(sc, name) search for name and add it to the List.
      %
      %   listAdd(sc, RA,DEC, {name}) add RA/DEC (in deg) to the List. 
      %   The name of the object can optionally be given.
      if nargin > 1
        if ischar(RA)
          self.selected = findobj(self, RA);
        elseif isnumeric(RA) && isnumeric(DEC)
          if nargin < 4 || isempty(name)
            name = sprintf('RA=%.2f DEC=%.2f', RA, DEC); 
          end
          self.selected = struct('RA', RA, 'DEC', DEC, 'NAME', name);
        elseif isfield(self.selected, 'RA') && isfield(self.selected, 'DEC')
          self.selected = RA;
        end
      end
      if ~isempty(self.selected) && isfield(self.selected, 'RA') && isfield(self.selected, 'DEC')
        if isempty(self.list)
          self.list = self.selected;
        else
          self.list(end+1) = self.selected;
        end
        disp([ mfilename ': Add to list: ' self.selected.NAME ]);
      end
      if ishandle(self.figure)
        plot(self, 1);
      end
      l = self.selected;
    end % listAdd
    
    function listClear(self)
      % listClear Clear the List.
      %   listClear(sc)
      self.list = [];
      if ishandle(self.figure)
        plot(self, 1);
      end
    end % listClear
    
    function l=listShow(self)
      % listShow Show the current List in a Dialogue window.
      %   listShow(sc)
      %   The Dialogue allows to select objects and remove them.
      
      ListString = {};
      for index=1:numel(self.list)
        ListString{end+1} = self.list(index).NAME;
      end
      disp([ mfilename ': Current List with Period ' num2str(self.list_period) ' [s] between items' ])
      disp(char(ListString))
      if ~isempty(ListString)
        [selection, ok] = listdlg('PromptString', ...
          { [ 'Current List with Period ' num2str(self.list_period) ' [s]' ],...
            'Select items to remove from list and Click REMOVE' }, ...
          'ListSize', [ 480 240 ], ...
          'ListString', ListString, ...
          'Name', 'SkyChart: Current List', ...
          'OKString', 'Remove Selection','CancelString', 'OK');
        if ok == 1 && ~isempty(selection)
          self.list(selection) = [];
        end
      end
      if ishandle(self.figure)
        plot(self, 1);
      end
      l = self.list;
    end % listShow
    
    function listRun(self)
      % listRun Start to execute a list of GOTO's.
      %   listRun(sc)
      self.list_start = true;
    end % listRun
    
    function listPeriod(self, dt)
      % listPeriod Open dialogue to change the List period (in [s]) .
      %   listPeriod(sc) open the dialogue.
      %
      %   listPeriod(sc, dt) set the List period to dt [s]
      if nargin > 1
        self.list_period = dt;
      else
        prompt = {'{\color{blue}Enter List period in [seconds]} (time between each GOTO. This is e.g. observation time)'};
        name = 'SkyChart: Set List Period';
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        answer=inputdlg(prompt,name, 1, {num2str(self.list_period)}, options);
        if isempty(answer), 
          return;
        else
          self.list_period = str2double(answer{1}); 
        end
      end
    end % listPeriod
    
    function l=listGrid(self, RA, DEC, n, da)
      % listGrid Build a grid of observations around an object for e.g. stitching.  
      %   listGrid(sc) build a 3x3 grid around selection with step 0.75 deg.
      %
      %   listGrid(sc, RA, DEC, N, da) build a N x N grid around RA/DEC with angular step da.
      %
      %   listGrid(sc, name   , N, da) build a N x N grid around named object.
      %
      %   The grid size can be given as N = [nDEC nRA] to specify a non-square grid
      %   as well as similarly for the angular step da = [dDEC dRA]
      %
      %   The angular step should be e.g. the field of view (FOV) in order to 
      %   build a panorama / stitch view.
      %   When using a focal length F with a camera sensor size S, the FOV is:
      %     FOV = S/F*57.3 [deg], where S and F should e.g. be in [mm]
      %
      %   With a 1200 mm focal length and an APS-C sensor 23.5x15.6, the FOV is:
      %     FOV = 0.74 and 1.12 [deg]
      %   With a 400 mm focal length and similar sensor:
      %     FOV = 2.23 and 3.36 [deg]
      
      l = [];
      
      if nargin < 2, RA =[]; end
      if nargin < 3, DEC=[]; end
      if nargin < 4, n  =[]; end
      if nargin < 5, da =[]; end
      name = '';
      
      % input as a name
      if ischar(RA) RA = findobj(self, RA); end
      
      % input as a struct (from findobj)
      if isstruct(RA) && isfield(RA, 'RA') && isfield(RA, 'DEC')
        if nargin >= 4, da=n;  n=[]; end
        if nargin >= 3, n=DEC; end
        selected = RA;
        if isempty(selected), return; end
        RA = selected.RA;
        DEC= selected.DEC;
        name = selected.NAME;
      end
      
      if isempty(RA)
        RA = self.selected.RA;
        name = self.selected.NAME;
      end
      if isempty(DEC)
        DEC= self.selected.DEC;
      end
      if isempty(n), n=3;     end
      if isempty(da) da=0.75; end
      
      % build the grid
      if all(isfinite(n)) && all(isfinite(da))
        if isscalar(n),  n  = [n  n]; end
        if isscalar(da), da = [da da]; end
        n = round(n);
        for dec = DEC+da(1)*((0:(n(1)-1))-(n(1)-1)/2)
          for ra = RA+da(2)*((0:(n(2)-1))-(n(2)-1)/2)
            l = [ l listAdd(self, ra, dec, ...
              sprintf('RA=%.2f DEC=%.2f %s', ra, dec, name)) ];
          end
        end
      end
    end % listGrid
    
    function l=grid(self, varargin)
      % GRID Build a grid around current selection.
      %   GRID(s, ...) is equivalent to listGrid.
      l=listGrid(self, varargin{:});
    end
   
  end % methods
  
end % skychart

% ------------------------------------------------------------------------------

function catalogs = getcatalogs
  % load catalogs: objects, stars
  persistent c
  
  catalogs = [];
  
  if isempty(c)
    disp([ mfilename ': Welcome ! Loading Catalogs:' ]);
    
    c = load(mfilename);
    
    % create planet catalog with empty coordinates
    c.planets = struct('Description','Planets - http://wise-obs.tau.ac.il/~eran/matlab.html','RA',1:9);
    
    % display available catalogs
    for f=fieldnames(c)'
      name = f{1};
      if ~isempty(c.(name))
        num  = numel(c.(name).RA);
        if isfield(c.(name), 'Description')
          desc = c.(name).Description;
        else desc = ''; end
        disp([ mfilename ': ' name ' with ' num2str(num) ' entries.' ]);
        disp([ '  ' desc ])
      end
    end
  end
  catalogs = c;
end % getcatalogs

function TimerCallback(src, evnt)
  % TimerCallback: executed regularly (5 sec)
  
  sc = get(src, 'UserData');
  
  if isvalid(sc), 
    % update: compute and plot
    sc.update_counter = sc.update_counter + 1;
    % do we need full update ?
    if sc.update_counter > sc.update_period % every 10 min
      sc.update_counter = 0;
      compute(sc,'force');
      if ~isempty(sc.figure) && ishandle(sc.figure) && ishandle(sc.axes)
        plot(sc,1);
      end
    else
      % fast update
      compute(sc);
      if ~isempty(sc.figure) && ishandle(sc.figure) && ishandle(sc.axes)
        plot(sc);
      end
    end
    
    % look if we are running a list
    if sc.list_start
      if isscalar(sc.list_start) || etime(clock, sc.list_start) > sc.list_period
        % time elapsed between items for GOTO
        % we reset the timer, and goto the first item, then clear it
        if ~isempty(sc.list)
          sc.list_start = clock;
          sc.selected   = sc.list(1);
          goto(sc);        % go there
          sc.list(1) = []; % remove that element from the list
        else
          % all list elements have been used. Cancel Run.
          sc.list_start = 0;
        end
      end
    end 
  else delete(src); end
  
end % TimerCallback

