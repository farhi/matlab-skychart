classdef skychart < handle
  % SKYCHART: a class to plot a sky chart with stars/objects
  %
  % This class comptes and plots the sky seen at given location and time. About
  % 43000 stars and 13000 deep sky objects are considered, as well as the Sun, the 
  % Moon and 7 planets. The actual number of rendered objects depends on the zomm 
  % level in the sky chart.
  %
  % You may zoom the plot using the Zoom tool (in the Toolbar). You may as well 
  % use the drag tool to move the visible area. Right-click shows a contextual
  % menu with the under-lying object properties (coordinates, type, ...).
  %
  % To use this code, type
  %
  % >> sc = skychart
  %
  % displays the view at the current UTC and location.
  %
  % Methods:
  %   skychart:   create the view
  %   date:       set/get the date (UTC)
  %   load:       load the catalogs. Done at start.
  %   getplace:   get the current GPS location from the network.
  %   plot:       plot/replot the view.
  %
  % You may force a re-computation and replot of the sky view with:
  %
  % >> compute(sc, 'force')
  % >> plot(sc, 'force)
  %
  % Credits:
  % E. Ofek     MAAT            GPL3 2004
  %   http://weizmann.ac.il/home/eofek/matlab/ 
  % F. Glineur  parse_json      BSD  2009
  %   https://fr.mathworks.com/matlabcentral/fileexchange/23393--another--json-parser
  % P. Tenenbaum Local Time to UTC, BSD 2008
  %   https://fr.mathworks.com/matlabcentral/fileexchange/22295-local-time-to-utc
  % Stars (~46000) data base  from http://astrosci.scimuze.com/stellar_data.htm
  % Deep sky objects (~13000) from http://klima-luft.de/steinicke/ngcic/ngcic_e.htm

  
  properties
  
    catalogs  = {};       % contains data bases as a cell of struct
    utc       = [];       % UTC
    place     = [];       % [ long lat in deg ]
    julianday = 0;
    update_time = 0;      % local time of last computation
    update_period = 120;  % in seconds
    figure    = 0;
    handles   = [];       % handles of the view to delete when updating the plot
    telescope = [];
    visibility= [];       % a polygon area in the X,Y coordinates (defined by user)
    xlim      = [0 0];
    ylim      = [0 0];
    
    % catalogs is a struct array of single catalog entries.
    % Each named catalog entry has fields:
    %  RA:    Right Ascension in deg
    %  DEC:   Declinaison in deg
    %  DIST:  Distance Mpc for objects, pc for stars
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
    function sc = skychart(telescope)
      
      % load catalogs
      load(sc);
      
      % populate with starting stuff
      getplace(sc);
      % connect to telesqcope (StarBook)
      if nargin
        sc.telescope = telescope;
      end
      
      % update all and compute Alt/Az, X/Y coordinates
      compute(sc);
      
      % plot the sky view
      plot(sc);
    end
    
    function load(self)
      % load catalogs: objects, stars
      disp([ mfilename ': Welcome ! Loading Catalogs:' ]);
      self.catalogs = load(mfilename);
      
      % create planet catalog with empty coordinates
      self.catalogs.planets = struct('Description','Planets - http://wise-obs.tau.ac.il/~eran/matlab.html','RA',1:9);
      
      % display available catalogs
      for f=fieldnames(self.catalogs)'
        name = f{1};
        if ~isempty(self.catalogs.(name))
          num  = numel(self.catalogs.(name).RA);
          if isfield(self.catalogs.(name), 'Description')
            desc = self.catalogs.(name).Description;
          else desc = ''; end
          disp([ mfilename ': ' name ' with ' num2str(num) ' entries.' ]);
          disp([ '  ' desc ])
        end
      end
    end
   
    function getplace(self, sb)
      % getplace(sc): get the location
      % getplace(sc, [lon lat]): set location in deg
      % getplace(sc, starbook): set location from StarBook
      if nargin > 1 && isa(sb, 'starbook')
        p = sb.place;
        e = double(p{2})+double(p{3})/60;
        if p{1} == 'W', e=-e; end
        n = double(p{5})+double(p{6})/60;
        if p{4} == 'S', n=-n; end
        self.place = [ n e ];
        disp([ mfilename ': Using location from Vixen StarBook [long lat]=' mat2str(self.place) ]);
      elseif nargin > 1 && isnumeric(sb) && numel(sb) == 2
        self.place = sb;
      else % use ip-api to get the location from the IP
        try
          % could also use: https://api.ipdata.co/
          ip = urlread('http://ip-api.com/json');
          ip = parse_json(ip);  % into struct (private)
          self.place = [ ip.lon ip.lat ];
          disp([ mfilename ': Your are located in ' ip.city ' ' ip.country ' [long lat]=' num2str(self.place) ]);
        catch
          self.place = [ 45.26 5.45 ];
        end
      end
    end
   
    function d=date(self, utc)
      % date(sc):      get the current UTC
      % date(sc, UTC): sets the date as given UTC
      if nargin > 1 && ~isempty(utc)
        self.utc = utc;
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
    end
    
    function compute(self, utc)
      % compute(sc):      compute and update all catalogs
      % compute(sc, utc): the same, but for a given UTC
      
      if nargin < 2, utc = []; end
    
      if strcmp(utc,'force') force = true; utc = [];
      else force = false; end
      self.date(utc); % set UTC

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
        [catalog.X, catalog.Y]    = pr_stereographic_polar( ...
          catalog.Az+90, catalog.Alt);
          
        % special case for constellations: compute lines for patterns
        if strcmp(f{1}, 'constellations')
          catalog = compute_constellations(catalog, ...
            self.julianday, self.place);
        end
          
        % update catalog
        catalog.julianday    = self.julianday;
        self.catalogs.(f{1}) = catalog;
        
      end
      
      self.update_time = self.julianday;
    end
    
    function h = plot(self, force)
      % plot(sc): plot the sky chart
      
      if nargin < 2, force = false; end
      if ~ishandle(self.figure), force=true; end
      
      % create or get current figure.
      [self.figure, xl, yl, new] = plot_frame(self.date, self);
      h  = self.figure;

      % only plot if the figure was closed, or zoom/visible area has changed
      if ~isempty(self.figure) && (isempty(force) || strcmp(force,'force') || ~force) ...
        && all(xl == self.xlim) && all(yl == self.ylim)
        return
      end

      delete(self.handles(ishandle(self.handles)));
      hold on
      % plot constellations and restore current xlim/ylim
      self.handles = plot_constellations(self.catalogs.constellations);
      self.xlim = xl;
      self.ylim = yl;

      % plot catalogs, restricting to magnitude and xlim/ylim
      [handles, self.catalogs] = ...
                       plot_catalogs(self.catalogs, self.xlim, self.ylim);
      self.handles = [ self.handles handles ];
      if new, plot_legend; end
    end
   
  end % methods
  
end % skychart

