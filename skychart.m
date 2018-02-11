classdef skychart < handle

  % SKYCHART: a class to plot a sky chart with stars/objects
  %
  % Credits:
  % E. Ofek     MAAT            GPL3 2004
  %   http://weizmann.ac.il/home/eofek/matlab/ 
  % F. Glineur  parse_json      BSD  2009
  %   https://fr.mathworks.com/matlabcentral/fileexchange/23393--another--json-parser
  
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
    
    % catalogs is a struct of single catalog entries.
    % Each named catalog entry has fields:
    %  RA:    Right Ascension in deg
    %  DEC:   Declinaison in deg
    %  DIST:  Distance Mpc for objects, pc for stars
    %  MAG:   Visual magnitude
    %  TYPE:  star (O B A F G K M), planet, galaxy, open cluster, ...
    %  NAME:  Usual name; Can be a cellstr or string
    %  SIZE:  Visual size, X or [X Y], 0 for stars
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
      self.catalogs.planets = struct('Description','Planets','RA',[]);
      
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
      self.julianday = skychart_julianday(datenum(self.utc));
     
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
      if ~isempty(self.figure) && (isempty(force) || ~force) ...
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








% ------------------------------------------------------------------------------
% private functions
% ------------------------------------------------------------------------------

function JD = skychart_julianday(d)
  % depends: julday, convertdms from MAAT
  if ischar(d)
    d = datenum(d); % make it a number
  end
  Date = str2num(datestr(d,'dd mm yyyy HH MM SS')); % vector
  Frac = convertdms(Date(:,4:6),'H','f');
  JD   = julday([Date(:,1:3), Frac]);
end % skychart_julianday

function [Az Alt] = radec2altaz(RA,DEC, julianday, place)
  % radec2altaz: convert RA/DEC coordinates to horizontal frame
  %
  % [az,alt] = radec2altaz(ra, dec, julianday, place): compute the Alt Az 
  %   (horizontal coordinates) from the RA DEC ones, all in [deg].
  %
  % DEC is the      'height' respective to North pole
  % Altitude is the 'height' respective to Zenith
  %
  % uses: horiz_coo, refraction (MAAT by E. Ofek)
  %       URL : http://wise-obs.tau.ac.il/~eran/matlab.html
  
  d2r = pi/180;
  r2d = 1/d2r;

  % taken from MAAT: horiz_coo and refraction (take [rad])
  % Eran O. Ofek http://weizmann.ac.il/home/eofek/matlab/
  HorizCoo          = horiz_coo([ RA(:) DEC(:) ]*d2r, ...
    julianday, place*d2r, 'h');
  HorizCoo(:,2)     = HorizCoo(:,2) + refraction(HorizCoo(:,2));
  HorizCoo          = HorizCoo*r2d;
  Az                = HorizCoo(:,1);
  Alt               = 90-HorizCoo(:,2);

end % radec2altaz

function [X,Y]=pr_stereographic_polar(Az,ZenithDist);
%------------------------------------------------------------------------------
% pr_stereographic_polar function                                     AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              Stereographic polar projection.
%              This projection preserves angles.
% Input  : - Vector of Azimuth, in deg.
%          - Vector of Zenith-distance, in deg.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                  November 2004  
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------
  d2r   = pi/180;
  Az    = Az*d2r; ZenithDist=ZenithDist*d2r;
  X     = cos(Az).*tan(0.5.*ZenithDist);
  Y     = sin(Az).*tan(0.5.*ZenithDist);
end % pr_stereographic_polar

function constellations = compute_constellations(constellations, julianday, place)
  % compute_constellations: compute Alt-Az coordinates for constellations

  % compute Alt-Az of constellation lines (pairs)
  [constellations.Lines_Az1, constellations.Lines_Alt1] = radec2altaz( ...
    constellations.Lines_RA1, constellations.Lines_DEC1, ...
    julianday, place);
  [constellations.Lines_Az2, constellations.Lines_Alt2] = radec2altaz( ...
    constellations.Lines_RA2, constellations.Lines_DEC2, ...
    julianday, place);

  [constellations.X1, constellations.Y1]     = pr_stereographic_polar( ...
          constellations.Lines_Az1+90, constellations.Lines_Alt1);
  [constellations.X2, constellations.Y2]     = pr_stereographic_polar( ...
          constellations.Lines_Az2+90, constellations.Lines_Alt2);
          
end % compute_constellations

function handles = plot_constellations(constellations)
  % plot_constellations: plot constellation lines and labels
  
  % identify the visible Constellations and plot
  X1 = constellations.X1;
  X2 = constellations.X2;
  Y1 = constellations.Y1;
  Y2 = constellations.Y2;
  v  = (X1.*X1+Y1.*Y1 < 1 | X2.*X2+Y2.*Y2 <1);
  X1 = X1(v); Y1 = Y1(v); 
  X2 = X2(v); Y2 = Y2(v);
  
  X = nan*ones(numel(X1)*3, 1); Y=X;
  X(1:3:(end-2)) = X1; X(2:3:(end-1)) = X2; 
  Y(1:3:(end-2)) = Y1; Y(2:3:(end-1)) = Y2; 
  handles = line(X,Y, 'Color','g','LineWidth',1);
  set(handles, 'Tag', 'SkyChart_Constellations');

  %--- Plot Constellation Names ---
  for Icn=1:1:length(constellations.X)
     X = constellations.X(Icn);
     Y = constellations.Y(Icn);
     if X*X+Y*Y > 1, continue; end
     Htext   = text(X, Y, constellations.Name{Icn}, 'FontSize',8,'Color','g');
     handles = [ handles Htext ];
  end

end % plot_constellations

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
    uimenu(m, 'Label', 'Update All',  'Separator','on', ...
      'Callback', @MenuCallback, 'Accelerator','u');
    % set callback for mouse click on e.g. objects
    set(gca, 'ButtonDownFcn',        @ButtonDownCallback);
  else
    set(0,'CurrentFigure', h); % select but not raise
    set(h, 'Name', [ 'SkyChart: ' datestr(Date) ' (UTC)' ]);
    new = false;
  end
  x = xlim(gca);
  y = ylim(gca);
  set(gca, 'Tag', 'SkyChart_Axes');

end % plot_frame

function planets = compute_planets(planets, julianday, place)

  % magnitude of objects -> size
  MagFun   = [-30 12;...
             -5 10;...
              0 8.5;...
              1 6.5;...
              2 4;...
              3 3;...
              4 2;...
              5 1;...
              9 1];
  Labels = {'Sun','Moon','Mercury','Venus','Mars', ...
            'Jupiter','Saturn','Uranus','Neptune'};
  r2d = 180/pi;
  
  
  if ~isfield(planets, 'DEC')
    planets.Description = 'Planets - http://wise-obs.tau.ac.il/~eran/matlab.html';
    planets.RA = [];
    planets.DEC = [];
    planets.MAG = [];
    planets.DIST = [];
    planets.SIZE = [];
    planets.TYPE = {};
    planets.NAME = {};
  end
    

  for index=1:numel(Labels)
    RA = ''; Coo=[];
    switch Labels{index}
    case 'Sun'
      [RA, Dec]  = suncoo(julianday,'a');
      Mag = -26; Dist = 149597870700; %km
    case 'Moon'
      [RA, Dec] = mooncool(julianday,place,'b');
      Mag = -5; Dist = 384400; % km
    case 'Mercury'
      [Coo,Dist,Ang,Mag] = planet_lowephem(julianday, ...
        'Mercury','Earth','SphericEq','date');
    case 'Venus'
      [Coo,Dist,Ang,Mag] = planet_lowephem(julianday, ...
        'Venus','Earth','SphericEq','date');
    case 'Mars'
      [Coo,Dist,Ang,Mag] = planet_lowephem(julianday, ...
        'Mars','Earth','SphericEq','date');
    case 'Jupiter'
      [Coo,Dist,Ang,Mag] = planet_lowephem(julianday, ...
        'Jupiter','Earth','SphericEq','date');
    case 'Saturn'
      [Coo,Dist,Ang,Mag] = planet_lowephem(julianday, ...
        'Saturn','Earth','SphericEq','date');
    case 'Uranus'
      [Coo,Dist,Ang,Mag] = planet_lowephem(julianday, ...
        'Uranus','Earth','SphericEq','date');
    case 'Neptune'
      [Coo,Dist,Ang,Mag] = planet_lowephem(julianday, ...
        'Neptune','Earth','SphericEq','date');
    end
    if isempty(RA) && numel(Coo) > 2
      RA = Coo(1); Dec = Coo(2);
    end
    
    planets.NAME{index} = Labels{index};
    planets.RA(index)   = RA*r2d;
    planets.DEC(index)  = Dec*r2d;
    planets.MAG(index)  = Mag(1);
    planets.DIST(index) = Dist(1);
    planets.SIZE(index) = 0;
    planets.TYPE{index} = 'planet';
    
    % updates are only done for the Sun and the Moon. Skip others which do not 'move'
    if numel(planets.RA) == numel(Labels) && index == 2
      break;
    end
  end
  
end % compute_planets

function [handles, catalogs] = plot_catalogs(catalogs, xl, yl)
  % plot_catalogs(catalogs, xl, yl): plot all visible objects from the catalog
  %
  handles = [];
  
  % we shall scale the symbols as well to the figure size
  p = get(gcf, 'Position'); p=mean(p(3:4));
  factor = p/1024;
  
  for f=fieldnames(catalogs)'
    h = [];

    catalog = catalogs.(f{1});
    if ~isfield(catalog, 'MAG'), continue; end
    
    % limit magnitude to show, so that we only have 1000 objects max
    mag_max = abs(2.5*log10(abs(4/diff(xl)/diff(yl))))+5;
    
    visible = (catalog.X.^2+catalog.Y.^2 < 1 ...
            & min(xl) < catalog.X & catalog.X < max(xl) ...
            & min(yl) < catalog.Y & catalog.Y < max(yl) ...
            & catalog.MAG(:) & catalog.MAG(:) < mag_max);

    catalog.visible = visible; % store it as visible state
    
    mag = catalog.MAG(visible);
    x   = catalog.X(visible);
    y   = catalog.Y(visible);
    sz  = catalog.SIZE(visible);
    typ = catalog.TYPE(visible);
    
    SZ       = 8-mag+min(mag); % from 10 down ...
    SZ(SZ<1) = 1;
    SZ=SZ.^(2*factor);

    % when out user-inpolygon, color is scaled down (grayed_out) / 2
    
    % scatter plots: o,ly the 'circle' can have non-uniform size
    switch f{1}
    case 'planets'
      % for planets: show name
      % For planets,   coloured disk (scatter) with size=Magnitude^2 colour 
      h = scatter(x,y, markersize(mag).^(2*factor), 'r', 'filled');
      
    case 'stars'
      % for star, when proper name, show it
      % For stars, use coloured hexagon (scatter) with size=Magnitude^2 colour from TYPE
      h = scatter(x,y, SZ, colour(typ), 'filled');
      
    case 'deep_sky_objects'
      % for DSO, when proper name, show it
      % For dso,   use circle        (scatter) with thickness=Magnitude, and given size
      h = scatter(x,y, markersize(mag).^(2*factor), colour(typ), 'o');
      
    otherwise
      h = scatter(x,y, markersize(mag), 'w');
    end
    set(h, 'Tag', [ 'SkyChart_' f{1} ], 'ButtonDownFcn', @ButtonDownCallback, ...
      'UserData', f{1});
    handles = [ handles h ];
    
    catalogs.(f{1}) = catalog;
  end
  
end % plot_catalogs

function m = markersize(mag)
  % markersize(mag): compute the marker size from the magnitude
  
  m = log(13-mag)*3+1;
  m(mag>=13) = 1;
  m = ceil(abs(m));
end % markersize

function c = colour(typ)
  % colour: determine the colour of objects for scatter3
  c = ones(numel(typ),3);  % initialise to white
  
  tokens = { 'star O',  [ 0 0   1 ]; ...
             'star B',  [ 0 0.5 1 ]; ...
             'star A',  [ 0 1   1 ]; ...
             'star F',  [ 0 1   0 ]; ...
             'star G',  [ 1 1   0 ]; ...
             'star K',  [ 1 0.5 0 ]; ...
             'star M',  [ 1 0   0 ]; ...
             'DSO C',   [ 1 0 0 ]; ...
             'DSO D',   [ 1 0 0 ]; ...
             'DSO E',   [ 1 0 0 ]; ...
             'DSO I',   [ 1 0 0 ]; ...
             'DSO P',   [ 1 0 0 ]; ...
             'DSO G',   [ 1 0 0 ]; ...
             'DSO S',   [ 1 0 0 ]; ...
             'DSO G',   [ 1 0 0 ];
             'DSO OCL', [ 0 0 1 ]; ...
             'DSO GCL', [ 0 1 1 ]; ...
             'DSO DN',  [ 0 1 0 ]; ...
             'DSO EN',  [ 0 1 0 ]; ...
             'DSO RN',  [ 0 1 0 ]; ...
             'DSO PN',  [ 0 1 0 ] };

  for index=1:size(tokens, 1)
    tok = tokens{index, 1};
    col = tokens{index, 2};
    ok  = strncmp(typ, tok, numel(tok));
    c(ok,1) = col(1); c(ok,2) = col(2); c(ok,3) = col(3);
  end
end % colour

function legend_h = plot_legend
  % plot_legend: create a legend for object categories
  types = { ...
    'Constellations', 'g.-'; ...
    'Planets',        'r0 '; ...
    'Stars',          'c0 '; ...
    'Galaxies',       'ro '; ...
    'Clusters',       'bo '; ...
    'Nebulae',        'go '};
    
  % we create fake plots and the corresponding legend
  handles = [];
  for index=1:size(types,1)
    this = types{index,2};
    h = plot(1,1, 'Color', this(1));
    if this(2) == '0'
      set(h, 'Marker', 'o', 'MarkerFaceColor', this(1));
    else set(h, 'Marker', this(2)); end
    if this(3) ~= ' '
      set(h, 'LineStyle', this(3)); 
    else set(h, 'LineStyle', 'none'); end
    handles = [ handles h ];
    set(h, 'DisplayName', types{index,1});
  end
  % populate the legend
  legend_h = legend(handles);
  set(legend_h, 'TextColor','w','Color','none','Tag','SkyChart_legend');
 
  
end % plot_legend

function MenuCallback(src, evnt)
  % MenuCallback: execute callback from menu.
  %   the action depends on the src Label (uimenu)
  sc = get(gcf,'UserData');
  
  switch lower(strtok(get(src, 'Label')))
  case 'update'
    compute(sc,'force');
    plot(sc, 1);
  case 'find'
    % find an object from its name 
  end
end % MenuCallback

function ButtonDownCallback(src, evnt)
  % ButtonDownCallback: callback when user clicks on the StarBook image

  % where the mouse click is
  xy = get(gca, 'CurrentPoint'); 
  x = xy(1,1); y = xy(1,2);
  
  % get the SkyChart object handle
  self=get(gcf, 'UserData');
  
  % search for closest object in all catalogs
  found.dist=inf; found.catalog=''; ; found.index=[];
  for f=fieldnames(self.catalogs)'
    catalog = self.catalogs.(f{1});
    if ~isfield(catalog, 'X') || ~isfield(catalog, 'MAG'), continue; end
    % compute distance to closest X/Y object in that catalog
    dist = (catalog.X - x).^2 + (catalog.Y - y).^2;
    if isfield(catalog, 'visible')
      dist(~catalog.visible) = inf;
    end
    [dist_min, dist_index] = min(dist(:));
    % check if that guess is closer than previous guess
    if dist_min < found.dist
      found.dist    = dist_min;
      found.catalog = f{1};
      found.index   = dist_index;
    end
  end
  % display our findings
  catalog = self.catalogs.(found.catalog);
  index   = found.index;
  fprintf(1, '%s: RA=%f DEC=%f MAG=%f TYPE=%s NAME=%s\n', ...
    found.catalog, catalog.RA(index), catalog.DEC(index), ...
    catalog.MAG(index), catalog.TYPE{index}, catalog.NAME{index});
  
end

