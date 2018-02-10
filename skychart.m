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
    telescope = [];
    visibility=[];        % a polygon area in the X,Y coordinates (defined by user)
    xlim      =[-1 1];
    ylim      =[-1 1];
    mlim      = 8;        % magnitude limit
    
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
      % sc.figure = plot(sc);
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
        disp([ mfilename ': Using location from Vixen StarBook [long lat]=' num2str(self.place) ]);
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
      % compute(sc): compute and update all catalogs
      
      if nargin < 2, utc = []; end
    
      self.date(utc); % get UTC
      
      % compute RA/DEC for planets
      % self.catalogs.planets = compute_planets(self.julianday, self.place);
      
      % catalogs -> store new stereographic polar coordinates and update time.
      for f=fieldnames(self.catalogs)'
        % update all valid catalogs
        catalog = self.catalogs.(f{1});
        if isempty(catalog) || ~isfield(catalog, 'RA') || isempty(catalog.RA)
          continue;
        end
        
        % do we need to compute/update ?
        if isfield(catalog,'julianday') && ...
          abs(self.julianday - catalog.julianday) < self.update_period/3600/24
          continue; % no need to update
        end
        
        % compute the horizontal coordinates (Alt-Az)
        [catalog.Az, catalog.Alt] = radec2altaz(catalog.RA, catalog.DEC, ...
          self.julianday, self.place);
        % compute the stereographic polar coordinates
        d2r = pi/180;
        [catalog.X, catalog.Y]    = pr_stereographic_polar( ...
          (catalog.Az(:)+90)*d2r, (90-catalog.Alt(:))*d2r);
          
        % special case for constellations: compute lines for patterns
        if strcmp(f{1}, 'constellations')
          catalog = compute_constellations(catalog, ...
            self.julianday, self.place);
        end
          
        % update catalog
        catalog.julianday = self.julianday;
        self.catalogs.(f{1}) = catalog;
        
        
      end
      
      self.update_time = self.julianday;
    end
    
    function plot(self)
      plot_frame(self.date);
      plot_constellations(self.catalogs.constellations);
    end
   
  end % methods
  
end % skychart

% ------------------------------------------------------------------------------
% private functions
% ------------------------------------------------------------------------------

function JD = skychart_julianday(d)
  % depends: julday, convertdms
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
  % uses: horiz_coo, refraction (MAAT by Ofek)
  
  d2r = pi/180;
  r2d = 1/d2r;

  % taken from MAAT: horiz_coo and refraction (take [rad])
  % Eran O. Ofek http://weizmann.ac.il/home/eofek/matlab/
  HorizCoo          = horiz_coo([ RA(:) DEC(:) ]*d2r, ...
    julianday, place*d2r, 'h');
  HorizCoo(:,2)     = HorizCoo(:,2) + refraction(HorizCoo(:,2));
  HorizCoo          = HorizCoo*r2d;
  Az                = HorizCoo(:,1);
  Alt               = HorizCoo(:,2); % > 0 when visible

end % radec2altaz

function [X,Y]=pr_stereographic_polar(Az,ZenithDist);
%------------------------------------------------------------------------------
% pr_stereographic_polar function                                     AstroMap
% Description: Project coordinates (longitude and latitude) using the
%              Stereographic polar projection.
%              This projection preservs angles.
% Input  : - Vector of Azimuth, in radians.
%          - Vector of Zenith-distance, in radians.
% Output : - Vector of X position
%          - Vector of Y position 
% Tested : Matlab 5.3
%     By : Eran O. Ofek                  November 2004  
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------
  X     = cos(Az).*tan(0.5.*ZenithDist);
  Y     = sin(Az).*tan(0.5.*ZenithDist);
end % pr_stereographic_polar

function constellations = compute_constellations(constellations, julianday, place)

  % compute Alt-Az of constellation lines (pairs)
  [constellations.Lines_Az1, constellations.Lines_Alt1] = radec2altaz( ...
    constellations.Lines_RA1, constellations.Lines_DEC1, ...
    julianday, place);
  [constellations.Lines_Az2, constellations.Lines_Alt2] = radec2altaz( ...
    constellations.Lines_RA2, constellations.Lines_DEC2, ...
    julianday, place);
  
  d2r = pi/180;
  [constellations.X1, constellations.Y1]     = pr_stereographic_polar( ...
          (constellations.Lines_Az1+90)*d2r, (90-constellations.Lines_Alt1)*d2r);
  [constellations.X2, constellations.Y2]     = pr_stereographic_polar( ...
          (constellations.Lines_Az2+90)*d2r, (90-constellations.Lines_Alt2)*d2r);
          
end % compute_constellations


function constellations = plot_constellations(constellations)

  
  % identify the visible Constellations and plot
  v            = find(constellations.Lines_Alt1>0 ...
                    | constellations.Lines_Alt2>0);
  
  X1 = constellations.X1(v);
  X2 = constellations.X2(v);
  Y1 = constellations.Y1(v);
  Y2 = constellations.Y2(v);
  hold on
  
  X = nan*ones(numel(X1)*3, 1); Y=X;
  X(1:3:(end-2)) = X1; X(2:3:(end-1)) = X2; 
  Y(1:3:(end-2)) = Y1; Y(2:3:(end-1)) = Y2; 
  line(X,Y, 'Color','g','LineWidth',1);

  %--- Plot Constellation Names ---
  for Icn=1:1:length(constellations.X)
     Htext = text(constellations.X(Icn),constellations.Y(Icn), ...
       constellations.Name{Icn});
     set(Htext,'FontSize',8,'Color','g');
  end

end % plot_constellations

function h = plot_frame(Date)

  h = findobj('Tag','SkyChart');
  if isempty(h)
    h = figure('Tag','SkyChart', 'Name', [ 'SkyChart: ' datestr(Date) ' (UTC)' ]);
      
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
  else figure(h);
  end

end % plot_frame
