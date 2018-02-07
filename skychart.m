classdef skychart < handle

  % SKYCHART: a class to plot a sky chart with stars/objects
  %
  % Credits:
  % E. Ofek     MAAT            GPL3 2004
  %   http://weizmann.ac.il/home/eofek/matlab/ 
  % F. Glineur  parse_json      BSD  2009
  %   https://fr.mathworks.com/matlabcentral/fileexchange/23393--another--json-parser
  
  properties
  
    catalogs  = {};   % contains data bases as a cell of struct
    date      = [];   % UTC
    place     = [];   % [ long lat in deg ]
    julianday = 0;
    
  end % properties
  
  methods
    function skychart
    end

    function plot
    end

    function update
    end
   
    function getplace(self, sb)
      % getplace(sc): get the location
      % getplace(sc, [lon lat]): set location in deg
      % getplace(sc, starbook): set location from StarBook
      if nargin && isa(sb, 'starbook')
        p = sb.place;
        e = double(p{2})+double(p{3})/60;
        if p{1} == 'W', e=-e; end
        n = double(p{5})+double(p{6})/60;
        if p{4} == 'S', n=-n; end
        self.place = [ n e ];
        disp([ mfilename ': Using location from Vixen StarBook [long lat]=' num2str(self.place) ]);
      elseif isnumeric(sb) && numel(sb) == 2
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
   
   function gettime(self, sb)
     if nargin && isa(sb, 'starbook')
     elseif nargin
       self.date = datestr(sb);
     else
       self.date = local_time_to_utc(now); % in private
       % could also use http://www.convert-unix-time.com/api?timestamp=now
     end
   end
   
   function getstatus
   end
   
   
  end % methods
  
end % skychart
  
  
  if isempty(constellations)
    % load: constellations stars deepskyobjects
    load skychart
  end
  
  %check/convert input
  
  %if Place not given:
  %   get it  from: 'starbook.getplace'
  %   get it  from: http://ip-api.com/json
  %   lat/long in degrees
  
  % separate computation from display
  % store alt-az coordinates, time of computation,
  %   and visible state in database struct
  
  Date        = now;
  Place       = [35 32 0]*pi/180;
  
  % plot all stuff
  JD          = skychart_julianday(Date);
  
  h = skychart_frame(datestr(Date)); hold on;  % get figure handle
  
  skychart_constellations(constellations, JD,Place);  % lines and names
  
  skychart_planets(JD,Place);                  % filled symbols and names

  % stars
  options.color   = 'w';
  options.filled  = 'filled';
  options.MAG_max = 10;
  skychart_stars(stars,          JD,Place, options);    % filled symbols 
  
  % deep sky objects
  options.color   = 'b';
  options.filled  = 'open';
  options.MAG_max = 10;
  skychart_dso(deepskyobjects, JD,Place, options);    % open symbols
  
end % skychart
  
  
% ==============================================================================
  
function h = skychart_frame(Date)

  h = findobj('Tag','SkyChart');
  if isempty(h)
    h = figure('Tag','SkyChart', 'Name', [ 'SkyChart: ' datestr(Date) ]);
      
    %--- Horizon ---
    Theta = (0:5:360)';
    X     = cosd(Theta);
    Y     = sind(Theta);
    plot(X,Y,'LineWidth',3); title(Date); axis tight

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

end % skychart_frame
  
function skychart_constellations(constellations, JD,GeodPos)

  AltLimit = 0;
  %--- Plot Constellations Lines ---
  CL_HorizCoo1 = horiz_coo([constellations.Lines_RA1, constellations.Lines_DEC1], ...
                           JD,GeodPos,'h');
  CL_HorizCoo2 = horiz_coo([constellations.Lines_RA2, constellations.Lines_DEC2], ...
                           JD,GeodPos,'h');
  % add refraction
  CL_HorizCoo1(:,2) = CL_HorizCoo1(:,2) + refraction(CL_HorizCoo1(:,2));
  CL_HorizCoo2(:,2) = CL_HorizCoo2(:,2) + refraction(CL_HorizCoo2(:,2));
  I            = find(CL_HorizCoo1(:,2)>AltLimit ...
                    | CL_HorizCoo2(:,2)>AltLimit);
  CL_HorizCat1 = [CL_HorizCoo1(I,1)+pi./2, pi./2-CL_HorizCoo1(I,2)];
  CL_HorizCat2 = [CL_HorizCoo2(I,1)+pi./2, pi./2-CL_HorizCoo2(I,2)];
  [X1,Y1]      = pr_stereographic_polar(CL_HorizCat1(:,1),CL_HorizCat1(:,2));
  [X2,Y2]      = pr_stereographic_polar(CL_HorizCat2(:,1),CL_HorizCat2(:,2));

  N_CL = length(I);
  for Icl=1:1:N_CL,
     Hcl = plot([X1(Icl);X2(Icl)],[Y1(Icl);Y2(Icl)],'-');
     set(Hcl,'Color','g','LineWidth',1);
  end

  %--- Plot Constellation Names ---
  [X,Y, RestCol, visible] = coo2xy([constellations.RA, constellations.DEC, ...
                          [1:1:length(constellations.RA)]' ], JD,GeodPos);
  N_CN   = length(X);
  for Icn=1:1:N_CN,
     Htext = text(X(Icn),Y(Icn), constellations.Name{RestCol(Icn)});
     set(Htext,'FontSize',8,'Color','g')
  end

end % skychart_constellations
  
function skychart_planets(JD,GeodPos)
  PlCoo    = zeros(0,3);
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
  for index=Labels
    PlCoo = '';
    switch index{1}
    case 'Sun'
      [SunRA, SunDec]  = suncoo(JD,'a');
      PlCoo            = [SunRA, SunDec, ones(length(SunRA),1).*-26];
    case 'Moon'
      [MoonRA,MoonDec] = mooncool(JD,GeodPos,'b');
      PlCoo            = [MoonRA, MoonDec, ones(length(MoonRA),1).*-5]; 
    case 'Mercury'
      [Coo,Dist,Ang,Mag] = planet_lowephem(JD,'Mercury','Earth','SphericEq','date');
    case 'Venus'
      [Coo,Dist,Ang,Mag] = planet_lowephem(JD,'Venus','Earth','SphericEq','date');
    case 'Mars'
      [Coo,Dist,Ang,Mag] = planet_lowephem(JD,'Mars','Earth','SphericEq','date');
    case 'Jupiter'
      [Coo,Dist,Ang,Mag] = planet_lowephem(JD,'Jupiter','Earth','SphericEq','date');
    case 'Saturn'
      [Coo,Dist,Ang,Mag] = planet_lowephem(JD,'Saturn','Earth','SphericEq','date');
    case 'Uranus'
      [Coo,Dist,Ang,Mag] = planet_lowephem(JD,'Uranus','Earth','SphericEq','date');
    case 'Neptune'
      [Coo,Dist,Ang,Mag] = planet_lowephem(JD,'Neptune','Earth','SphericEq','date');
    end
    if isempty(PlCoo)
      PlCoo              = [Coo(:,1:2), Mag(:,1)];
    end
    %--- convert to horizontal coordinates in the map epoch! ---
    [X,Y, RestCol] = coo2xy(PlCoo,JD,GeodPos);

    % plot planets
    if ~isempty(X)
      MagMarkSize = interp1(MagFun(:,1),MagFun(:,2),RestCol(:,1),'linear');
      Hpl         = plot(X,Y,'LineStyle','-','Marker','o', ...
        'MarkerFaceColor','r','MarkerEdgeColor','r','Color','r',...
        'MarkerSize',mean(MagMarkSize));
      text(X,Y,0, index{1},'Color','r');
    end
  end

end % skychart_planets

function skychart_stars(RA_DEC_MAG, JD,Place, options)

  if nargin < 4, options = ''; end
  if ~isfield(options, 'color'),  options.color  = 'w';    end
  if ~isfield(options, 'filled'), options.filled = 'open'; end
  if ~isfield(options, 'MAG_max'),options.MAG_max= 6;end
  
  RA  = RA_DEC_MAG.RA *2*pi/24;
  DEC = RA_DEC_MAG.DEC*pi/180;
  MAG = RA_DEC_MAG.MAG;
  
  [X,Y, ~, visible] = coo2xy([RA(:) DEC(:)], JD,Place);

  % make a symbol size with magnitude
  % mag= MAG_max      -> 0
  % mag =min(MAG)     -> 10
  MAG = MAG(visible);
  SZ       = 10-MAG+min(MAG); % from 10 down ...
  SZ(SZ<1) = 1;
  toplot = find(MAG<options.MAG_max);
  SZ=SZ.^2;
  
  if strcmp(options.filled, 'filled')
    scatter(X(toplot),Y(toplot), SZ(toplot), options.color, 'filled');
  else
    scatter(X(toplot),Y(toplot), SZ(toplot), options.color);
  end
  
end % skychart_stars

function skychart_dso(RA_DEC_MAG, JD,Place, options)

  if nargin < 4, options = ''; end
  if ~isfield(options, 'color'),  options.color  = 'w';    end
  if ~isfield(options, 'filled'), options.filled = 'open'; end
  
  RA  = RA_DEC_MAG.RA;
  DEC = RA_DEC_MAG.DEC;
  MAG = RA_DEC_MAG.MAG;
  SZ  = RA_DEC_MAG.SIZE;
  
  [X,Y, ~, visible] = coo2xy([RA(:) DEC(:)], JD,Place);
  
  if strcmp(options.filled, 'filled')
    scatter(X,Y, SZ, options.color, 'filled');
  else
    scatter(X,Y, SZ, options.color);
  end
  
end % skychart_dso

function skychart_objects(RA,DEC,SZ, JD,Place, color)
  skychart_stars(RA,DEC, -SZ, JD,Place, 'c', 'circle');
end
% ==============================================================================

function JD = skychart_julianday(d)
  % depends: julday, convertdms
  if ischar(d)
    d = datenum(d); % make it a number
  end
  Date = str2num(datestr(d,'dd mm yyyy HH MM SS')); % vector
  Frac = convertdms(Date(:,4:6),'H','f');
  JD   = julday([Date(:,1:3), Frac]);
end % skychart_julianday
  
function [X,Y, RestCol, Iabove] = coo2xy(CooData,JD,GeodPos);
  %--------------------------------------------------------
  % Given [RA, Dec, additional columns]
  % calculate the x,y, position in sterographic polar projection
  % and return [X, Y, Rest of colums] for objects found
  % above the horizon
  % Eran O. Ofek http://weizmann.ac.il/home/eofek/matlab/
  %--------------------------------------------------------
  AltLimit          = 0;
  HorizCoo          = horiz_coo(CooData(:,1:2),JD,GeodPos,'h');
  HorizCoo(:,2)     = HorizCoo(:,2) + refraction(HorizCoo(:,2));
  Iabove            = find(HorizCoo(:,2)>AltLimit);
  HorizCoo          = HorizCoo(Iabove,:);
  if (size(CooData,2)>2),
     RestCol        = CooData(Iabove,3:end);
  else
     RestCol        = [];
  end
  N_data            = size(HorizCoo,1);
  HorizCat          = [HorizCoo(:,1)+pi./2, pi./2-HorizCoo(:,2), ones(N_data,2)];
  [X,Y]             = pr_stereographic_polar(HorizCat(:,1), HorizCat(:,2));
end % coo2xy

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

