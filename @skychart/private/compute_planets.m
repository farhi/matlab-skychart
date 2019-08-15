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

  au = 149597870;     % km
  for index=1:numel(Labels)
    RA = ''; Coo=[];
    switch Labels{index}
    case 'Sun'
      [RA, Dec]  = suncoo(julianday,'a');
      Mag = -26; Dist = 149597870.7/au; %au
    case 'Moon'
      [RA, Dec] = mooncool(julianday,place,'b');
      Mag = -5; Dist = 384400/au; % au
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
