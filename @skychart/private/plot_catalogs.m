function self = plot_catalogs(self)
  % plot_catalogs: plot all visible objects from the catalog
  %
  
  % we shall scale the symbols as well to the figure size
  p = get(self.figure, 'Position'); p=mean(p(3:4));
  sz_fact = p/1024;
  
  for f=fieldnames(self.catalogs)'
    h = [];

    catalog = self.catalogs.(f{1});
    if ~isfield(catalog, 'MAG'), continue; end
    
    % limit magnitude to show, so that we only have 1000 objects max
    mag_max = abs(2.5*log10(abs(4/diff(self.xlim)/diff(self.ylim))))+5;
    
    visible = (catalog.X.^2+catalog.Y.^2 < 1 ...
            & min(self.xlim) < catalog.X & catalog.X < max(self.xlim) ...
            & min(self.ylim) < catalog.Y & catalog.Y < max(self.ylim) ...
            & catalog.MAG(:) < mag_max);

    catalog.visible = visible; % store it as visible state
    
    mag = catalog.MAG(visible);
    x   = catalog.X(visible);
    y   = catalog.Y(visible);
    sz  = catalog.SIZE(visible);
    typ = catalog.TYPE(visible);
    nam = catalog.NAME(visible);
    
    % remove previous plot and labels
    h = findall(self.figure, 'Tag', [ 'SkyChart_' f{1} ]);
    if ~isempty(h) delete(h); end
    h = findall(self.figure, 'Tag', [ 'SkyChart_Labels_' f{1} ]);
    if ~isempty(h) delete(h); end
    
    % scatter plots: o,ly the 'circle' can have non-uniform size
    switch f{1}
    case 'planets'
      % for planets: show name
      % For planets,   coloured disk (scatter) with size=Magnitude^2 colour 
      h = scatter(x,y, markersize(mag).^(2*sz_fact), 'r', 'filled');
      
    case 'stars'
      % for star, when proper name, show it
      % For stars, use coloured (scatter) with size=Magnitude^2 colour from TYPE
      SZ       = 8-mag+min(mag); % from 10 down ...
      SZ(SZ<1) = 1;
      h = scatter(x,y, SZ.^2, colour(typ,mag), 'filled');
      
    case 'deep_sky_objects'
      % for DSO, when proper name, show it
      % For dso,   use circle        (scatter) with thickness=Magnitude, and given size
      
      % the size of the screen is max(diff(xl) diff(yl)) in 0-1/pixels
      SZ1 = markersize(mag); SZ2=zeros(size(sz));
      index = find(sz>0);
      if ~isempty(index)
        % in [deg] then in 0-1 then in pixels
        SZ2(index) = sz(index)/60/90*p/max([ diff(self.xlim) diff(self.ylim) ]);  
      end
      SZ = max(SZ1(:).^(2*sz_fact), SZ2(:).^2);
      h = scatter(x,y, SZ, colour(typ,mag), 'o');
      mag(mag == 0) = 20; % disable show-names for the objects showing regions (mag=0)
      
    otherwise
      h = scatter(x,y, markersize(mag), 'w');
    end

    % add mouse handling to get closest object and set menu labels
    set(h, 'Tag', [ 'SkyChart_' f{1} ], ...
      'ButtonDownFcn', @ButtonDownCallback, ...
      'UserData', f{1});
      
    uicm = uicontextmenu('Parent', self.figure);
    uimenu(uicm, 'Label', 'NAME', 'Tag', [ 'SkyChart_MenuNAME_' f{1} ]);
    uimenu(uicm, 'Label', 'RADEC', 'Tag', [ 'SkyChart_MenuRADEC_' f{1} ]);
    uimenu(uicm, 'Label', 'ALTAZ', 'Tag', [ 'SkyChart_MenuALTAZ_' f{1} ]);
    uimenu(uicm, 'Label', 'MAGTYPE', 'Tag', [ 'SkyChart_MenuMAGTYPE_' f{1} ]);
    uimenu(uicm, 'Label', 'Send Scope here Now', 'Callback', @MenuCallback, ...
      'Tag','SkyChart_MenuGOTO', 'Separator','on');  
    uimenu(uicm, 'Label', 'Add to List', 'Callback', @MenuCallback);  
    set(h, 'UIContextMenu', uicm);
    
    % display the brightest objects name
    [~, index] = sort(mag(:)); % sort magnitude in ascending order
    h = [];
    if numel(index) > 20, index=index(1:20); end
    for n=index(:)'
      name = nam{n};
      t = text(x(n)+diff(self.xlim)*.01, y(n), strtok(name,';'));
      h = [ h t ];
      set(t, 'Color','w');
    end
    set(h, 'Tag', [ 'SkyChart_Labels_' f{1} ]);
    
    self.catalogs.(f{1}) = catalog;
  end
  
end % plot_catalogs

function ButtonDownCallback(src, evnt)
  % ButtonDownCallback: callback when user clicks on an object

  % get the SkyChart object handle
  h = findall(0, 'Tag','SkyChart_Axes');
  if numel(h) > 1, delete(h(2:end)); h=h(1); end
  if isempty(h) || ~ishandle(h)
    return
  end
  self=get(h, 'UserData');
  
  % where the mouse click is
  xy = get(self.axes, 'CurrentPoint'); 
  x = xy(1,1); y = xy(1,2);
  
  % search for closest object in the corresponding catalogs
  found.dist=inf; found.catalog=''; ; found.index=[];
  % we try with the catalog for the given clicked handle, else try all
  f = get(src, 'UserData');
  if any(strcmp(f, fieldnames(self.catalogs))), catalogs = { f };
  else catalogs = fieldnames(self.catalogs); end
  % search in catalogs for closest object
  for f=catalogs(:)'
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
      found.RA      = catalog.RA(found.index);
      found.DEC     = catalog.DEC(found.index);
      found.Alt     = catalog.Alt(found.index);
      found.Az      = catalog.Az(found.index);
      found.X       = catalog.X(found.index);
      found.Y       = catalog.Y(found.index);
      found.MAG     = catalog.MAG(found.index);
      found.TYPE    = catalog.TYPE{found.index};
      found.NAME    = catalog.NAME{found.index};
    end
  end
  % display our findings
  labels = { ...
    'RADEC',   sprintf('RA=%f DEC=%f', found.RA, found.DEC);
    'ALTAZ',   sprintf('Az=%f Alt=%f', found.Az, found.Alt);
    'MAGTYPE', sprintf('%s MAG=%f TYPE=%s', found.catalog, found.MAG, found.TYPE);
    'NAME',    found.NAME };

  for index=1:size(labels,1)
    uicm = findall(self.figure, 'Tag', ...
      [ 'SkyChart_Menu' labels{index,1} '_' found.catalog ]);
    if isempty(uicm)
      disp([ mfilename ': not found ' 'SkyChart_Menu' labels{index,1} '_' found.catalog ])
    end
    set(uicm, 'Label', labels{index,2});
  end
  self.selected = found;
  disp([ 'skychart: Selected ' labels{1,2} ' "' found.NAME '"' ]);
  title(self.axes, [ 'Selected ' found.NAME '"' ]);
  
end % ButtonDownCallback

function MenuCallback(src, evnt)
  % MenuCallback: callback from UIContext Menu (right click on object/star)
  
  % get the last selected target
  h = findall(0, 'Tag','SkyChart_Axes');
  if numel(h) > 1, delete(h(2:end)); h=h(1); end
  if isempty(h) || ~ishandle(h)
    return
  end
  self=get(h, 'UserData');
  
  switch lower(get(src, 'Label'))
  case 'send scope here now'
    goto(self);
  case 'add to list'
    listAdd(self);
  end
  
end % MenuCallback
