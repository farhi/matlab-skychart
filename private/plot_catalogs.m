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
            & catalog.MAG(:) < mag_max);

    catalog.visible = visible; % store it as visible state
    
    mag = catalog.MAG(visible);
    x   = catalog.X(visible);
    y   = catalog.Y(visible);
    sz  = catalog.SIZE(visible);
    typ = catalog.TYPE(visible);
    
    % when out user-inpolygon, color is scaled down (grayed_out) / 2
    
    % scatter plots: o,ly the 'circle' can have non-uniform size
    switch f{1}
    case 'planets'
      % for planets: show name
      % For planets,   coloured disk (scatter) with size=Magnitude^2 colour 
      h = scatter(x,y, markersize(mag).^(2*factor), 'r', 'filled');
      
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
        SZ2(index) = sz(index)/60/90*p/max([ diff(xl) diff(yl) ]);  
      end
      SZ = max(SZ1, SZ2).^(2*factor);
      h = scatter(x,y, SZ, colour(typ,mag), 'o');
      
    otherwise
      h = scatter(x,y, markersize(mag), 'w');
    end
    set(h, 'Tag', [ 'SkyChart_' f{1} ], ...
      'ButtonDownFcn', @ButtonDownCallback, ...
      'UserData', f{1});
    uicm = uicontextmenu;
    uimenu(uicm, 'Label', 'NAME', 'Tag', [ 'SkyChart_MenuNAME_' f{1} ]);
    uimenu(uicm, 'Label', 'RADEC', 'Tag', [ 'SkyChart_MenuRADEC_' f{1} ]);
    uimenu(uicm, 'Label', 'ALTAZ', 'Tag', [ 'SkyChart_MenuALTAZ_' f{1} ]);
    uimenu(uicm, 'Label', 'MAGTYPE', 'Tag', [ 'SkyChart_MenuMAGTYPE_' f{1} ]);
    % uimenu(uicm, 'Label', 'Properties', 'Separator', 'on');
    % uimenu(uicm, 'Label', 'GOTO Now');
    % uimenu(uicm, 'Label', 'Add to Selection');
    set(h, 'UIContextMenu', uicm);
    handles = [ handles h ];
    
    catalogs.(f{1}) = catalog;
  end
  
end % plot_catalogs
