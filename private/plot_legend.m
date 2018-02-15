function legend_h = plot_legend(fig)
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
  set(0, 'CurrentFigure', fig); % select but not raise
  hold on
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
  set(legend_h, 'TextColor','w','Color','none','Tag','SkyChart_legend', ...
    'XColor','y', 'YColor','y');
 
  
end % plot_legend
