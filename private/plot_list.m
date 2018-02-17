function plot_list(self)
  % plot_list: plot a marker on list/selected object
  
  % clean up previous selection
  delete(findobj(self.figure, 'Tag','SkyChart_Selection'));
  
  % get current selection
  toplot = [ self.selected self.list ];
  if isempty(toplot), return; end
  
  RA = [ toplot.RA ];
  DEC= [ toplot.DEC ];
  
  if ~isempty(RA)
    % compute Alt-Az and stereographic polar coords
    [Az, Alt] = radec2altaz(RA, DEC, self.julianday, self.place);
    [X, Y]    = pr_stereographic_polar(Az+90, Alt);
    % first remove any previous pointer
    
    % the plot the pointer at selection location
    plot(X,Y, 'wx', 'MarkerSize', 15, 'Tag','SkyChart_Selection'); 
  end
end % plot_telescope
