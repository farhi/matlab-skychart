function plot_list(self)
  % plot_list: plot a marker on list/selected object
  
  % clean up previous selection
  delete(findobj(self.figure, 'Tag','SkyChart_Selection'));
  
  % get current selection
  RA = [];
  DEC= [];
  for index=0:numel(self.list)
    if index==0
      if ~isempty(self.selected)
        RA(end+1)  = self.selected.RA;
        DEC(end+1) = self.selected.DEC;
      end
    else
      RA(end+1)  = self.list(index).RA;
      DEC(end+1) = self.list(index).DEC;
    end
  end
  if isempty(RA), return; end
  
  if ~isempty(RA)
    % compute Alt-Az and stereographic polar coords
    [Az, Alt] = radec2altaz(RA, DEC, self.julianday, self.place);
    [X, Y]    = pr_stereographic_polar(Az+90, Alt);
    % first remove any previous pointer
    
    % the plot the pointer at selection location
    plot(X,Y, 'wx', 'MarkerSize', 15, 'Tag','SkyChart_Selection'); 
  end
end % plot_telescope
