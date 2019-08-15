function plot_telescope(self)
  % plot_telescope: plot a marker at the telescope location
  
  % FOV for focal length and camera sensor:
  % FOV=sensor size / focal length in radians -> one for H and V
  % https://www.dpreview.com/forums/post/54555442
  %
  % pointer could be a sqare symbol with actual mean dimension (enevr know in which orientation)
  
  % first remove any previous pointer
    delete(findobj(self.figure, 'Tag','SkyChart_Pointer1'));
    delete(findobj(self.figure, 'Tag','SkyChart_Pointer2'));
    
  % get the scope location
  RA = []; DEC=[];
  if ~isempty(self.telescope) && isvalid(self.telescope)
    try
      getstatus(self.telescope);
      % get the current scope location
      RA = (self.telescope.ra.h   +self.telescope.ra.min/60)*15; % h -> deg
      DEC=  self.telescope.dec.deg+self.telescope.dec.min/60;
    catch ME
      disp(getReport(ME))
    end
  end
  
  % plot
  if ~isempty(RA)
    % compute Alt-Az and stereographic polar coords
    [Az, Alt] = radec2altaz(RA, DEC, self.julianday, self.place);
    [X, Y]    = pr_stereographic_polar(Az+90, Alt);
    
    % the plot the pointer at scope location (cross + circle), 0.5 deg
    plot(X,Y, 'ro', 'MarkerSize', 20, 'Tag','SkyChart_Pointer1'); 
    plot(X,Y, 'r+', 'MarkerSize', 20, 'Tag','SkyChart_Pointer2');
  end
end % plot_telescope
