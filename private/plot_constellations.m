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
