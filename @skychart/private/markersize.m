function m = markersize(mag)
  % markersize(mag): compute the marker size from the magnitude
  
  m = log(13-mag)*3+1;
  m(mag>=13) = 1;
  m = ceil(abs(m));
end % markersize
