function c = colour(typ, mag)
  % colour: determine the colour of objects for scatter3
  c = ones(numel(typ),3);  % initialise to white
  
  tokens = { 'star O',  [ 0 0   1 ]; ...
             'star B',  [ 0 0.5 1 ]; ...
             'star A',  [ 0 1   1 ]; ...
             'star F',  [ 0 1   0 ]; ...
             'star G',  [ 1 1   0 ]; ...
             'star K',  [ 1 0.5 0 ]; ...
             'star M',  [ 1 0   0 ]; ...
             'DSO C',   [ 1 0 0 ]; ...
             'DSO D',   [ 1 0 0 ]; ...
             'DSO E',   [ 1 0 0 ]; ...
             'DSO I',   [ 1 0 0 ]; ...
             'DSO P',   [ 1 0 0 ]; ...
             'DSO G',   [ 1 0 0 ]; ...
             'DSO S',   [ 1 0 0 ]; ...
             'DSO OCL', [ 0 0 1 ]; ...
             'DSO GCL', [ 0 1 1 ]; ...
             'DSO DN',  [ 0 1 0 ]; ...
             'DSO EN',  [ 0 1 0 ]; ...
             'DSO RN',  [ 0 1 0 ]; ...
             'DSO PN',  [ 0 1 0 ]; ...
             'DSO *',   [ 0 0 0 ]; ...
             'DSO NF',  [ 0 0 0 ]};
  mag1 = 1-(mag-min(mag))*.1;
  mag1(mag1 < .5 | mag==0) = .5;

  for index=1:size(tokens, 1)
    tok = tokens{index, 1};
    col = tokens{index, 2};
    ok  = strncmp(typ, tok, numel(tok));
    c(ok,1) = col(1).*mag1(ok); 
    c(ok,2) = col(2).*mag1(ok); 
    c(ok,3) = col(3).*mag1(ok);
  end
end % colour
