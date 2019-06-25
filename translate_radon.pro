pro translate_radon,input,output,shiftpix
  ;translates a radon transform image along the theta axis. It wraps/flips at the edges
  ;shift is in pixel units, not theta (could try fixing this later, but interpolation would be needed)
  sz=size(input)
  sel=where(input ne -1e-10)
  input_indices = array_indices(input,sel) ;2xn_element array
  new_indices = input_indices
  new_indices[0,*]=new_indices[0,*]+shiftpix
  sel=where(new_indices[0,*] gt sz[1]-1,count)
  if count ne 0 then begin
    new_indices[0,sel]=new_indices[0,sel]-sz[1]
    new_indices[1,sel]=(sz[2]-1)-new_indices[1,sel]
  endif
  
  sel=where(new_indices[0,*] lt 0,count)
  if count ne 0 then begin
      new_indices[0,sel]=new_indices[0,sel]+sz[1]
      new_indices[1,sel]=(sz[2]-1)-new_indices[1,sel]
  endif
  output=input*0
  output[new_indices[0,*],new_indices[1,*]] = input[input_indices[0,*],input_indices[1,*]]
end
