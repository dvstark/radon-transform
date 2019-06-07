pro unwrap_pa,rho_arr,theta_arr,srtarr=srt

  srt = sort(rho_arr)
  rho_arr = rho_arr[srt]
  theta_arr = theta_arr[srt]
  
  ind0 = (where(rho_arr eq 0))[0]
  theta_new = theta_arr
  
                                ;upper part
  theta_last = theta_arr[ind0]
  for m=ind0+1,n_elements(theta_arr)-1,1 do begin
     if theta_new[m] ne -1 then begin
        dtheta = theta_new[m]-theta_last
                                ;print,theta_last/!dtor,theta_new[m]/!dtor,dtheta/!dtor
        if dtheta gt !pi/2. then theta_new[m] = theta_new[m]-!pi
        if dtheta lt -!pi/2. then theta_new[m] = theta_new[m]+!pi
        theta_last = theta_new[m]
                                ;print,theta_new[m]/!dtor
     endif
  endfor
                                ;lower part
  theta_last = theta_arr[ind0]
  for m=ind0-1,0,-1 do begin
     if theta_new[m] ne -1 then begin
        dtheta = theta_new[m]-theta_last
        ;print,''
                               ; print,rho_arr[m],theta_last/!dtor,theta_new[m]/!dtor,dtheta/!dtor
        if dtheta gt !pi/2. then theta_new[m] = theta_new[m]-!pi
        if dtheta lt -!pi/2. then theta_new[m] = theta_new[m]+!pi
        theta_last = theta_new[m]
                                ;lprint,theta_new[m]/!dtor
     endif
  endfor

  theta_arr = theta_new
;  return,theta_new

end
