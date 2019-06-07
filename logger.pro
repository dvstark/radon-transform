pro logger,logfile,message
  openw,1,logfile,/append
  printf,1,''
  printf,1,message
  printf,1,''
  close,1
end
