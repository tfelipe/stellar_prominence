


pro write_profile_hazel,file_name,wavelength,ii,qq,uu,vv,ii_n,qq_n,uu_n,vv_n



  nwavelength=n_elements(wavelength)
  
  FORMAT='(F10.8,5X,F10.8,5X,F10.8,5X,F10.8,5X,F10.8,5X,F10.8,5X,F10.8,5X,F10.8,5X,F10.8)' 
  
  OPENW,1,file_name, WIDTH=150
  
  PRINTF,1,nwavelength
  
  for i=0,nwavelength-1 do begin
  
    PRINTF,1,wavelength[i],ii[i],qq[i],uu[i],vv[i],ii_n[i],qq_n[i],uu_n[i],vv_n[i]
    
  endfor
    
  CLOSE,1    
    
    
end