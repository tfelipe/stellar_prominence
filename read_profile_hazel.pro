


pro read_profile_hazel,file_name,wavelength,ii,qq,uu,vv

  
  OPENR,1,file_name
  
  READF,1,nlines
 
 
  wavelength=dblarr(nlines)
  ii=dblarr(nlines)  
  qq=dblarr(nlines)  
  uu=dblarr(nlines)    
  vv=dblarr(nlines)    
 
 
 
  for i=0,nlines-1 do begin
  
    READF,1,lambdai,ii_i,qq_i,uu_i,vv_i
    
    
    wavelength[i]=lambdai
    ii[i]=ii_i
    qq[i]=qq_i
    uu[i]=uu_i
    vv[i]=vv_i
    
    
  endfor
    
  CLOSE,1    
    
    
end