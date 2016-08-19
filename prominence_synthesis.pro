



pro prominence_synthesis



sintesis=1	      ;switch on synthesis (1) 	
nt=11		      ;number of time steps
periodo=0.514d0       ;star period
time_end=0.514d0      ;temporal span of simulation (days)
t0=periodo/4.d0       ;time of prominence pass at long=0 


inclination=60.d0     ;star axis inclination [Degrees] 
latitude=30.d0        ;latitude of the prominence [Degrees] 
R_star_in=0.98d0      ;star radius [Solar radius] 
height_in=3.0d0       ;prominence height from the surface of the star [Solar radius]
		       

cloud_dlong_in=13.d0	;Size of the cloud in longitude [Degrees]
cloud_dlat_in=13.d0    	;Size of the cloud in latitude [Degrees]

nlat=5.d0	 	;Number of points for sampling the cloud in latitude
nlong=5.d0	  	;Number of points for sampling the cloud in longitude


b_prominence=10.        ;Prominence magnetic field strength [G]
theta_prominence=90.	;Prominence magnetic field inclination in local reference frame [Degrees]
chi_prominence=0.	;Prominence magnetic field azimuth in local reference frame (Degrees)

outname='ABDor_He10830_B10Gchi0final'        ;Header output file

boundaryStokes='4.098e-5, 0.0, 0.0, 0.0'   ;Boundary conditions for the Stokes vector used in the solution of the radiative transfer equation (/HAZEL/IDL_routines/solar_field.pro)


multiplet='10830'        		;5876 or 10830

lambda_end=30.		;Wavelength range and sampling 	
lambda_start=-30.
nl=800




Temperature=10000.	;Temperature [k] 
m_He=4.002602		;He atomic mass
uu=1.660538921e-27	;atomic mass unit [kg]
Kb=1.38064852e-23	;Boltzmann constant [IS]


v_therm=sqrt(2.*kb*Temperature/(uu*m_He))*1.e-3	;Thermal Doppler velocity




;Output IDL save file
outname_sav=outname+'_i'+strcompress(string(inclination),/remove_all)+'_lat'+strcompress(string(latitude),/remove_all)+'_Rstar'+strcompress(string(R_star_in),/remove_all)+'_per'+strcompress(string(periodo),/remove_all)+'_height'+strcompress(string(height_in),/remove_all)+'_dlong'+strcompress(string(round(cloud_dlong_in)),/remove_all)+'_dlat'+strcompress(string(round(cloud_dlat_in)),/remove_all)+'.sav'


;Rotational velocity
omega=2.d0*!pi/periodo  		;radian/day
omega_s=2.d0*!pi/(periodo*24.*3600.)   	;radian/s


;height

UA=149597870700.d0*1.e-6   ;Mm
arcsecMm=tan(1./60.d0/60.d0*!pi/180.d0)*UA   ;1 arcsec [Mm] 

R_sun=696.d0  ;Mm
R_star=R_sun*R_star_in

height_Mm=(height_in)*R_star        ;prominence height [Mm]
height=height_Mm/arcsecMm           ;prominence height [arcsec at 1 Astronomical Unit] 


        
      
        
;Cloud size in radians        
cloud_dlong=cloud_dlong_in*!pi/180.d0  ;rad
cloud_dlat=cloud_dlat_in*!pi/180.d0    ;rad
        
        
;Some arrays        
time=findgen(nt)/(nt-1)*time_end    
      
x_cloud=dblarr(nlong,nlat)
y_cloud=dblarr(nlong,nlat)

long_cloud=x_cloud
lat_cloud=x_cloud   
   
        
;Write HAZEL control file
        
file_name='sintesis.ini'        

file_obs='"input.prof"'
file_inv_profiles='"output.prof"'
file_inv_parameters='"output.parameters"'
file_error_inv='"output.errors"'        
      
Wavelength_axis=strcompress(string(lambda_end),/remove_all)+','+strcompress(string(lambda_start),/remove_all)+','+strcompress(string(nl),/remove_all)   
      
        
@parameters_control_file




 
;Define stellar disk 
 
nx=300
ny=300

pixel_sample=dblarr(nx,ny)
Stokes_cloud=dblarr(nlong,nlat,4,nl)

xx=(findgen(nx)/(nx-1)*2-1.)*height_Mm*1.3
yy=(findgen(ny)/(ny-1)*2-1.)*height_Mm*1.3


rad=dblarr(nx,ny)

for ix=0,nx-1 do for iy=0,ny-1 do rad[ix,iy]=sqrt(xx[ix]^2+yy[iy]^2)

pixel_sample[where(rad le R_star)]=1




;Array with the central position of the cloud for each time step
xpath=dblarr(nt)
ypath=xpath


;Radian
  latitude=latitude*!pi/180.   
  inclination=inclination*!pi/180.

  
Stokes_cloud_star_av=dblarr(4,nl,nt)  
  

;Temporal loop

for it=0,nt-1 do begin

  print,'it=', it
  
  
  longitude=omega*(time[it]-t0)  ;Radian
    
      
;Local reference frame to inclined local reference frame

  latitude_los=asin(sin(latitude)*cos(!pi/2-inclination)-cos(latitude)*sin(!pi/2-inclination)*sin(longitude+!pi/2))
  longitude_los=acos(cos(longitude+!pi/2)*cos(latitude)/cos(latitude_los))-!pi/2.    

  if longitude gt !pi/2.d0 then longitude_los=!pi-longitude_los     
     
      
  
  xpos_center=sin(!pi/2.d0-latitude_los)*sin(longitude_los)*height_Mm
  ypos_center=height_Mm*cos(!pi/2.d0-latitude_los) 

    
  xpath[it]=xpos_center
  ypath[it]=ypos_center  
  

  ;Sampling in longitude loop
  
  for ilong=0,nlong-1 do begin
    
    
    longitude_ilong=longitude+cloud_dlong/(nlong-1)*ilong-cloud_dlong/2.d0    
    if longitude_ilong ge 3.d0*!pi/2.d0 then longitude_ilong=longitude_ilong-2.d0*!pi
    if longitude_ilong lt -!pi/2.d0 then longitude_ilong=longitude_ilong+2.d0*!pi    

  ;Sampling in latitude loop  
  
    for ilat=0,nlat-1 do begin
  
      latitude_ilat=latitude+cloud_dlat/(nlat-1)*ilat-cloud_dlat/2. 

  
      latitude_los=asin(sin(latitude_ilat)*cos(!pi/2.d0-inclination)-cos(latitude_ilat)*sin(!pi/2.d0-inclination)*sin(longitude_ilong+!pi/2.d0))
      longitude_los=acos(cos(longitude_ilong+!pi/2.d0)*cos(latitude_ilat)/cos(latitude_los))-!pi/2.d0    
      if longitude_ilong gt !pi/2. then longitude_los=!pi-longitude_los    
      


  
  
;Position on the plane of sky  
      xpos=sin(!pi/2.-latitude_los)*sin(longitude_los)*height_Mm
      ypos=height_Mm*cos(!pi/2.-latitude_los) 

    
      x_cloud[ilong,ilat]=xpos
      y_cloud[ilong,ilat]=ypos   
 
      R_distance=sqrt(xpos^2+ypos^2)
      
;Boundary condition for incident illumination:

;0 if prominence is out of the stellar disk)
 
      if (R_distance gt R_star) then boundary_condition='0.0, 0.0, 0.0, 0.0'
      
;Solar radiation field if prominence is on the stellar disk

      if (R_distance le R_star) then boundary_condition= boundaryStokes

 
      long_cloud[ilong,ilat]=longitude_los
      lat_cloud[ilong,ilat]=latitude_los    
 
    
      kk=min(abs(xx-xpos),xpos_px)  
      kk=min(abs(yy-ypos),ypos_px)
  



;Angle between local vertical and LOS
      theta=acos(cos(longitude_los)*sin(!pi/2.-latitude_los))

      
;Angle between x_sl^i and plane formed by LOS and local vertical
chi=atan(sin(latitude_los)*cos(longitude_los)/sin(longitude_los))

             
      
;Position on the plane of sky    
      xpos=sin(!pi/2.-latitude_los)*sin(longitude_los)*height_Mm
      ypos=height_Mm*cos(!pi/2.-latitude_los)       
      
      
      
;gamma (defines the arbitrary positive reference direction for Stokes Q)        
      
     gamma=atan(xpos/ypos) 
  
      
           

;radian to degrees    
    theta=theta*180./!pi
    chi=chi*180./!pi  
    gamma=gamma*180./!pi    
         
      
      
      
;Rotational velocity in the LOS (vmac1_real)

  r_centro=height_Mm*1.e6*cos(latitude_ilat)   ;m

  vlineal=omega_s*r_centro   ;m/s      
     
            
  ux_rot=cos(longitude_ilong)
  uy_rot=-sin(latitude_ilat)*sin(longitude_ilong)           
  uz_rot=cos(latitude_ilat)*sin(longitude_ilong)         
      
;Rotation matrix                     
  Rxx=cos(!pi/2.d0-inclination)+ux_rot^2*(1.d0-cos(!pi/2.d0-inclination))
  Ryx=uy_rot*ux_rot*(1.d0-cos(!pi/2.d0-inclination))+uz_rot*sin(!pi/2.d0-inclination)
  Rzx=uz_rot*ux_rot*(1.d0-cos(!pi/2.d0-inclination))-uy_rot*sin(!pi/2.d0-inclination)

  vlineal_x2=Rxx*vlineal
  vlineal_y2=Ryx*vlineal
  vlineal_z2=Rzx*vlineal

  v_los=-sin(longitude_ilong)*vlineal_x2-sin(latitude_ilat)*cos(longitude_ilong)*vlineal_y2+cos(latitude_ilat)*cos(longitude_ilong)*vlineal_z2

  vmac1_real=-v_los*1.e-3   ;positivo para redshift, km/s  
  
  
      
      
;Magnetic field in inclined local reference frame

if b_prominence ne 0 then begin


  bx=b_prominence*sin(theta_prominence*!pi/180.)*cos(chi_prominence*!pi/180.)
  by=b_prominence*sin(theta_prominence*!pi/180.)*sin(chi_prominence*!pi/180.)
  bz=b_prominence*cos(theta_prominence*!pi/180.d0)


  bh=sqrt(bx^2+by^2)     
      
;Rotation matrix    
  Ryy=cos(!pi/2.d0-inclination)+uy_rot^2*(1.d0-cos(!pi/2.d0-inclination))
  Rzz=cos(!pi/2.d0-inclination)+uz_rot^2*(1.d0-cos(!pi/2.d0-inclination))

  Rxy=uy_rot*ux_rot*(1.d0-cos(!pi/2.d0-inclination))-uz_rot*sin(!pi/2.d0-inclination)
  Rxz=uz_rot*ux_rot*(1.d0-cos(!pi/2.d0-inclination))+uy_rot*sin(!pi/2.d0-inclination)      
      
  Rzy=uz_rot*uy_rot*(1.d0-cos(!pi/2.d0-inclination))+ux_rot*sin(!pi/2.d0-inclination)      
  Ryz=uz_rot*uy_rot*(1.d0-cos(!pi/2.d0-inclination))-ux_rot*sin(!pi/2.d0-inclination)         


  bx_rot=Rxx*bx+Rxy*by+Rxz*bz      
  by_rot=Ryx*bx+Ryy*by+Ryz*bz           
  bz_rot=Rzx*bx+Rzy*by+Rzz*bz           
     
  if abs(bx_rot/b_prominence) le 1.e-7 then bx_rot=0.    
  if abs(by_rot/b_prominence) le 1.e-7 then by_rot=0.       
     
 
 
  bx_rot0=-sin(longitude_ilong)*bx_rot-sin(latitude_ilat)*cos(longitude_ilong)*by_rot+cos(latitude_ilat)*cos(longitude_ilong)*bz_rot
  by_rot0=cos(longitude_ilong)*bx_rot-sin(latitude_ilat)*sin(longitude_ilong)*by_rot+cos(latitude_ilat)*sin(longitude_ilong)*bz_rot
  bz_rot0=cos(latitude_ilat)*by_rot+sin(latitude_ilat)*bz_rot


;Correct ambiguity
  if inclination lt 45.*!pi/180. and longitude_ilong gt !pi/2. then longitude_los=!pi-longitude_los   


  bx2=-sin(longitude_los)*bx_rot0+cos(longitude_los)*by_rot0
  by2=-sin(latitude_los)*cos(longitude_los)*bx_rot0-sin(latitude_los)*sin(longitude_los)*by_rot0+cos(latitude_los)*bz_rot0
  bz2=cos(latitude_los)*cos(longitude_los)*bx_rot0+cos(latitude_los)*sin(longitude_los)*by_rot0+sin(latitude_los)*bz_rot0

  if abs(bx2/b_prominence) le 1.e-7 then bx2=0.    
  if abs(by2/b_prominence) le 1.e-7 then by2=0.       
  if abs(bz2/b_prominence) le 1.e-7 then bz2=0.       
 

  bh2=sqrt(bx2^2+by2^2)     

  theta_prominence_rot=atan(bh2/bz2)*180./!pi
  chi_prominence_rot=atan(by2/bx2)*180./!pi       
    
  if bx2 lt 0 then chi_prominence_rot=chi_prominence_rot+180.
  if bz2 lt 0 then theta_prominence_rot=theta_prominence_rot+180.
    
    
       
  thetaB1=strcompress(string(theta_prominence_rot),/remove_all)
  chiB1=strcompress(string(chi_prominence_rot),/remove_all)
      
endif       
      

;Write HAZEL control file

      LOS_angles=strcompress(string(theta),/remove_all)+','+strcompress(string(chi),/remove_all)+','+strcompress(string(gamma),/remove_all)

      vmac1=strcompress(string(vmac1_real),/remove_all)
      
      write_controlfile_hazel,file_name,input_model,file_obs,file_inv_profiles,file_inv_parameters,file_error_inv,action,verbose,linear_solver,stopping_volume,synthesis_mode,stimulated_emission,magnetic_field,Paschen_Back,atomic_polarization,magneto_optical,stimulated_emission_RT,multiplet,LOS_angles,Wavelength_axis,number_slabs,boundary_condition,aa,height,beta,ff,BB1,thetaB1,chiB1,vdopp1,tau1,vmac1,BB2,thetaB2,chiB2,vdopp2,tau2,vmac2,aa_ranges,beta_ranges,ff_ranges,BB1_ranges,thetaB1_ranges,chiB1_ranges,vdopp1_ranges,tau1_ranges,vmac1_ranges,BB2_ranges,thetaB2_ranges,chiB2_ranges,vdopp2_ranges,tau2_ranges,vmac2_ranges,iterationsLM,ncycles,inversion_modes,aa_cycles,beta_cycles,ff_cycles,BB1_cycles,thetaB1_cycles,chiB1_cycles,vdopp1_cycles,tau1_cycles,vmac1_cycles,BB2_cycles,thetaB2_cycles,chiB2_cycles,vdopp2_cycles,tau2_cycles,vmac2_cycles,Weights_I,Weights_Q,Weights_U,Weights_V


        
;Sintesis 

;If prominence point is not covered by the stellar disk, do synthesis

if ((sintesis eq 1) and (long_cloud[ilong,ilat] le !pi/2.) and (long_cloud[ilong,ilat] ge -!pi/2.)) or ((sintesis eq 1) and (sqrt(x_cloud[ilong,ilat]^2+y_cloud[ilong,ilat]^2) ge R_star)) then begin


      spawn,'run.py sintesis.ini > log.txt'

             
             
;Read profiles
  read_profile_hazel,'output.prof',lambda_sint,i_sint,q_sint,u_sint,v_sint
          
    

      Stokes_cloud[ilong,ilat,0,*]=i_sint
      Stokes_cloud[ilong,ilat,1,*]=q_sint     
      Stokes_cloud[ilong,ilat,2,*]=u_sint
      Stokes_cloud[ilong,ilat,3,*]=v_sint      
  
   
endif   




endfor
endfor



;Average profiles (weighted by the area of each point)


    mitadx=x_cloud*0    
    mitadx[1:nlong-1,*]=abs(x_cloud[1:nlong-1,*]-x_cloud[0:nlong-2,*])/2.

    ladox=x_cloud*0.
    ladox[0:nlong-2,*]=mitadx[0:nlong-2,*]+mitadx[1:nlong-1,*]
    
    ladox[0,*]=2.*mitadx[1,*]
    ladox[nlong-1,*]=2.*mitadx[nlong-1,*]    



    mitady=y_cloud*0    
    mitady[*,1:nlat-1]=abs(y_cloud[*,1:nlat-1]-y_cloud[*,0:nlat-2])/2.

    ladoy=y_cloud*0.
    ladoy[*,0:nlat-2]=mitady[*,0:nlat-2]+mitady[*,1:nlat-1]
    
    ladoy[*,0]=2.*mitady[*,1]
    ladoy[*,nlat-1]=2.*mitady[*,nlat-1]    





;Area of the cloud covering the stellar disk  
  
    R_distance=sqrt(x_cloud^2+y_cloud^2)    
    
    area_eclipse=0.
    for ilong=0,nlong-1 do for ilat=0,nlat-1 do begin
    
      if ((R_distance[ilong,ilat] le R_star) and (long_cloud[ilong,ilat] lt !pi/2.) and (long_cloud[ilong,ilat] gt -!pi/2.)) then area_eclipse=area_eclipse+ladox[ilong,ilat]*ladoy[ilong,ilat]
    
    endfor


    area=0.
    
    for ilong=0,nlong-1 do for ilat=0,nlat-1 do begin   
      if ((long_cloud[ilong,ilat] le !pi/2.) and (long_cloud[ilong,ilat] ge -!pi/2.)) or (R_distance[ilong,ilat] ge R_star) then area=area+ladoy[ilong,ilat]*ladox[ilong,ilat]
    endfor
      
    
    print,'area=',area/!pi/R_star^2,' A_star'
    
    Stokes_cloud_av=dblarr(4,nl)
    for ilong=0,nlong-1 do for ilat=0,nlat-1 do begin
      
      if ((long_cloud[ilong,ilat] le !pi/2.) and (long_cloud[ilong,ilat] ge -!pi/2.)) or (R_distance[ilong,ilat] ge R_star) then Stokes_cloud_av=Stokes_cloud_av+Stokes_cloud[ilong,ilat,*,*]*ladox[ilong,ilat]*ladoy[ilong,ilat]
    
    endfor
    
    
    if area eq 0 then area=1   
    
    Stokes_cloud_av=reform(Stokes_cloud_av)/area
  
        
    print,'area_eclipse=',area_eclipse/!pi/R_star^2,' A_star'    
    
    

    
    continuo=dblarr(4,nl)
    continuo[0,*]=1.
    

      
    Stokes_cloud_star_av[*,*,it]=(Stokes_cloud_av*area+(!pi*R_star^2-area_eclipse)*continuo)/height_Mm^2
    
            
        
    
;plots    

  !p.multi=[0,3,2]  

  tvframe,pixel_sample,xrange=[xx[0],xx[nx-1]],yrange=[yy[0],yy[ny-1]],/aspect
  oplot,[0,0],[R_star*sin(inclination),R_star*sin(inclination)],psym=4,color=0  
  oplot,[xpos_center,xpos_center],[ypos_center,ypos_center],psym=2,color=200
  oplot,xpath[0:it],ypath[0:it],color=100,line=1
  
  
  x_cloud_plot=x_cloud
  tapado=where(long_cloud gt !pi/2  and R_distance le R_star or long_cloud lt -!pi/2 and R_distance le R_star)
  if tapado[0] ne -1 then x_cloud_plot[tapado]=10.*height_Mm
  
  y_cloud_plot=y_cloud
  if tapado[0] ne -1 then y_cloud_plot[tapado]=10.*height_Mm
  
  oplot,x_cloud_plot,y_cloud_plot,psym=4,color=200
  
  
  

 
 !x.charsize=2
 !y.charsize=2   
 
    plot,lambda_sint,Stokes_cloud_star_av[0,*,it],title='I',xtitle='!4Dk!3 (!3'+string(197B)+'!X)',thick=2.5,yrange=[min(Stokes_cloud_star_av[0,*,it])-0.1,max(abs(Stokes_cloud_star_av[0,*,it]))+0.05],/yst
            
    plot,lambda_sint,Stokes_cloud_star_av[1,*,it],title='Q',xtitle='!4Dk!3 (!3'+string(197B)+'!X)',thick=2.5,yrange=[-max(abs(Stokes_cloud_star_av[1,*,it])),max(abs(Stokes_cloud_star_av[1,*,it]))]
    
    
    !p.multi=[2,3,2]
    
    plot,lambda_sint,Stokes_cloud_star_av[2,*,it],title='U',xtitle='!4Dk!3 (!3'+string(197B)+'!X)',thick=2.5,yrange=[-max(abs(Stokes_cloud_star_av[2,*,it])),max(abs(Stokes_cloud_star_av[2,*,it]))]

    plot,lambda_sint,Stokes_cloud_star_av[3,*,it],title='V',xtitle='!4Dk!3 (!3'+string(197B)+'!X)',thick=2.5,yrange=[-max(abs(Stokes_cloud_star_av[3,*,it])),max(abs(Stokes_cloud_star_av[3,*,it]))]
    
    
    
    

endfor


;Normalize to continuum intensity
Stokes_cloud_star_av=Stokes_cloud_star_av/Stokes_cloud_star_av[0,0,0]




save,filename=outname_sav,Stokes_cloud_star_av,lambda_sint,time,nt,periodo,inclination,latitude,R_star_in,height_in,cloud_dlong_in,cloud_dlat_in,nlat,nlong,t0,time_end,b_prominence,theta_prominence,chi_prominence










stop
end