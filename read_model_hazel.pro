


pro read_model_hazel,file_name,BB,thetaB,chiB,v_doppler,beta,tau,macro,damping,height


BB=0.d0
thetaB=0.d0
chiB=0.d0
beta=0.d0
tau=0.d0
macro=0.d0
damping=0.d0
height=0.d0
v_doppler=0.d0

kk='a'
number=0.d0
  
  OPENR,1,file_name
  
    READF,1,kk

    
    
    READF,1,kk
    READF,1,kk 
        
    READF,1,kk
            
    READF,1,kk
       
    READF,1,kk
     
    READF,1,kk
     
    READF,1,kk
         
    READF,1,kk
            
    READF,1,kk
           
    READF,1,kk
     
    READF,1,kk
     
    READF,1,kk
         
    READF,1,kk
            
    READF,1,kk
           
    READF,1,kk
      
    READF,1,kk
            
    READF,1,kk
           
    READF,1,kk
        
    READF,1,kk
           
    READF,1,kk
      
    READF,1,kk
            
    READF,1,kk
           

    
    READF,1,BB,thetaB,chiB
    READF,1,kk     
    READF,1,kk     

    READF,1,height
    
    READF,1,kk     
    READF,1,kk 
    
    READF,1,tau    
    
    READF,1,kk     
    READF,1,kk 
    READF,1,beta    
    
    READF,1,kk      
    READF,1,kk
    READF,1,kk
    READF,1,kk    
    READF,1,kk       
    READF,1,kk  
    READF,1,kk
    READF,1,kk
    READF,1,kk    
    READF,1,kk       
    READF,1,kk      
    READF,1,kk  
    READF,1,kk
        
    READF,1,kk
        
    READF,1,kk
        
    READF,1,kk     
        
    READF,1,kk 
    
    
    READF,1,number,macro,damping
    
    READF,1,kk
    READF,1,kk   
    READF,1,v_doppler
        
    
  CLOSE,1    
    
    
end