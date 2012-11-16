#INLINE F90_RATES
    REAL(kind=dp) FUNCTION extarr(type, var1,var2,var3,var4,var5, &
                      var6, var7,var8,var9,var10, var11)
         
        use dsmacc_global, only: lat, lon, jday, szas, svj_tj, bs,cs,ds
        use constants,  ONLY:zenang, seval
    
        CHARACTER(LEN=*) type
        REAL var1, var2, var3, var4, var5, var6, &
             var7, var8, var9, var10, var11, var12
        REAL(kind=dp) k0,ki
        REAL*8 tmp(19), tmp2(19),b(19), c(19), d(19), time2, lat2, theta
        INTEGER n,jl, i
        include 'tuv_old/params'
    
        extarr=0.0d0
        IF (type .eq. '0.00') then 
             ! Standard arhenius equation. No designation number so becomes '0.00'
             extarr=arr(var2, var3, -var4)
        ENDIF    
    
        IF (type .eq. 'P') then 
    ! standard pressure dependent reacions
            k0=arr(var2, var3, -var4)
            ki=arr(var9, var10, -var11)
            extarr=(k0*m/(1+k0*m/ki))*var6 ** ( &
                (1e0+(log10(k0*m/ki)) ** 2.0e0) ** -1.0e0)
        ENDIF    
    
        IF (type .eq. 'hv') then 
    ! photolysis rates
    ! this is copied from the constants.f90 bit of the mcm
            
    
        Time2=mod(Time/(60.*60.), 24.)
        IF (TIME2 .LT. 0) TIME2=TIME2+24.
            LAT2=LAT
            THETA=ZENANG(int(jday),Time2,LAT2)*180./PI
            
            IF (theta .le. 90) then     
                n=19 
        
                jl=-1
                IF (var1 .eq. 36.1) jl=2 
                IF (var1 .eq. 35.1) jl=4
                
        
                IF (jl .ne. -1) then 
                do i=1,n    
                    tmp(i)=szas(i)
                    tmp2(i)=svj_tj(i,jl)
        
                    b(i)=bs(i,jl)
                    c(i)=cs(i,jl)
                    d(i)=ds(i,jl)
                enddo
        
                    extarr=seval(n, theta, tmp, tmp2, b,c,d)
            ELSE
                    extarr=0.e0
            ENDIF 
                WRITE(OUTPUT_UNIT,*) 'hv',jl, var1, extarr
            ENDIF
                
        ENDIF            
    END FUNCTION extarr
#ENDINLINE
