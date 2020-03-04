      subroutine wrtcdf(var,nbdim,dims)
      
# include "netcdf.inc"

      include 'mpif.h'


      !local
      integer::ierr
      integer::iproc,i,j,k,kint,l,Nin,kcount,mpierr
      integer,save:: nbtime,timeid
      integer,save::dummyid
      integer,dimension(40),save:: dimid,varid
      character*64 :: name
      integer,dimension(MPI_STATUS_SIZE):: istat

      select case(iflag)
      case(0) 		! initialization (open netcdf file)
      if (mynode.eq.0) then

         
! Create the netCDF file. The NF_NETCDF4 flag causes a
! HDF5/netCDF-4 file to be created. The comm and info parameters
! cause parallel I/O to be enabled.

         specname = 'diag_spec.nc' 
         call check(nf_create(trim(specname), NF_CLOBBER, ncidspec)  ) 

         call check(nf_def_dim(ncidspec, 'kk',Nax,dimid(1)) )
         call check(nf_def_dim(ncidspec, 'level',Ninterp,dimid(2)) )
         nbtime = int((ntimes-ntstart)/nwrtspec)+1                  
         call check(nf_def_dim(ncidspec, 'time',NF_UNLIMITED,dimid(3)) )
! put global attributes

         call check(nf_put_att_int(ncidspec,NF_GLOBAL,'nx',NF_INT,1,(LLm)))
         call check(nf_put_att_int(ncidspec,NF_GLOBAL,'ny',NF_INT,1,(MMm)))        
         call check(nf_put_att_double(ncidspec,NF_GLOBAL,'dx',
     &              NF_DOUBLE,1,1./pm(1,1)))  
         call check(nf_put_att_double(ncidspec,NF_GLOBAL,'dy',
     &              NF_DOUBLE,1,1./pn(1,1)))
         call check(nf_put_att_double(ncidspec,NF_GLOBAL,'dz',
     &              NF_DOUBLE,1,Hz(1,1,1)))

! define variables

         call check(nf_def_var(ncidspec,'time' ,
     &                    nf_double,1,(/dimid(3)/),timeid))
         do l=1,nspec
            name=name_cospec(l)            
            call check(nf_def_var(ncidspec,trim(name),NF_DOUBLE,3,  
     &                       (/dimid(1),dimid(2),dimid(3)/),varid(l)))                        
         enddo        
    
         call check(nf_enddef(ncidspec) )
          
      endif 
              
      case(1) ! write spec

! each processor write it's own cospec in the parallel netcdf file if needed
          
          nrecspec=nrecspec+1
          
          cospec(:,:,10) = 0.5 * cospec(:,:,10) 
          cospec(:,:,18) = 0.5 * cospec(:,:,18) 
                   
     
          do l=11,nspec             
             do k=1,nint
                ktot=mynode*nint+k               
                if (ktot <= Ninterp) then             
                     cospec(:,k,l) = - g/drhordz_interp(ktot)/rho0*cospec(:,k,l)
                endif
             enddo
          enddo

! cospec divided because of averaging, duplicate domain EW, duplicate domain NS          

          cospec = cospec / float(nwrtspec/navgspec) / coefewper / coefnsper     

          if (mynode==0)  call check(nf_put_vara_double(ncidspec,timeid,  
     &                     (/nrecspec/),(/1/),time*sec2day) ) ! time in days      
          
          
          ! reads/writes to happen together on every processor.
          do l=1,nspec

             if (mynode==0) then
                   Nin=1
                   kcount=nint 
                   
                   call check(nf_put_vara_double(ncidspec, varid(l),  
     &                  (/1,Nin,nrecspec/),  
     &                  (/Nax,kcount,1/),
     &                  real(cospec(1:Nax,1:kcount,l))) )  
             else
                    Nin=mynode*nint+1
                    if (Nin <= Ninterp) then                    
!                      call MPI_SEND(cospec,Nax*nint*nspec,
                       call MPI_SEND(cospec(:,:,l),Nax*nint,                
     &                       MPI_DOUBLE_COMPLEX,0,mynode,  
     &                       MPI_COMM_WORLD,mpierr)  
                    endif
             endif                     
               
             if (mynode==0) then
                  do iproc=1,NNODES-1
                     Nin=iproc*nint+1   
                     if (Nin <= Ninterp) then                    
!                       call MPI_IRECV(cospec,  
!    &                      Nax*nint*nspec, 
                        call MPI_IRECV(cospec(:,:,l),  
     &                      Nax*nint,  
     &                      MPI_DOUBLE_COMPLEX,iproc,iproc,  
     &                      MPI_COMM_WORLD,reqid(iproc),mpierr)
                        call MPI_WAIT(reqid(iproc),istat,mpierr)
                        kcount = min(nint,Ninterp-Nin+1)                         
                        call check(nf_put_vara_double(ncidspec, varid(l),  
     &                     (/1,Nin,nrecspec/),  
     &                     (/Nax,kcount,1/),
     &                     real(cospec(1:Nax,1:kcount,l)))  )     		   
                     endif
                  enddo
               endif
            enddo
                              
            if (nrecspec==nbtime) then
               ierr = nf_CLOSE(ncidspec)
            else
               ierr = nf_sync(ncidspec)
            endif
            cospec = 0. 
                 
      end select

      
      end subroutine wrt_diags_spec  
