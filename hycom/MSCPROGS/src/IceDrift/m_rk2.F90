module m_rk2

contains

   ! Subroutine advances drift over obe day using second order RK
   subroutine rk2(u,v,scpx,scpy,nx,ny,x,y,drnx,drny,delt,undef)
   implicit none

   integer,                      intent(in)    :: drnx,drny,nx,ny
   real, dimension(  nx,  ny,2), intent(in)    :: u,v
   real, dimension(  nx,  ny)  , intent(in)    :: scpx,scpy
   real, dimension(drnx,drny)  , intent(inout) :: x,y
   real,                         intent(in)    :: delt
   real,                         intent(in)    :: undef


   real, dimension(drnx,drny)  :: drscpx, drscpy,  &
      urk1, vrk1, urk2, vrk2, urk, vrk, xpred, ypred

   integer :: rkdt
   integer :: istep2, irk
   integer :: i,j,ip,jp
   real :: wx,wy,a1,a2,a3,a4,wt


      ! Cycle time steps of this day
      do istep2=1,1/(delt/86400.)

         !print *,'istep2 ',istep2

         xpred=x
         ypred=y

         !print *,xpred,ypred

         ! Runge kutta time step
         do irk=1,2

            !Pivot points
            do j=1,drny
            do i=1,drnx

               if (x(i,j)/=undef) then
                  ! pivot points in MODEL
                  ip=floor(xpred(i,j))
                  jp=floor(ypred(i,j))

                  !print *,ip,jp

                  ! bilinear coefficients - spatial weight
                  wx=x(i,j)-ip
                  wy=y(i,j)-jp
                  a1=(1-wx)*(1-wy);
                  a2=wx*(1-wy);
                  a3=wx*wy;
                  a4=(1-wx)*wy;

                  ! spatial weighted for t=1,2
                  urk1(i,j)=a1*u(ip  ,jp  ,1) + &
                            a2*u(ip+1,jp  ,1) + &
                            a3*u(ip+1,jp+1,1) + &
                            a4*u(ip  ,jp+1,1) 
                  urk2(i,j)=a1*u(ip  ,jp  ,2) + &
                            a2*u(ip+1,jp  ,2) + &
                            a3*u(ip+1,jp+1,2) + &
                            a4*u(ip  ,jp+1,2) 

                  vrk1(i,j)=a1*v(ip  ,jp  ,1) + &
                            a2*v(ip+1,jp  ,1) + &
                            a3*v(ip+1,jp+1,1) + &
                            a4*v(ip  ,jp+1,1) 
                  vrk2(i,j)=a1*v(ip  ,jp  ,2) + &
                            a2*v(ip+1,jp  ,2) + &
                            a3*v(ip+1,jp+1,2) + &
                            a4*v(ip  ,jp+1,2) 

                  drscpx(i,j)=a1*scpx(ip  ,jp  ) + &
                              a2*scpx(ip+1,jp  ) + &
                              a3*scpx(ip+1,jp+1) + &
                              a4*scpx(ip  ,jp+1) 
                  drscpy(i,j)=a1*scpy(ip  ,jp  ) + &
                              a2*scpy(ip+1,jp  ) + &
                              a3*scpy(ip+1,jp+1) + &
                              a4*scpy(ip  ,jp+1) 
               else

                  urk1(i,j)=0.
                  vrk1(i,j)=0.
                  urk2(i,j)=0.
                  vrk2(i,j)=0.
                  drscpx(i,j)=1.
                  drscpy(i,j)=1.
               end if
            end do
            end do

            ! Temporal weight
            if (irk/=1) then
               rkdt=delt
            else
               rkdt=delt/2
            end if
            wt=((istep2-1)*delt)/86400.



            ! Time weighted
            urk=(1.-wt)*urk1 + wt*urk2
            vrk=(1.-wt)*vrk1 + wt*vrk2

            if (irk==1) then
               xpred=x + urk*rkdt/drscpx;
               ypred=y + vrk*rkdt/drscpy;
            else
               x=x+ urk*rkdt/drscpx;
               y=y+ vrk*rkdt/drscpy;
            end if
         enddo ! rk steps
      enddo ! istep

   end subroutine rk2
end module m_rk2
