c      Subroutine to compute the galactic velocity relative to the box
c      center in Hidodynamical simulations (HART)

       subroutine findvelc (x,y,z,vxt,vyt,vzt,npart,vxc,vyc,vzc)

       include 'parameters.h'

       dimension vxt(npartm),vyt(npartm),vzt(npartm),npart(5)
       dimension x(npartm),y(npartm),z(npartm)

      vxc=0.
      vyc=0.
      vzc=0.

      num=0

      do i=npart(2)+1,npart(2)+1+lspecies(2)
        if(dsqrt(x(i)**2.0d0+y(i)**2.0d0+z(i)**2.0d0).lt.100.d0) then
        vxc=vxc+vxt(i)
        vyc=vyc+vyt(i)
        vzc=vzc+vzt(i)
        num=num+1
        endif
      enddo

        vxc=vxc/dble(num)
        vyc=vyc/dble(num)
        vzc=vzc/dble(num)

      write(*,*)'Velocity of the galaxy with respect to the box center',
     +vxc,vyc,vzc


       return
       end
