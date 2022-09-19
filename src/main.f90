! file: main.f90

      module arvo_main
      use arvo_env, only: stderr, wp
      implicit none(type, external)

      private
      public :: arvo

      contains
! call arvo(n_atoms, coordinates, radii, probe_radius, V, S, stat, errmsg_f)
      subroutine arvo(ns,sc,sr,pr,volume,area,stat,errmsg)
      integer, intent(in) :: ns
!f2py intent(in) ns
      real(wp), intent(in) :: sc(3, ns), sr(ns), pr
!f2py intent(in) sc, sr, pr
!f2py depend(ns) sc, sr
      real(wp), intent(out) ::  area, volume
!f2py intent(out) area, volume
      integer, intent(out) :: stat
!f2py intent(out) stat
      character(len=:), allocatable, intent(out) :: errmsg
!f2py intent(out)
      real(wp) :: av(2)
      integer :: neighbors_indices(100000)
      real(wp), allocatable, dimension(:,:) :: spheres
      integer, allocatable, dimension(:) :: neighbors_number,index_start
      integer :: i
!     Computing surface area and volume of the overlapping spheres
      real(wp), parameter :: pi=3.14159265358979323846264d0, sa=0.324d0 
! "Random" sin value sa
      integer, parameter :: ks=1000,kl=300,ka=5000,ki=30000

!                ks - maximal sphere number
!                kl - maximal neighbors number of one sphere
!                ka - maximal angles or arcs number
!                ki - maximal neighbors relations number = cca.
!                     spheres number * maximal neighbors number
!        eps_north_pole - accuracy level in the function North_Pole_test
!        eps_deltat - accuracy level in the subroutine circles_intersect
!        eps_angle - accuracy level in the subroutine delete_equal (angl
!          
!      dimension spheres(ns,4),neighbors_number(ns),index_start(ns),&
!     &          neighbors_indices(ki),av(2)

      allocate(spheres(ks,4))
      allocate(neighbors_number(ks))
      allocate(index_start(ks))
! Validate input
      if (size(sc, dim=2) /= size(sr)) then
      errmsg = "Mismatch in dimension between coordinates and radii."
      stat = 1
      return
      endif

      if (size(sr) > ks) then
      errmsg = "That many spheres can not be processed by arvo."
      stat = 1
      return
      endif

      if (size(sc, dim=2).EQ.1 ) then
      volume = (4.0_wp/3.0_wp) * pi * sr(1)**3 
      area = 4.0_wp * pi * sr(1)**2
      stat = 0 
      return
      endif

!     spheres(i,1)=xi 
!     spheres(i,2)=yi     - ith sphere center coordinates
!     spheres(i,3)=zi 
!     spheres(i,4)=ri     - ith sphere radius 
!     neighbors_number, index_start, neighbors_indices description
!            is given in the    subroutine imake_neighbors
      do i=1,ns
          spheres(i,1) = sc(1,i)
          spheres(i,2) = sc(2,i)
          spheres(i,3) = sc(3,i)
          spheres(i,4) = sr(i) + pr
      enddo
!     ns - spheres number

!     Study the neighborhood relations
      call imake_neighbors(1,ns,spheres,neighbors_number,index_start,&
     &          neighbors_indices,ks,kl,ns,ki)

!   If some North Pole is close to the other atoms surface molecules rot
      do while (North_Pole_test(1,ns,spheres,neighbors_number,&
     &          index_start,neighbors_indices,ks,ki).EQ.0) 
           call ispheres_rotation(spheres,ks,ns,sa) ! random molecule rot
      enddo

!     Computation of area and volume as a sum of surface integrals
      volume=0d0
      area=0d0
      do i=1,ns
          call areavolume(i,spheres,neighbors_number,index_start,&
     &          neighbors_indices,ks,kl,ka,ki,av)
          volume=volume+av(1)
          area=area+av(2)
      enddo
     
      stat = 0
      return
      end subroutine arvo

      subroutine imake_neighbors(i1,i2,spheres,neighbors_number,&
     &           index_start,neighbors_indices,ks,kl,ns,ki)
!  
!       Determination of neighbors for all atoms. We construct next stru
!            neighbors_number(i)=neighbors number for ith atom
!          index_start(i)=start of neighbors indices for ith atom in arr
!                                 neighbors_indices
!          neighbors_indices - array of neighbors indices for each atom
!   neighbors_indices(index_start(i)):neighbors_ind(index_start(i)+neigh
!   
!          For example: 1. atom has neighbors with indices 2, 4, 7
!                       2. atom has neighbors with indices 1, 3
!                       3. atom has neighbors with indices 2, 4
!                       4. atom has neighbors with indices 1, 3
!                       5. atom is subset of some atom 
!                       6. atom has no neighbors 
!                       7. atom has neighbors with index 1
!          then we have
!                 neighbors_number=(3,2,2,2,-1,0,1)
!               index_start=(1,4,6,8,10,10,10,11)
!               neighbors_indices(2,4,7,1,3,2,4,1,3,1)
!                 
      real(wp), dimension(:,:), intent(in) :: spheres
      integer, intent(in) :: i1, i2, ks, kl, ns, ki
      integer, dimension(:), intent(inout) :: &
     &      neighbors_number, index_start
      integer, dimension(:), intent(inout) :: neighbors_indices
      integer :: ind(kl)
      integer :: i,j
        index_start(i1)=1
        do i=i1,i2
           neighbors_number(i)=ineighbors(i,spheres,ind,ks,kl,ns)      
           if (neighbors_number(i).LE.0) then 
!             sphere is subset ot there are no neighbors  
              index_start(i+1)=index_start(i)
           else  ! there are neighbors
              index_start(i+1)=index_start(i)+neighbors_number(i)
              do j=1,neighbors_number(i)
                 neighbors_indices(index_start(i)+j-1)=ind(j)
              enddo              
           endif
        enddo

      return
      end subroutine imake_neighbors




      integer function ineighbors(i,spheres,ind,ks,kl,ns)
!     
!     Function 
!  
!       If ith sphere is a subset of other sphere, index_number(i)=-1 an
!       radius in matrix spheres to -radius !!!
!       If some other sphere is subset of ith sphere, than we change its
!  
      integer, intent(in) :: i, ks, kl, ns
      real(wp), dimension(:,:), intent(in) :: spheres
      integer, intent(inout) :: ind(kl)
      real(wp) :: dd, xi, yi, zi, ri, rk
      integer :: k, neighbors_num

      neighbors_num=0
!     i-th sphere data
      xi=spheres(i,1)
      yi=spheres(i,2)
      zi=spheres(i,3)
      ri=spheres(i,4)
        do k=1,ns
           if (k .NE. i) then
!           first simple test 
              if(dabs(xi-spheres(k,1)).LT.ri+spheres(k,4)) then
                 dd=dsqrt((xi-spheres(k,1))**2+(yi-spheres(k,2))**2+&
     &                (zi-spheres(k,3))**2)
                 rk=spheres(k,4)
                 if (dd.LT.ri+rk) then
                    if (dd+ri.LE.rk) then
!                    ith sphere is inside of other sphere 
                       neighbors_num=-1
                     exit
                    elseif (dd+rk.GT.ri) then  
!                    kth sphere is neighbor 
                       neighbors_num=neighbors_num+1
                       ind(neighbors_num)=k
                    endif
                 endif
              endif
           endif
        enddo
      ineighbors=neighbors_num
      return
      end function ineighbors


      integer function North_Pole_test(i1,i2,spheres,neighbors_number,&
     &                    index_start,neighbors_indices,ks,ki)
!          Here we check, that North Pole of no sphere lies on other 
!          neighbor sphere
!  
!          dmin - square of minimal distance of the North Pole to neighb
!  

      real(wp), dimension(:,:), intent(in) :: spheres
      integer, intent(in) :: i1, i2, ks, ki
      integer, allocatable, dimension(:), intent(inout) :: &
     &      neighbors_number, index_start
      integer, dimension(:), intent(inout) :: neighbors_indices
      integer :: i,k,ink,npt
      real(wp) :: eps_north_pole, dmin, d

!       Test precision - MAY BE CHANGED
      eps_north_pole=1d-4
      dmin=10000d0
      do i=i1,i2  
         do k=1,neighbors_number(i)
              ink=neighbors_indices(index_start(i)+k-1) ! kth neighbor i
              d=dabs(dsqrt((spheres(i,1)-spheres(ink,1))**2&
     &          +(spheres(i,2)-spheres(ink,2))**2&
     &          +(spheres(i,3)+spheres(i,4)-spheres(ink,3))**2)&
     &          -spheres(ink,4))
              if (d.LT.dmin) then 
                dmin=d 
              endif
          enddo
      enddo


!     minimal distance = dmin   
      if (dmin.LT.eps_north_pole) then
         npt=0 ! Bad news!
      else
         npt=1 ! O.K.
      endif
         
      North_Pole_test=npt   
           
      return
      end function North_Pole_test


      subroutine ispheres_rotation(spheres,ks,ns,sa)
!       Random rotation of molecule about the y-axis
!       after bad North Pole test.
!          Some North Pole is near other spheres surface
      real(wp), dimension(:,:), intent(inout) :: spheres
      integer, intent(in) :: ks, ns
      real(wp), intent(in) :: sa
      real(wp) :: x, z, ca
      integer :: i, j
      
      ca=dsqrt(1d0-sa*sa)     
      do i=1,ns
         x=spheres(i,1)
         z=spheres(i,3)
         spheres(i,1)=ca*x-sa*z
         spheres(i,3)=sa*x+ca*z       
      enddo
      return
      end subroutine ispheres_rotation


      subroutine areavolume(i,spheres,neighbors_number,index_start,&
     &                        neighbors_indices,ks,kl,ka,ki,av)
!  
!     Function computes i-th part of the whole volume -
!      - the volume of domain inside i-th and outside of
!        all other spheres
!     
      real(wp), dimension(:,:), intent(in) :: spheres
      integer, allocatable, dimension(:), intent(in) :: &
     &      neighbors_number, index_start
      integer, dimension(:), intent(in) :: neighbors_indices
      integer, intent(in) :: ks, kl, ki, ka
      real(wp), intent(out) :: av(2)
      real(wp) ::avi(2)
      real(wp), dimension(:,:), allocatable :: circles, arcs,&
     & sphere_local
      integer :: ind(kl)
      integer :: i,j,nls,npos,narcs
      real(wp) :: z1,r1
!  
!       circles, arcs, sphere_local are described below
!  
      real(wp), parameter :: pi=3.14159265358979323846264d0
      allocate(circles(kl,4))
      allocate(arcs(ka,3))
      allocate(sphere_local(kl,4))
!  
!      Determination of i-th sphere's neighbors (row indices in matrix s
!  
      if (neighbors_number(i).LT.0) then
!  
!     ith sphere is subset of other sphere
!           sphere(i,4) will be done negative
!  
         av(1)=0d0
         av(2)=0d0
        elseif (neighbors_number(i).EQ.0) then 
!  
!       there are no neighbors (nls - number of local spheres = ith sphe
!  
           av(1)=4d0*pi*spheres(i,4)**3/3.d0
         av(2)=4d0*pi*spheres(i,4)**2
      else 
!  
!        there are neighbors
!  
         nls=neighbors_number(i)+1
         ind(1)=i
         do j=1,(nls-1)
              ind(j+1)=neighbors_indices(index_start(i)+j-1)
         enddo

!          we will work only with ith and neighbors spheres             
           call local_spheres(spheres,ind,sphere_local,nls,ks,kl)
           av(1)=0d0
         av(2)=0d0

         call make_ts_circles(sphere_local,circles,kl,nls)
           narcs=circles_to_arcs(circles,arcs,kl,nls,ka)

         npos=0
         do j=1,(nls-1)
              if (circles(j,4).GT.0) then 
               npos=npos+1
            endif
         enddo

         z1=sphere_local(1,3)
         r1=sphere_local(1,4)
         if (npos.GT.0) then   
!             there exists positive oriented circle 
              call avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi) 
            av(1)=av(1)+avi(1)
            av(2)=av(2)+avi(2)
           else 
!             all circles are negative oriented - we compute complement
              call avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi) 
            av(1)=av(1)+avi(1)+4d0*pi*sphere_local(1,4)**3/3d0 
            av(2)=av(2)+avi(2)+4d0*pi*sphere_local(1,4)**2 
           endif
        endif

      return
      end subroutine areavolume




      subroutine local_spheres(spheres,ind,sphere_local,nls,ks,kl)
!  
!     Take sphere_local out of the main array spheres
!  
      integer, intent(in) :: nls, ks, kl
      real(wp), intent(in), dimension(:,:) :: spheres
      real(wp), intent(inout), dimension(:,:) :: sphere_local
      integer, intent(in) :: ind(kl)
      integer :: i, j

      do i=1,nls
          do j=1,4
              sphere_local(i,j)=spheres(ind(i),j)
          enddo 
        enddo 

      return
      end subroutine local_spheres





      subroutine make_ts_circles(sphere_local,circles,kl,nls)
!  
!       Preparing circles structure for 1st sphere in array   circles
!       according to the paper Hayrjan, Dzurina, Plavka, Busa
!       
!       circles(i,1)=ti 
!       circles(i,2)=si     - ith circle's center coordinates
!       circles(i,3)=ri    - ith circle's radius 
!       circles(i,4)=+1/-1 - circle orientation 
!   
      integer, intent(in) :: kl, nls
      real(wp), intent(in), dimension(:,:) :: sphere_local
      real(wp), intent(inout), dimension(:,:) :: circles
      integer :: i,k
      real(wp) :: r1, dx, dy, a, b, c, d
      
      r1=sphere_local(1,4)
      do k=1,(nls-1)
          dx=sphere_local(1,1)-sphere_local(k+1,1)
          dy=sphere_local(1,2)-sphere_local(k+1,2)    
          a=dx*dx+dy*dy+(sphere_local(1,3)+r1-sphere_local(k+1,3))**2-&
     &       sphere_local(k+1,4)**2
          b=8d0*r1*r1*dx
          c=8d0*r1*r1*dy
          d=4d0*r1*r1*(dx*dx+dy*dy+(sphere_local(1,3)-r1-&
     &       sphere_local(k+1,3))**2-sphere_local(k+1,4)**2)
          circles(k,1)=-b/(2d0*a)
          circles(k,2)=-c/(2d0*a)       
          circles(k,3)=dsqrt((b*b+c*c-4d0*a*d)/(4d0*a*a))   
          if (a.GT.0) then
              circles(k,4)=-1
          else 
              circles(k,4)=1
          endif
      enddo
      return    
      end subroutine make_ts_circles




      integer function circles_to_arcs(circles,arcs,kl,nls,ka)
!  
!       Computing integration arcs
!       
!       arcs(i,1)=ci     - corresponding circle index
!       arcs(i,2)=sigma     - starting arc angle 
!       arcs(i,3)=delta     - oriented arc angle
!  
!       Arcs (with their orientation) are parts of circles, which
!       bounds are circles intersection points. If the center of
!     arc lies inside all other positive and outside all other
!       negative circles, then we will put it inside arcs structure
!  
      integer, intent(in) :: kl, nls, ka
      real(wp), intent(in), dimension(:,:) :: circles
      real(wp), intent(inout), dimension(:,:) :: arcs
      real(wp), allocatable, dimension(:,:) :: arcsnew
      integer :: number_arc, i, j, k, nna
      real(wp), parameter :: pi=3.14159265358979323846264d0

      allocate(arcsnew(ka,3))

      number_arc=0
      if (nls.EQ.2) then
!         we have only 1 circle
            number_arc=1
            arcs(1,1)=1
            arcs(1,2)=0d0
            arcs(1,3)=2d0*pi*circles(1,4)
      else 
!         more than 1 circle
          do i=1,(nls-1)
              nna=new_arcs(i,circles,arcsnew,kl,ka,nls)
              if (nna.GT.0) then
                  do j=1,nna
                      do k=1,3
                          arcs(number_arc+j,k)=arcsnew(j,k)
                      enddo                   
                  enddo
                  number_arc=number_arc+nna
              endif    
          enddo
      endif

      deallocate(arcsnew)
      circles_to_arcs=number_arc
      return
      end function circles_to_arcs




      integer function new_arcs(i,circles,arcsnew,kl,ka,nls)
!  
!     Function prepares arcs, which are part of i-th circle
!     in circle structure circles.
!     Interesting are these arcs, which are inside other positive
!     circles or outside other negative circles
!   
!     Matrix arcsnew in each row has elements
!  
!       arcsnew(i,1)=ic - ic is the index of arc-circle in circle 
!       arcsnew(i,2)=sigma - sigma is the starting angle of arc
!       arcsnew(i,3)=delta - delta is oriented arc angle
!  
      integer, intent(in) :: i, kl, ka, nls
      real(wp), intent(in), dimension(:,:) :: circles
      real(wp), intent(inout), dimension(:,:) :: arcsnew
      real(wp) :: angles(ka)
      real(wp), parameter :: pi=3.14159265358979323846264d0
      integer :: num_arc, num_angle, number_cond, j, jj, na 
      real(wp) :: ti, si, ri, t, s, r, d, a1, a2, b1, b2

      num_arc=0
      num_angle=0

      ti=circles(i,1)
      si=circles(i,2)
      ri=circles(i,3)
      do j=1,(nls-1) 
!        composition of angles vector, consisting of intersection points
          if (j .NE. i) then
              t=circles(j,1)
              s=circles(j,2)
              r=circles(j,3)
              d=dsqrt((ti-t)**2+(si-s)**2)
              if ( (d.LT.r+ri) .AND. (dabs(r-ri).LT.d) ) then
!                  2 intersection points exist
                  call circles_intersection(i,j,circles,kl,a1,a2,b1,b2)
                    angles(num_angle+1)=a1
                    angles(num_angle+2)=a2
                    num_angle=num_angle+2            
              endif
          endif
      enddo
      if (num_angle .EQ. 0) then
!         there are no double intersections of i-th circles with others
          number_cond=0 
!         if i-th circle is inside of all other positive and outside of 
!           all other negative circles, it will be new arc
          do j=1,(nls-1)
             if (j.NE.i) then
               number_cond=number_cond+circle_in_circle(i,j,circles,kl)
             endif
          enddo
          if (number_cond.EQ.(nls-2)) then
!             all conditions hold
              num_arc=1
              arcsnew(1,1)=i
              arcsnew(1,2)=0d0
              arcsnew(1,3)=2d0*pi*circles(i,4)
          endif
      else 
!         there are double intersection points
          if (circles(i,4).GT.0) then
              call mysort(angles,ka,num_angle)
          else
              call mydsort(angles,ka,num_angle)
          endif
          na=idelete_equal(angles,ka,num_angle)
          num_angle=na
          do j=1,(na-1)
              number_cond=0
              do jj=1,(nls-1)
                  if (jj.NE.i) then
                      t=ti+ri*dcos((angles(j)+angles(j+1))/2d0)
                      s=si+ri*dsin((angles(j)+angles(j+1))/2d0)
              number_cond=number_cond+point_in_circle(t,s,jj,circles,&
     &        kl)
                  endif
              enddo
              if (number_cond.EQ.(nls-2)) then
!                 all conditions hold
                  num_arc=num_arc+1  
                  arcsnew(num_arc,1)=i
                  arcsnew(num_arc,2)=angles(j)
                  arcsnew(num_arc,3)=angles(j+1)-angles(j)
              endif
          enddo
          number_cond=0
          do j=1,(nls-1)
              if (j.NE.i) then
                  t=ti+ri*dcos((angles(1)+2d0*pi+angles(na))/2d0)
                  s=si+ri*dsin((angles(1)+2d0*pi+angles(na))/2d0)
              number_cond=number_cond+point_in_circle(t,s,j,circles,kl)
              endif
          enddo
          if (number_cond.EQ.(nls-2)) then
!             all conditions hold
              num_arc=num_arc+1  
              arcsnew(num_arc,1)=i
              arcsnew(num_arc,2)=angles(na)
          arcsnew(num_arc,3)=angles(1)+circles(i,4)*2d0*pi-angles(na)
          endif
      endif

      new_arcs=num_arc
      return
      end function new_arcs
          


      subroutine circles_intersection(ic1,ic2,circles,kl,a1,a2,b1,b2)
!  
!     Function returns angles of two intersection points
!     of circles with indices ic1 and ic2 in circles structure circles
!     (we will use it ONLY IN CASE, WHEN 2 INTERSECTION POINTS EXIST!!!)
!  
!     a1 and a2 are corresponding angles with respect to the center of 1
!     b1 and b2 are corresponding angles with respect to the center of 2
!  
      real(wp), intent(in), dimension(:,:) :: circles
      integer, intent(in) :: ic1, ic2, kl
      real(wp), intent(inout) :: a1, a2, b1, b2
      real(wp) :: t1, s1, r1, t2, s2, r2, A, B, C, D
      real(wp), parameter :: pi=3.14159265358979323846264d0
      real(wp), parameter :: eps_deltat = 1d-12
!       (t,s) - circle center, r - circle radius
      t1=circles(ic1,1)
      s1=circles(ic1,2)
      r1=circles(ic1,3) 
      t2=circles(ic2,1)
      s2=circles(ic2,2)
      r2=circles(ic2,3) 
      if (dabs(t2-t1).LT.eps_deltat) then
!           t2 .EQ. t1
          B=((r1*r1-r2*r2)/(s2-s1)-(s2-s1))/2d0
          A=dsqrt(r2*r2-B*B)
          if (B.EQ.0) then
              b1=0d0
              b2=pi
          elseif (B.GT.0) then
              b1=datan(dabs(B/A))
              b2=pi-b1
          else
              b1=pi+datan(dabs(B/A))
              b2=3d0*pi-b1        
          endif
          B=B+s2-s1
          if (B.EQ.0) then
              a1=0d0
              a2=pi
          elseif (B.GT.0) then
              a1=datan(dabs(B/A))
              a2=pi-a1
          else
              a1=pi+datan(dabs(B/A))
              a2=3d0*pi-a1        
          endif    
      else 
!       t2 .NE. t1
          C=((r1*r1-r2*r2-(s2-s1)**2)/(t2-t1)-(t2-t1))/2d0
          D=(s1-s2)/(t2-t1)
          B=(-C*D+dsqrt((D*D+1d0)*r2*r2-C*C))/(D*D+1d0)
          A=C+D*B
          if (A.EQ.0) then
              if (B.GT.0) then
                  b1=pi/2d0
              else 
                  b1=-pi/2d0
              endif
          elseif (A.GT.0) then
              b1=datan(B/A)
          else
              b1=pi+datan(B/A)
          endif
          B=B+s2-s1
          A=A+t2-t1
          if (A.EQ.0) then 
              if (B.GT.0) then
                  a1=pi/2d0
              else 
                  a1=-pi/2d0
              endif
          elseif (A.GT.0) then
              a1=datan(B/A)
          else
              a1=pi+datan(B/A)
          endif
          B=(-C*D-dsqrt((D*D+1d0)*r2*r2-C*C))/(D*D+1d0)
          A=C+D*B
          if (A.EQ.0) then
              if (B.GT.0) then
                  b2=pi/2d0
              else 
                  b2=-pi/2d0
              endif
          elseif (A.GT.0) then
              b2=datan(B/A)
          else
              b2=pi+datan(B/A)
          endif
          B=B+s2-s1
          A=A+t2-t1
          if (A.EQ.0) then
              if (B.GT.0) then
                  a2=pi/2d0
              else 
                  a2=-pi/2d0
              endif
          elseif (A.GT.0) then
              a2=datan(B/A)
          else
              a2=pi+datan(B/A)
          endif
      endif
      if (a1.LT.0) a1=a1+2d0*pi 
      if (a2.LT.0) a2=a2+2d0*pi
      if (b1.LT.0) b1=b1+2d0*pi
      if (b2.LT.0) b2=b2+2d0*pi

      return
      end subroutine circles_intersection
      
          
          


      integer function circle_in_circle(i,k,circles,kl)
!  
!        1  if i-th circle is inside k-th positive circle or 
!                            outside k-th negative circle
!      0  - otherwise
!  
!        WE KNOW, THAT CIRCLES HAVE LESS THAN 2 INTERSECTION POINTS !!!
!  
      integer, intent(in) :: i, k, kl
      real(wp), intent(in), dimension(:,:) :: circles
      real(wp) :: d

      d=dsqrt((circles(i,1)+circles(i,3)-circles(k,1))**2+&
     &        (circles(i,2)-circles(k,2))**2)
      if (d.LT.circles(k,3)) then
          if (circles(k,4).GT.0) then
              circle_in_circle=1
          else
              circle_in_circle=0
          endif    
      elseif (d.GT.circles(k,3)) then 
          if (circles(k,4).GT.0) then
              circle_in_circle=0
          else
              circle_in_circle=1
          endif
      else 
!           d=circles(k,3) - right point on k-th circle - touching of ci
          d=dsqrt((circles(i,1)-circles(k,1))**2+&
     &   (circles(i,2)-circles(k,2))**2)
          if (d.LT.circles(k,3)) then
              if (circles(k,4).GT.0) then
                  circle_in_circle=1
              else
                  circle_in_circle=0
              endif
          else
              if (circles(k,4).GT.0) then
                  circle_in_circle=0
              else
                  circle_in_circle=1
              endif
          endif
      endif

      return
      end function circle_in_circle



      integer function point_in_circle(t,s,k,circles,kl)
!  
!  
!       1  if point (t,s) is inside k-th positive circle 
!                        or outside k-th negative circle
!       0  - otherwise
!  
!     WE KNOW, THAT POINT IS NOT ON THE CIRCLE !!!
!  
      integer, intent(in) :: k, kl
      real(wp), intent(in), dimension(:,:) :: circles
      real(wp), intent(in) :: t, s 
      real(wp) :: d

      d=dsqrt((t-circles(k,1))**2+(s-circles(k,2))**2)
      if (d.LT.circles(k,3)) then
          if (circles(k,4).GT.0) then
              point_in_circle=1
          else
              point_in_circle=0
          endif
      else
          if (circles(k,4).GT.0) then
              point_in_circle=0
          else
              point_in_circle=1
          endif
      endif

      return
      end function point_in_circle




      subroutine mysort(angles,ka,num_angle)
!  
!     Sorting array angles in increasing order
!     num_angle is the angles array length
!  
      real(wp), intent(inout), dimension(:) :: angles
      integer, intent(in) :: num_angle, ka
      real(wp) :: amax
      integer :: i, j, ii, k
      
      do i=1,(num_angle-1)
          ii=i
          amax=angles(i)
          do j=i+1,num_angle
              if (amax.GT.angles(j)) then
                  ii=j
                  amax=angles(j)      
              endif
          enddo
          if (ii.NE.i) then
              angles(ii)=angles(i)
              angles(i)=amax
          endif
      enddo


      return
      end subroutine mysort 




      subroutine mydsort(angles,ka,num_angle)
!  
!     Sorting array angles in decreasing order
!     num_angle is the angles array length
!  
      real(wp), intent(inout), dimension(:) :: angles
      integer, intent(in) :: num_angle, ka
      real(wp) :: amin
      integer :: i, j, ii, k
      
      do i=1,(num_angle-1)
          ii=i
          amin=angles(i)
          do j=i+1,num_angle
              if (amin.LT.angles(j)) then
                  ii=j
                  amin=angles(j)      
              endif
          enddo
          if (ii.NE.i) then
              angles(ii)=angles(i)
              angles(i)=amin
          endif
      enddo


      return
      end subroutine mydsort 



      integer function idelete_equal(angles,ka,num_angle)
!  
!     Deletion of "equal" (to some precision eps_angle)
!     angles in sorted vector angles
!  
      integer, intent(in) :: ka
      integer, intent(in) :: num_angle 
      real(wp), intent(inout), dimension(:) :: angles
      real(wp) :: anglesnew(ka)
      real(wp) :: angle
      integer :: i, m, delete_equal
      real(wp), parameter :: eps_angle=1d-12

      m=1 
      angle=angles(1)
      anglesnew(1)=angle
      do i=2,num_angle
          if (dabs(angles(i)-angle).GT.eps_angle) then
              angle=angles(i)
              m=m+1
              anglesnew(m)=angle
          endif
      enddo
      delete_equal=m
      do i=1,m
          angles(i)=anglesnew(i)
      enddo
      idelete_equal = delete_equal
      return
      end function idelete_equal


      subroutine avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi)
!  
!     Computing integrals over arcs given in arc structure
!       according to paper by Hayrian, Dzurina, Plavka, Busa
!       
      integer, intent(in) :: kl, ka, narcs
      real(wp), intent(in), dimension(:,:) :: circles, arcs
      real(wp), intent(in) :: r1, z1
      real(wp), intent(inout) :: avi(2)
      real(wp), parameter :: pi=3.14159265358979323846264d0
      real(wp), parameter :: eps_two_pi=1d-12
      integer :: i, k
      real(wp) :: t, s, r, A, B, C, rr
      real(wp) :: vIone, vItwo, vIthree, vJone, vJtwo, vJthree
      real(wp) :: delta_vint, delta_aint
      real(wp) :: sb, cb, sa, ca, al, be

      avi(1)=0d0
      avi(2)=0d0

      do k=1,narcs 
!         cycle over all arcs
          t=circles(int(arcs(k,1)),1)
          s=circles(int(arcs(k,1)),2)
          r=circles(int(arcs(k,1)),3)
          A=(4d0*r1*r1+t*t+s*s+r*r)/2d0
          B=t*r
          C=s*r
          S=dsqrt(A*A-B*B-C*C) 
          rr=r*r-A
          if (dabs(dabs(arcs(k,3))-2d0*pi).LT.eps_two_pi) then
!             full circle arc
              vIone=2d0*pi/S
              vItwo=2d0*pi*A/(S**3)
              vIthree=pi*(2d0*A*A+B*B+C*C)/(S**5)
              vJone=pi+rr/2d0*vIone
              vJtwo=(vIone+rr*vItwo)/4d0
              vJthree=(vItwo+rr*vIthree)/8d0
              delta_vint=(128d0*vJthree*r1**7+8d0*vJtwo*r1**5+&
     &         2d0*vJone*r1**3)/3d0-8d0*r1**4*vJtwo*(z1+r1)
              delta_aint=2d0*vJone*r1**2
              if (arcs(k,3).LT.0) then
                 delta_vint=-delta_vint
                 delta_aint=-delta_aint
                endif
              avi(1)=avi(1)+delta_vint    
              avi(2)=avi(2)+delta_aint    
          else
!           integration over arcs
              if (arcs(k,3).LT.0) then 
                  al=arcs(k,2)+arcs(k,3)
                  be=arcs(k,2) 
              else
                  be=arcs(k,2)+arcs(k,3)
                  al=arcs(k,2) 
              endif 
              vIone=2d0*(pi/2d0-datan((A*dcos((be-al)/2d0)+&
     &             B*dcos((al+be)/2d0)+C*dsin((al+be)/2d0))/&
     &            (S*dsin((be-al)/2d0))))/S
              sb=dsin(be)
              cb=dcos(be)
              sa=dsin(al)
              ca=dcos(al)
              vItwo=(fract(A,B,C,sb,cb,1)-fract(A,B,C,sa,ca,1)+&
     &             A*vIone)/(S*S)
              vIthree=(fract(A,B,C,sb,cb,2)-fract(A,B,C,sa,ca,2)+&
     &           (fract(A,B,C,sb,cb,1)-fract(A,B,C,sa,ca,1))/A+&
     &                  (2d0*A*A+B*B+C*C)*vItwo/A)/(2d0*S*S)
              vJone=((be-al)+rr*vIone)/2d0
              vJtwo=(vIone+rr*vItwo)/4d0
              vJthree=(vItwo+rr*vIthree)/8d0
              delta_vint=(128d0*vJthree*r1**7+8d0*vJtwo*r1**5+&
     &           2d0*vJone*r1**3)/3d0-8d0*r1**4*vJtwo*(z1+r1)
              delta_aint=2d0*vJone*r1**2
              if (arcs(k,3).LT.0) then
                 delta_vint=-delta_vint
                 delta_aint=-delta_aint
                endif
              avi(1)=avi(1)+delta_vint    
              avi(2)=avi(2)+delta_aint    
          endif
      enddo

      return
      end subroutine avintegral



      real(wp) function fract(A,B,C,sinphi,cosphi,k)
!  
!     Fraction evaluation for integral
!  
      real(wp), intent(in) :: A, B, C, sinphi, cosphi
      integer, intent(in) :: k

        fract=(-B*sinphi+C*cosphi)/(A+B*cosphi+C*sinphi)**k

      return
      end function fract


      end module arvo_main
!   end file: main.f90
