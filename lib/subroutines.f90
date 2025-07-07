!
!   how to compile?:
!   python -m numpy.f2py -c -m mylib subroutines.f90
!
subroutine print_coord(pos,box)

    implicit none
    real(kind=8), intent(in)    :: pos(:,:)
    real(kind=8), intent(in)    :: box(:)

    integer :: jpart

    print*, "# of particles:", size(pos,1)
    print*, "Box vector (a,b,c,alpha,beta,gamma):", box

    do jpart = 1, size(pos,1)
       print*, pos(jpart,1), pos(jpart,2), pos(jpart,3)
    enddo
     
end subroutine

subroutine rdf_self(pos,box,bin,histo) 

    implicit none
    real(kind=8), intent(in)    :: pos(:,:)
    real(kind=8), intent(in)    :: box(:)
    real(8)     , intent(in)    :: bin

    integer :: jpart, kpart, num_part, ibin
    real(8) :: dx, dy, dz, d, rho, vol, addval

    real(kind=8), intent(out), dimension(4000) :: histo
    real(8) :: r,pi=4d0*atan(1d0)

    histo(:) = 0d0
    num_part = size(pos,1)
    vol = box(1)*box(2)*box(3)
    rho = dble(num_part)/vol
    addval = 2d0/rho/bin/dble(num_part)

    do jpart = 1, num_part-1

       do kpart = jpart+1, num_part
         dx = pos(jpart,1) - pos(kpart,1)
         dy = pos(jpart,2) - pos(kpart,2)
         dz = pos(jpart,3) - pos(kpart,3)

         ! pbc 
         dx = dx - box(1)*dble( nint(dx/box(1)) )
         dy = dy - box(2)*dble( nint(dy/box(2)) )
         dz = dz - box(3)*dble( nint(dz/box(3)) )
 
         d = dsqrt( dx*dx + dy*dy + dz*dz )
         ibin = int(d/bin) + 1
         if ( ibin .le. 4000 ) then
         histo(ibin) = histo(ibin) + addval
         endif

       enddo
    enddo

    do ibin=1,4000
       r = (dble(ibin)-0.5d0)*bin
       histo(ibin) = histo(ibin)/(4d0*pi*r*r)
       !print*, r, gr
    enddo

end subroutine

subroutine rdf(pos1, pos2, box,bin,histo) 

    implicit none
    real(kind=8), intent(in)    :: pos1(:,:), pos2(:,:)
    real(kind=8), intent(in)    :: box(:)
    real(8)     , intent(in)    :: bin

    integer :: jpart, kpart, num_part1, num_part2, ibin
    real(8) :: dx, dy, dz, d, rho, vol, addval

    real(kind=8), intent(out), dimension(4000) :: histo
    real(8) :: r, pi=4d0*atan(1d0)

    histo(:) = 0d0
    num_part1 = size(pos1,1)
    num_part2 = size(pos2,1)
    vol = box(1)*box(2)*box(3)
    rho = dble(num_part2)/vol
    addval = 1d0/rho/bin/dble(num_part1)

    do jpart = 1, num_part1

       do kpart = 1, num_part2
         dx = pos1(jpart,1) - pos2(kpart,1)
         dy = pos1(jpart,2) - pos2(kpart,2)
         dz = pos1(jpart,3) - pos2(kpart,3)

         ! pbc 
         dx = dx - box(1)*dble( nint(dx/box(1)) )
         dy = dy - box(2)*dble( nint(dy/box(2)) )
         dz = dz - box(3)*dble( nint(dz/box(3)) )
 
         d = dsqrt( dx*dx + dy*dy + dz*dz )
         ibin = int(d/bin) + 1
         if ( ibin .le. 4000 ) then
         histo(ibin) = histo(ibin) + addval
         endif

       enddo
    enddo

    do ibin=1,4000
       r = dble(ibin)*bin
       histo(ibin) = histo(ibin)/(4d0*pi*r*r)
       !print*, r, gr
    enddo

end subroutine

subroutine center(pos,weights,num_residues,center_pos)
    implicit none
    real(kind=8), intent(in)    :: pos(:,:), weights(:)
    integer, intent(in)         :: num_residues

    real(kind=8), intent(out), dimension(num_residues,3) :: center_pos
    integer :: num_parts, size_resi

    integer :: iresi, jpart, idx
    real(8) :: sumval
    real(8), allocatable, dimension(:) :: weights_norm

    num_parts = size(pos,1)
    size_resi = size(weights,1)

    if ( mod(num_parts, size_resi) .ne. 0 ) then
       print*, "Number of particles and Size of residue unmatched!"
       print*,num_parts,size_resi
       print*,mod(num_parts,size_resi)
       stop
    endif
    if( num_residues .ne. num_parts/size_resi ) then
       print*, "Number of residues not matched!"
       stop
    endif

    ! Initialize array
    center_pos(:,:) = 0d0

    ! Normalize weights
    sumval = sum( weights )
    if ( sumval .eq. 0 ) then
       print*, "Sum of weights is zero!"
       stop
    endif

    allocate(weights_norm(size_resi))
    do jpart = 1, size_resi
       weights_norm(jpart) = weights(jpart)/sumval
    enddo

    ! Calculate center of weights 
    idx = 1
    do iresi = 1, num_residues
       do jpart = 1, size_resi
          center_pos(iresi,1) = center_pos(iresi,1) + pos(idx, 1) * weights_norm(jpart)
          center_pos(iresi,2) = center_pos(iresi,2) + pos(idx, 2) * weights_norm(jpart)
          center_pos(iresi,3) = center_pos(iresi,3) + pos(idx, 3) * weights_norm(jpart)
          idx = idx + 1
       enddo
    enddo
 
end subroutine

subroutine rdf_zz(pos, box, weights, bin, histo) 

    implicit none
    real(kind=8), intent(in)    :: pos(:,:)
    real(kind=8), intent(in)    :: box(:), weights(:)
    real(8)     , intent(in)    :: bin

    integer :: jpart, kpart, num_part, ibin
    real(8) :: dx, dy, dz, d, rho, vol, addval

    real(kind=8), intent(out), dimension(4000) :: histo
    real(8) :: r, pi=4d0*atan(1d0)

    histo(:) = 0d0
    num_part = size(pos,1)
    vol = box(1)*box(2)*box(3)
    rho = dble(num_part)/vol
    addval = 2d0/rho/bin/dble(num_part)

    do jpart = 1, num_part-1

       do kpart = jpart+1, num_part
         dx = pos(jpart,1) - pos(kpart,1)
         dy = pos(jpart,2) - pos(kpart,2)
         dz = pos(jpart,3) - pos(kpart,3)

         ! pbc 
         dx = dx - box(1)*dble( nint(dx/box(1)) )
         dy = dy - box(2)*dble( nint(dy/box(2)) )
         dz = dz - box(3)*dble( nint(dz/box(3)) )
 
         d = dsqrt( dx*dx + dy*dy + dz*dz )
         ibin = int(d/bin) + 1
         if ( ibin .le. 4000 ) then
         histo(ibin) = histo(ibin) + weights(jpart)*weights(kpart)*addval
         endif

       enddo
    enddo

    do ibin=1,4000
       r = (dble(ibin)-0.5d0)*bin
       histo(ibin) = histo(ibin)/(4d0*pi*r*r)
       !print*, r, gr
    enddo

end subroutine


