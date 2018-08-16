! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE partlist

  USE shared_data
#ifdef PHOTONS
  USE random_generator
#endif
  USE utilities

  IMPLICIT NONE

  SAVE

  INTEGER :: nvar

  TYPE pointer_item
    TYPE(particle), POINTER :: part
    TYPE(pointer_item), POINTER :: next
  END TYPE pointer_item

  TYPE pointer_list
    TYPE(pointer_item), POINTER :: head, tail
  END TYPE pointer_list

  REAL(num), DIMENSION(:), ALLOCATABLE :: packed_particle_data

CONTAINS

  SUBROUTINE setup_partlists

    nvar = 3 + c_ndims
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    nvar = nvar+1
#endif
#ifdef DELTAF_METHOD
    nvar = nvar+1
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    nvar = nvar+2
#endif
#ifdef PARTICLE_DEBUG
    nvar = nvar+2
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    nvar = nvar+1
#endif
#ifdef COLLISIONS_TEST
    nvar = nvar+1
#endif
#ifdef PHOTONS
    nvar = nvar+2
#ifdef TRIDENT_PHOTONS
    nvar = nvar+1
#endif
#endif
    ALLOCATE(packed_particle_data(nvar))

  END SUBROUTINE setup_partlists



  SUBROUTINE deallocate_partlists

    INTEGER :: stat

    IF (ALLOCATED(packed_particle_data)) &
        DEALLOCATE(packed_particle_data, STAT=stat)

  END SUBROUTINE deallocate_partlists



  SUBROUTINE create_empty_partlist(partlist, use_store_in)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    LOGICAL, INTENT(IN), OPTIONAL :: use_store_in
    LOGICAL :: use_store

    IF(.NOT. PRESENT(use_store_in)) THEN
      use_store = .FALSE.
    ELSE
      use_store = use_store_in
    ENDIF

    NULLIFY(partlist%head)
    NULLIFY(partlist%tail)
    partlist%count = 0
    partlist%id_update = 0
    partlist%safe = .TRUE.
    partlist%use_store = use_store
    IF (use_store) &
        CALL create_particle_store(partlist, sublist_size, .FALSE., .TRUE.)

  END SUBROUTINE create_empty_partlist



  SUBROUTINE create_particle_store(partlist, n_els_min, &
      link_el_in, no_pad_store)


    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER(i8), INTENT(IN) :: n_els_min
    INTEGER(i8) :: actual_elements, i_part, els_allocated
    INTEGER(i8) :: i_sub, n_subs, last, link_to
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: els_to_allocate
    LOGICAL, INTENT(IN), OPTIONAL :: link_el_in, no_pad_store
    LOGICAL :: link_el

    !If sublists, need to deallocate any
    CALL destroy_store(partlist%store)

    !Allocate one sublist of correct size
    actual_elements = FLOOR(MAX(n_els_min, sublist_size) &
        *(1.0_num + sublist_slack), i8)
    IF(present(no_pad_store)) THEN
      IF(no_pad_store) THEN
        actual_elements = MAX(n_els_min, sublist_size)
      ENDIF
    ENDIF
    actual_elements = MAX(n_els_min, sublist_size)

    link_el = .TRUE.
    IF(present(link_el_in)) THEN
      link_el = link_el_in
    ENDIF

    !TODO Here goes logic on how to split up actual request
    !Not sure if we want to allocate a large chunk and then add
    !pieces, or all in pieces, or other
    !For now, all pieces

    n_subs = CEILING(FLOAT(actual_elements)/FLOAT(sublist_size))
    ALLOCATE(els_to_allocate(n_subs))
    els_to_allocate = sublist_size

    DO i_sub = 1, n_subs
      link_to = els_to_allocate(i_sub)
      IF(.NOT. link_el) link_to = 0
      IF(i_sub == n_subs) THEN
        last = n_els_min - SUM(els_to_allocate(1:n_subs-1))
        IF(link_to > 0) link_to = last
        !TODO This might be off by 1??
      ENDIF
      CALL create_linked_substore(partlist%store, els_to_allocate(i_sub),&
          link_to)
    ENDDO

    partlist%store%next_slot => partlist%store%tail%store(&
        partlist%store%tail%first_free_element)
    IF(partlist%store%tail%first_free_element > 1) THEN
      partlist%store%next_slot%prev &
          => partlist%store%tail%store(partlist%store%tail%first_free_element-1)
    ELSE IF(ASSOCIATED(partlist%store%tail%prev)) THEN
      !TODO this wont work properly in general, e.g. if prev is also empty
      !Do we ever need to consider this case?
      partlist%store%next_slot%prev &
          => partlist%store%tail%prev%store(&
            partlist%store%tail%prev%first_free_element-1)
    ELSE
      NULLIFY(partlist%store%next_slot%prev)
    ENDIF

    IF(link_to > 0) THEN
      !TODO actually the following MAY NOT be list tail, should do better
      partlist%tail => partlist%store%next_slot%prev
      partlist%head => partlist%store%head%head
    ELSE
      NULLIFY(partlist%head, partlist%tail)
    ENDIF
  END SUBROUTINE create_particle_store


  !Create and a new store segment object
  !in given store
  !Takes total_size, link_upto

  !Called in two circumstances
  !Creating a list at first
  !And adding a chunk to a list
  SUBROUTINE create_linked_substore(store, total_size, link_upto_in)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: substore, tmp1, tmp2
    INTEGER(i8), INTENT(IN) :: total_size, link_upto_in
    INTEGER(i8) :: link_upto
    TYPE(particle), POINTER :: current, prev
    INTEGER(i8) :: i_part

    !Nothing to do
    IF(total_size <= 0) RETURN

    link_upto = link_upto_in
    IF(link_upto_in > total_size) link_upto = total_size

    !Allocate substore object
    ALLOCATE(substore)
    !Set length, nullify links
    substore%length = total_size
    NULLIFY(substore%prev, substore%next)
    !Allocate backing memory
    ALLOCATE(substore%store(total_size))
    substore%head => substore%store(1)
    !If creating list, we can link it all up already
    !And setup the list
    !Then calling code just sets positions etc
    IF(link_upto > 0) THEN
      substore%first_free_element = link_upto + 1
    ELSE
      substore%first_free_element = 1
    ENDIF

    DO i_part = 1, total_size
      IF (i_part > 1 .AND. i_part < link_upto) THEN
        !Each particle slot up to request length should be linked up
        substore%store(i_part)%prev => substore%store(i_part-1)
        substore%store(i_part)%next => substore%store(i_part+1)
      ELSE
        !Nullify pointers
        NULLIFY(substore%store(i_part)%prev, &
            substore%store(i_part)%next)
      ENDIF
      !Set not-live state
      substore%store(i_part)%live = -1
    ENDDO

    IF(link_upto > 1) THEN
      !Correct links for 0th and n_elements-th particles
      substore%store(1)%next => substore%store(2)
      IF(total_size > 1 .AND. link_upto >= total_size) THEN
        substore%store(total_size)%prev => substore%store(total_size-1)
      ELSE IF(link_upto > 1) THEN
        substore%store(link_upto)%prev => substore%store(link_upto-1)
      ENDIF
    ENDIF

    IF(ASSOCIATED(store%tail)) THEN
      !Now link substore into store
      tmp1 => store%tail
      tmp2 => store%tail%next
      store%tail%next => substore
      tmp2 => store%tail%next

      substore%prev => store%tail
      tmp2 => substore%prev

      !And create link in partlist between prior last particle
      !and first particle of new sub
      IF(link_upto > 1) THEN
        IF(store%tail%first_free_element > 1) THEN
          prev => store%tail%store(store%tail%first_free_element-1)
        ELSE
          !Have a completely empty sub before this one
          !Probably means it's been compacted out?
       ENDIF
        current => substore%store(1)
        prev%next => current
        current%prev => prev
      ENDIF
      store%tail => substore
      tmp2 => store%tail
      NULLIFY(substore%next)
    ELSE
      store%tail => substore
    ENDIF
    IF(.NOT. ASSOCIATED(store%head)) store%head => substore
    store%n_subs = store%n_subs + 1
    store%total_length = store%total_length + total_size

  END SUBROUTINE create_linked_substore


  SUBROUTINE create_empty_substore(store, total_size)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: substore, tmp1, tmp2
    INTEGER(i8), INTENT(IN) :: total_size
    TYPE(particle), POINTER :: current, prev
    INTEGER(i8) :: i_part

    !Nothing to do
    IF(total_size <= 0) RETURN

    !Allocate substore object
    ALLOCATE(substore)
    !Set length, nullify links
    substore%length = total_size
    NULLIFY(substore%prev, substore%next)
    !Allocate backing memory
    ALLOCATE(substore%store(total_size))

    substore%first_free_element = 1

    DO i_part = 1, total_size
      !Nullify pointers
      NULLIFY(substore%store(i_part)%prev, &
          substore%store(i_part)%next)
      !Set not-live state
      substore%store(i_part)%live = 0
    ENDDO

    IF(ASSOCIATED(store%tail)) THEN
      !Now link substore into store
      tmp1 => store%tail
      tmp2 => store%tail%next
      store%tail%next => substore
      tmp2 => store%tail%next

      substore%prev => store%tail
      tmp2 => substore%prev

      store%tail => substore
      tmp2 => store%tail
      NULLIFY(substore%next)
    ELSE
      store%tail => substore
    ENDIF
    IF(.NOT. ASSOCIATED(store%head)) store%head => substore
    store%n_subs = store%n_subs + 1
    store%total_length = store%total_length + total_size

  END SUBROUTINE create_empty_substore



  SUBROUTINE destroy_store(store)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: section

    section => store%head
    IF (.NOT. ASSOCIATED(section)) RETURN
    DO WHILE(ASSOCIATED(section))
      DEALLOCATE(section%store)
      section => section%next
    ENDDO
    NULLIFY(store%head, store%tail)
    store%n_subs = 0
    store%total_length = 0

  END SUBROUTINE destroy_store



  SUBROUTINE test_store(species)

    TYPE(particle_species), INTENT(IN) :: species
    TYPE(particle), POINTER :: current
    TYPE(particle_sub_store), POINTER :: sub
    INTEGER(i8) :: counta, countb, countc, i, a_count, cell_x, cell_y, countd, j
    REAL(num) :: idx, idy, part_x, part_y
    INTEGER :: ierr

    !First check general integrity
    counta = 0
    countb = 0
    i = 0
    sub => species%attached_list%store%head
    DO WHILE(ASSOCIATED(sub))
      i = i + 1
      IF(ASSOCIATED(sub%store)) THEN
        counta = counta + 1
      ELSE
        WRITE(100+rank, *) 'Bad substore', i
      ENDIF
      IF(ASSOCIATED(species%attached_list%store%tail, TARGET = sub)) &
        countb = i
      sub => sub%next
    ENDDO
    WRITE(100+rank, *) 'Number of sublists ', &
        species%attached_list%store%n_subs, counta
    FLUSH(100+rank)
    IF (species%attached_list%use_store) THEN
      IF(species%attached_list%store%n_subs .GT. 0 .AND. &
          counta /= species%attached_list%store%n_subs) &
          CALL abort_with_trace(7)
    ENDIF

    sub => species%attached_list%store%head

    WRITE(100+rank, *) "Head matches ", &
        ASSOCIATED(sub%head, TARGET=species%attached_list%head)

    j = 1
    DO WHILE(ASSOCIATED(sub))
      DO i = 1, sub%length
        IF(ASSOCIATED(species%attached_list%tail, &
            TARGET=sub%store(i))) &
            WRITE(100+rank, *) 'Tail is at index ', i, 'in', j
      ENDDO
      sub => sub%next
      j = j + 1
    ENDDO
    sub => species%attached_list%store%head
    j = 1
    DO WHILE(ASSOCIATED(sub))
      DO i = 1, sub%length
        IF(ASSOCIATED(species%attached_list%store%next_slot, &
            TARGET=sub%store(i))) &
            WRITE(100+rank, *) 'Next slot is at index ', i, 'in', j
      ENDDO
      sub => sub%next
      j = j + 1
    ENDDO

    counta = 0
    current => species%attached_list%head
    DO WHILE (ASSOCIATED(current))
      counta = counta + 1
      current => current%next
    END DO

    WRITE(100+rank, *)  "Checking partlists"
    WRITE(100+rank, *) counta, species%attached_list%count
    FLUSH(100+rank)

    countb = 0
    current => species%attached_list%head
    DO i=1, species%attached_list%count
      IF (ASSOCIATED(current)) THEN
        countb = countb + 1
        current => current%next
      ELSE
        CONTINUE
      ENDIF
      IF (ASSOCIATED(current, species%attached_list%tail)) CONTINUE
    ENDDO
    WRITE(100+rank, *) countb, ASSOCIATED(current)
    FLUSH(100+rank)

    countc = 0
    a_count = 0
    sub => species%attached_list%store%head
    DO WHILE(ASSOCIATED(sub))
      DO i = 1, sub%length
        IF (sub%store(i)%live > 0) THEN
        countc = countc + 1
        IF (ASSOCIATED(sub%store(i)%next)) &
            a_count = a_count + 1
       ENDIF
      ENDDO
      sub => sub%next
    ENDDO
    WRITE(100+rank, *)  countc, species%attached_list%count, a_count+1
    FLUSH(100+rank)

    idx = 1.0_num / dx
    idy = 1.0_num / dy


    WRITE(100+rank, *) 'Checking all positions'
    FLUSH(100+rank)
    current => species%attached_list%head

    countd = 0
    DO WHILE (ASSOCIATED(current))
      part_x  = current%part_pos(1)
      part_y  = current%part_pos(2)
      cell_x = part_x * idx
      cell_y = part_y * idy
      IF( part_x .GT. x_max_local  .OR. part_x .LT. x_min_local) THEN
        WRITE(100+rank, *) 'Error, particle out of range, x', cell_x, part_x, &
            x_min_local, x_max_local, current%live, current%id
        countd = countd + 1
      ENDIF
      IF(part_y .GT. y_max_local  .OR. part_y .LT. y_min_local) THEN
        WRITE(100+rank, *) 'Error, particle out of range, y', cell_y
        countd = countd + 1
      ENDIF
      current => current%next

    END DO
    WRITE(100+rank, *) "Positions Done"
    WRITE(100+rank, *) "Checking Counts"

    IF (species%attached_list%use_store) THEN
      IF(counta /= countb .OR. countb /= countc) THEN
        CALL abort_with_trace(1)
      ENDIF
      IF(counta /= species%attached_list%count) &
        CALL abort_with_trace(2)
      IF(countc > 0 .AND. countb /= a_count+1) THEN
          CALL abort_with_trace(3)
      ENDIF
    ELSE
      IF(counta /= countb .OR. countb /= species%attached_list%count) &
          CALL abort_with_trace(4)
    ENDIF

    IF(countd /= 0) &
        CALL abort_with_trace(5)

    FLUSH(100+rank)

  END SUBROUTINE test_store


  SUBROUTINE test_list_positions(list, inside_is_err)

   TYPE(particle_list), INTENT(IN) :: list
   TYPE(particle), POINTER :: current
   INTEGER(i8) :: countd
   LOGICAL, INTENT(IN) :: inside_is_err
   REAL(num) :: part_x, part_y
   INTEGER :: ierr

   countd = 0
   current => list%head

   DO WHILE (ASSOCIATED(current))
      countd = countd + 1
      part_x  = current%part_pos(1)
      part_y  = current%part_pos(2)
      WRITE(100+rank, *) 'Particle ',countd, ' at ', part_x, part_y
      IF(inside_is_err .AND. part_x < x_max_local  .AND. part_x > x_min_local &
          .AND. part_y < y_max_local .AND. part_y > y_min_local) THEN
          WRITE(100+rank, *) 'Error, particle in domain', countd, current%live
          WRITE(100+rank, *) part_x, part_y
      ELSE IF(.NOT. inside_is_err .AND. (part_x >= x_max_local  .OR. part_x <= x_min_local &
          .AND. part_y >= y_max_local .OR. part_y <= y_min_local)) THEN
          WRITE(100+rank, *) 'Error, particle in domain', countd, current%live
          WRITE(100+rank, *) part_x, part_y
      ENDIF
      current => current%next
   END DO
   WRITE(100+rank, *) "Positions Done"

  END SUBROUTINE test_list_positions


  SUBROUTINE test_store_walk(store)

    TYPE(particle_store), INTENT(IN) :: store
    TYPE(particle), POINTER :: current, prev
    TYPE(particle_sub_store), POINTER :: sub
    INTEGER(i8) :: counta, countb, countc, i, a_count, cell_x, cell_y, countd, j
    INTEGER :: ierr

    WRITE(100+rank, *) 'Walking'
    !First check general integrity
    counta = 0
    countb = 0
    i = 0
    sub => store%head
    DO WHILE(ASSOCIATED(sub))
      i = i + 1
      IF(ASSOCIATED(sub%store)) THEN
        counta = counta + 1
      ELSE
        WRITE(100+rank, *) 'Bad substore', i
      ENDIF
      IF(ASSOCIATED(store%tail, TARGET = sub)) &
        countb = i
      sub => sub%next
    ENDDO
    WRITE(100+rank, *) 'Number of sublists ', &
        store%n_subs, counta
    FLUSH(100+rank)
    IF(store%n_subs > 0 .AND. counta /= store%n_subs) &
          CALL abort_with_trace(11)

    sub => store%head
    j = 1
    DO WHILE(ASSOCIATED(sub))
      DO i = 1, sub%length
        IF(ASSOCIATED(store%next_slot, &
            TARGET=sub%store(i))) &
            WRITE(100+rank, *) 'Next slot is at index ', i, 'in', j
      ENDDO
      sub => sub%next
      j = j + 1
    ENDDO

    counta = 0
    current => store%head%head
    DO WHILE (ASSOCIATED(current))
      counta = counta + 1
      current => current%next
    END DO

    WRITE(100+rank, *)  "Checking store"
    WRITE(100+rank, *) counta
    FLUSH(100+rank)

    countc = 0
    a_count = 0
    sub => store%head
    DO WHILE(ASSOCIATED(sub))
      DO i = 1, sub%length
        IF (sub%store(i)%live > 0) THEN
        countc = countc + 1
        IF (ASSOCIATED(sub%store(i)%next)) &
            a_count = a_count + 1
       ENDIF
      ENDDO
      sub => sub%next
    ENDDO
    WRITE(100+rank, *)  countc, a_count+1

    current => store%head%head
    NULLIFY(prev)
    i = 1
    WRITE(100+rank, *) "Checking prevs"
    DO WHILE (ASSOCIATED(current))
      IF(ASSOCIATED(prev) .AND. .NOT. ASSOCIATED(current%prev, TARGET=prev)) THEN
        WRITE(100+rank, *) "Bad prev in walk at ", i
        FLUSH(100+rank)
        CALL abort_with_trace(12)
      ENDIF
      prev => current
      current => current%next
      i = i + 1
    END DO


    WRITE(100+rank, *) 'Walk complete'
    FLUSH(100+rank)
    IF(countc /= 0 .AND. counta /= countc .OR. counta /= a_count+1) &
          CALL abort_with_trace(13)


  END SUBROUTINE test_store_walk



  SUBROUTINE create_unsafe_partlist(partlist, a_particle, n_elements)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: a_particle
    INTEGER(i8), INTENT(IN) :: n_elements
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart

    CALL create_empty_partlist(partlist)

    partlist%safe = .FALSE.
    current => a_particle
    ipart = 1
    DO WHILE (ASSOCIATED(current) .AND. ipart < n_elements)
      ipart = ipart+1
      current => current%next
    ENDDO
    partlist%head => a_particle
    partlist%tail => current
    partlist%count = ipart

  END SUBROUTINE create_unsafe_partlist



  SUBROUTINE create_unsafe_partlist_by_tail(partlist, head, tail)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: head, tail
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart

    CALL create_empty_partlist(partlist)

    partlist%safe = .FALSE.
    partlist%head => head
    partlist%tail => tail

    current => head
    ipart = 0
    DO WHILE (ASSOCIATED(current))
      ipart = ipart+1
      current => current%next
      IF (ASSOCIATED(current)) THEN
        IF (ASSOCIATED(current%prev, TARGET=tail)) EXIT
      ENDIF
    ENDDO

    partlist%count = ipart

  END SUBROUTINE create_unsafe_partlist_by_tail



  SUBROUTINE create_allocated_partlist(partlist, n_elements, use_store_in)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER(i8), INTENT(IN) :: n_elements
    LOGICAL, INTENT(IN), OPTIONAL :: use_store_in
    LOGICAL :: use_store
    TYPE(particle), POINTER :: new_particle
    INTEGER(i8) :: ipart

    IF(.NOT. PRESENT(use_store_in)) THEN
      use_store = .FALSE.
    ELSE
      use_store = use_store_in
    ENDIF

    IF (use_store) THEN
      CALL create_particle_store(partlist, n_elements)
    !  IF(n_elements > 0) THEN
    !    !TODO and if not condition??
    !    !TODO actually the following MAY NOT be list tail, should do better
    !    partlist%tail => partlist%store%next_slot%prev
    !    partlist%head => partlist%store%head%head
    !  ELSE
    !    NULLIFY(partlist%head, partlist%tail)
    !  ENDIF
      partlist%count = n_elements
    ELSE
      CALL create_empty_partlist(partlist)

      DO ipart = 0, n_elements-1
        CALL create_particle(new_particle)
        CALL add_particle_to_partlist(partlist, new_particle)
        NULLIFY(new_particle)
      ENDDO
    ENDIF

    partlist%use_store = use_store

  END SUBROUTINE create_allocated_partlist



  !Walk the particle store and regenerate a linked list
  !from the slots which hold live particles
  !to make list a valid linked list again
  SUBROUTINE relink_partlist(list, recount)

    TYPE(particle_list), INTENT(INOUT) :: list
    TYPE(particle), POINTER :: current, previous
    TYPE(particle_sub_store), POINTER :: sub
    INTEGER(i8) ::cnt, icurr
    LOGICAL, INTENT(IN) :: recount

    NULLIFY(previous)

    cnt = 0
    sub => list%store%head
    DO WHILE(ASSOCIATED(sub))
      DO icurr = 1, sub%length
      current => sub%store(icurr)
      IF (current%live > 0) THEN
        cnt = cnt + 1
        IF (ASSOCIATED(previous)) THEN
          previous%next => current
        ELSE
          list%head => current
        ENDIF
        current%prev => previous
        previous => current
      ENDIF
      ENDDO
      sub => sub%next
    ENDDO

    IF (ASSOCIATED(previous)) THEN
      NULLIFY(previous%next)
      list%tail => previous
    ENDIF
    list%store%head%head => list%head

    IF(recount) list%count = cnt

  END SUBROUTINE relink_partlist


  !Increment the position of the next slot in store
  !If this overflows, then make space
  !First try compacting the list
  !If that doesn't help, rellocate larger
  !THIS ROUTINE MAY INVALIDATE POINTERS!!
  !ANY CALLER MUST CHECK LIST IS UNCHANGED OR REPOINT!!
  !Once add multiple linked sub_stores and/or more sophisticated
  !trcking of free slots, this must get more clever
  SUBROUTINE increment_next_free_element(list)

    TYPE(particle_list), INTENT(INOUT) :: list
    INTEGER(i8) :: new_size, i

    IF(list%store%tail%first_free_element >= list%store%tail%length -1) THEN
      !First resort: compact store
      IF(list%count > 0) CALL compact_backing_store(list%store, list)
      IF(list%store%tail%first_free_element >= list%store%tail%length -1) THEN
        !Compacting insufficient, must grow
     WRITE(100+rank, *) "Adding new sub"
        CALL create_empty_substore(list%store, sublist_size)
      ENDIF
    ELSE
      !Do this only if we've not created anything new
      list%store%tail%first_free_element = list%store%tail%first_free_element + 1
    ENDIF
    !Do this always, as if we grew, we have a nice empty, valid substore
    list%store%next_slot => &
        list%store%tail%store(list%store%tail%first_free_element)

  END SUBROUTINE increment_next_free_element



  !Update the next free element, e.g. after store gets compacted
  SUBROUTINE update_store_endpoint(list, last_full_slot)

    TYPE(particle_list), INTENT(INOUT) :: list
    INTEGER(i8) :: last_full_slot

    list%store%tail%first_free_element = last_full_slot + 1
    list%store%next_slot => list%store%tail%store(last_full_slot + 1)

  END SUBROUTINE



  SUBROUTINE create_filled_partlist(partlist, data_in, n_elements)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    REAL(num), DIMENSION(:), INTENT(IN) :: data_in
    INTEGER(i8), INTENT(IN) :: n_elements
    TYPE(particle), POINTER :: new_particle
    INTEGER(i8) :: ipart, cpos = 0

    CALL create_empty_partlist(partlist)

    DO ipart = 0, n_elements-1
      ALLOCATE(new_particle)
      cpos = ipart*nvar+1
      CALL unpack_particle(data_in(cpos:cpos+nvar-1), new_particle)
#ifdef PARTICLE_DEBUG
      new_particle%processor = rank
#endif
      CALL add_particle_to_partlist(partlist, new_particle)
      NULLIFY(new_particle)
    ENDDO

  END SUBROUTINE create_filled_partlist



  FUNCTION test_partlist(partlist)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: current
    INTEGER :: test_partlist
    INTEGER(i8) :: test_ct

    test_partlist = 0
    test_ct = 0

    ! Empty list is OK
    IF (.NOT. ASSOCIATED(partlist%head) &
        .AND. .NOT. ASSOCIATED(partlist%tail)) THEN
      test_partlist = 0
      RETURN
    ENDIF

    ! List with head or tail but not both is broken
    IF (.NOT. ASSOCIATED(partlist%head) &
        .OR. .NOT. ASSOCIATED(partlist%tail)) THEN
      test_partlist = -1
      RETURN
    ENDIF

    ! Having head and tail elements which are not the end of a list are OK for
    ! unsafe partlists
    IF (ASSOCIATED(partlist%head%prev) .AND. partlist%safe) &
        test_partlist = IOR(test_partlist, 1)
    IF (ASSOCIATED(partlist%tail%next) .AND. partlist%safe) &
        test_partlist = IOR(test_partlist, 2)

    ! Since we don't KNOW that count is OK (that's what we're checking)
    ! Have to check both for end of list and for having reached the tail item
    current => partlist%head
    DO WHILE (ASSOCIATED(current))
      test_ct = test_ct+1
      current => current%next
      IF (ASSOCIATED(current)) THEN
        ! This tests if we've just jumped to the tail element
        ! Allows testing of unsafe partlists
        IF (ASSOCIATED(current%prev, TARGET=partlist%tail)) EXIT
      ENDIF
    ENDDO

    IF (test_ct /= partlist%count) test_partlist = IOR(test_partlist, 4)

  END FUNCTION test_partlist



  SUBROUTINE destroy_partlist(partlist)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: new_particle, next
    INTEGER(i8) :: ipart

    IF(partlist%use_store) THEN
      CALL destroy_store(partlist%store)
    ELSE
      ! Go through list and delete all the particles in the list
      new_particle => partlist%head
      ipart = 0
      DO WHILE (ipart < partlist%count)
        next => new_particle%next
        DEALLOCATE(new_particle)
        new_particle => next
        ipart = ipart+1
      ENDDO

      CALL create_empty_partlist(partlist)
    ENDIF

  END SUBROUTINE destroy_partlist



  SUBROUTINE copy_partlist(partlist1, partlist2)

    TYPE(particle_list), INTENT(INOUT) :: partlist1, partlist2

    partlist2%head => partlist1%head
    partlist2%tail => partlist1%tail
    partlist2%count = partlist1%count
    partlist2%id_update = partlist1%id_update

  END SUBROUTINE copy_partlist




  SUBROUTINE append_partlist(list, newlist, ignore_live_in)

    TYPE(particle_list), INTENT(INOUT) :: list, newlist
    LOGICAL, INTENT(IN), OPTIONAL :: ignore_live_in
    LOGICAL :: ignore_live

    IF(newlist%count == 0) RETURN
    IF (.NOT. list%safe .OR. .NOT. newlist%safe) THEN
      IF (rank == 0) &
          PRINT *, 'Unable to append partlists because one is not safe'
      RETURN
    ENDIF

    IF(PRESENT(ignore_live_in)) THEN
      ignore_live = ignore_live_in
    ELSE
      ignore_live = .FALSE.
    ENDIF

    !Do the appending
    !TODO this misses the case of a linked list and a store backed newlist
    IF (.NOT. (list%use_store .OR. newlist%use_store)) THEN
      IF (ASSOCIATED(list%tail)) THEN
        list%tail%next => newlist%head
      ELSE
        list%head => newlist%head
      ENDIF
      IF (ASSOCIATED(newlist%head)) newlist%head%prev => list%tail
      IF (ASSOCIATED(newlist%tail)) list%tail => newlist%tail
      list%count = list%count + newlist%count
      list%id_update = list%id_update + newlist%id_update

    ELSE
      CALL add_partlist_to_list_and_store(list, newlist, ignore_live)
    ENDIF
    !Clean up newlist
    CALL create_empty_partlist(newlist)

  END SUBROUTINE append_partlist



  !Take a list and append its content to list-with-store
  SUBROUTINE add_partlist_to_list_and_store(list, newlist, override_live)

    TYPE(particle_list), INTENT(INOUT) :: newlist, list
    LOGICAL, INTENT(IN) :: override_live !Override any live states in newlist
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: add_count, next_index !TODO both temporary diagnostics

    IF(newlist%count < 1) RETURN

    current => newlist%head
    NULLIFY(next)
    IF ( .NOT. ASSOCIATED(list%head) .AND. newlist%count > 0) THEN
      !List was empty, so hook up the head
      list%head => list%store%next_slot
    ENDIF

    add_count = 0
    DO WHILE(ASSOCIATED(current))
      !Only consider live particles, unless overrriding
      IF ( override_live .OR. current%live == 1) THEN

        !Diagnostic only...
        next_index = list%store%head%first_free_element
        !End diagnostic

        CALL create_particle_in_list(next, list)
        next%live = 1
        CALL copy_particle(current, next)
        add_count = add_count + 1
      ENDIF

      current => current%next
    END DO
    list%id_update = list%id_update + newlist%id_update


  END SUBROUTINE add_partlist_to_list_and_store



  SUBROUTINE add_particle_to_partlist(partlist, new_particle)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: new_particle
    INTEGER :: ierr

    ! Note that this will work even if you are using an unsafe particle list
    ! BE CAREFUL if doing so, it can cause unexpected behaviour

    ! if (!particle) return;
    IF (.NOT. ASSOCIATED(new_particle)) RETURN
    IF (partlist%use_store) CALL abort_code(8)
    NULLIFY(new_particle%next, new_particle%prev)

    ! Add particle count
    partlist%count = partlist%count + 1
    partlist%id_update = 1
    IF (.NOT. ASSOCIATED(partlist%tail)) THEN
      ! partlist is empty
      partlist%head => new_particle
      partlist%tail => new_particle
      RETURN
    ENDIF

    partlist%tail%next => new_particle
    new_particle%prev => partlist%tail
    NULLIFY(new_particle%next)
    partlist%tail => new_particle

  END SUBROUTINE add_particle_to_partlist




  SUBROUTINE remove_particle_from_partlist(partlist, a_particle, &
      fromstore_in, nocopy_in)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: a_particle, tmp_particle
    LOGICAL, INTENT(IN), OPTIONAL :: fromstore_in, nocopy_in
    LOGICAL :: fromstore, nocopy

    ! Note that this will work even if you are using an unsafe particle list
    ! BE CAREFUL if doing so, it can cause unexpected behaviour

    IF (PRESENT(fromstore_in)) THEN
      fromstore = fromstore_in
    ELSE
      fromstore = .TRUE.
    ENDIF
    IF (PRESENT(nocopy_in)) THEN
      nocopy = nocopy_in
    ELSE
      nocopy = .FALSE.
    ENDIF


    ! Check whether particle is head or tail of list and unlink
    IF (ASSOCIATED(partlist%head, TARGET=a_particle)) THEN
      partlist%head => a_particle%next
      IF (partlist%use_store) &
          partlist%store%head%head => a_particle%next
    ENDIF

    IF (ASSOCIATED(partlist%tail, TARGET=a_particle)) &
        partlist%tail => a_particle%prev

    ! Link particles on either side together
    IF (ASSOCIATED(a_particle%next)) a_particle%next%prev => a_particle%prev
    IF (ASSOCIATED(a_particle%prev)) a_particle%prev%next => a_particle%next

    NULLIFY(a_particle%next, a_particle%prev)
    IF (fromstore) a_particle%live = 0
    ! Decrement counter
    partlist%count = partlist%count-1

    IF (partlist%use_store .AND. .NOT. nocopy) THEN
      !If a_particle is in a store, make a copy
      !Then what comes back is a valid, FREE particle
      CALL create_particle(tmp_particle)
      CALL copy_particle(a_particle, tmp_particle)
      !Return a live particle
      tmp_particle%live = 1
      a_particle => tmp_particle
    ENDIF


  END SUBROUTINE remove_particle_from_partlist



  SUBROUTINE pack_particle(array, a_particle)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    TYPE(particle), POINTER :: a_particle
    INTEGER(i8) :: cpos

    cpos = 1
    array(cpos:cpos+c_ndims-1) = a_particle%part_pos
    cpos = cpos+c_ndims
    array(cpos:cpos+2) = a_particle%part_p
    cpos = cpos+3
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    array(cpos) = a_particle%weight
    cpos = cpos+1
#endif
#ifdef DELTAF_METHOD
    array(cpos) = a_particle%pvol
    cpos = cpos+1
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    array(cpos) = a_particle%charge
    array(cpos+1) = a_particle%mass
    cpos = cpos+2
#endif
#ifdef PARTICLE_DEBUG
    array(cpos) = REAL(a_particle%processor, num)
    array(cpos+1) = REAL(a_particle%processor_at_t0, num)
    cpos = cpos+2
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    array(cpos) = REAL(a_particle%id, num)
    cpos = cpos+1
#endif
#ifdef COLLISIONS_TEST
    array(cpos) = REAL(a_particle%coll_count, num)
    cpos = cpos+1
#endif
#ifdef PHOTONS
    array(cpos) = a_particle%optical_depth
    array(cpos+1) = a_particle%particle_energy
    cpos = cpos+2
#ifdef TRIDENT_PHOTONS
    array(cpos) = a_particle%optical_depth_tri
    cpos = cpos+1
#endif
#endif

  END SUBROUTINE pack_particle



  SUBROUTINE unpack_particle(array, a_particle)

    REAL(num), DIMENSION(:), INTENT(IN) :: array
    TYPE(particle), POINTER :: a_particle
    INTEGER(i8) :: cpos

    cpos = 1
    a_particle%part_pos = array(cpos:cpos+c_ndims-1)
    cpos = cpos+c_ndims
    a_particle%part_p = array(cpos:cpos+2)
    cpos = cpos+3
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    a_particle%weight = array(cpos)
    cpos = cpos+1
#endif
#ifdef DELTAF_METHOD
    a_particle%pvol = array(cpos)
    cpos = cpos+1
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    a_particle%charge = array(cpos)
    a_particle%mass = array(cpos+1)
    cpos = cpos+2
#endif
#ifdef PARTICLE_DEBUG
    a_particle%processor = rank
    a_particle%processor_at_t0 = NINT(array(cpos+1))
    cpos = cpos+2
#endif
#ifdef PARTICLE_ID4
    a_particle%id = NINT(array(cpos))
    cpos = cpos+1
#elif PARTICLE_ID
    a_particle%id = NINT(array(cpos),i8)
    cpos = cpos+1
#endif
#ifdef COLLISIONS_TEST
    a_particle%coll_count = NINT(array(cpos))
    cpos = cpos+1
#endif
#ifdef PHOTONS
    a_particle%optical_depth = array(cpos)
    a_particle%particle_energy = array(cpos+1)
    cpos = cpos+2
#ifdef TRIDENT_PHOTONS
    a_particle%optical_depth_tri = array(cpos)
    cpos = cpos+1
#endif
#endif
    a_particle%live = 1 !TODO actually send this info

  END SUBROUTINE unpack_particle



  SUBROUTINE init_particle(new_particle)

    TYPE(particle), POINTER :: new_particle

    new_particle%part_p = 0.0_num
    new_particle%part_pos = 0.0_num
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    new_particle%weight = 0.0_num
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    new_particle%charge = 0.0_num
    new_particle%mass = 0.0_num
#endif
#ifdef PARTICLE_DEBUG
    new_particle%processor = 0
    new_particle%processor_at_t0 = 0
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    new_particle%id = 0
#endif
#ifdef COLLISIONS_TEST
    new_particle%coll_count = 0
#endif
#ifdef PHOTONS
    ! This assigns an optical depth to newly created particle
    new_particle%particle_energy = 0.0_num
    new_particle%optical_depth = LOG(1.0_num / (1.0_num - random()))
#ifdef TRIDENT_PHOTONS
    new_particle%optical_depth_tri = LOG(1.0_num / (1.0_num - random()))
#endif
#endif
    new_particle%live = 1

  END SUBROUTINE init_particle



  SUBROUTINE copy_particle(src_particle, new_particle)

    TYPE(particle), POINTER :: src_particle, new_particle

    new_particle%part_p = src_particle%part_p
    new_particle%part_pos = src_particle%part_pos
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    new_particle%weight = src_particle%weight
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    new_particle%charge = src_particle%charge
    new_particle%mass = src_particle%mass
#endif
#ifdef PARTICLE_DEBUG
    new_particle%processor = src_particle%processor
    new_particle%processor_at_t0 = src_particle%processor_at_t0
#endif
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    new_particle%id = src_particle%id
#endif
#ifdef COLLISIONS_TEST
    new_particle%coll_count = src_particle%coll_count
#endif
#ifdef PHOTONS
    ! This assigns an optical depth to newly created particle
    new_particle%particle_energy = src_particle%particle_energy
    new_particle%optical_depth = src_particle%optical_depth
#ifdef TRIDENT_PHOTONS
    new_particle%optical_depth_tri = src_particle%optical_depth_tri
#endif
#endif

    new_particle%live = src_particle%live

  END SUBROUTINE copy_particle




  SUBROUTINE create_particle(new_particle)

    TYPE(particle), POINTER :: new_particle

    ALLOCATE(new_particle)
    CALL init_particle(new_particle)

  END SUBROUTINE create_particle



  SUBROUTINE create_particle_in_list(new_particle, list)

    TYPE(particle), POINTER, INTENT(INOUT) :: new_particle
    TYPE(particle_list), INTENT(INOUT) :: list

    IF (list%use_store) THEN
      new_particle => list%store%next_slot
      CALL init_particle(new_particle)
      new_particle%live = 1
      new_particle%prev => list%tail
      NULLIFY(new_particle%next)
      list%tail => new_particle
      IF(ASSOCIATED(new_particle%prev)) new_particle%prev%next => new_particle

      CALL increment_next_free_element(list)
      new_particle => list%tail !In case reallocation occurred
    ELSE
      CALL create_particle(new_particle)
      new_particle%prev => list%tail
      NULLIFY(new_particle%next)
      list%tail => new_particle
      IF(ASSOCIATED(new_particle%prev)) new_particle%prev%next => new_particle
    ENDIF
    IF (.NOT. ASSOCIATED(list%head)) list%head => new_particle
    list%count = list%count + 1

  END SUBROUTINE create_particle_in_list



  SUBROUTINE display_particle(a_particle)

    TYPE(particle), POINTER :: a_particle

    PRINT *, 'Position', a_particle%part_pos
    PRINT *, 'Momentum', a_particle%part_p

  END SUBROUTINE display_particle



  FUNCTION compare_particles(part1, part2)

    TYPE(particle), POINTER :: part1, part2
    LOGICAL :: compare_particles

    compare_particles = .TRUE.
    IF (MAXVAL(ABS(part1%part_pos-part2%part_pos)) > c_tiny) &
        compare_particles = .FALSE.
    IF (MAXVAL(ABS(part1%part_p - part2%part_p)) > c_tiny) &
        compare_particles = .FALSE.

#ifndef PER_SPECIES_WEIGHT
    IF (ABS(part1%weight - part2%weight) > c_tiny) &
        compare_particles = .FALSE.
#endif

#ifdef PER_PARTICLE_CHARGE_MASS
    IF (ABS(part1%charge - part2%charge) > c_tiny) &
        compare_particles = .FALSE.
    IF (ABS(part1%mass - part2%mass) > c_tiny) &
        compare_particles = .FALSE.
#endif

    IF (.NOT. compare_particles) THEN
      CALL display_particle(part1)
      CALL display_particle(part2)
    ENDIF

  END FUNCTION compare_particles



  FUNCTION test_packed_particles(partlist, array, npart_in_data)

    TYPE(particle_list), INTENT(IN) :: partlist
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER(i8), INTENT(IN) :: npart_in_data
    TYPE(particle), POINTER :: current
    TYPE(particle), POINTER :: a_particle
    LOGICAL :: test_packed_particles
    INTEGER(i8) :: ipart

    test_packed_particles = .FALSE.

    IF (npart_in_data * nvar /= SIZE(array)) THEN
      PRINT *, 'Size of data array does not match specified on', rank, &
          npart_in_data, SIZE(array)
      RETURN
    ENDIF
    IF (partlist%count /= npart_in_data) THEN
      PRINT *, 'Size of data array does not match partlist on', rank
      RETURN
    ENDIF

    ALLOCATE(a_particle)

    current => partlist%head
    DO ipart = 0, npart_in_data-1
      CALL unpack_particle(array(ipart*nvar+1:(ipart+1)*nvar), a_particle)
      IF (.NOT. compare_particles(a_particle, current)) THEN
        PRINT *, 'BAD PARTICLE ', ipart, 'on', rank
        RETURN
      ENDIF
      current => current%next
    ENDDO

    DEALLOCATE(a_particle)

    test_packed_particles = .TRUE.

  END FUNCTION test_packed_particles



  SUBROUTINE partlist_send_nocount(partlist, dest)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: dest
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    INTEGER :: ipart, nsend, cpos
    TYPE(particle), POINTER :: current

    nsend = INT(partlist%count) * nvar
    ALLOCATE(array(nsend))
    array = 0.0_num

    current => partlist%head
    ipart = 0
    cpos = 0
    DO WHILE (ipart < partlist%count)
      cpos = ipart * nvar + 1
      CALL pack_particle(array(cpos:cpos+nvar-1), current)
      ipart = ipart + 1
      current => current%next
    ENDDO

    CALL MPI_SEND(array, nsend, mpireal, dest, tag, comm, errcode)

    DEALLOCATE(array)

  END SUBROUTINE partlist_send_nocount



  SUBROUTINE partlist_send(partlist, dest)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: dest
    INTEGER(i8) :: send_buf(2)

    send_buf(1) = partlist%count
    send_buf(2) = partlist%id_update

    CALL MPI_SEND(send_buf, 2, MPI_INTEGER8, dest, tag, comm, errcode)

    CALL partlist_send_nocount(partlist, dest)

  END SUBROUTINE partlist_send



  SUBROUTINE partlist_recv_nocount(partlist, src, count)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: src
    INTEGER(i8), INTENT(IN) :: count
    INTEGER :: nrecv
    REAL(num), DIMENSION(:), ALLOCATABLE :: array

    CALL create_empty_partlist(partlist)

    nrecv = INT(count) * nvar
    ALLOCATE(array(nrecv))
    array = 0.0_num

    CALL MPI_RECV(array, nrecv, mpireal, src, tag, comm, status, errcode)
    CALL create_filled_partlist(partlist, array, count)

    DEALLOCATE(array)

  END SUBROUTINE partlist_recv_nocount



  SUBROUTINE partlist_recv(partlist, src)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: src
    INTEGER(i8) :: count, recv_buf(2)

    recv_buf = 0
    CALL MPI_RECV(recv_buf, 2, MPI_INTEGER8, src, tag, comm, status, errcode)
    count = recv_buf(1)
    partlist%id_update = partlist%id_update + INT(recv_buf(2))

    CALL partlist_recv_nocount(partlist, src, count)

  END SUBROUTINE partlist_recv



  SUBROUTINE partlist_sendrecv(partlist_send, partlist_recv, dest, src)

    TYPE(particle_list), INTENT(INOUT) :: partlist_send, partlist_recv
    INTEGER, INTENT(IN) :: dest, src
    REAL(num), DIMENSION(:), ALLOCATABLE :: data_send, data_recv
    INTEGER(i8) :: cpos = 0, ipart = 0
    INTEGER(i8) :: npart_recv, send_buf(2), recv_buf(2)
    INTEGER :: nsend, nrecv
    TYPE(particle), POINTER :: current

    ! This subroutine doesn't try to use memory efficient buffering, it sends
    ! all the particles at once. This should work for boundary calls, but
    ! don't try it for any other reason

    recv_buf = 0
    send_buf(1) = partlist_send%count
    send_buf(2) = partlist_send%id_update

    CALL MPI_SENDRECV(send_buf, 2, MPI_INTEGER8, dest, tag, recv_buf, 2, &
        MPI_INTEGER8, src, tag, comm, status, errcode)

    npart_recv = recv_buf(1)
    nsend = INT(send_buf(1)) * nvar
    nrecv = INT(npart_recv) * nvar
    partlist_recv%id_update = partlist_recv%id_update + INT(recv_buf(2))

    ! Copy the data for the particles into a buffer
    ALLOCATE(data_send(nsend))
    ALLOCATE(data_recv(nrecv))

    ! Pack particles to send into buffer
    current => partlist_send%head
    ipart = 0
    DO WHILE (ipart < partlist_send%count)
      cpos = ipart * nvar + 1
      CALL pack_particle(packed_particle_data, current)
      data_send(cpos:cpos+nvar-1) = packed_particle_data
      ipart = ipart + 1
      current => current%next
    ENDDO

    ! No longer need the sending partlist, so destroy it to save some memory
    !Sendlist no longer owns the memory
    CALL destroy_partlist(partlist_send)

    ! Actual MPI commands
    CALL MPI_SENDRECV(data_send, nsend, mpireal, dest, tag, &
        data_recv, nrecv, mpireal, src, tag, comm, status, errcode)

    DEALLOCATE(data_send)
    CALL create_filled_partlist(partlist_recv, data_recv, npart_recv)
    DEALLOCATE(data_recv)

  END SUBROUTINE partlist_sendrecv



  SUBROUTINE add_particle_to_list(part, list)

    TYPE(particle), POINTER :: part
    TYPE(pointer_list) :: list
    TYPE(pointer_item), POINTER :: item

    ALLOCATE(item)
    item%part => part
    NULLIFY(item%next)

    list%tail%next => item
    list%tail => item

  END SUBROUTINE add_particle_to_list



  SUBROUTINE generate_particle_ids(partlist)

    USE constants

    TYPE(particle_list) :: partlist
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    INTEGER(i8), ALLOCATABLE :: nid_all(:)
    INTEGER(i8) :: nid, part_id
    INTEGER :: i, id_update
    TYPE(particle), POINTER :: current
    TYPE(pointer_list) :: idlist
    TYPE(pointer_item), POINTER :: idcurrent, idnext

    id_update = partlist%id_update

    CALL MPI_ALLREDUCE(id_update, partlist%id_update, 1, MPI_INTEGER, &
        MPI_MAX, comm, errcode)

    IF (partlist%id_update == 0) RETURN

    ALLOCATE(idlist%head)
    idlist%tail => idlist%head
    NULLIFY(idlist%head%next)
    NULLIFY(idlist%head%part)

    ! Scan through particle list and identify particles which need
    ! an ID to be assigned.
    nid = 0
    current => partlist%head
    DO WHILE(ASSOCIATED(current))
      IF (current%id == 0) THEN
        nid = nid + 1
        CALL add_particle_to_list(current, idlist)
      ENDIF
      current => current%next
    ENDDO

    ALLOCATE(nid_all(nproc))

    CALL MPI_ALLGATHER(nid, 1, MPI_INTEGER8, nid_all, 1, MPI_INTEGER8, &
        comm, errcode)

    ! Count number of particles on ranks zero to rank-1
    nid = 0
    DO i = 1, rank
      nid = nid + nid_all(i)
    ENDDO
    part_id = particles_max_id + nid

    ! Count remaining particles
    DO i = rank+1, nproc
      nid = nid + nid_all(i)
    ENDDO

    particles_max_id = particles_max_id + nid

    DEALLOCATE(nid_all)

    ! Number each particle with a unique id
    idcurrent => idlist%head%next
    DO WHILE(ASSOCIATED(idcurrent))
      part_id = part_id + 1
#if PARTICLE_ID
      idcurrent%part%id = part_id
#else
      idcurrent%part%id = INT(part_id,i4)
#endif
      idnext => idcurrent%next
      DEALLOCATE(idcurrent)
      idcurrent => idnext
    ENDDO

    DEALLOCATE(idlist%head)

    partlist%id_update = 0
#endif

  END SUBROUTINE generate_particle_ids


  !TODO Currently works with only a single store
  SUBROUTINE compact_backing_store(store, list)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: section
    TYPE(particle_list), INTENT(INOUT) :: list
    TYPE(particle), POINTER :: original, new
    INTEGER(i8) :: i, last_placed, count, moved, j

    !TODO Should we compact between stores?
    !TODO remove moved diagnostic
    section => store%head
    last_placed = 1
    count = 0
    moved = 0

    !TODO is there any point compacting the segments or should we
    ! just remove them when they become empty?
    !We can compact the last chunk, but if we do the above this might be pointless too

    DO j = 1, store%n_subs
      !We're doing before we update first_free_element, so it may contain particle
      last_placed = 1

      DO i = 1, section%first_free_element
        IF( i > section%length) EXIT
        IF(section%store(i)%live > 0) THEN
          count = count + 1
        ENDIF
        IF(section%store(i)%live > 0 .AND. i - last_placed > 1) THEN
          moved = moved + 1
          original => section%store(i)
          new => section%store(last_placed + 1)
          CALL copy_particle(original, new)
          last_placed = last_placed + 1
        ELSE IF(section%store(i)%live > 0) THEN
          last_placed = i
        ENDIF
      ENDDO

      DO i = last_placed+1 , section%length
        new => section%store(i)
        new%live = 0
      ENDDO
      section%head => section%store(1)
      section => section%next
    ENDDO
    CALL relink_partlist(list, .FALSE.)
    !CALL update_store_endpoint(list, last_placed-1)
    CALL update_store_endpoint(list, last_placed)


  END SUBROUTINE compact_backing_store


  !Completely relink all particle lists
  !Slow, for test purposes
  SUBROUTINE relink_partlists

    INTEGER(i8) :: ispecies
    TYPE(particle_list), POINTER :: list

    DO ispecies = 1, n_species

      list => species_list(ispecies)%attached_list

      CALL relink_partlist(list, .TRUE.)

    ENDDO


  END SUBROUTINE

END MODULE partlist
