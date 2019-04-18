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
  USE particle_id_hash_mod
#ifdef PHOTONS
  USE random_generator
#endif

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

  SUBROUTINE set_partlist_size

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
#ifdef WORK_DONE_INTEGRATED
    nvar = nvar+6
#endif
    ! Persistent IDs
    IF (any_persistent_subset) nvar = nvar+1

  END SUBROUTINE set_partlist_size



  SUBROUTINE setup_partlists

    LOGICAL :: old_any_persistent_subset

    old_any_persistent_subset = any_persistent_subset
    any_persistent_subset = .TRUE.

    CALL set_partlist_size

    any_persistent_subset = old_any_persistent_subset

    ALLOCATE(packed_particle_data(nvar))

  END SUBROUTINE setup_partlists



  SUBROUTINE deallocate_partlists

    INTEGER :: stat

    IF (ALLOCATED(packed_particle_data)) &
        DEALLOCATE(packed_particle_data, STAT=stat)

  END SUBROUTINE deallocate_partlists



  SUBROUTINE create_empty_partlist(partlist, use_store_in, holds_copies)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    LOGICAL, INTENT(IN), OPTIONAL :: use_store_in, holds_copies
    LOGICAL :: use_store

    IF(.NOT. PRESENT(use_store_in)) THEN
      use_store = .FALSE.
    ELSE
      use_store = use_store_in
    END IF

    NULLIFY(partlist%head)
    NULLIFY(partlist%tail)
    partlist%count = 0
    partlist%id_update = 0
    partlist%safe = .TRUE.
    partlist%use_store = use_store
    partlist%locked_store = .FALSE.
    IF (use_store) &
        CALL create_particle_store(partlist, sublist_size, .FALSE., .TRUE.)

    IF (PRESENT(holds_copies)) THEN
      partlist%holds_copies = holds_copies
    ELSE
      partlist%holds_copies = .FALSE.
    END IF

  END SUBROUTINE create_empty_partlist



  SUBROUTINE create_particle_store(partlist, n_els_min, &
      link_el_in, no_pad_store, live_state)


    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER(i8), INTENT(IN) :: n_els_min
    INTEGER(i8) :: actual_elements
    INTEGER(i8) :: i_sub, n_subs, last, link_to
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: els_to_allocate
    LOGICAL, INTENT(IN), OPTIONAL :: link_el_in, no_pad_store
    INTEGER(i4), INTENT(IN), OPTIONAL :: live_state
    LOGICAL :: link_el

    !If sublists, need to deallocate any
    CALL destroy_store(partlist%store)

    ! Make sure store is at least one sublist long
    ! and at least one extra element exists (for next_slot)
    actual_elements = MAX(n_els_min + 1, sublist_size)

    IF (PRESENT(no_pad_store)) THEN
      IF(no_pad_store) THEN
        actual_elements = n_els_min + 1
      END IF
    END IF

    link_el = .TRUE.
    IF (PRESENT(link_el_in)) THEN
      link_el = link_el_in
    END IF

    ! This is currently splitting up the memory into
    ! uniform chunks of size sublist_size
    ! This is optimal for memory usage.
    ! Could alternately allocate a large chunk, say half total
    ! and then add pieces

    n_subs = CEILING(REAL(actual_elements, num)/REAL(sublist_size, num))
    ALLOCATE(els_to_allocate(n_subs))
    els_to_allocate = sublist_size

    ! DO NOT set any element of  els_to_allocate to zero, or one!
    DO i_sub = 1, n_subs
      link_to = els_to_allocate(i_sub)
      IF(.NOT. link_el) link_to = 0
      IF(i_sub == n_subs) THEN
        last = n_els_min - SUM(els_to_allocate(1:n_subs-1))
        IF(link_to > 0) link_to = last
      END IF
      CALL create_linked_substore(partlist%store, els_to_allocate(i_sub),&
          link_to, live_state=live_state)
    END DO

    ! Next slot is the first unlinked element. This is either the
    ! first element of the whole thing (if unlinked), or the first element of
    ! the last chunk of the store if linked
    IF (link_el) THEN
      partlist%store%next_slot => partlist%store%tail%store(&
          partlist%store%tail%first_free_element)
    ELSE
      partlist%store%next_slot => partlist%store%head%store(1)
    END IF
    ! First free element should link to last filled element
    ! If we linked the store, then the last element is
    ! At slot number actual_elements, which is in either
    ! last substore, or at the very end of the one before that
    ! If we didn't, no elements are filled
    IF(link_el) THEN
      partlist%head => partlist%store%head%store(1)
      IF(partlist%store%tail%first_free_element > 1) THEN
        partlist%tail &
            => partlist%store%tail%store(&
              partlist%store%tail%first_free_element-1)
      ELSE IF(ASSOCIATED(partlist%store%tail%prev)) THEN
        partlist%tail &
            => partlist%store%tail%prev%store(&
              partlist%store%tail%prev%first_free_element-1)
      ELSE
        !This should never happen
        NULLIFY(partlist%tail)
      END IF
    ELSE
      NULLIFY(partlist%head, partlist%tail)
    END IF


  END SUBROUTINE create_particle_store


  ! Create and a new store segment object in given store
  ! Takes total_size of new chunk, and how many 'live'
  ! linked particles it should contain
  ! Optionally, if using to add a new linked chunk to an existing
  ! store, supply link_to_particle - the particle the new chunk
  ! should link onto (usually list%tail)

  SUBROUTINE create_linked_substore(store, total_size, link_upto_in, &
      link_to_particle, live_state)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: substore
    INTEGER(i8), INTENT(IN) :: total_size, link_upto_in
    TYPE(particle), POINTER, INTENT(IN), OPTIONAL :: link_to_particle
    INTEGER(i4), INTENT(IN), OPTIONAL :: live_state
    INTEGER(i8) :: link_upto
    INTEGER(i4) :: live_state_set
    TYPE(particle), POINTER :: current, prev
    INTEGER(i8) :: i_part

    !Nothing to do
    IF(total_size <= 0) RETURN

    link_upto = link_upto_in
    IF(link_upto > total_size) link_upto = total_size

    ! Default to pseudo-live - allocated but not populated. Can override
    live_state_set = -1
    IF(PRESENT(live_state)) THEN
      live_state_set = live_state
    END IF

    !Allocate substore object
    ALLOCATE(substore)
    !Set length, nullify links
    substore%length = total_size
    NULLIFY(substore%prev, substore%next)
    !Allocate backing memory
    ALLOCATE(substore%store(total_size))
    !If creating list, we can link it all up already
    !And setup the list
    !Then calling code just sets positions etc
    IF(link_upto > 0) THEN
      substore%first_free_element = link_upto + 1
    ELSE
      substore%first_free_element = 1
    END IF

    DO i_part = 1, total_size
      IF (i_part > 1 .AND. i_part < link_upto) THEN
        !Each particle slot up to request length should be linked up
        substore%store(i_part)%prev => substore%store(i_part-1)
        substore%store(i_part)%next => substore%store(i_part+1)
        current => substore%store(i_part)
        CALL init_particle(current)
      ELSE IF ((i_part == 1 .AND. i_part < link_upto) &
          .OR. i_part == link_upto) THEN
        current => substore%store(i_part)
        CALL init_particle(current)
        NULLIFY(substore%store(i_part)%prev, &
            substore%store(i_part)%next)
      ELSE
        !Nullify pointers
        NULLIFY(substore%store(i_part)%prev, &
            substore%store(i_part)%next)
      END IF
      !Set not-live state
      substore%store(i_part)%live = live_state_set
    END DO

    ! Correct links for 0th and n_elements-th particles
    IF(link_upto > 1) THEN
      substore%store(1)%next => substore%store(2)
      ! Link_upto was capped to at most total_size
      substore%store(link_upto)%prev => substore%store(link_upto-1)
    END IF

    IF(ASSOCIATED(store%tail)) THEN
      !Now link substore into store
      store%tail%next => substore
      substore%prev => store%tail

      !And create link in partlist between prior last particle
      !and first particle of new sub
      IF(link_upto > 1) THEN
        IF(PRESENT(link_to_particle)) THEN
          prev => link_to_particle
        ELSE IF(store%tail%first_free_element > 1) THEN
          ! Assume previous store is contiguously filled
          prev => store%tail%store(store%tail%first_free_element-1)
        END IF
        current => substore%store(1)
        prev%next => current
        current%prev => prev
      END IF
      store%tail => substore
      NULLIFY(substore%next)
    ELSE
      store%tail => substore
    END IF
    IF(.NOT. ASSOCIATED(store%head)) store%head => substore
    store%n_subs = store%n_subs + 1
    store%total_length = store%total_length + total_size

  END SUBROUTINE create_linked_substore


  SUBROUTINE create_empty_substore(store, total_size)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: substore
    INTEGER(i8), INTENT(IN) :: total_size
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
    END DO

    IF(ASSOCIATED(store%tail)) THEN
      !Now link substore into store
      store%tail%next => substore
      substore%prev => store%tail
      store%tail => substore
      NULLIFY(substore%next)
    ELSE
      store%tail => substore
    END IF
    IF(.NOT. ASSOCIATED(store%head)) store%head => substore
    store%n_subs = store%n_subs + 1
    store%total_length = store%total_length + total_size

  END SUBROUTINE create_empty_substore



  SUBROUTINE destroy_store(store)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: section, prev_section

    section => store%head
    IF (.NOT. ASSOCIATED(section)) RETURN
    DO WHILE(ASSOCIATED(section))
      DEALLOCATE(section%store)
      prev_section => section
      section => section%next
      DEALLOCATE(prev_section)
    END DO
    NULLIFY(store%head, store%tail)
    store%n_subs = 0
    store%total_length = 0

  END SUBROUTINE destroy_store



  FUNCTION test_store(list, skip_tail)

    TYPE(particle_list), INTENT(IN) :: list
    TYPE(particle), POINTER :: current, prev
    TYPE(particle_sub_store), POINTER :: sub
    INTEGER(i8) :: counta, countb, countc, i, a_count, countd, j, b_pos
    REAL(num) :: part_x
    LOGICAL, INTENT(IN), OPTIONAL :: skip_tail ! Ignore position of tail particle - use if calling during compact
    LOGICAL :: stop_check
    INTEGER(i8) :: test_store

    test_store = 0
    IF(.NOT. list%use_store) THEN
      !Is not a store, probably an error
      test_store = 1
      RETURN
    END IF
    !First check general integrity
    counta = 0
    countb = 0
    i = 0
    sub => list%store%head
    DO WHILE(ASSOCIATED(sub))
      i = i + 1
      IF(ASSOCIATED(sub%store)) THEN
        counta = counta + 1
      ELSE
        WRITE(100+rank, *) 'Bad substore', i
      END IF
      IF(ASSOCIATED(list%store%tail, TARGET = sub)) &
        countb = i
      sub => sub%next
    END DO
    WRITE(100+rank, *) 'Number of sublists ', &
        list%store%n_subs, counta
    FLUSH(100+rank)
    IF (list%use_store) THEN
      IF(list%store%n_subs .GT. 0 .AND. &
          counta /= list%store%n_subs) &
          test_store = test_store + 2
    END IF

    sub => list%store%head

    j = 1
    DO WHILE(ASSOCIATED(sub))
      DO i = 1, sub%length
        IF(ASSOCIATED(list%tail, &
            TARGET=sub%store(i))) &
            WRITE(100+rank, *) 'Tail is at index ', i, 'in', j
      END DO
      sub => sub%next
      j = j + 1
    END DO
    sub => list%store%head
    j = 1
    DO WHILE(ASSOCIATED(sub))
      DO i = 1, sub%length
        IF(ASSOCIATED(list%store%next_slot, &
            TARGET=sub%store(i))) &
            WRITE(100+rank, *) 'Next slot is at index ', i, 'in', j
      END DO
      sub => sub%next
      j = j + 1
    END DO

    counta = 0
    current => list%head
    DO WHILE (ASSOCIATED(current))
      counta = counta + 1
      current => current%next
    END DO

    WRITE(100+rank, *)  "Checking partlists"
    WRITE(100+rank, *) counta, list%count
    FLUSH(100+rank)

    countb = 0
    current => list%head
    DO i=1, list%count
      IF (ASSOCIATED(current)) THEN
        countb = countb + 1
        current => current%next
      ELSE
        CONTINUE
      END IF
      IF (ASSOCIATED(current, list%tail)) CONTINUE
    END DO
    WRITE(100+rank, *) countb, ASSOCIATED(current)
    FLUSH(100+rank)

    countc = 0
    a_count = 0
    sub => list%store%head
    DO WHILE(ASSOCIATED(sub))
      DO i = 1, sub%length
        IF (sub%store(i)%live > 0) THEN
        countc = countc + 1
        IF (ASSOCIATED(sub%store(i)%next)) &
            a_count = a_count + 1
       END IF
      END DO
      sub => sub%next
    END DO
    WRITE(100+rank, *)  countc, list%count, a_count+1
    FLUSH(100+rank)

    current => list%head
    NULLIFY(prev)
    i = 1
    WRITE(100+rank, *) "Checking prevs"
    DO WHILE (ASSOCIATED(current))
      IF(ASSOCIATED(prev) .AND. &
          .NOT. ASSOCIATED(current%prev, TARGET=prev)) THEN
        WRITE(100+rank, *) "Bad prev in walk at ", i
        FLUSH(100+rank)
        test_store = test_store + 4
      END IF
      prev => current
      current => current%next
      i = i + 1
    END DO


    WRITE(100+rank, *) 'Checking all positions'
    FLUSH(100+rank)
    current => list%head

    countd = 0
    b_pos = 1
    IF(PRESENT(skip_tail)) THEN
      stop_check = skip_tail
    ELSE
      stop_check = .FALSE.
    END IF
    DO WHILE (ASSOCIATED(current))
      IF (stop_check .AND. ASSOCIATED(current, TARGET=list%tail)) EXIT
      part_x  = current%part_pos
      IF( part_x .GT. x_max_local  .OR. part_x .LT. x_min_local) THEN
        WRITE(100+rank, *) 'Error, particle out of range, x', part_x, b_pos
        countd = countd + 1
      END IF
      current => current%next
      b_pos = b_pos + 1
    END DO
    WRITE(100+rank, *) "Positions Done"
    WRITE(100+rank, *) "Checking Counts"
    FLUSH(100+rank)
    IF (list%use_store) THEN
      IF(counta /= countb .OR. countb /= countc) THEN
        test_store = test_store + 8
      END IF
      IF(counta /= list%count) &
          test_store = test_store + 16
      IF(countc > 0 .AND. countb /= a_count+1) THEN
        test_store = test_store + 32
      END IF
    ELSE
      IF(counta /= countb .OR. countb /= list%count) &
          test_store = test_store + 64
    END IF

    IF(countd /= 0) &
        test_store = test_store + 128

    FLUSH(100+rank)

  END FUNCTION test_store



  SUBROUTINE create_unsafe_partlist(partlist, a_particle, n_elements, &
      holds_copies)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: a_particle
    INTEGER(i8), INTENT(IN) :: n_elements
    LOGICAL, INTENT(IN), OPTIONAL :: holds_copies
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart

    CALL create_empty_partlist(partlist, holds_copies=holds_copies)

    partlist%safe = .FALSE.
    current => a_particle
    ipart = 1
    DO WHILE (ASSOCIATED(current) .AND. ipart < n_elements)
      ipart = ipart+1
      current => current%next
    END DO
    partlist%head => a_particle
    partlist%tail => current
    partlist%count = ipart

  END SUBROUTINE create_unsafe_partlist



  SUBROUTINE create_unsafe_partlist_by_tail(partlist, head, tail, holds_copies)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: head, tail
    LOGICAL, INTENT(IN), OPTIONAL :: holds_copies
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart

    CALL create_empty_partlist(partlist, holds_copies=holds_copies)

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
      END IF
    END DO

    partlist%count = ipart

  END SUBROUTINE create_unsafe_partlist_by_tail



  SUBROUTINE create_allocated_partlist(partlist, n_elements, use_store_in, &
      holds_copies, make_live)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER(i8), INTENT(IN) :: n_elements
    LOGICAL, INTENT(IN), OPTIONAL :: use_store_in, holds_copies, make_live
    LOGICAL :: use_store
    TYPE(particle), POINTER :: new_particle
    INTEGER(i4) :: live_state
    INTEGER(i8) :: ipart

    IF(.NOT. PRESENT(use_store_in)) THEN
      use_store = .FALSE.
    ELSE
      use_store = use_store_in
    END IF

    ! Default to pseudo-live - allocated but not populated. Can override
    live_state = -1
    IF(PRESENT(make_live)) THEN
      IF (make_live) live_state = 1
    END IF

    IF (use_store) THEN
      CALL create_particle_store(partlist, n_elements, live_state=live_state)
      partlist%count = n_elements
    ELSE
      CALL create_empty_partlist(partlist, holds_copies=holds_copies)

      DO ipart = 0, n_elements-1
        CALL create_particle(new_particle)
        new_particle%live = live_state
        CALL add_particle_to_partlist(partlist, new_particle)
        NULLIFY(new_particle)
      END DO
    END IF

    partlist%use_store = use_store

  END SUBROUTINE create_allocated_partlist



  !Walk the particle store and regenerate a linked list
  !from the slots which hold live particles
  !to make list a valid linked list again

  !If particles were 'removed' without removing from store
  ! recount will restore the correct count!
  SUBROUTINE relink_partlist(list, recount)

    TYPE(particle_list), INTENT(INOUT) :: list
    TYPE(particle), POINTER :: current, previous
    TYPE(particle_sub_store), POINTER :: sub
    INTEGER(i8) :: cnt, icurr
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
          END IF
          current%prev => previous
          previous => current
        END IF
      END DO
      sub => sub%next
    END DO

    IF (ASSOCIATED(previous)) THEN
      NULLIFY(previous%next)
      list%tail => previous
    ELSE
      ! Can only happen if there is nothing in list
      NULLIFY(list%tail)
      IF (store_debug) PRINT*, "ERROR, empty list ", cnt
    END IF

    IF(recount) list%count = cnt

  END SUBROUTINE relink_partlist


  !Increment the position of the next slot in store
  ! This is meant to be done AFTER you have filled the
  ! next_slot and want it to be updated to a new value
  !If this overflows, then make space
  !First try compacting the list
  !If that doesn't help, allocate a new sub-store chunk
  !THIS ROUTINE MAY INVALIDATE POINTERS!!
  !ANY CALLER MUST CHECK LIST IS UNCHANGED OR REPOINT!!
  ! Next_slot (==list%tail) when this routine is entered
  ! will remain the list%tail particle
  SUBROUTINE increment_next_free_element(list)

    TYPE(particle_list), INTENT(INOUT) :: list

    IF(list%store%tail%first_free_element >= list%store%tail%length) THEN
      ! Any path here will leave first_free_element incremented
      !First resort: compact store
      IF(list%count > 0 .AND..NOT. list%locked_store .AND. &
          REAL(list%count)/REAL(list%store%total_length) < fill_factor) THEN
        ! Remove empty space. Mostly this does very little
        ! but it doesn't cost much and helps if chunks are small
        CALL update_first_frees(list)
        CALL remove_empty_subs(list)
        ! If that wasn't enough, compact completely
        IF (REAL(list%count)/REAL(list%store%total_length) < fill_factor) THEN
          IF (fold_compact) THEN
            CALL fold_compact_backing_store(list%store, list)
          ELSE
            CALL compact_backing_store(list%store, list)
          END IF
        END IF
      END IF
      IF(list%store%tail%first_free_element >= list%store%tail%length) THEN
        ! Compacting not possible or insufficient - have to add subs
        CALL create_empty_substore(list%store, sublist_size)
        list%store%next_slot => &
            list%store%tail%store(list%store%tail%first_free_element)
      END IF
    ELSE
      !Do this only if we've not created anything new
      list%store%tail%first_free_element = &
          list%store%tail%first_free_element + 1
      list%store%next_slot => &
          list%store%tail%store(list%store%tail%first_free_element)

   END IF

  END SUBROUTINE increment_next_free_element



  SUBROUTINE create_filled_partlist(partlist, data_in, n_elements, holds_copies)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    REAL(num), DIMENSION(:), INTENT(IN) :: data_in
    INTEGER(i8), INTENT(IN) :: n_elements
    LOGICAL, INTENT(IN), OPTIONAL :: holds_copies
    TYPE(particle), POINTER :: new_particle
    INTEGER(i8) :: ipart, cpos = 0

    CALL set_partlist_size
    CALL create_empty_partlist(partlist, holds_copies=holds_copies)

    DO ipart = 0, n_elements-1
      ALLOCATE(new_particle)
      cpos = ipart*nvar+1
      CALL unpack_particle(data_in(cpos:cpos+nvar-1), new_particle)
#ifdef PARTICLE_DEBUG
      new_particle%processor = rank
#endif
      CALL add_particle_to_partlist(partlist, new_particle)
      NULLIFY(new_particle)
    END DO

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
    END IF

    ! List with head or tail but not both is broken
    IF (.NOT. ASSOCIATED(partlist%head) &
        .OR. .NOT. ASSOCIATED(partlist%tail)) THEN
      test_partlist = -1
      RETURN
    END IF

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
      END IF
    END DO

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
        ! A partlist that holds copies or an unsafe partlist should not cause
        ! unlinking
        CALL destroy_particle(new_particle, &
            partlist%holds_copies .OR. .NOT.partlist%safe)
        new_particle => next
        ipart = ipart+1
      END DO

      CALL create_empty_partlist(partlist)
    END IF


  END SUBROUTINE destroy_partlist



  SUBROUTINE copy_partlist(partlist1, partlist2)

    TYPE(particle_list), INTENT(INOUT) :: partlist1, partlist2

    partlist2%head => partlist1%head
    partlist2%tail => partlist1%tail
    partlist2%count = partlist1%count
    partlist2%id_update = partlist1%id_update
    partlist2%holds_copies = partlist1%holds_copies

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
    END IF

    IF(PRESENT(ignore_live_in)) THEN
      ignore_live = ignore_live_in
    ELSE
      ignore_live = .FALSE.
    END IF

    !Do the appending
    IF (.NOT. (list%use_store .OR. newlist%use_store)) THEN
      IF (ASSOCIATED(list%tail)) THEN
        list%tail%next => newlist%head
      ELSE
        list%head => newlist%head
      END IF
      IF (ASSOCIATED(newlist%head)) newlist%head%prev => list%tail
      IF (ASSOCIATED(newlist%tail)) list%tail => newlist%tail
      list%count = list%count + newlist%count
      list%id_update = list%id_update + newlist%id_update
    ELSE IF (.NOT. list%use_store .AND. newlist%use_store) THEN
      !This is an error and should never arise.
      IF (rank == 0) &
          PRINT *, 'Unable to append partlists'
      RETURN
    ELSE
      CALL add_partlist_to_list_and_store(list, newlist, ignore_live)
    END IF
    !Clean up newlist
    CALL create_empty_partlist(newlist)


  END SUBROUTINE append_partlist


  !Take a list and append its content to list-with-store
  SUBROUTINE add_partlist_to_list_and_store(list, newlist, override_live)

    TYPE(particle_list), INTENT(INOUT) :: newlist, list
    LOGICAL, INTENT(IN) :: override_live !Override any live states in newlist
    TYPE(particle), POINTER :: current, next

    IF (newlist%count < 1) RETURN !Nothing to append, make no change

    current => newlist%head
    NULLIFY(next)

    DO WHILE(ASSOCIATED(current))
      !Only consider live particles, unless overrriding
      IF (override_live .OR. current%live == 1) THEN
        CALL create_particle_in_list(next, list, .TRUE.)

        CALL copy_particle(current, next)
        next%live = 1 ! Required if over-riding, does nothing else
      END IF

      current => current%next
    END DO
    list%id_update = list%id_update + newlist%id_update

  END SUBROUTINE add_partlist_to_list_and_store




  SUBROUTINE add_particle_to_partlist(partlist, new_particle)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: new_particle

    ! Note that this will work even if you are using an unsafe particle list
    ! BE CAREFUL if doing so, it can cause unexpected behaviour
    ! DO NOT USE with store-backed lists

    ! if (!particle) return;
    IF (.NOT. ASSOCIATED(new_particle)) RETURN
    NULLIFY(new_particle%next, new_particle%prev)

    ! Add particle count
    partlist%count = partlist%count + 1
    partlist%id_update = 1
    IF (.NOT. ASSOCIATED(partlist%tail)) THEN
      ! partlist is empty
      partlist%head => new_particle
      partlist%tail => new_particle
      RETURN
    END IF

    partlist%tail%next => new_particle
    new_particle%prev => partlist%tail
    NULLIFY(new_particle%next)
    partlist%tail => new_particle

  END SUBROUTINE add_particle_to_partlist


  SUBROUTINE remove_particle_from_partlist(partlist, a_particle, &
      destroy)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: a_particle, tmp_particle
    LOGICAL, INTENT(IN), OPTIONAL :: destroy
    LOGICAL :: destroy_in
    ! Remove a particle from a partlist completely. If the list is store-
    ! backed, then a copy will be returned in place of a_particle - unless
    ! destroy_particle is true, when the particle is completely deleted and
    ! NULL pointer returned

    ! Note that this will work even if you are using an unsafe particle list
    ! BE CAREFUL if doing so, it can cause unexpected behaviour
    IF( .NOT. ASSOCIATED(a_particle)) RETURN
    IF( a_particle%live /= 1) RETURN

    IF (PRESENT(destroy)) THEN
      destroy_in=destroy
    ELSE
      destroy_in = .FALSE.
    END IF

    CALL unlink_particle_from_partlist(partlist, a_particle)

    IF (partlist%use_store) THEN
      ! Setting live effectively removes particle from store -
      ! it will be overwritten on next compact etc
      a_particle%live = 0
    END IF

    IF (partlist%use_store .AND. .NOT. destroy_in) THEN
      !If a_particle is in a store, make a copy
      !Then what comes back is a valid, FREE particle
      CALL create_particle(tmp_particle, .TRUE.)
      CALL copy_particle(a_particle, tmp_particle)
      !Return a live particle
      tmp_particle%live = 1
      a_particle => tmp_particle
    ELSE IF (partlist%use_store) THEN
      !Don't need to actually destroy, just hand back NULL
      NULLIFY(a_particle)
    ELSE IF (destroy_in) THEN
      CALL destroy_particle(a_particle)
      NULLIFY(a_particle)
    END IF


  END SUBROUTINE remove_particle_from_partlist



  SUBROUTINE unlink_particle_from_partlist(partlist, a_particle)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: a_particle
    ! Unlink a particle from a partlist. It remains
    ! a valid particle in the store, but not part of the list it backs:
    ! this means a relink will ADD it back!

    ! Note that this will work even if you are using an unsafe particle list
    ! BE CAREFUL if doing so, it can cause unexpected behaviour
    IF( .NOT. ASSOCIATED(a_particle)) RETURN
    IF( a_particle%live /= 1) RETURN

    ! Check whether particle is head or tail of list and unlink
    IF (ASSOCIATED(partlist%head, TARGET=a_particle)) THEN
      partlist%head => a_particle%next
    END IF

    IF (ASSOCIATED(partlist%tail, TARGET=a_particle)) THEN
      partlist%tail => a_particle%prev
    END IF

    ! Link particles on either side together
    IF (ASSOCIATED(a_particle%next)) a_particle%next%prev => a_particle%prev
    IF (ASSOCIATED(a_particle%prev)) a_particle%prev%next => a_particle%next

    NULLIFY(a_particle%next, a_particle%prev)

    ! Decrement counter
    partlist%count = partlist%count-1


  END SUBROUTINE unlink_particle_from_partlist


  SUBROUTINE pack_particle(array, a_particle)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    TYPE(particle), POINTER :: a_particle
    INTEGER(i8) :: cpos, temp_i8

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
#ifdef PARTICLE_ID4
    temp_i8 = INT(a_particle%id, i8)
    array(cpos) = TRANSFER(temp_i8, 1.0_num)
    cpos = cpos+1
#elif PARTICLE_ID
    array(cpos) = TRANSFER(a_particle%id, 1.0_num)
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
#ifdef WORK_DONE_INTEGRATED
    array(cpos) = a_particle%work_x
    array(cpos+1) = a_particle%work_y
    array(cpos+2) = a_particle%work_z
    array(cpos+3) = a_particle%work_x_total
    array(cpos+4) = a_particle%work_y_total
    array(cpos+5) = a_particle%work_z_total
    cpos = cpos+6
#endif
    IF (any_persistent_subset) THEN
      temp_i8 = id_registry%map(a_particle)
      array(cpos) = TRANSFER(temp_i8, 1.0_num)
      cpos = cpos+1
    END IF

  END SUBROUTINE pack_particle



  SUBROUTINE unpack_particle(array, a_particle)

    REAL(num), DIMENSION(:), INTENT(IN) :: array
    TYPE(particle), POINTER :: a_particle
    INTEGER(i8) :: cpos, temp_i8

    cpos = 1
    a_particle%part_pos = array(cpos)
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
    temp_i8 = TRANSFER(array(cpos), temp_i8)
    a_particle%id = INT(temp_i8)
    cpos = cpos+1
#elif PARTICLE_ID
    a_particle%id = TRANSFER(array(cpos), a_particle%id)
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
    a_particle%live = 1 !Only live particles ever sent
#ifdef WORK_DONE_INTEGRATED
    a_particle%work_x = array(cpos)
    a_particle%work_y = array(cpos+1)
    a_particle%work_z = array(cpos+2)
    a_particle%work_x_total = array(cpos+3)
    a_particle%work_y_total = array(cpos+4)
    a_particle%work_z_total = array(cpos+5)
    cpos = cpos+6
#endif
    IF (any_persistent_subset) THEN
      CALL id_registry%add_with_map(a_particle, TRANSFER(array(cpos), temp_i8))
      cpos = cpos+1
    END IF

  END SUBROUTINE unpack_particle



#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
  FUNCTION generate_id()

    INTEGER(idkind) :: generate_id

    highest_id = highest_id + 1_idkind
    generate_id = cpu_id + highest_id

  END FUNCTION generate_id
#endif



  SUBROUTINE init_particle(new_particle, no_gen_id)

    TYPE(particle), POINTER :: new_particle
    LOGICAL, INTENT(IN), OPTIONAL :: no_gen_id

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
    IF (PRESENT(no_gen_id)) THEN
      IF (.NOT. no_gen_id) THEN
        new_particle%id = generate_id()
      END IF
    ELSE
      new_particle%id = generate_id()
    END IF
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




  SUBROUTINE create_particle(new_particle, no_gen_id)

    TYPE(particle), POINTER :: new_particle
    LOGICAL, INTENT(IN), OPTIONAL :: no_gen_id

    ALLOCATE(new_particle)
    CALL init_particle(new_particle, no_gen_id)

  END SUBROUTINE create_particle



  SUBROUTINE create_particle_in_list(new_particle, list, no_gen_id)

    TYPE(particle), POINTER, INTENT(INOUT) :: new_particle
    TYPE(particle_list), INTENT(INOUT) :: list
    LOGICAL, INTENT(IN), OPTIONAL :: no_gen_id

    IF (list%use_store) THEN
      new_particle => list%store%next_slot
      CALL init_particle(new_particle, no_gen_id)
      new_particle%prev => list%tail
      NULLIFY(new_particle%next)
      list%tail => new_particle
      IF(ASSOCIATED(new_particle%prev)) THEN
        new_particle%prev%next => new_particle
      ELSE
        IF (list%count == 0) THEN
          list%head => new_particle
        ELSE IF (.NOT. ASSOCIATED(list%head, TARGET=new_particle)) THEN
          IF (store_debug) PRINT*, &
               rank, 'Error creating in list - bad prev or head'
        END IF
      END IF
      list%count = list%count + 1

      CALL increment_next_free_element(list)
      ! If compact occured on line above, our slot can have moved, so reset
      new_particle => list%tail
    ELSE
      CALL create_particle(new_particle, no_gen_id)
      new_particle%prev => list%tail
      NULLIFY(new_particle%next)
      list%tail => new_particle
      IF(ASSOCIATED(new_particle%prev)) new_particle%prev%next => new_particle
      list%count = list%count + 1

    END IF
    IF (.NOT. ASSOCIATED(list%head)) list%head => new_particle

  END SUBROUTINE create_particle_in_list



  SUBROUTINE destroy_particle(part, is_copy)

    ! Routine to delete a particle. This routine is only safe to use on
    ! a particle that is not in a partlist
    TYPE(particle), POINTER :: part
    LOGICAL, INTENT(IN), OPTIONAL :: is_copy

    IF (any_persistent_subset) THEN
      IF (PRESENT(is_copy)) THEN
        IF (.NOT. is_copy) CALL id_registry%delete_all(part)
      ELSE
        CALL id_registry%delete_all(part)
      END IF
    END IF

    DEALLOCATE(part)

  END SUBROUTINE destroy_particle



  SUBROUTINE display_particle(a_particle)

    TYPE(particle), POINTER :: a_particle

    PRINT *, 'Position', a_particle%part_pos
    PRINT *, 'Momentum', a_particle%part_p

  END SUBROUTINE display_particle



  FUNCTION compare_particles(part1, part2)

    TYPE(particle), POINTER :: part1, part2
    LOGICAL :: compare_particles

    compare_particles = .TRUE.
    IF (ABS(part1%part_pos-part2%part_pos) > c_tiny) &
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
    END IF

  END FUNCTION compare_particles



  FUNCTION test_packed_particles(partlist, array, npart_in_data)

    TYPE(particle_list), INTENT(IN) :: partlist
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER(i8), INTENT(IN) :: npart_in_data
    TYPE(particle), POINTER :: current
    TYPE(particle), POINTER :: a_particle
    LOGICAL :: test_packed_particles
    INTEGER(i8) :: ipart

    CALL set_partlist_size

    test_packed_particles = .FALSE.

    IF (npart_in_data * nvar /= SIZE(array)) THEN
      PRINT *, 'Size of data array does not match specified on', rank, &
          npart_in_data, SIZE(array)
      RETURN
    END IF
    IF (partlist%count /= npart_in_data) THEN
      PRINT *, 'Size of data array does not match partlist on', rank
      RETURN
    END IF

    ALLOCATE(a_particle)

    current => partlist%head
    DO ipart = 0, npart_in_data-1
      CALL unpack_particle(array(ipart*nvar+1:(ipart+1)*nvar), a_particle)
      IF (.NOT. compare_particles(a_particle, current)) THEN
        PRINT *, 'BAD PARTICLE ', ipart, 'on', rank
        RETURN
      END IF
      current => current%next
    END DO

    DEALLOCATE(a_particle)  ! DO NOT REPLACE WITH CALL TO destroy_particle

    test_packed_particles = .TRUE.

  END FUNCTION test_packed_particles



  SUBROUTINE partlist_send_nocount(partlist, dest)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: dest
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    INTEGER :: ipart, nsend, cpos
    TYPE(particle), POINTER :: current

    CALL set_partlist_size

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
    END DO

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

    CALL set_partlist_size
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

    CALL set_partlist_size

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
      data_send(cpos:cpos+nvar-1) = packed_particle_data(1:nvar)
      ipart = ipart + 1
      current => current%next
    END DO

    ! No longer need the sending partlist, so destroy it to save some memory
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



  SUBROUTINE update_particle_count

    ! This routine ensures that the particle count for the species_list
    ! objects is accurate. This makes some things easier, but increases
    ! communication
    INTEGER :: ispecies
    LOGICAL, SAVE :: update = .TRUE.

    IF (.NOT.update) RETURN

    DO ispecies = 1, n_species
      CALL MPI_ALLREDUCE(species_list(ispecies)%attached_list%count, &
          species_list(ispecies)%count, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
      species_list(ispecies)%count_update_step = step
    END DO

    update = use_particle_count_update

  END SUBROUTINE update_particle_count



  !The following goes through the backing store as an array, packing
  ! particles into a contiguous chunk. Any empty substores are skipped
  ! over and then removed. In some cases this can save some copying

  SUBROUTINE compact_backing_store(store, list)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: src_section, dest_section
    TYPE(particle_list), INTENT(INOUT) :: list
    TYPE(particle), POINTER :: original, place_into, prev_placed
    INTEGER(i8) :: i, next_place_index, j

    IF (store_debug) THEN
      PRINT*, 'Compacting backing store on ',  rank
    END IF

    src_section => store%head
    dest_section => store%head
    next_place_index = 1

    place_into => store%head%store(1)
    NULLIFY(prev_placed)
    outer : DO j = 1, store%n_subs
      DO i = 1, src_section%first_free_element
        IF (i > src_section%length) EXIT
        ! Check if current particle is live and not already in right place
        place_into%prev => prev_placed
        IF (src_section%store(i)%live > 0 .AND. &
            .NOT. ASSOCIATED(place_into, TARGET=src_section%store(i))) THEN
          original => src_section%store(i)
          CALL copy_particle(original, place_into)
          original%live = 0
          next_place_index = next_place_index + 1
          prev_placed => place_into
        ELSE IF (src_section%store(i)%live > 0) THEN
          !Otherwise we're leaving it where it is
          next_place_index = i + 1
          prev_placed => place_into
        END IF
        IF (next_place_index > dest_section%length) THEN
          !Move dest_section to next, before advancing place_into
          !Update dest_section count
          dest_section%first_free_element = dest_section%length + 1

          !Increment dest_section to next chunk
          !There must be one if we still have live particles
          dest_section => dest_section%next
          next_place_index = 1
        END IF
        IF (ASSOCIATED(dest_section)) THEN
          place_into => dest_section%store(next_place_index)
        ELSE
          EXIT outer
        END IF
        IF(ASSOCIATED(prev_placed)) prev_placed%next => place_into
      END DO
      src_section => src_section%next
    END DO outer

    !Set first_free for any remaining destination sections
    IF(ASSOCIATED(dest_section)) THEN
      dest_section%first_free_element = next_place_index
      DO WHILE(ASSOCIATED(dest_section))
        dest_section => dest_section%next
        IF(ASSOCIATED(dest_section)) dest_section%first_free_element = 1
      END DO
    END IF

    CALL remove_empty_subs(list)
    CALL set_next_slot(list)

    list%head => store%head%store(1)
    list%tail => prev_placed
    NULLIFY(list%head%prev, list%tail%next)


  END SUBROUTINE compact_backing_store



  !This does the same as above, but in a slightly simpler way
  ! However, this one can't be adapted to skip empty sections
  SUBROUTINE alt_compact_backing_store(store, list)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: dest_section
    TYPE(particle_list), INTENT(INOUT) :: list
    TYPE(particle), POINTER :: original, place_into
    INTEGER(i8) :: i, j

    IF (store_debug) THEN
      PRINT*, 'Compacting backing store on ',  rank
    END IF

    dest_section => store%head
    original => list%head
    place_into => store%head%store(1)
    outer : DO j = 1, store%n_subs
      DO i = 1, dest_section%length
        IF (i > dest_section%length) EXIT
        ! Check if current particle is live and not already in right place
        IF (dest_section%store(i)%live > 0) THEN
          original => original%next
        ELSE
          place_into => dest_section%store(i)
          CALL copy_particle(original, place_into)
          original%live = 0
          original => original%next
        END IF
        IF (.NOT. ASSOCIATED(original)) THEN
          IF (dest_section%store(i)%live == 1) THEN
            dest_section%first_free_element = i+1
          ELSE
            dest_section%first_free_element = i
          END IF
          EXIT outer
        END IF
      END DO
      dest_section%first_free_element = dest_section%length + 1
      dest_section => dest_section%next
    END DO outer

    !Set first_free for any remaining destination sections
    IF(ASSOCIATED(dest_section)) THEN
      DO WHILE(ASSOCIATED(dest_section))
        dest_section => dest_section%next
        IF(ASSOCIATED(dest_section)) dest_section%first_free_element = 1
      END DO
    END IF

    CALL remove_empty_subs(list)
    CALL set_next_slot(list)

    CALL relink_partlist(list, .FALSE.)

  END SUBROUTINE alt_compact_backing_store



  !The following goes through the backing store as an array, packing
  ! particles into a contiguous chunk. This is done by wrapping
  ! particles from the tail of the list into empty spaces.
  ! THIS DOES NOT preserve list ordering!
  ! However, for various reasons it DOES preserve the tail particle

  SUBROUTINE fold_compact_backing_store(store, list)

    TYPE(particle_store), INTENT(INOUT) :: store
    TYPE(particle_sub_store), POINTER :: dest_section
    TYPE(particle_list), INTENT(INOUT) :: list
    TYPE(particle), POINTER :: original, place_into
    INTEGER(i8) :: i, d_count, offset
    LOGICAL :: patched_tail

    IF (store_debug) THEN
      PRINT*, 'Folding backing store on ',  rank
    END IF

   ! THIS should only be called if compact is needed
   ! Assume there is some space unused
    dest_section => store%head

    ! List can only have no tail if empty or corrupt
    IF (.NOT. ASSOCIATED(list%tail) .AND. store_debug) THEN
      PRINT*, "Error - list has no tail, cannot fold"
    END IF

    ! Below means only one element left in list, the head/tail
    ! Folding makes no sense in that case
    ! But that element may be in unhelpful position
    ! So we do want to move it back to the first slot
    IF (.NOT. ASSOCIATED(list%tail%prev)) THEN
      place_into => store%head%store(1)
      CALL copy_particle(list%tail, place_into)
      list%tail%live = 0
      NULLIFY(place_into%prev, place_into%next)
      list%head => place_into
      list%tail => list%head
      store%head%first_free_element = 2
      store%next_slot => store%head%store(2)
      RETURN
    END IF

    patched_tail = .FALSE.
    original => list%tail%prev
    offset = 0
    d_count = 0
    outer: DO WHILE (offset <= list%count)
      DO i = 1, dest_section%length
        ! At worst we copy every particle
        ! offset tracks the total of either live, or filled by copying
        ! We're done when this is 1 below the list count as the
        ! tail is handled separately
        offset = offset + 1
        IF (offset >= list%count) THEN
          dest_section%first_free_element = i+1
          place_into => dest_section%store(i)
          CALL copy_particle(list%tail, place_into)
          list%tail%live = 0
          patched_tail = .TRUE.
          EXIT outer
        END IF
        ! Find next empty slot
        IF (dest_section%store(i)%live == 1) CYCLE
        ! Got a slot, copy it and move on src
        place_into => dest_section%store(i)
        CALL copy_particle(original, place_into)
        original%live = 0
        original => original%prev
      END DO
      ! For all except the last filled segment, below is correct
      dest_section%first_free_element = dest_section%length+1
      dest_section => dest_section%next
    END DO outer

    IF (.NOT. patched_tail) THEN
      ! Must have ended exactly on a section boundary
      ! hence exited from WHILE cond, not EXIT outer
      ! That cannot have been the last segment
      dest_section => dest_section%next
      place_into => dest_section%store(1)
      CALL copy_particle(list%tail, place_into)
      list%tail%live = 0
      dest_section%first_free_element = 2
    END IF

    !Set first_free for any remaining sections
    DO WHILE(ASSOCIATED(dest_section))
      dest_section => dest_section%next
      IF(ASSOCIATED(dest_section)) THEN
        dest_section%first_free_element = 1
      END IF
    END DO

    CALL remove_empty_subs(list)
    CALL set_next_slot(list)

    CALL relink_partlist(list, .FALSE.)

  END SUBROUTINE fold_compact_backing_store



  !Remove empty sublists, updating store accordingly
  SUBROUTINE remove_empty_subs(list)

    TYPE(particle_list), INTENT(INOUT) :: list
    TYPE(particle_sub_store), POINTER :: current, next
    INTEGER(i8) :: length, st

    current => list%store%head
    length = list%store%total_length
    st = 0
    !Do nothing if only one sub
    IF (.NOT. ASSOCIATED(current%next)) RETURN

    DO WHILE(ASSOCIATED(current) .AND. list%store%n_subs > 1)
      next => current%next
      st = st + 1
      IF(current%first_free_element == 1) THEN
        !Delete empty segments
        IF (store_debug) THEN
          !Check store is truly empty!
          IF (count_live(current) > 0) THEN
            PRINT*, 'Severe error in particle stores on rank ', rank, &
              ' segment ', st, ' is not empty'
          END IF
        END IF

        IF(ASSOCIATED(current%prev)) THEN
          current%prev%next => next
        ELSE
          !Is head of list
          list%store%head => next
        END IF

        IF(ASSOCIATED(next)) THEN
          next%prev => current%prev
        ELSE
          !Is tail of list
          list%store%tail => current%prev
        END IF

        list%store%total_length = list%store%total_length - current%length
        DEALLOCATE(current)
        list%store%n_subs = list%store%n_subs - 1
      END IF
      current => next
    END DO

  END SUBROUTINE remove_empty_subs


  ! Update indices of first free element in each sub
  SUBROUTINE update_first_frees(list)

    TYPE(particle_list), INTENT(INOUT) :: list
    TYPE(particle_sub_store), POINTER :: current
    INTEGER(i8) :: i, old, start

    current => list%store%head
    !Don't touch tail-most store
    DO WHILE(ASSOCIATED(current) .AND. ASSOCIATED(current%next))
      old = current%first_free_element
      IF (current%first_free_element > 1) THEN
        current%first_free_element = 1
        ! Check if we can move it down: start at old value, hunt
        ! for last actual live particle
        start = MIN(old, current%length)
        DO i = start, 1, -1
          IF (current%store(i)%live == 1) THEN
            current%first_free_element = i+1
            EXIT
          END IF
        END DO
      END IF
      current => current%next
    END DO

  END SUBROUTINE update_first_frees



  FUNCTION count_live(sub_store)

    INTEGER(KIND=i8) :: count_live, i
    TYPE(particle_sub_store), POINTER :: sub_store

    count_live = 0
    DO i = 1, sub_store%length
      IF (sub_store%store(i)%live == 1) count_live = count_live + 1
    END DO

  END FUNCTION count_live



  SUBROUTINE set_next_slot(list)

    TYPE(particle_list), INTENT(INOUT) :: list

    IF(list%store%tail%first_free_element <= list%store%tail%length) THEN
      list%store%next_slot => &
          list%store%tail%store(list%store%tail%first_free_element)
    ELSE
      CALL create_empty_substore(list%store, sublist_size)
      list%store%next_slot => &
          list%store%tail%store(list%store%tail%first_free_element)
    END IF

  END SUBROUTINE set_next_slot



END MODULE partlist
