************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
#ifdef _DMRG_
!***********************************************************************
!                                                                      *
!  Purpose:  read (and check) DMRG input in QCMaquis style             *
!                                                                      *
!***********************************************************************
      subroutine qcmaquis_rdinp(luinput,switch,nr_lines)

      use qcmaquis_interface_cfg
      use qcmaquis_interface_utility_routines, only:
     &    lower_to_upper, find_qcmaquis_keyword

      implicit none

      integer, intent(in)    :: luinput
      integer, intent(in)    :: switch
      integer, intent(inout) :: nr_lines

      character(len=500)     :: line, line2
      integer                :: i, io, j

      ! Leon: flag for the compression keyword parsing
      logical                :: compr_flag
      ! Leon: the evil common block that contains the MPSCompressM
      ! option for the MPS compression
#include "nevptp.fh"

      compr_flag = .false.

      select case(switch)
        !> find maximum number of input lines
        case(1)
          nr_lines = 0
          do
            line(1:500) = ' '
            read(luinput,*,IOSTAT=io) line
            if (is_iostat_end(io)) exit
            if (io>0) stop 'problem reading QCMaquis input'
            call lower_to_upper(line(1:5))
            if(line(1:5) == 'ENDRG' .or. line(1:5) == 'ENDDM') exit
            nr_lines = nr_lines + 1
          end do
          dmrg_input%nr_qcmaquis_input_lines = nr_lines
        !> read input lines
        case(2)
          if(nr_lines == 0) return
          allocate(dmrg_input%qcmaquis_input(nr_lines))
          dmrg_input%qcmaquis_input = ' '

          do i = 1, nr_lines
            read(luinput,*,IOSTAT=io) line
            if (is_iostat_end(io)) exit
            if (io>0) then
              write (6,*) "Problem reading QCMaquis input"
              call Quit_OnUserError()
            end if
            ! Leon: handle a special case where we input the MPS compression
            ! parameter in the RGINPUT block but it does not get passed to
            ! QCMaquis at all

            dmrg_input%qcmaquis_input(i) = trim(line)
            if (.not.compr_flag) then
              line2=adjustl(trim(line))
              call lower_to_upper(line2(1:5))
              if (line2(1:5) == 'COMPR') then
                ! process the compression keyword
                ! set the compression flag to true and read the next line.
                compr_flag = .true.
              end if
            else
              read(line,*,IOSTAT=io) MPSCompressM
              if (io /= 0) then
                write (6,*) "Problem reading QCMaquis input"
                call Quit_OnUserError()
              end if
              compr_flag = .false.
            end if
          end do

        !> sanity check input for ALL mandatory keywords
        case(3)
          do j = 1,2
            line(1:500) = ' '
            select case(j)
              case(1)
                line = 'NSWEEPS'
              case(2)
                line = 'MAX_BOND_DIMENSION'
            end select
            i = 0
            call find_qcmaquis_keyword(dmrg_input%qcmaquis_input,
     &                                 nr_lines,
     &                                 line,
     &                                 i
     &                                )
            if(i <= 0)then
              !> check for keyword sweep_bond_dimension which is an alternative
              if(trim(line) == 'MAX_BOND_DIMENSION')then
                line2(1:500) = ' '
                line2        = 'SWEEP_BOND_DIMENSIONS'
                i = 0
                call find_qcmaquis_keyword(dmrg_input%qcmaquis_input,
     &                                     nr_lines,
     &                                     line2,
     &                                     i
     &                                    )
                if(i > 0) cycle
              end if
              Call WarningMessage(2,'Error in input preprocessing.')
              write(6,*)' qcmaquis_rdinp: mandatory keyword ',
     &        trim(line),' missing in QCMaquis DMRG input section'
              nr_lines = -1; return
            end if
          end do
        case default
          write(6,*) ' QCMaquis input reader - you should have never'//
     &             ' reached this spot...'
          call Quit_OnUserError()
      end select

      end subroutine qcmaquis_rdinp
#elif defined (NAGFOR)
c Some compilers do not like empty files
      subroutine empty_qcmaquis_rdinp()
      end subroutine empty_qcmaquis_rdinp
#endif
