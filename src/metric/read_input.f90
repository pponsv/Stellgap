module read_input

   use constants
   use read_boozer_mod

   implicit none

contains

   subroutine read_args_and_boozer
      integer :: ierr
      character(len=128) :: boozer_extension

      ! Check the number of command line arguments
      if (iargc() .ne. 1) then
         print *, ' MUST ENTER FILE SUFFIX ON COMMAND LINE'
         stop
      end if

      ! Read the file suffix from the command line
      call getarg(1, boozer_extension)

      ! Read the boozer file
      call read_boozer_file(boozer_extension, ierr)

      ! Check for errors
      if (ierr .ne. 0) then
         print *, ' ERROR READING BOOZER FILE'
         stop
      end if

   end subroutine read_args_and_boozer

end module read_input
