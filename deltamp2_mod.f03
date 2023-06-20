      module test_mod
!
!     This module supports the program scfEnergyTerms.
!
!     -H. P. Hratchian, 2020.
!
!
!     USE Connections
!
      use mqc_general
      use mqc_molecule
      use mqc_gaussian
      use mqc_algebra2
      use mqc_algebra
      use iso_fortran_env
      use OMP_LIB
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
!
!
!     Module Procedures
!
      CONTAINS
!
!
!PROCEDURE commandLineArgs
      subroutine commandLineArgs(iPrint,matrixFilename,doN8,doSlowN5,  &
        doRegularN5,fail)
!
!     This subroutine is used to process the command line arguments.
!
      integer,intent(OUT)::iPrint
      character(len=512),intent(OUT)::matrixFilename
      logical,intent(OUT)::doN8,doSlowN5,doRegularN5,fail
!
      integer::nCommands,nMatrixFilenames,i
      character(len=512)::tmpString,lowercase
!
!     Format statements.
!
 9000 format(1x,'Unknown command line switch: ',A,'.')
 9100 format(1x,'More than 1 matrix filename found!')
!
!     Set defaults.
!
      iPrint = 0
      doN8 = .false.
      doSlowN5 = .true.
      doRegularN5 = .true.
      fail = .false.
      nMatrixFilenames = 0
!
!     Determine the number of command line arguments. Then, loop through the
!     list of arguments to set options.
!
      nCommands = command_argument_count()
      do i = 1,nCommands
        call get_command_argument(i,tmpString)
        write(*,*) i,TRIM(tmpString)
        if(tmpString(1:1).eq.'-') then
          call String_Change_Case(tmpString(2:),'l',lowercase)
          select case(TRIM(lowercase))
          case('don8')
            doN8 = .true.
          case('skipn8')
            doN8 = .false.
          case('doslown5')
            doSlowN5 = .true.
          case('skipslown5')
            doSlowN5 = .false.
          case('doregularn5')
            doRegularN5 = .true.
          case('skipregularn5')
            doRegularN5 = .false.
          case('debug')
            iPrint = 1
          case default
            fail = .true.
            write(iOut,9000) TRIM(tmpString)
            return
          endSelect
        else
          nMatrixFilenames = nMatrixFilenames + 1
          if(nMatrixFilenames.gt.1) then
            fail = .true.
            write(iOut,9100)
            return
          endIf
          matrixFilename = tmpString
        endIf
      endDo
!
      return
      end subroutine commandLineArgs

!
!
!PROCEDURE integralTransformationN8
      subroutine integralTransformationN8sameSpin(nBasis,nBasisUse,C,aoInts,moInts)
!
!     This subroutine carries out AO-->MO integral transformations for a same-spin case.
!
!
!     H. P. Hratchian, 2022
!
      implicit none
      integer(kind=int64),intent(in)::nBasis,nBasisUse
      real(kind=real64),dimension(nBasis,nBasisUse)::C
      real(kind=real64),dimension(nBasis,nBasis,nBasis,nBasis)::aoInts
      real(kind=real64),dimension(nBasisUse,nBasisUse,nBasisUse,nBasisUse)::moInts
      integer(kind=int64)::p,q,r,s,mu,nu,lambda,sigma
!
!     Do the transformation with the straightforward O(N^8) algorithm.
!
      do p = 1,nBasisUse
        do q = 1,nBasisUse
          do r = 1,nBasisUse
            do s = 1,nBasisUse
              do mu = 1,nBasis
                do nu = 1,nBasis
                  do lambda = 1,nBasis
                    do sigma = 1,nBasis
                      moInts(p,q,r,s) = moInts(p,q,r,s)  &
                        + C(mu,p)*C(nu,q)*C(lambda,r)*C(sigma,s)*aoInts(mu,nu,lambda,sigma)
                    endDo
                  endDo
                endDo
              endDo
            endDo
          endDo
        endDo
      endDo
!
      return
      end subroutine integralTransformationN8sameSpin
!
!
      end module test_mod
