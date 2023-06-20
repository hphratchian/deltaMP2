      module integraltransformation_mod
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
      subroutine commandLineArgs(iPrint,nOMP,matrixFilename,doN8,  &
        doSlowN5,doRegularN5,useBLAS,fail)
!
!     This subroutine is used to process the command line arguments.
!
      integer,intent(OUT)::iPrint,nOMP
      character(len=512),intent(OUT)::matrixFilename
      logical,intent(OUT)::doN8,doSlowN5,doRegularN5,useBLAS,fail
!
      integer::nCommands,nMatrixFilenames,i
      character(len=512)::tmpString,lowercase
      logical::getNProc
!
!     Format statements.
!
 9000 format(1x,'Unknown command line switch: ',A,'.')
 9100 format(1x,'More than 1 matrix filename found!')
!
!     Set defaults.
!
      iPrint = 0
      nOMP = 1
      doN8 = .false.
      doSlowN5 = .false.
      doRegularN5 = .false.
      useBLAS = .true.
      fail = .false.
      nMatrixFilenames = 0
      getNProc = .false.
!
!     Determine the number of command line arguments. Then, loop through the
!     list of arguments to set options.
!
      nCommands = command_argument_count()
      do i = 1,nCommands
        call get_command_argument(i,tmpString)
        if(getNProc) then
          read(tmpString,'(i)') nOMP
          getNProc = .false.
          cycle
        endIf
        if(tmpString(1:1).eq.'-') then
          call String_Change_Case(tmpString(2:),'l',lowercase)
          select case(TRIM(lowercase))
          case('debug')
            iPrint = 1
          case('nproc')
            getNProc = .true.
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
          case('useblas')
            useBLAS = .true.
          case('usematmul')
            useBLAS = .false.
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
      end module integraltransformation_mod
