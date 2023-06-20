INCLUDE 'integraltransformation_mod.f03'
      Program integraltransformation
!
!     This program reads AO integrals from a Gaussian matrix file and times
!     AO-to-MO integral transformations.
!
!     -H. P. Hratchian, 2022.
!
!
!     USE Connections
!
      use integraltransformation_mod
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64)::nCommands,iPrint=0,nOMP,nAtoms,nBasis,  &
        nBasisUse,nElectrons,nElectronsAlpha,nElectronsBeta
      integer(kind=int64)::mu,nu,lambda,sigma,p,q,r,s,pq,rs,pqrs,iCount,i,j,a,b
      real(kind=real64)::timeStart,timeEnd,time0,time1,tmpReal,  &
        deltaIJAB,numerator,E2AA,E2BB,E2AB,E2BA
      real(kind=real64),dimension(:),allocatable::moEnergiesAlpha,moEnergiesBeta
      real(kind=real64),dimension(:,:),allocatable::CAlpha,CBeta,  &
        tmpMatrix1,tmpMatrix2
      real(kind=real64),dimension(:,:,:,:),allocatable::aoInts,moInts,  &
        partialInts1,partialInts2
      type(MQC_Variable)::ERIs,mqcTmpArray
      character(len=512)::tmpString,matrixFilename
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      logical::fail=.false.,doN8=.false.,doSlowN5=.true.,  &
        doRegularN5=.true.,useBLAS=.true.
!
!     Format Statements
!
 1000 Format(1x,'Enter Test Program integralTransformation.')
 1010 Format(1x,'Matrix File: ',A,/)
 1020 Format(1x,'Use ',I3,' shared memory processors.')
 1100 Format(1x,'nAtoms    =',I4,6x,'nBasis  =',I4,6x,'nBasisUse=',I4,/,  &
             1x,'nElectrons=',I4,6x,'nElAlpha=',I4,6x,'nElBeta  =',I4)
 1500 Format(/,1x,'Carrying out O(N^8) transformation.')
 1510 Format(/,1x,'Skipping O(N^8) transformation.')
 2000 Format(1x,'<',I3,',',I3,' || ',I3,',',I3,' > ... pq=',I3,'  rs=',I3,'  pqrs=',I3)
 3000 Format(/,1x,'E(2)-SS = ',f15.10,' a.u.',4x,'E(2)-OS = ',f15.10,' a.u.')
 5000 Format(1x,'Time (',A,'): ',f8.1,' s.')
 8999 Format(/,1x,'END OF PROGRAM integralTransformation.')
!
!
      call cpu_time(timeStart)
      write(IOut,1000)
      call mqc_version_print(iOut)
      if(.not.mqc_version_check(newerThanMajor=22,newerThanMinor=12,newerThanRevision=1))  &
        call mqc_error('MQCPack version is too old.')
!
!     Open the Gaussian matrix file and load the number of atomic centers.

      nCommands = command_argument_count()
      if(nCommands.eq.0)  &
        call mqc_error('No command line arguments provided. The input Gaussian matrix file name is required.')
      call get_command_argument(1,tmpString)
      call commandLineArgs(iPrint,nOMP,matrixFilename,doN8,doSlowN5,  &
        doRegularN5,useBLAS,fail)
      call omp_set_num_threads(nOMP)
      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      write(iOut,1020) nOMP
      nAtoms = GMatrixFile%getVal('nAtoms')
      nBasis = Int(GMatrixFile%getVal('nbasis'))
      nBasisUse = Int(GMatrixFile%getVal('nbasisuse'))
      nElectrons = Int(GMatrixFile%getVal('nelectrons'))
      nElectronsAlpha = Int(GMatrixFile%getVal('nAlpha'))
      nElectronsBeta = Int(GMatrixFile%getVal('nBeta'))
      write(IOut,1100) nAtoms,nBasis,nBasisUse,nElectrons,  &
        nElectronsAlpha,nElectronsBeta
!
!     Load the orbital eigenvalues.
!
      call GMatrixFile%getArray('ALPHA ORBITAL ENERGIES',mqcVarOut=mqcTmpArray)
      call mqcTmpArray%print(header='Alpha MO Energies')
      flush(iOut)
      moEnergiesAlpha = mqcTmpArray
      if(GMatrixFile%isUnrestricted()) then
        call mqc_error('UHF/UKS NYI.')
      else
        moEnergiesBeta = moEnergiesAlpha
      endIf
!
!     Load the MO coefficients.
!
      write(*,*)' Hrant -  isUnrestricted: ',GMatrixFile%isUnrestricted()
      call GMatrixFile%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=mqcTmpArray)
      CAlpha = mqcTmpArray
      if(GMatrixFile%isUnrestricted()) then
        call GMatrixFile%getArray('BETA  MO COEFFICIENTS',mqcVarOut=mqcTmpArray)
        CBeta = mqcTmpArray
      else
        CBeta = CAlpha
      endIf
!
!     Read in and report out the (AO) ERIs.
!
      call GMatrixFile%getArray('REGULAR 2E INTEGRALS',mqcVarOut=ERIs)
      if(iPrint.ge.2) call ERIs%print(IOut,' ERIs=')
!
!     Allocate space for AO and MO integrals. Then fill the intrinsic array of
!     AO integrals.
!
      allocate(aoInts(nBasis,nBasis,nBasis,nBasis),  &
        moInts(nBasisUse,nBasisUse,nBasisUse,nBasisUse))

      write(*,*)' Hrant - FLAG B'
      flush(iOut)

!hph+
!      aoInts = ERIs
      call dpReshape4(nBasis,nBasis,nBasis,nBasis,ERIs%realArray,aoInts)
!hph-

      if(iPrint.ge.2) call mqc_print_rank4Tensor_array_real(iOut,  &
        aoInts,header='Intrinsic AO Integrals')

      write(*,*)' Hrant - FLAG C'
      flush(iOut)

!
!     Do N^8 AO --> MO transformation.
!
      if(doN8) then
        write(iOut,1500)
        call cpu_time(time0)
        moInts = float(0)
        call cpu_time(time1)
        write(iOut,5000) 'Init moInts',time1-time0
        flush(iOut)
!
        call cpu_time(time0)
        do p = 1,nBasisUse
          do q = 1,nBasisUse
            do r = 1,nBasisUse
              do s = 1,nBasisUse
                do mu = 1,nBasis
                  do nu = 1,nBasis
                    do lambda = 1,nBasis
                      do sigma = 1,nBasis
                        moInts(p,q,r,s) = moInts(p,q,r,s)  &
                          + CAlpha(mu,p)*CAlpha(nu,q)*CAlpha(lambda,r)*CAlpha(sigma,s)*aoInts(mu,nu,lambda,sigma)
                      endDo
                    endDo
                  endDo
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'O(N^8) Transformation',time1-time0
      else
        write(iOut,1510)
      endIf
      flush(iOut)
!
!     Do quarter transformations now.
!
      call cpu_time(time0)
      moInts = float(0)
      Allocate(partialInts1(nBasisUse,nBasis,nBasis,nBasis))
      partialInts1 = float(0)
      call cpu_time(time1)
      write(iOut,5000) 'Init moInts and partialInts1',time1-time0
      flush(iOut)
!
      if(doSlowN5) then
        call cpu_time(time0)
        do p = 1,nBasisUse
          do mu = 1,nBasis
            do nu = 1,nBasis
              do lambda = 1,nBasis
                do sigma = 1,nBasis
                  partialInts1(p,nu,lambda,sigma) = partialInts1(p,nu,lambda,sigma)  &
                    + CAlpha(mu,p)*aoInts(mu,nu,lambda,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 1a',time1-time0
        flush(iOut)
        partialInts1 = float(0)
      endIf
!
      if(doRegularN5) then
        call cpu_time(time0)
        do nu = 1,nBasis
          do lambda = 1,nBasis
            do sigma = 1,nBasis
              do p = 1,nBasisUse
                tmpReal = float(0)
                do mu = 1,nBasis
                  partialInts1(p,nu,lambda,sigma) = partialInts1(p,nu,lambda,sigma)  &
                    + CAlpha(mu,p)*aoInts(mu,nu,lambda,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 1b',time1-time0
        flush(iOut)
      endIf
!
      Allocate(partialInts2(nBasisUse,nBasisUse,nBasis,nBasis))
      if(doSlowN5) then
        call cpu_time(time0)
        do p = 1,nBasisUse
          do q = 1,nBasisUse
            do nu = 1,nBasis
              do lambda = 1,nBasis
                do sigma = 1,nBasis
                  partialInts2(p,q,lambda,sigma) = partialInts2(p,q,lambda,sigma)  &
                    + CAlpha(nu,q)*partialInts1(p,nu,lambda,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 2a',time1-time0
        flush(iOut)
        partialInts2 = float(0)
      endIf
!
      if(doRegularN5) then
        call cpu_time(time0)
        do lambda = 1,nBasis
          do sigma = 1,nBasis
            do q = 1,nBasisUse
              do p = 1,nBasisUse
                do nu = 1,nBasis
                  partialInts2(p,q,lambda,sigma) = partialInts2(p,q,lambda,sigma)  &
                    + CAlpha(nu,q)*partialInts1(p,nu,lambda,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 2b',time1-time0
        flush(iOut)
      endIf
!
      DeAllocate(partialInts1)
      Allocate(partialInts1(nBasisUse,nBasisUse,nBasisUse,nBasis))
      partialInts1 = float(0)
!
      if(doSlowN5) then
        call cpu_time(time0)
        do p = 1,nBasisUse
          do q = 1,nBasisUse
            do r = 1,nBasisUse
              do lambda = 1,nBasis
                do sigma = 1,nBasis
                  partialInts1(p,q,r,sigma) = partialInts1(p,q,r,sigma)  &
                    + CAlpha(lambda,r)*partialInts2(p,q,lambda,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 3a',time1-time0
        flush(iOut)
        partialInts1 = float(0)
      endIf
!
      if(doRegularN5) then
        call cpu_time(time0)
        do p = 1,nBasisUse
          do q = 1,nBasisUse
                do sigma = 1,nBasis
            do r = 1,nBasisUse
              do lambda = 1,nBasis
                  partialInts1(p,q,r,sigma) = partialInts1(p,q,r,sigma)  &
                    + CAlpha(lambda,r)*partialInts2(p,q,lambda,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 3b',time1-time0
        flush(iOut)
      endIf
!
      DeAllocate(partialInts2)
      moInts = float(0)
!
      if(doRegularN5.or.doSlowN5) then
        call cpu_time(time0)
        do p = 1,nBasisUse
          do q = 1,nBasisUse
            do r = 1,nBasisUse
              do s = 1,nBasisUse
                do sigma = 1,nBasis
                  moInts(p,q,r,s) = moInts(p,q,r,s)  &
                    + CAlpha(sigma,s)*partialInts1(p,q,r,sigma)
                endDo
              endDo
            endDo
          endDo
        endDo
        call cpu_time(time1)
        write(iOut,5000) 'Quarter Transformation 4',time1-time0
        flush(iOut)
        DeAllocate(partialInts1)
        if(iPrint.ge.1) call mqc_print_rank4Tensor_array_real(iOut,  &
          moInts,header='Transformed MO Integrals 2')
      endIf
!
!     Load up AA MO ERIs from the matrixfile and print them out to ensure the
!     explicit AO-->MO transformations above gave the right answers.
!
      call GMatrixFile%getArray('AA MO 2E INTEGRALS',mqcVarOut=mqcTmpArray)
      if(iPrint.ge.1) call mqcTmpArray%print(header='MO Ints from Matrix File')
!
!     Now, try some matrix multiplication based approach to the quarter
!     transformations using the intrinsic Transform and MatMul functions.
!
      Allocate(tmpMatrix1(nBasisUse,nBasis*nBasis*nBasis))
      Allocate(tmpMatrix2(nBasis,nBasis*nBasis*nBasisUse))
      call cpu_time(time0)
      if(useBLAS) then
        call dgemm('t','n',nBasisUse,nBasis*nBasis*nBasis,  &
          nBasis,float(1),CAlpha,nBasis,aoInts,nBasis,float(0),  &
          tmpMatrix1,nBasisUse)
      else
        tmpMatrix1 = MatMul(Transpose(CAlpha),  &
          Reshape(aoInts,(/nBasis,nBasis*nBasis*nBasis/)))
      endIf
      call cpu_time(time1)
      write(iOut,5000) 'Quarter Transformation 1 w/MatMul',time1-time0
      flush(iOut)
!
      call cpu_time(time0)
      tmpMatrix2 = Reshape(Transpose(tmpMatrix1),  &
        (/nBasis,nBasis*nBasis*nBasisUse/))
      call cpu_time(time1)
      write(iOut,5000) 'tmpMatrix1 transpose 1',time1-time0
      flush(iOut)
      call cpu_time(time0)
      DeAllocate(tmpMatrix1)
      Allocate(tmpMatrix1(nBasisUse,nBasis*nBasis*nBasisUse))
      if(useBLAS) then
        call dgemm('t','n',nBasisUse,nBasis*nBasis*nBasisUse,  &
          nBasis,float(1),CAlpha,nBasis,tmpMatrix2,nBasis,float(0),  &
          tmpMatrix1,nBasisUse)
      else
        tmpMatrix1 = MatMul(Transpose(CAlpha),tmpMatrix2)
      endIf
      call cpu_time(time1)
      write(iOut,5000) 'Quarter Transformation 2 w/MatMul',time1-time0
      flush(iOut)
!
      call cpu_time(time0)
      DeAllocate(tmpMatrix2)
      Allocate(tmpMatrix2(nBasis,nBasis*nBasisUse*nBasisUse))
      tmpMatrix2 = Reshape(Transpose(tmpMatrix1),  &
        (/nBasis,nBasis*nBasisUse*nBasisUse/))
      call cpu_time(time1)
      write(iOut,5000) 'tmpMatrix1 transpose 2',time1-time0
      flush(iOut)
      call cpu_time(time0)
      DeAllocate(tmpMatrix1)
      Allocate(tmpMatrix1(nBasisUse,nBasis*nBasisUse*nBasisUse))
      if(useBLAS) then
        call dgemm('t','n',nBasisUse,nBasis*nBasisUse*nBasisUse,  &
          nBasis,float(1),CAlpha,nBasis,tmpMatrix2,nBasis,float(0),  &
          tmpMatrix1,nBasisUse)
      else
        tmpMatrix1 = MatMul(Transpose(CAlpha),tmpMatrix2)
      endIf
      call cpu_time(time1)
      write(iOut,5000) 'Quarter Transformation 3 w/MatMul',time1-time0
      flush(iOut)
!
      call cpu_time(time0)
      DeAllocate(tmpMatrix2)
      Allocate(tmpMatrix2(nBasis,nBasisUse*nBasisUse*nBasisUse))
      tmpMatrix2 = Reshape(Transpose(tmpMatrix1),  &
        (/nBasis,nBasisUse*nBasisUse*nBasisUse/))
      call cpu_time(time1)
      write(iOut,5000) 'tmpMatrix1 transpose 3',time1-time0
      flush(iOut)
      call cpu_time(time0)
      DeAllocate(tmpMatrix1)
      Allocate(tmpMatrix1(nBasisUse,nBasisUse*nBasisUse*nBasisUse))
      if(useBLAS) then
        call dgemm('t','n',nBasisUse,nBasisUse*nBasisUse*nBasisUse,  &
          nBasis,float(1),CAlpha,nBasis,tmpMatrix2,nBasis,float(0),  &
          tmpMatrix1,nBasisUse)
      else
        tmpMatrix1 = MatMul(Transpose(CAlpha),tmpMatrix2)
      endIf
      call cpu_time(time1)
      write(iOut,5000) 'Quarter Transformation 4 w/MatMul',time1-time0
      flush(iOut)
!
      call cpu_time(time0)
      moInts = Reshape(Transpose(tmpMatrix1),  &
        (/nBasisUse,nBasisUse,nBasisUse,nBasisUse/))
      DeAllocate(tmpMatrix1,tmpMatrix2)
      call cpu_time(time1)
      write(iOut,5000) 'Finalize moInts w/MatMul',time1-time0
      flush(iOut)

!hph+
!!
!!     Try this again...
!!
!      Allocate(tmpMatrix1(nBasisUse,nBasis*nBasis*nBasis))
!      Allocate(tmpMatrix2(nBasis,nBasis*nBasis*nBasisUse))
!      call cpu_time(time0)
!      call dgemm('t','n',nBasisUse,nBasis*nBasis*nBasis,  &
!        nBasis,float(1),CAlpha,nBasis,aoInts,nBasis,float(0),  &
!        tmpMatrix1,nBasisUse)
!      call cpu_time(time1)
!      write(iOut,5000) 'Quarter Transformation 1',time1-time0
!      flush(iOut)
!!
!      call cpu_time(time0)
!      tmpMatrix2 = Reshape(Transpose(tmpMatrix1),  &
!        (/nBasis,nBasis*nBasis*nBasisUse/))
!      call cpu_time(time1)
!      write(iOut,5000) 'tmpMatrix1 transpose 1',time1-time0
!      flush(iOut)
!      call cpu_time(time0)
!      DeAllocate(tmpMatrix1)
!      Allocate(tmpMatrix1(nBasisUse,nBasis*nBasis*nBasisUse))
!      if(useBLAS) then
!        call dgemm('t','n',nBasisUse,nBasis*nBasis*nBasisUse,  &
!          nBasis,float(1),CAlpha,nBasis,tmpMatrix2,nBasis,float(0),  &
!          tmpMatrix1,nBasisUse)
!      else
!        tmpMatrix1 = MatMul(Transpose(CAlpha),tmpMatrix2)
!      endIf
!      call cpu_time(time1)
!      write(iOut,5000) 'Quarter Transformation 2 w/MatMul',time1-time0
!      flush(iOut)
!!
!      call cpu_time(time0)
!      DeAllocate(tmpMatrix2)
!      Allocate(tmpMatrix2(nBasis,nBasis*nBasisUse*nBasisUse))
!      tmpMatrix2 = Reshape(Transpose(tmpMatrix1),  &
!        (/nBasis,nBasis*nBasisUse*nBasisUse/))
!      call cpu_time(time1)
!      write(iOut,5000) 'tmpMatrix1 transpose 2',time1-time0
!      flush(iOut)
!      call cpu_time(time0)
!      DeAllocate(tmpMatrix1)
!      Allocate(tmpMatrix1(nBasisUse,nBasis*nBasisUse*nBasisUse))
!      if(useBLAS) then
!        call dgemm('t','n',nBasisUse,nBasis*nBasisUse*nBasisUse,  &
!          nBasis,float(1),CAlpha,nBasis,tmpMatrix2,nBasis,float(0),  &
!          tmpMatrix1,nBasisUse)
!      else
!        tmpMatrix1 = MatMul(Transpose(CAlpha),tmpMatrix2)
!      endIf
!      call cpu_time(time1)
!      write(iOut,5000) 'Quarter Transformation 3 w/MatMul',time1-time0
!      flush(iOut)
!!
!      call cpu_time(time0)
!      DeAllocate(tmpMatrix2)
!      Allocate(tmpMatrix2(nBasis,nBasisUse*nBasisUse*nBasisUse))
!      tmpMatrix2 = Reshape(Transpose(tmpMatrix1),  &
!        (/nBasis,nBasisUse*nBasisUse*nBasisUse/))
!      call cpu_time(time1)
!      write(iOut,5000) 'tmpMatrix1 transpose 3',time1-time0
!      flush(iOut)
!      call cpu_time(time0)
!      DeAllocate(tmpMatrix1)
!      Allocate(tmpMatrix1(nBasisUse,nBasisUse*nBasisUse*nBasisUse))
!      if(useBLAS) then
!        call dgemm('t','n',nBasisUse,nBasisUse*nBasisUse*nBasisUse,  &
!          nBasis,float(1),CAlpha,nBasis,tmpMatrix2,nBasis,float(0),  &
!          tmpMatrix1,nBasisUse)
!      else
!        tmpMatrix1 = MatMul(Transpose(CAlpha),tmpMatrix2)
!      endIf
!      call cpu_time(time1)
!      write(iOut,5000) 'Quarter Transformation 4 w/MatMul',time1-time0
!      flush(iOut)
!!
!      call cpu_time(time0)
!      moInts = Reshape(Transpose(tmpMatrix1),  &
!        (/nBasisUse,nBasisUse,nBasisUse,nBasisUse/))
!      DeAllocate(tmpMatrix1,tmpMatrix2)
!      call cpu_time(time1)
!      write(iOut,5000) 'Finalize moInts w/MatMul',time1-time0
!      flush(iOut)
!
!
!!
!!     Test the timing of N^2,N^2 transpose.
!!
!      call cpu_time(time0)
!      Allocate(tmpMatrix1(nBasisUse*nBasisUse,nBasisUse*nBasisUse))
!      tmpMatrix1 = Reshape(moInts,  &
!        (/nBasisUse*nBasisUse,nBasisUse*nBasisUse/))
!      call cpu_time(time1)
!      write(iOut,5000) 'Time to allocate tmpMatrix and reshape mo ERIs.',time1-time0
!      call cpu_time(time0)
!      Allocate(tmpMatrix2(nBasisUse*nBasisUse,nBasisUse*nBasisUse))
!      tmpMatrix2 = Transpose(tmpMatrix1)
!      call cpu_time(time1)
!      write(iOut,5000) 'Time for ERI N^2 x N^2 transpose.',time1-time0
!hph-

      
      
!
!     Evaluate the E(2) AA, BB, AB, and BA contribution.
!
      write(*,*)
      write(*,*)' Same Spin E2...'
      E2AA = float(0)
      call cpu_time(time0)
      do i = 1,nElectronsAlpha
        do j = 1,nElectronsAlpha
          do a = nElectronsAlpha+1,nBasisUse
            do b = nElectronsAlpha+1,nBasisUse
              deltaIJAB = moEnergiesAlpha(i) + moEnergiesAlpha(j)  &
                - moEnergiesAlpha(a) - moEnergiesAlpha(b)
              numerator = moInts(i,a,j,b) - moInts(i,b,j,a)
              numerator = numerator*numerator
              if(iPrint.ge.1) write(*,*)' num, denom = ',numerator,deltaIJAB
              E2AA = E2AA + numerator/(float(4)*deltaIJAB)
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Same Spin E2',time1-time0
      E2BB = E2AA
      write(*,*)
      write(*,*)' Opposite Spin E2...'
      E2AB = float(0)
      call cpu_time(time0)
      do i = 1,nElectronsAlpha
        do j = 1,nElectronsBeta
          do a = nElectronsAlpha+1,nBasisUse
            do b = nElectronsBeta+1,nBasisUse
              deltaIJAB = moEnergiesAlpha(i) + moEnergiesBeta(j)  &
                - moEnergiesAlpha(a) - moEnergiesBeta(b)
              numerator = moInts(i,a,j,b)*moInts(i,a,j,b)
              if(iPrint.ge.1) write(*,*)' num, denom = ',numerator,deltaIJAB
              E2AB = E2AB + numerator/deltaIJAB
            endDo
          endDo
        endDo
      endDo
      call cpu_time(time1)
      write(iOut,5000) 'Opposite Spin E2',time1-time0
      E2BA = E2AB
      write(iOut,3000) E2AA,E2AB
      write(*,*)' Total E2 = ',E2AA+E2BB+E2AB
!
  999 Continue
      call cpu_time(timeEnd)
      write(iOut,5000) 'TOTAL JOB TIME',timeEnd-timeStart
      write(iOut,8999)
      end program integraltransformation


      subroutine dpReshape4(N1,N2,N3,N4,arrayIn,r4ArrayOut)
!
      use iso_fortran_env
      implicit none
      integer(kind=int64)::N1,N2,N3,N4
      real(kind=real64),dimension(N1,N2,N3,N4)::arrayIn,r4ArrayOut
!
      r4ArrayOut = arrayIn
!
      return
      end subroutine dpReshape4
