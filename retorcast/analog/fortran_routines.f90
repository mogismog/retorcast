! FILE: fortran_routines.f90
!  f2py -c -m fortran_routines fortran_routines.f90
!  This file will have a collection of fortran routines that can be interfaced to python

! ======================

Subroutine rank_analog(trainField,fcstField,svrField,&
     trainnum,members,n_members,iNum,jNum,&
     startLat,endLat,startLon,endLon,allLats,allLons,window,outProbs)
INTEGER, INTENT(IN) :: iNum,jNum,startLat,endLat,startLon,endLon,window,n_members,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
integer, intent(in), dimension(n_members) :: members
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField,svrField
real, intent(out), dimension(n_members,iNum,jNum) :: outProbs

! --- Now, some other variables we'll need
integer :: i,j,np,nmem,startLatIdx,startLonIdx,endLonIdx,endLatIdx,n_grdpts
real, DIMENSION(:,:,:),allocatable :: ranks
real, DIMENSION(:,:,:),allocatable :: allData
integer, DIMENSION(:),allocatable :: finalRanks,dummyIdx,dummy_finalRanks,dummy_ranks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: dummy_array,rankDiffs

!f2py intent(in) iNum,jNum,trainnum,startLat,endLat,startLon,endLon,window,allLats,allLons,trainField,svrField
!f2py intent(in) n_members,members,fcstField
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(n_members) members
!f2py depend(iNum,jNum) fcstField
!f2py depend(trainnum,iNum,jNum) trainField,svrField
!f2py intent(out) outProbs
!f2py depend(n_members,iNum,jNum) outProbs

outProbs(:,:,:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. startLat) then
      startLatIdx = i
   elseif (allLats(i) .eq. endLat) then
      endLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. startLon) then
      startLonIdx = j
   elseif (allLons(j) .eq. endLon) then
      endLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+1)*((window*2)+1)

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(ranks(trainnum+1,iNum,jNum)) ! --- Array of rankings
allocate(dummy_ranks(trainnum+1)) ! --- temporary ranking array
allocate(dummy_array(trainnum+1)) ! --- temporary data array
allocate(trainData(trainnum+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(rankDiffs(trainnum)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum)) ! --- rankings of rankDiffs array
allocate(finalRanks(trainnum)) ! --- rankings of rankDiffs array
allocate(dummyIdx(trainnum)) ! --- indices of training dates

! --- We need to put the training and fcst data in one array
allData(:trainnum,:,:) = trainField(:,:,:) ! --- Training data
allData(trainnum+1,:,:) = fcstField(:,:) ! --- Forecast data

! --- Ranking each grid point in this loop
do j = 1,jNum
   do i = 1,iNum
      dummy_array(:) = allData(:,i,j)
      dummy_ranks(:)=0
      call real_rank(dummy_array,size(dummy_array),dummy_ranks)
      call ties(dummy_ranks,dummy_array,size(dummy_array),ranks(:,i,j))
   end do
end do


! --- Now time to find ranking differences at each grid point given a window of grid points
do j = startLonIdx,endLonIdx
   do i = startLatIdx,endLatIdx

      rankDiffs(:) = 999999
      ! --- First, extract data
      trainData(:,:) = reshape(ranks(:,i-window:i+window,j-window:j+window),(/trainnum+1,n_grdpts/))

      ! --- Find sum of absolute value of rank differences
      do np = 1,trainnum
         rankDiffs(np) = sum(abs(trainData(np,:)-trainData(trainnum+1,:)),dim=1)
      end do

      ! --- Ranking the ranked diffs to make it easier to get probabilities
      call real_rank(rankDiffs,size(rankDiffs),dummy_finalRanks)
      do np = 1,size(dummy_finalRanks)
         finalRanks(dummy_finalRanks(np)) = np
      end do

      !print*,rankDiffs(minloc(finalRanks(:dumbuffer))),rankDiffs(maxloc(finalRanks(:dumbuffer)))
      ! --- Now we find the probabilities
      do nmem = 1,n_members
         outProbs(nmem,i,j) = (sum( svrField(:,i,j),dim=1,&
              mask= (finalRanks(:) .le. members(nmem) )) &
              / real(members(nmem))) * 100.
      end do
   end do
end do

! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(ranks) ! --- Array of rankings

deallocate(dummy_ranks) ! --- temporary ranking array

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

deallocate(finalRanks) ! --- rankings of rankDiffs array

deallocate(dummyIdx) ! --- indices of training dates

deallocate(dummy_array) ! --- temporary data array

return

end Subroutine rank_analog

Subroutine rank_analog_pct(trainField,fcstField,svrField,&
     pctField,pctFcst,trainnum,members,n_members,iNum,jNum,&
     startLat,endLat,startLon,endLon,allLats,allLons,window,outProbs)
INTEGER, INTENT(IN) :: iNum,jNum,startLat,endLat,startLon,endLon,window,n_members,trainnum
real, intent(IN), dimension(iNum) :: allLats
real, intent(IN), dimension(jNum) :: allLons
integer, intent(in), dimension(n_members) :: members
real, INTENT(IN), DIMENSION(iNum,jNum) :: fcstField,pctFcst
real, INTENT(IN), DIMENSION(trainnum,iNum,jNum) :: trainField,svrField,pctField
real, intent(out), dimension(n_members,iNum,jNum) :: outProbs

! --- Now, some other variables we'll need
integer :: i,j,np,nmem,startLatIdx,startLonIdx,endLonIdx,endLatIdx,n_grdpts,extraDts,ii,jj,extrabuff,dumbuffer
real, DIMENSION(:,:),allocatable :: ranks
real, DIMENSION(:,:,:),allocatable :: allData,allPct
integer, DIMENSION(:),allocatable :: finalRanks,dummyIdx,dummy_finalRanks,dummy_ranks
real, dimension(:,:),allocatable :: trainData
real, dimension(:), allocatable :: dummy_array,rankDiffs,trainSvr

!f2py intent(in) iNum,jNum,trainnum,startLat,endLat,startLon,endLon,window,allLats,allLons,trainField,svrField
!f2py intent(in) n_members,members,fcstField,pctField,pctFcst
!f2py depend(iNum) allLats
!f2py depend(jNum) allLons
!f2py depend(n_members) members
!f2py depend(iNum,jNum) fcstField,pctFcst
!f2py depend(trainnum,iNum,jNum) trainField,svrField,pctField
!f2py intent(out) outProbs
!f2py depend(n_members,iNum,jNum) outProbs

outProbs(:,:,:) = -9999.9

! --- Get latitude indices to iterate over
do i = 1,iNum
   if (allLats(i) .eq. startLat) then
      startLatIdx = i
   elseif (allLats(i) .eq. endLat) then
      endLatIdx = i
   end if
end do

! --- Get longitude indices to iterate over
do j = 1,jNum
   if (allLons(j) .eq. startLon) then
      startLonIdx = j
   elseif (allLons(j) .eq. endLon) then
      endLonIdx = j
   end if
end do

! --- Start big coefficient loop
n_grdpts = ((window*2)+1)*((window*2)+1)
extraDts = 1000 ! --- Maximum number of additional dates found to be near same percentile

! --- Let's allocate some arrays
allocate(allData(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(allPct(trainnum+1,iNum,jNum)) ! --- train data + forecast data
allocate(ranks(trainnum+extraDts+1,n_grdpts)) ! --- Array of rankings
allocate(dummy_ranks(trainnum+extraDts+1)) ! --- temporary ranking array
allocate(dummy_array(trainnum+extraDts+1)) ! --- temporary data array
allocate(trainData(trainnum+extraDts+1,((window*2)+1)*((window*2)+1))) ! --- 2-D training data
allocate(trainSvr(trainnum+extraDts)) ! --- 2-D training data
allocate(rankDiffs(trainnum+extraDts)) ! --- mean absolute difference of ranks
allocate(dummy_finalRanks(trainnum+extraDts)) ! --- rankings of rankDiffs array
allocate(finalRanks(trainnum+extraDts)) ! --- rankings of rankDiffs array
allocate(dummyIdx(trainnum)) ! --- indices of training dates

! --- We need to put the training and fcst data in one array
allData(2:trainnum+1,:,:) = trainField(:,:,:) ! --- Training data
allData(1,:,:) = fcstField(:,:) ! --- Forecast data
allPct(1,:,:) = pctFcst(:,:)
allPct(2:trainnum+1,:,:) = pctField(:,:,:) ! --- Training data

! --- Ranking each grid point in this loop
!do j = 1,jNum
!   do i = 1,iNum
!      dummy_array(:) = allData(:,i,j)
!      dummy_ranks(:)=0
!      call real_rank(dummy_array,size(dummy_array),dummy_ranks)
!      call ties(dummy_ranks,dummy_array,size(dummy_array),ranks(:,i,j))
!   end do
!end do

dumbuffer = trainnum
! --- Now time to find ranking differences at each grid point given a window of grid points
do j = startLonIdx,endLonIdx
   do i = startLatIdx,endLatIdx
    
     rankDiffs(:) = 999999
     extrabuff=0
     ! --- First, extract data
     !print*,"extracting training data"
     trainData(:dumbuffer+1,:) = reshape(allData(:dumbuffer+1,i-window:i+window,j-window:j+window),(/dumbuffer+1,n_grdpts/))
     trainSvr(:dumbuffer) = svrField(:,i,j)

     ! --- Now, let's find similar percentiles around the fcst domain and collect their
     !print*,"finding similar percentile fcsts..."
     ! --- fcst data and tornado reports.

     if (allPct(1,i,j) .ge. .99) then
       dateloop: do np = 2,dumbuffer
           do jj = startLonIdx+10,endLonIdx-20
              do ii = startLatIdx+6,endLatIdx-10
                 if (extrabuff .ge. extraDts) exit dateloop
                 if (ii.le.i+1 .and. ii.ge.i-1 .and. jj.le.j+1 .and. jj.ge.j-1) cycle
                 !if (allPct(np,ii,jj) .ge. (allPct(1,i,j)-.01) .and. &
                 !     allPct(np,ii,jj) .le. (allPct(1,i,j)+.01)) then
                 if ((allPct(np,ii,jj) .ge. allPct(1,i,j)) .and. (allData(np,ii,jj).gt.0)) then
                 !if (allPct(np,ii,jj) .ge. .99) then
                    extrabuff = extrabuff+1
                    !print*,"Putting in more training data",np,ii,jj,dumbuffer+extrabuff+1
                    trainData(extrabuff+dumbuffer+1,:) = reshape(allData(np,&
                         ii-window:ii+window,jj-window:jj+window),&
                         (/n_grdpts/))
                     trainSvr(extrabuff+dumbuffer) = svrField(np-1,ii,jj)
                  end if
               end do
            end do
         end do dateloop
      end if
     ! --- Rank the newly attached samples

     do ii = 1,n_grdpts
        call real_rank(trainData(:extrabuff+dumbuffer+1,ii),(extrabuff+dumbuffer+1),dummy_ranks(:extrabuff+dumbuffer+1))
        call ties(dummy_ranks(:extrabuff+dumbuffer+1),trainData(:extrabuff+dumbuffer+1,ii), &
             extrabuff+dumbuffer+1,ranks(:extrabuff+dumbuffer+1,ii))
     end do

         !print*,"finding differences..."
     ! --- Take MAE of absolute value of rank differences
     do np = 2,dumbuffer+extrabuff+1
        rankDiffs(np-1) = sum(abs(ranks(np,:)-ranks(1,:)),dim=1)
     end do

     ! --- Ranking the ranked diffs to make it easier to get probabilities
     call real_rank(rankDiffs(:dumbuffer+extrabuff),size(rankDiffs(:dumbuffer+extrabuff))&
          ,dummy_finalRanks(:dumbuffer+extrabuff))
     finalRanks(:) = 999999
     do np = 1,size(dummy_finalRanks(:dumbuffer+extrabuff))
        finalRanks(dummy_finalRanks(np)) = np
     end do


     !print*,rankDiffs(minloc(finalRanks(:dumbuffer))),rankDiffs(maxloc(finalRanks(:dumbuffer)))
     ! --- Now we find the probabilities
     do nmem = 1,n_members
        outProbs(nmem,i,j) = (sum( trainSvr(:dumbuffer+extrabuff),dim=1,&
             mask= (finalRanks(:dumbuffer+extrabuff) .le. members(nmem) )) &
             / real(members(nmem))) * 100.
     end do
   end do
end do

! --- Time to deallocate those arrays
deallocate(allData) ! --- train data + forecast data

deallocate(ranks) ! --- Array of rankings

deallocate(dummy_ranks) ! --- temporary ranking array

deallocate(trainData) ! --- 2-D training data

deallocate(rankDiffs) ! --- mean absolute difference of ranks

deallocate(dummy_finalRanks) ! --- rankings of rankDiffs array

deallocate(finalRanks) ! --- rankings of rankDiffs array

deallocate(dummyIdx) ! --- indices of training dates

deallocate(dummy_array) ! --- temporary data array

return

end Subroutine rank_analog_pct


Subroutine real_rank (XDONT, NOBS, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! _________________________________________________________
      integer, intent(in) :: NOBS
      Real, Dimension (NOBS), Intent (In) :: XDONT
      Integer, Dimension (NOBS), Intent (Out) :: IRNGT
! __________________________________________________________
      Real :: XVALA, XVALB
!
      Integer, Dimension (NOBS) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
     ! print*, "Here!"
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine real_rank

Subroutine int_rank (XDONT, n_obs, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Integer, intent(in) :: n_obs
      Integer, Dimension (n_obs), Intent (In)  :: XDONT
      Integer, Dimension (n_obs), Intent (Out) :: IRNGT
! __________________________________________________________
      Integer :: XVALA, XVALB
!
      Integer, Dimension (n_obs) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine int_rank

Subroutine slow_real_rank (XVALT, nobs, IRNGT)
!   Ranks array XVALT into index array IRNGT, using merge-sort
! __________________________________________________________
!   This version is not optimized for performance, and is thus
!   not as difficult to read as some other ones.
!   Michel Olagnon - April 2000
! __________________________________________________________
! _________________________________________________________
      integer, intent(in) :: nobs
      Real, Dimension (nobs), Intent (In) :: XVALT
      Integer, Dimension (nobs), Intent (Out) :: IRNGT
! __________________________________________________________
!
      Integer, Dimension (:), Allocatable :: JWRKT
      Integer :: LMTNA, LMTNC
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      If (NVAL <= 0) Then
         Return
      End If
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XVALT(IIND-1) < XVALT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      Allocate (JWRKT(1:NVAL))
      LMTNC = 2
      LMTNA = 2
!
!  Iteration. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
         IWRK = 0
!
!   Loop on merges of A and B into C
!
         Do
            IINDA = IWRKF
            IWRKD = IWRKF + 1
            IWRKF = IINDA + LMTNC
            JINDA = IINDA + LMTNA
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDB = JINDA
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B (no need to do anything)
!
            If (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(JINDA+1))) Then
               IWRK = IWRKF
               Cycle
            End If
!
!  One steps in the C subset, that we create in the final rank array
!
            Do
               If (IWRK >= IWRKF) Then
!
!  Make a copy of the rank array for next iteration
!
                  IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
                  Exit
               End If
!
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (IINDA < JINDA) Then
                  If (IINDB < IWRKF) Then
                     If (XVALT(IRNGT(IINDA+1)) > XVALT(IRNGT(IINDB+1))) &
                    & Then
                        IINDB = IINDB + 1
                        JWRKT (IWRK) = IRNGT (IINDB)
                     Else
                        IINDA = IINDA + 1
                        JWRKT (IWRK) = IRNGT (IINDA)
                     End If
                  Else
!
!  Only A still with unprocessed values
!
                     IINDA = IINDA + 1
                     JWRKT (IWRK) = IRNGT (IINDA)
                  End If
               Else
!
!  Only B still with unprocessed values
!
                  IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
                  IWRK = IWRKF
                  Exit
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!  Clean up
!
      Deallocate (JWRKT)
      Return
!
End Subroutine slow_real_rank

Subroutine cdf_loop(quantile_array,indices,pct,lats,lons,final_array)
integer, intent(in) :: pct,lats,lons
integer, intent(in), dimension(lats,lons) :: indices
real, intent(in), dimension(pct,lats,lons) :: quantile_array
real, intent(out), dimension(lats,lons) :: final_array
integer :: i,j,idx

!f2py intent(in) pct,lats,lons,indices,quantile_array
!f2py depend(lats,lons) indices
!f2py depend(pct,lats,lons) quantile_array
!f2py intent(out) final_array
!f2py depend(lats,lons) final_array

!final_array(:,:) = 0.0
do i = 1,lats
   do j = 1,lons
      idx = (indices(i,j))+1 !python index counting -> fortran index counting
      final_array(i,j) = quantile_array(idx,i,j)
   end do
end do
return
end Subroutine cdf_loop

Subroutine pct_loop(indices,percentages,pct,lats,lons,final_array)
integer, intent(in) :: pct,lats,lons
real, intent(in), dimension(pct) :: percentages
integer, intent(in), dimension(lats,lons) :: indices
real, intent(out), dimension(lats,lons) :: final_array
integer :: i,j,idx

!f2py intent(in) pct,lats,lons,indices,percentages
!f2py depend(pct) percentages
!f2py depend(lats,lons) indices
!f2py intent(out) final_array
!f2py depend(lats,lons) final_array

!final_array(:,:) = 0.0
do i = 1,lats
   do j = 1,lons
      idx = (indices(i,j))+1 !python index counting -> fortran index counting
      final_array(i,j) = percentages(idx)
   end do
end do
return
end Subroutine pct_loop

subroutine ties(ranks,vals,n_ranks,outranks)
! --- Subroutine designed to account for ties in rankings, something
! --- this other ranking program I have doesn't do (ugh).

integer, intent(in) :: n_ranks
integer, intent(in) :: ranks(n_ranks)
real, intent(in) :: vals(n_ranks)
real, intent(out) :: outranks(n_ranks)
integer :: idx,begin_idx,idx2,intcount,intcount2
real :: curr_val,counter,rank_sum

counter=0.
intcount = 0
curr_val = 0.
rank_sum=0.
begin_idx=1
outranks(:)=0.
do idx = 1,n_ranks
   !write(71,*),"HERE:",vals(ranks(idx)),curr_val,rank_sum,counter,begin_idx,idx,n_ranks
   if (vals(ranks(idx)) .eq. curr_val) then
      counter=counter+1.
      intcount = intcount+1
      rank_sum=rank_sum+real(idx)
      !write(83,*),"RANK_SUM:",rank_sum,idx,real(idx)
   else if (vals(ranks(idx)) .ne. curr_val) then
      if (counter .gt.0) then
         intcount2=1
         idx2 = (idx-intcount)
         do while (intcount2 .le. intcount)
            outranks(ranks(idx2)) = (rank_sum/counter)
            idx2=idx2+1
            intcount2=intcount2+1
         end do
         !print*,rank_sum,counter,(rank_sum/counter),outranks(begin_idx),outranks(idx-1)
         !write(83,*),"ranksum & counter",rank_sum,counter,curr_val,vals(ranks(idx-1))
      else if (idx.eq.1) then
         counter=0.
         intcount = 0
         curr_val = vals(ranks(idx))
         rank_sum=0.
         begin_idx=idx
      else if (counter .eq. 0 .and. (idx.ne.n_ranks) .and. (idx.ne.1)) then
         outranks(ranks(idx-1)) = idx-1
         !write(83,*),"ranking...",idx,outranks((ranks(idx-1)))
      else if ((idx.eq.n_ranks) .and. (counter .eq. 0)) then
         outranks(ranks(idx-1)) = idx-1
         outranks(ranks(idx)) = idx
         !print*,"last!",outranks(ranks(idx)),outranks(ranks(idx-1)),idx,n_ranks
      else if ((idx.eq.n_ranks) .and. (counter .gt. 0)) then
         rank_sum = rank_sum+real(idx)
         intcount2=1
         idx2 = (idx-intcount)
         do while (intcount2 .le. intcount)
            outranks(ranks(idx2)) = (rank_sum/counter)
            idx2=idx2+1
            intcount2=intcount2+1
         end do
         outranks(ranks(idx)) = idx

      end if
      counter=0.
      intcount = 0
      curr_val = vals(ranks(idx))
      rank_sum=0.
      begin_idx=idx
   end if
end do

return

end subroutine ties
