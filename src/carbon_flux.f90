PROGRAM CARBON

    ! ************************************************************************
    !
    !
    !         CO2 FLUXES PER LATITUDE BAND
    !
    !
    ! ************************************************************************

    ! Version 1.0: written by Aitor Aldama-Campino, Autumn 2018

        
    USE netcdf

    IMPLICIT NONE

    INTEGER, PARAMETER     :: IX = 1440, JYD = 1, JYU = 1021, JY = 1021, KZ = 75
    INTEGER, DIMENSION(4)  :: start2D, count2D, start3D, count3D
    INTEGER                :: ncid, varid, ierr, ii, jj, jp, kk, nindex, mm

    REAL, DIMENSION(IX,JYD:JY)    :: sou, inp, atl, glo, dx, dy, co2,lat
    REAL, DIMENSION(IX,JYD:JY,9)  :: mask
    REAL, DIMENSION(IX,JYD:JY,KZ) :: ctot
    REAL, DIMENSION(JY)           :: latt
    REAL, DIMENSION(120,JY)       :: co2lat

    REAL, DIMENSION(120,9):: co2flux

    LOGICAL :: file_exists

    CHARACTER (LEN=400)                :: filepath, commandline
    CHARACTER (LEN=20), DIMENSION(120) :: label1, label2, label3

    ! ***************************************************************

    ! Path to the data input directory
    filepath = '/User/Input/'

    ! READING REGIONS AND MESH
    ! ===============================================================


    start2D = [ 1,JYD,1,1]
    count2D = [IX,JY-JYD+1,1,1]

    ! Basins
    ierr=NF90_OPEN(trim(filepath)//'regions.nc',NF90_NOWRITE,ncid)
    IF(ierr.ne.0) STOP 'No region file'
    ierr=NF90_INQ_VARID(ncid,'tmasksou',varid) ! the main data fields
    IF(ierr.ne.0) STOP 1
    ierr=NF90_GET_VAR(ncid,varid,sou,start2d,count2d)
    IF(ierr.ne.0) STOP 2
    ierr=NF90_INQ_VARID(ncid,'tmaskinp',varid) ! the main data fields
    IF(ierr.ne.0) STOP 3
    ierr=NF90_GET_VAR(ncid,varid,inp,start2d,count2d)
    IF(ierr.ne.0) STOP 4
     ierr=NF90_INQ_VARID(ncid,'tmaskatl',varid) ! the main data fields
    IF(ierr.ne.0) STOP 5
    ierr=NF90_GET_VAR(ncid,varid,atl,start2d,count2d)
    IF(ierr.ne.0) STOP 6
    ierr=NF90_CLOSE(ncid)

    ! Dx and Dz
    ierr=NF90_OPEN(trim(filepath)//'mesh_mask.nc',NF90_NOWRITE,ncid)
    IF(ierr.ne.0) STOP 'No mesh mask file'
    ierr=NF90_INQ_VARID(ncid,'e1t',varid) ! the main data fields
    IF(ierr.ne.0) STOP 7
    ierr=NF90_GET_VAR(ncid,varid,dx,start2d,count2d)
    IF(ierr.ne.0) STOP 8
    ierr=NF90_INQ_VARID(ncid,'e2t',varid) ! the main data fields
    IF(ierr.ne.0) STOP 9
    ierr=NF90_GET_VAR(ncid,varid,dy,start2d,count2d)
    IF(ierr.ne.0) STOP 10
    ierr=NF90_INQ_VARID(ncid,'nav_lat',varid) ! the main data fields
    IF(ierr.ne.0) STOP 11
    ierr=NF90_GET_VAR(ncid,varid,latt,[401,1],[1,JY])
    IF(ierr.ne.0) STOP 12
    ierr=NF90_CLOSE(ncid)

    glo = atl + inp + sou

    DO ii=1,IX
        lat(ii,:) = latt
    END DO

    PRINT*, '- Reading mesh and regions done'

    mask = 0
    mask(:,:302,1) = 1
    mask(:,303:425,2) = 1
    mask(:885,:,2) = 0
    mask(1241:,:,2) = 0
    mask(:,303:425,3) = 1
    mask(886:1240,:,3) = 0
    mask(:,426:572,4) = 1
    mask(:,426:572,5) = 1
    mask(:,573:721,6) = 1
    mask(:,573:721,7) = 1
    mask(:,722:,8) = 1

    PRINT*, '- Masking mask'

    ! ***************************************************************

    ! DATE
    ! ===============================================================
    kk = 1
    DO ii = 2000,2009
        DO jj = 1,12
            WRITE(label1(kk),'(i4,A5,i4,A4)') ii,'0101_',ii,'1230'
            IF (jj<10) THEN
                WRITE(label2(kk),'(i4,A1,i1)') ii,'0',jj
            ELSE
                WRITE(label2(kk),'(i4,i02)') ii,jj
            ENDIF
            WRITE(label3(kk),'(i4)') ii

            kk = kk + 1
        END DO
    END DO

    ! ***************************************************************

    ! BIG LOOP
    ! ===============================================================

    start2D = [ 1,JYD,1,1]
    count2D = [IX,JY-JYD+1,1,1]

    co2flux = 0.
    co2lat  = 0.

    DO nindex = 1, SIZE(label1)

        PRINT*, '- Time step ::', nindex

        ! Tracer file
        ! Read the CO_2 fluxes
        ierr=NF90_OPEN(trim(filepath)//"trim_ORCA025-ROAM_1m_"&
            //TRIM(label1(nindex))//"_diad_T_"//TRIM(label2(nindex))//'-'//TRIM(label2(nindex))//".nc",NF90_NOWRITE,ncid)
        IF(ierr.ne.0) STOP 'No C2 file'
        ierr=NF90_INQ_VARID(ncid,'CO2FLUX',varid) ! the main data fields
        IF(ierr.ne.0) STOP 18
        ierr=NF90_GET_VAR(ncid,varid,co2,start2d,count2d)
        IF(ierr.ne.0) STOP 19
        ierr=NF90_CLOSE(ncid)

        WHERE(co2>1e19) co2 = 0.

        co2flux(nindex,1) = SUM(SUM(co2*dx*dy*mask(:,:,1)*(sou+inp+atl),1),1)
        co2flux(nindex,2) = SUM(SUM(co2*dx*dy*mask(:,:,2)*(sou+atl),1),1)
        co2flux(nindex,3) = SUM(SUM(co2*dx*dy*mask(:,:,3)*(sou+inp),1),1)
        co2flux(nindex,4) = SUM(SUM(co2*dx*dy*mask(:,:,4)*(atl),1),1)
        co2flux(nindex,5) = SUM(SUM(co2*dx*dy*mask(:,:,5)*(inp),1),1)
        co2flux(nindex,6) = SUM(SUM(co2*dx*dy*mask(:,:,6)*(atl),1),1)
        co2flux(nindex,7) = SUM(SUM(co2*dx*dy*mask(:,:,7)*(inp),1),1)
        co2flux(nindex,8) = SUM(SUM(co2*dx*dy*mask(:,:,8)*(atl+inp),1),1)
        co2flux(nindex,9) = SUM(SUM(co2*dx*dy*(glo),1),1)
        co2lat(nindex,:) = SUM(co2*dx*dy*(glo),1)

    END DO

    co2flux = co2flux*365.25*12.01*1e-18
    co2lat = co2lat*365.25*12.01*1e-18

    OPEN(UNIT=11,FILE=TRIM(filepath)//"result_co2flux.bin",STATUS='REPLACE', FORM='UNFORMATTED')
    WRITE(11) co2flux
    CLOSE(11)

    OPEN(UNIT=11,FILE=TRIM(filepath)//"result_co2lat.bin",STATUS='REPLACE', FORM='UNFORMATTED')
    WRITE(11) co2lat
    CLOSE(11)

END PROGRAM CARBON
