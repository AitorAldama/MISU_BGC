PROGRAM CARBON

    ! ************************************************************************
    !
    !
    !         CARBON ACCUMULATION PER LATITUDE BAND
    !
    !
    ! ************************************************************************

    ! Version 1.0: written by Aitor Aldama-Campino, Autumn 2018

        
    USE netcdf

    IMPLICIT NONE

    INTEGER, PARAMETER     :: IX = 1440, JYD = 1, JYU = 1021, JY = 1021, KZ = 75
    INTEGER, DIMENSION(4)  :: start2D, count2D, start3D, count3D
    INTEGER                :: ncid, varid, ierr, ii, jj, jp, kk, nindex, mm

    REAL, DIMENSION(IX,JYD:JY)    :: dx, dy, kmt
    REAL, DIMENSION(IX,JYD:JY,KZ) :: dz, ctoti,ctotf, tmask
    REAL, DIMENSION(120,JYD:JYU)  :: accu
    REAL, PARAMETER           :: dt = 1., tmin = -350., trend = 9.11518581e-03

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

    start3D = [ 1, JYD,1,1]
    count3D = [IX, JY-JYD+1,KZ,1]


    ! Dx and Dz / Bathimetry
    ierr=NF90_OPEN(trim(filepath)//'mesh_mask.nc',NF90_NOWRITE,ncid)
    IF(ierr.ne.0) STOP 'No mesh mask file'
    ierr=NF90_INQ_VARID(ncid,'e1t',varid) ! the main data fields
    IF(ierr.ne.0) STOP 7
    ierr=NF90_GET_VAR(ncid,varid,dx,start2d,count2d)
    IF(ierr.ne.0) STOP 8
    ierr=NF90_INQ_VARID(ncid,'e2t',varid) ! the main data fields
    IF(ierr.ne.0) STOP 7
    ierr=NF90_GET_VAR(ncid,varid,dy,start2d,count2d)
    IF(ierr.ne.0) STOP 8
    ierr=NF90_INQ_VARID(ncid,'e3t',varid) ! the main data fields
    IF(ierr.ne.0) STOP 9
    ierr=NF90_GET_VAR(ncid,varid,dz,start3d,count3d)
    IF(ierr.ne.0) STOP 10
    ierr=NF90_INQ_VARID(ncid,'tmask',varid) ! the main data fields
    IF(ierr.ne.0) STOP 11
    ierr=NF90_GET_VAR(ncid,varid,tmask,start3d,count3d)
    IF(ierr.ne.0) STOP 12
    ierr=NF90_INQ_VARID(ncid,'bathy',varid) ! the main data fields
    IF(ierr.ne.0) STOP 13
    ierr=NF90_GET_VAR(ncid,varid,kmt,start2d,count2d)
    ierr=NF90_CLOSE(ncid)


    PRINT*, '- Reading mesh and regions done'

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

    start3D = [ 1, JYD,1,1]
    count3D = [IX, JY-JYD+1,KZ,1]

    accu = 0.
    DO nindex = 1,SIZE(label1)

        IF (nindex == 1) THEN
            PRINT*, '- Time step ::', nindex

            ! Tracer file
            ! Read the DIC
            ierr=NF90_OPEN(trim(filepath)//"/ORCA025-ROAM_1m_" &
                //TRIM(label1(nindex))//"_ptrc_T_"//TRIM(label2(nindex))//'-'//TRIM(label2(nindex))//".nc",NF90_NOWRITE,ncid)
            IF(ierr.ne.0) STOP 'No C1 file'
            ierr=NF90_INQ_VARID(ncid,'DIC',varid) ! the main data fields
            IF(ierr.ne.0) STOP 16
            ierr=NF90_GET_VAR(ncid,varid,ctoti,start3d,count3d)
            IF(ierr.ne.0) STOP 17
            ierr=NF90_CLOSE(ncid)
            WHERE (ctoti > 1e6) ctoti = 0.

            CALL SYSTEM("rm /Users/aitoraldama/Desktop/dumBC1/C1.nc")
        END IF

        PRINT*, '- Time step ::', nindex

        ! Tracer file
        ! Read the DIC
        ierr=NF90_OPEN(trim(filepath)//"/ORCA025-ROAM_1m_" &
            //TRIM(label1(nindex))//"_ptrc_T_"//TRIM(label2(nindex))//'-'//TRIM(label2(nindex))//".nc",NF90_NOWRITE,ncid)
        IF(ierr.ne.0) STOP 'No C1 file'
        ierr=NF90_INQ_VARID(ncid,'DIC',varid) ! the main data fields
        IF(ierr.ne.0) STOP 16
        ierr=NF90_GET_VAR(ncid,varid,ctotf,start3d,count3d)
        IF(ierr.ne.0) STOP 17
        ierr=NF90_CLOSE(ncid)
        WHERE (ctotf > 1e6) ctotf = 0.

        DO jj = 1, JY
            DO ii = 1, IX
                DO kk = 1, KZ
                    accu(nindex,jj) = accu(nindex,jj) + dz(ii,jj,kk)*dx(ii,jj)*dy(ii,jj)*&
                        &(ctotf(ii,jj,kk)-ctoti(ii,jj,kk))*tmask(ii,jj,kk)
                END DO
            END DO
        END DO


    END DO

    ! Save file
    OPEN(UNIT=11,FILE=TRIM(filepath)//"result_accu.bin", STATUS='REPLACE', FORM='UNFORMATTED')
    WRITE(11) accu
    CLOSE(11)


END PROGRAM CARBON
