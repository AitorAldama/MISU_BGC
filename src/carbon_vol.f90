PROGRAM CARBON

    ! ************************************************************************
    !
    !
    !         VOLUME DISTRIBUTION IN CARBON-LATITUDE SPACE
    !
    !
    ! ************************************************************************

    ! Version 1.0: written by Aitor Aldama-Campino, Autumn 2018

    USE netcdf

    IMPLICIT NONE

    INTEGER, PARAMETER     :: IX = 1440, JYD = 1, JYU = 1021, JY = 1021, KZ = 75
    INTEGER, DIMENSION(4)  :: start2D, count2D, start3D, count3D
    INTEGER                :: ncid, varid, ierr, ii, jj, jp, kk, nindex, mm

    REAL, DIMENSION(IX,JYD:JY)    :: v2,v1,sou, inp, atl, dx, dy
    REAL, DIMENSION(IX,JYD:JY,KZ) :: dz, vflux, ctot, dx3, vmask, csat
    REAL, DIMENSION(JYD:JY,3251,0:3,2)      :: vol

    CHARACTER (LEN=400)                :: filepath, commandline
    CHARACTER (LEN=20), DIMENSION(120) :: label1, label2, label3
    REAL, PARAMETER           :: dt = 1., tmin = -350.

    ! ***************************************************************

    ! Path to the data input directory
    filepath = '/User/Input/'

    ! READING REGIONS AND MESH
    ! ===============================================================

    start2D = [ 1,JYD,1,1]
    count2D = [IX,JY-JYD+1,1,1]

    start3D = [ 1, JYD,1,1]
    count3D = [IX, JY-JYD+1,KZ,1]

    ! Regions

    ierr=NF90_OPEN(trim(filepath)//'regions.nc',NF90_NOWRITE,ncid)
    IF(ierr.ne.0) STOP 'No region file'
    ierr=NF90_INQ_VARID(ncid,'vmasksou',varid) ! the main data fields
    IF(ierr.ne.0) STOP 1
    ierr=NF90_GET_VAR(ncid,varid,sou,start2d,count2d)
    IF(ierr.ne.0) STOP 2
    ierr=NF90_INQ_VARID(ncid,'vmaskinp',varid) ! the main data fields
    IF(ierr.ne.0) STOP 3
    ierr=NF90_GET_VAR(ncid,varid,inp,start2d,count2d)
    IF(ierr.ne.0) STOP 4
    ierr=NF90_INQ_VARID(ncid,'vmaskatl',varid) ! the main data fields
    IF(ierr.ne.0) STOP 5
    ierr=NF90_GET_VAR(ncid,varid,atl,start2d,count2d)
    IF(ierr.ne.0) STOP 6
    ierr=NF90_CLOSE(ncid)

    ! Dx and Dz

    ierr=NF90_OPEN(trim(filepath)//'mesh_mask.nc',NF90_NOWRITE,ncid)
    IF(ierr.ne.0) STOP 'No mesh mask file'
    ierr=NF90_INQ_VARID(ncid,'e1v',varid) ! the main data fields
    IF(ierr.ne.0) STOP 7
    ierr=NF90_GET_VAR(ncid,varid,dx,start2d,count2d)
    IF(ierr.ne.0) STOP 8
    ierr=NF90_INQ_VARID(ncid,'e2v',varid) ! the main data fields
    IF(ierr.ne.0) STOP 7
    ierr=NF90_GET_VAR(ncid,varid,dy,start2d,count2d)
    IF(ierr.ne.0) STOP 8
    ierr=NF90_INQ_VARID(ncid,'e3v',varid) ! the main data fields
    IF(ierr.ne.0) STOP 9
    ierr=NF90_GET_VAR(ncid,varid,dz,start3d,count3d)
    IF(ierr.ne.0) STOP 10
    ierr=NF90_INQ_VARID(ncid,'vmask',varid) ! the main data fields
    IF(ierr.ne.0) STOP 11
    ierr=NF90_GET_VAR(ncid,varid,vmask,start3d,count3d)
    IF(ierr.ne.0) STOP 12
    ierr=NF90_CLOSE(ncid)

    PRINT*, '- Reading mesh and regions done'


    ! ***************************************************************

    ! BIG LOOP
    ! ===============================================================

    start3D = [ 1, JYD,1,1]
    count3D = [IX, JY-JYD+1,KZ,1]

    vol = 0.

    ! Ctot
    ierr=NF90_OPEN(TRIM(filepath)//'ORCA025_clim_ptrc_T.nc',NF90_NOWRITE,ncid)
    IF(ierr.ne.0) STOP 15
    ierr=NF90_INQ_VARID(ncid,'DIC',varid) ! the main data fields
    IF(ierr.ne.0) STOP 16
    ierr=NF90_GET_VAR(ncid,varid,ctot,start3d,count3d)
    IF(ierr.ne.0) STOP 17
    ierr=NF90_CLOSE(ncid)
    WHERE (ctot > 1e6) ctot = 0.

    ! Csat
    ierr=NF90_INQ_VARID(ncid,'CSAT_TOT',varid) ! the main data fields
    IF(ierr.ne.0) STOP 19
    ierr=NF90_GET_VAR(ncid,varid,csat,start3d,count3d)
    IF(ierr.ne.0) STOP 20
    ierr=NF90_CLOSE(ncid)
    WHERE (csat > 1e6) csat = 0.



    DO ii = 1,IX
        DO jj = JYD,JYU
            jp = jj +1
            if (jj == 1021) jp = 1021
            IF (vmask(ii,jj,1) == 1) THEN
                DO kk = 1,75

                    ! C total
                    mm = NINT((0.5*(ctot(ii,jj,kk)+ctot(ii,jp,kk))-tmin)/dt) + 1

                    mm = MIN(mm, 3251)
                    mm = MAX(mm,    1)

                    vol(jj,mm,0,1) = vol(jj,mm,0,1) + dx(ii,jj)*dy(ii,jj)*dz(ii,jj,kk)*vmask(ii,jj,kk)


                    IF (sou(ii,jj)==1) THEN
                        vol(jj,mm,1,1) = vol(jj,mm,1,1) + dx(ii,jj)*dy(ii,jj)*dz(ii,jj,kk)*vmask(ii,jj,kk)
                    ELSEIF (atl(ii,jj)==1) THEN
                        vol(jj,mm,2,1) = vol(jj,mm,2,1) + dx(ii,jj)*dy(ii,jj)*dz(ii,jj,kk)*vmask(ii,jj,kk)
                    ELSEIF (inp(ii,jj)==1) THEN
                        vol(jj,mm,3,1) = vol(jj,mm,3,1) + dx(ii,jj)*dy(ii,jj)*dz(ii,jj,kk)*vmask(ii,jj,kk)
                    ENDIF


                    ! C saturation
                    mm = NINT((0.5*(csat(ii,jj,kk)+csat(ii,jp,kk))-tmin)/dt) + 1

                    mm = MIN(mm, 3251)
                    mm = MAX(mm,    1)

                    vol(jj,mm,0,2) = vol(jj,mm,0,2) + dx(ii,jj)*dy(ii,jj)*dz(ii,jj,kk)*vmask(ii,jj,kk)

                    IF (sou(ii,jj)==1) THEN
                        vol(jj,mm,1,2) = vol(jj,mm,1,2) + dx(ii,jj)*dy(ii,jj)*dz(ii,jj,kk)*vmask(ii,jj,kk)
                    ELSEIF (atl(ii,jj)==1) THEN
                        vol(jj,mm,2,2) = vol(jj,mm,2,2) + dx(ii,jj)*dy(ii,jj)*dz(ii,jj,kk)*vmask(ii,jj,kk)
                    ELSEIF (inp(ii,jj)==1) THEN
                        vol(jj,mm,3,2) = vol(jj,mm,3,2) + dx(ii,jj)*dy(ii,jj)*dz(ii,jj,kk)*vmask(ii,jj,kk)

                    ENDIF

                END DO
            END IF

        END DO
    END DO


    ! Save file
    OPEN(UNIT=11,FILE=TRIM(filepath)//"result_vol.bin",STATUS='REPLACE', FORM='UNFORMATTED')
    WRITE(11) vol
    CLOSE(11)


END PROGRAM CARBON
