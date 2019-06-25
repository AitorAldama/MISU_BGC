PROGRAM CARBON

    ! ************************************************************************
    !
    !
    !         CARBON-LATITUDE STREAMFUNCTION
    !
    !
    ! ************************************************************************

    ! Version 1.0: written by Aitor Aldama-Campino, Autumn 2018

    USE netcdf

    IMPLICIT NONE

    INTEGER, PARAMETER     :: IX = 1440, JYD = 1, JYU = 1021, JY = 1021, KZ = 75
    INTEGER, DIMENSION(4)  :: start2D, count2D, start3D, count3D
    INTEGER                :: ncid, varid, ierr, ii, jj, jp, kk, nindex, mm

    REAL, DIMENSION(IX,JYD:JY)    :: sou, inp, atl, dx, vmask, kmt, kmv
    REAL, DIMENSION(IX,JYD:JY,KZ) :: dz, vflux, ctot, csat, ccarb, csoft, cdis, cbio
    REAL, DIMENSION(JYD:JYU,3251,0:3,4) :: fluxb
    REAL, DIMENSION(JYD:JYU,0:3,4)      :: trp
    REAL, PARAMETER           :: dt = 1., tmin = -350.

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
    ierr=NF90_INQ_VARID(ncid,'e3v',varid) ! the main data fields
    IF(ierr.ne.0) STOP 9
    ierr=NF90_GET_VAR(ncid,varid,dz,start3d,count3d)
    IF(ierr.ne.0) STOP 10
    ierr=NF90_INQ_VARID(ncid,'vmask',varid) ! the main data fields
    IF(ierr.ne.0) STOP 11
    ierr=NF90_GET_VAR(ncid,varid,vmask,start2d,count2d)
    IF(ierr.ne.0) STOP 12
    ierr=NF90_INQ_VARID(ncid,'bathy',varid) ! the main data fields
    IF(ierr.ne.0) STOP 13
    ierr=NF90_GET_VAR(ncid,varid,kmt,start2d,count2d)
    ierr=NF90_CLOSE(ncid)

    DO jj = JYD,JYU
        DO ii = 1,IX
            jp = jj + 1
            IF (jp == JY+1) jp = JY

            kmv(ii,jj) = MIN(kmt(ii,jj), kmt(ii, jp))
        END DO
    END DO

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

    trp   = 0.
    fluxb = 0.

    DO nindex = 1, SIZE(label1)

        PRINT*, '- Time step ::', nindex

        ! Velocity file
        ierr=NF90_OPEN(trim(filepath)//"ORCA025-ROAM_1m_" &
            //TRIM(label1(nindex))//"_grid_V_"//TRIM(label2(nindex))//'-'//TRIM(label2(nindex))//".nc",NF90_NOWRITE,ncid)
        IF(ierr.ne.0) STOP 'No V file'
        ierr=NF90_INQ_VARID(ncid,'vomecrty',varid) ! the main data fields
        IF(ierr.ne.0) STOP 14
        ierr=NF90_GET_VAR(ncid,varid,vflux,start3d,count3d)
        IF(ierr.ne.0) STOP 15
        ierr=NF90_CLOSE(ncid)
        WHERE (vflux > 1e6) vflux = 0.

        ! Tracer file
        ierr=NF90_OPEN(trim(filepath)//"ORCA025-ROAM_1m_" &
                //TRIM(label1(nindex))//"_ptrc_T_"//TRIM(label2(nindex))//'-'//TRIM(label2(nindex))//".nc",NF90_NOWRITE,ncid)
        IF(ierr.ne.0) STOP 'No C1 file'
        ierr=NF90_INQ_VARID(ncid,'DIC',varid) ! the main data fields
        IF(ierr.ne.0) STOP 16
        ierr=NF90_GET_VAR(ncid,varid,ctot,start3d,count3d)
        IF(ierr.ne.0) STOP 17
        ierr=NF90_CLOSE(ncid)
        WHERE (ctot > 1e6) ctot = 0.

        ierr=NF90_OPEN(trim(filepath)//"ORCA25_"&
            //TRIM(label2(nindex))//"_DIC_analysis.nc",NF90_NOWRITE,ncid)
        IF(ierr.ne.0) STOP 'No C2 file'
        ierr=NF90_INQ_VARID(ncid,'CSAT_TOT',varid) ! the main data fields
        IF(ierr.ne.0) STOP 18
        ierr=NF90_GET_VAR(ncid,varid,csat,start3d,count3d)
        IF(ierr.ne.0) STOP 19
        ierr=NF90_INQ_VARID(ncid,'CSOFT',varid) ! the main data fields
        IF(ierr.ne.0) STOP 20
        ierr=NF90_GET_VAR(ncid,varid,csoft,start3d,count3d)
        IF(ierr.ne.0) STOP 21
        ierr=NF90_INQ_VARID(ncid,'CCARB',varid) ! the main data fields
        IF(ierr.ne.0) STOP 22
        ierr=NF90_GET_VAR(ncid,varid,ccarb,start3d,count3d)
        IF(ierr.ne.0) STOP 23
        ierr=NF90_CLOSE(ncid)

        WHERE (csat > 1e6) csat = 0.
        WHERE (csoft > 1e6) csoft = 0.
        WHERE (ccarb > 1e6) ccarb = 0.

        cbio = csoft + ccarb
        cdis = ctot - (cbio+csat)


        ! Computation of fluxes in carbon-latitude space and meridional transports
        DO ii = 1,IX
            DO jj = JYD,JYU

                IF (vmask(ii,jj) == 1) THEN
                    DO kk = 1,INT(kmv(ii,jj))

                        ! C total
                        mm = NINT((0.5*(ctot(ii,jj,kk)+ctot(ii,jj+1,kk))-tmin)/dt) + 1

                        mm = MIN(mm, 3251)
                        mm = MAX(mm,    1)

                        fluxb(jj,mm,0,1) = fluxb(jj,mm,0,1) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        trp(jj,0,1) = trp(jj,0,1) + 0.5*(ctot(ii,jj,kk)+ctot(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6

                        IF (sou(ii,jj)==1) THEN
                            fluxb(jj,mm,1,1) = fluxb(jj,mm,1,1) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,1,1) = trp(jj,1,1) + 0.5*(ctot(ii,jj,kk)+ctot(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        ELSEIF (atl(ii,jj)==1) THEN
                            fluxb(jj,mm,2,1) = fluxb(jj,mm,2,1) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,2,1) = trp(jj,2,1) + 0.5*(ctot(ii,jj,kk)+ctot(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        ELSEIF (inp(ii,jj)==1) THEN
                            fluxb(jj,mm,3,1) = fluxb(jj,mm,3,1) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,3,1) = trp(jj,3,1) + 0.5*(ctot(ii,jj,kk)+ctot(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        END IF

                        ! C saturation
                        mm = NINT((0.5*(csat(ii,jj,kk)+csat(ii,jj+1,kk))-tmin)/dt) + 1

                        mm = MIN(mm, 3251)
                        mm = MAX(mm,    1)

                        fluxb(jj,mm,0,2) = fluxb(jj,mm,0,2) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        trp(jj,0,2) = trp(jj,0,2) + 0.5*(csat(ii,jj,kk)+csat(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        IF (sou(ii,jj)==1) THEN
                            fluxb(jj,mm,1,2) = fluxb(jj,mm,1,2) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,1,2) = trp(jj,1,2) + 0.5*(csat(ii,jj,kk)+csat(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        ELSEIF (atl(ii,jj)==1) THEN
                            fluxb(jj,mm,2,2) = fluxb(jj,mm,2,2) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,2,2) = trp(jj,2,2) + 0.5*(csat(ii,jj,kk)+csat(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        ELSEIF (inp(ii,jj)==1) THEN
                            fluxb(jj,mm,3,2) = fluxb(jj,mm,3,2) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,3,2) = trp(jj,3,2) + 0.5*(csat(ii,jj,kk)+csat(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        END IF

                        ! C biology
                        mm = NINT((0.5*(cbio(ii,jj,kk)+cbio(ii,jj+1,kk))-tmin)/dt) + 1

                        mm = MIN(mm, 3251)
                        mm = MAX(mm,    1)

                        fluxb(jj,mm,0,3) = fluxb(jj,mm,0,3) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        trp(jj,0,3) = trp(jj,0,3) + 0.5*(cbio(ii,jj,kk)+cbio(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        IF (sou(ii,jj)==1) THEN
                            fluxb(jj,mm,1,3) = fluxb(jj,mm,1,3) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,1,3) = trp(jj,1,3) + 0.5*(cbio(ii,jj,kk)+cbio(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        ELSEIF (atl(ii,jj)==1) THEN
                            fluxb(jj,mm,2,3) = fluxb(jj,mm,2,3) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,2,3) = trp(jj,2,3) + 0.5*(cbio(ii,jj,kk)+cbio(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        ELSEIF (inp(ii,jj)==1) THEN
                            fluxb(jj,mm,3,3) = fluxb(jj,mm,3,3) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,3,3) = trp(jj,3,3) + 0.5*(cbio(ii,jj,kk)+cbio(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        END IF

                        ! C disequilibrium
                        mm = NINT((0.5*(cdis(ii,jj,kk)+cdis(ii,jj+1,kk)) -tmin)/dt) + 1

                        mm = MIN(mm, 3251)
                        mm = MAX(mm,    1)

                        fluxb(jj,mm,0,4) = fluxb(jj,mm,0,4) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        trp(jj,0,4) = trp(jj,0,4) + 0.5*(cdis(ii,jj,kk)+cdis(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        IF (sou(ii,jj)==1) THEN
                            fluxb(jj,mm,1,4) = fluxb(jj,mm,1,4) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,1,4) = trp(jj,1,4) + 0.5*(cdis(ii,jj,kk)+cdis(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        ELSEIF (atl(ii,jj)==1) THEN
                            fluxb(jj,mm,2,4) = fluxb(jj,mm,2,4) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,2,4) = trp(jj,2,4) + 0.5*(cdis(ii,jj,kk)+cdis(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        ELSEIF (inp(ii,jj)==1) THEN
                            fluxb(jj,mm,3,4) = fluxb(jj,mm,3,4) + vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                            trp(jj,3,4) = trp(jj,3,4) + 0.5*(cdis(ii,jj,kk)+cdis(ii,jj+1,kk))&
                                        *vflux(ii,jj,kk)*dx(ii,jj)*dz(ii,jj,kk)*1e-6
                        END IF

                    END DO
                END IF

            END DO
        END DO
    END DO

    ! Compute time mean
    fluxb = fluxb/SIZE(label1)
    trp   = trp/SIZE(label1)

    ! Save file
    OPEN(UNIT=11,FILE=TRIM(filepath)//"result_fluxestrp.bin", STATUS='REPLACE', FORM='UNFORMATTED')
    WRITE(11) fluxb
    WRITE(11) trp
    CLOSE(11)


END PROGRAM CARBON
