module giss_orbpar_subs
    
contains
    
subroutine GISS_orbpars(year_type, year, eccen, obliq_deg, perih_deg, precc)

    implicit none
    
    ! past ages are negative, e.g. 21 ka = 21,000 cal yr BP = -21000, and 1950 CE = 0 cal yr BP = 0 here
    ! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)
    
    character(2), intent(in)    :: year_type    ! year type (AD/BC, CE/BCE, BP (1950 CE))
    real(8), intent(in)         :: year        ! Year (negative, e.g. 1900 CE = -50 BP)
    real(8), intent(out)        :: eccen        ! eccentricity of orbital ellipse 
    real(8), intent(out)        :: obliq_deg    ! obliquity (degrees)
    real(8), intent(out)        :: perih_deg    ! longitude of perihelion (degrees)
    real(8), intent(out)        :: precc        ! climatological precession parameter = eccen * sin(omegvp)
 
    real(8)                     :: obliq        ! obliquity (latitude of Tropic of Cancer) (radians)
    real(8)                     :: omegvp       ! longitude of perihelion = spatial angle from vernal equinox to perihelion (radians)

    real(8)                     :: YearCE       ! YearCE -- for consistency with old code (input to ORBPAR() is YearCE)
    real(8)                     :: YearBP       ! YearBP
    real(8)                     :: pi           ! pi
    
    pi=4.0d0*datan(1.0d0)

    ! NOTE:  Year CE/AD = 0 is assumed to exist, and is equivalent to 1950 BP (-1950)
    ! subroutine orbpar() expects real-valued YearCE, but converts to YearBP for calculations
    select case (year_type)
    case ('CE', 'AD', 'ce', 'ad')
        YearCE = year
        YearBP = year - 1950.0d0
    case ('BP', 'bp')
        YearCE = year + 1950.0d0
        YearBP = year
    case default
        stop 'year_type'
    end select
    
    call ORBPAR(yearCE, eccen, obliq, omegvp)
    
    obliq_deg=obliq/((2.0d0*pi)/360.0d0)
    precc=eccen*dsin(omegvp)
    perih_deg=omegvp/((2.0d0*pi)/360.0d0)

end subroutine GISS_orbpars
    
subroutine orbpar(year, eccen, obliq, omegvp)   ! year should be YearCE

! NOTE:  Year CE/AD = 0 is assumed to not exist
    
! .f90 version of the ORBPAR() subroutine contained in
!  https://data.giss.nasa.gov/ar5/SOLAR/ORBPAR.FOR downloaded 2017-09-04 17:17
!  also retrievable from https://web.archive.org/web/20150920211936/http://data.giss.nasa.gov/ar5/solar.html
! 
!  ORBPAR calculates the three orbital parameters as a function of
!  YEAR.  The source of these calculations is: Andre L. Berger,
!  1978, "Long-Term Variations of Daily Insolation and Quaternary
!  Climatic Changes", JAS, v.35, p.2362.  Also useful is: Andre L.
!  Berger, May 1978, "A Simple Algorithm to Compute Long Term
!  Variations of Daily Insolation", published by the Institut
!  d'Astronomie et de Geophysique, Universite Catholique de Louvain,
!  Louvain-la Neuve, No. 18.
! 
!  Tables and equations refer to the first reference (JAS).  The
!  corresponding table or equation in the second reference is
!  enclosed in parentheses.  The coefficients used in this
!  subroutine are slightly more precise than those used in either
!  of the references.  The generated orbital parameters are precise
!  within plus or minus 1000000 years from present.
! 
!  Input:  YEAR   = years A.D. are positive, B.C. are negative
!  Output: ECCEN  = eccentricity of orbital ellipse
!          OBLIQ  = latitude of Tropic of Cancer in radians
!          OMEGVP = longitude of perihelion 
!                 = spatial angle from vernal equinox to perihelion
!                   in radians with sun as angle vertex
! 
    implicit none
    
    real(8), intent(in)     :: year     ! Year CE 
    real(8), intent(out)    :: eccen, obliq, omegvp
    
    real(8)                 :: table1(3,47),table4(3,19),table5(3,78)
    real(8)                 :: pi, twopi, piz180
    real(8)                 :: ym1950, sumc, arg, obliqd, esinpi, ecospi, pie, fsinfd, psi
    
    integer(4)              :: i
    
! Table 1 (2).  Obliquity relative to mean ecliptic of date: OBLIQD
    data table1 / &
        -2462.2214466d0, 31.609974d0, 251.9025d0, &
         -857.3232075d0, 32.620504d0, 280.8325d0, &
         -629.3231835d0, 24.172203d0, 128.3057d0, &
         -414.2804924d0, 31.983787d0, 292.7252d0, &
         -311.7632587d0, 44.828336d0,  15.3747d0, &
          308.9408604d0, 30.973257d0, 263.7951d0, &
         -162.5533601d0, 43.668246d0, 308.4258d0, &
         -116.1077911d0, 32.246691d0, 240.0099d0, &
          101.1189923d0, 30.599444d0, 222.9725d0, &
          -67.6856209d0, 42.681324d0, 268.7809d0, &
           24.9079067d0, 43.836462d0, 316.7998d0, &
           22.5811241d0, 47.439436d0, 319.6024d0, &
          -21.1648355d0, 63.219948d0, 143.8050d0, &
          -15.6549876d0, 64.230478d0, 172.7351d0, &
           15.3936813d0,  1.010530d0,  28.9300d0, &
           14.6660938d0,  7.437771d0, 123.5968d0, &
          -11.7273029d0, 55.782177d0,  20.2082d0, &
           10.2742696d0,  0.373813d0,  40.8226d0, &
            6.4914588d0, 13.218362d0, 123.4722d0, &
            5.8539148d0, 62.583231d0, 155.6977d0, &
           -5.4872205d0, 63.593761d0, 184.6277d0, &
           -5.4290191d0, 76.438310d0, 267.2772d0, &
            5.1609570d0, 45.815258d0,  55.0196d0, &
            5.0786314d0,  8.448301d0, 152.5268d0, &
           -4.0735782d0, 56.792707d0,  49.1382d0, &
            3.7227167d0, 49.747842d0, 204.6609d0, &
            3.3971932d0, 12.058272d0,  56.5233d0, &
           -2.8347004d0, 75.278220d0, 200.3284d0, &
           -2.6550721d0, 65.241008d0, 201.6651d0, &
           -2.5717867d0, 64.604291d0, 213.5577d0, &
           -2.4712188d0,  1.647247d0,  17.0374d0, &
            2.4625410d0,  7.811584d0, 164.4194d0, &
            2.2464112d0, 12.207832d0,  94.5422d0, &
           -2.0755511d0, 63.856665d0, 131.9124d0, &
           -1.9713669d0, 56.155990d0,  61.0309d0, &
           -1.8813061d0, 77.448840d0, 296.2073d0, &
           -1.8468785d0,  6.801054d0, 135.4894d0, &
            1.8186742d0, 62.209418d0, 114.8750d0, &
            1.7601888d0, 20.656133d0, 247.0691d0, &
           -1.5428851d0, 48.344406d0, 256.6114d0, &
            1.4738838d0, 55.145460d0,  32.1008d0, &
           -1.4593669d0, 69.000539d0, 143.6804d0, &
            1.4192259d0, 11.071350d0,  16.8784d0, &
           -1.1818980d0, 74.291298d0, 160.6835d0, &
            1.1756474d0, 11.047742d0,  27.5932d0, &
           -1.1316126d0,  0.636717d0, 348.1074d0, &
            1.0896928d0, 12.844549d0,  82.6496d0/
    
! Table 4 (1).  Fundamental elements of the ecliptic: ECCEN sin(pi)
    data table4/ &
          .01860798d0,  4.207205d0,  28.620089d0, &
          .01627522d0,  7.346091d0, 193.788772d0, &
         -.01300660d0, 17.857263d0, 308.307024d0, &
          .00988829d0, 17.220546d0, 320.199637d0, &
         -.00336700d0, 16.846733d0, 279.376984d0, &
          .00333077d0,  5.199079d0,  87.195000d0, &
         -.00235400d0, 18.231076d0, 349.129677d0, &
          .00140015d0, 26.216758d0, 128.443387d0, &
          .00100700d0,  6.359169d0, 154.143880d0, &
          .00085700d0, 16.210016d0, 291.269597d0, &
          .00064990d0,  3.065181d0, 114.860583d0, &
          .00059900d0, 16.583829d0, 332.092251d0, &
          .00037800d0, 18.493980d0, 296.414411d0, &
         -.00033700d0,  6.190953d0, 145.769910d0, &
          .00027600d0, 18.867793d0, 337.237063d0, &
          .00018200d0, 17.425567d0, 152.092288d0, &
         -.00017400d0,  6.186001d0, 126.839891d0, &
         -.00012400d0, 18.417441d0, 210.667199d0, &
          .00001250d0,  0.667863d0,  72.108838d0/
    
!  Table 5 (3).  General precession in longitude: psi
    DATA TABLE5/ &
         7391.0225890d0, 31.609974d0, 251.9025d0, &
         2555.1526947d0, 32.620504d0, 280.8325d0, &
         2022.7629188d0, 24.172203d0, 128.3057d0, &
        -1973.6517951d0,  0.636717d0, 348.1074d0, &
         1240.2321818d0, 31.983787d0, 292.7252d0, &
          953.8679112d0,  3.138886d0, 165.1686d0, &
         -931.7537108d0, 30.973257d0, 263.7951d0, &
          872.3795383d0, 44.828336d0,  15.3747d0, &
          606.3544732d0,  0.991874d0,  58.5749d0, &
         -496.0274038d0,  0.373813d0,  40.8226d0, &
          456.9608039d0, 43.668246d0, 308.4258d0, &
          346.9462320d0, 32.246691d0, 240.0099d0, &
         -305.8412902d0, 30.599444d0, 222.9725d0, &
          249.6173246d0,  2.147012d0, 106.5937d0, &
         -199.1027200d0, 10.511172d0, 114.5182d0, &
          191.0560889d0, 42.681324d0, 268.7809d0, &
         -175.2936572d0, 13.650058d0, 279.6869d0, &
          165.9068833d0,  0.986922d0,  39.6448d0, &
          161.1285917d0,  9.874455d0, 126.4108d0, &
          139.7878093d0, 13.013341d0, 291.5795d0, &
         -133.5228399d0,  0.262904d0, 307.2848d0, &
          117.0673811d0,  0.004952d0,  18.9300d0, &
          104.6907281d0,  1.142024d0, 273.7596d0, &
           95.3227476d0, 63.219948d0, 143.8050d0, &
           86.7824524d0,  0.205021d0, 191.8927d0, &
           86.0857729d0,  2.151964d0, 125.5237d0, &
           70.5893698d0, 64.230478d0, 172.7351d0, &
          -69.9719343d0, 43.836462d0, 316.7998d0, &
          -62.5817473d0, 47.439436d0, 319.6024d0, &
           61.5450059d0,  1.384343d0,  69.7526d0, &
          -57.9364011d0,  7.437771d0, 123.5968d0, &
           57.1899832d0, 18.829299d0, 217.6432d0, &
          -57.0236109d0,  9.500642d0,  85.5882d0, &
          -54.2119253d0,  0.431696d0, 156.2147d0, &
           53.2834147d0,  1.160090d0,  66.9489d0, &
           52.1223575d0, 55.782177d0,  20.2082d0, &
          -49.0059908d0, 12.639528d0, 250.7568d0, &
          -48.3118757d0,  1.155138d0,  48.0188d0, &
          -45.4191685d0,  0.168216d0,   8.3739d0, &
          -42.2357920d0,  1.647247d0,  17.0374d0, &
          -34.7971099d0, 10.884985d0, 155.3409d0, &
           34.4623613d0,  5.610937d0,  94.1709d0, &
          -33.8356643d0, 12.658184d0, 221.1120d0, &
           33.6689362d0,  1.010530d0,  28.9300d0, &
          -31.2521586d0,  1.983748d0, 117.1498d0, &
          -30.8798701d0, 14.023871d0, 320.5095d0, &
           28.4640769d0,  0.560178d0, 262.3602d0, &
          -27.1960802d0,  1.273434d0, 336.2148d0, &
           27.0860736d0, 12.021467d0, 233.0046d0, &
          -26.3437456d0, 62.583231d0, 155.6977d0, &
           24.7253740d0, 63.593761d0, 184.6277d0, &
           24.6732126d0, 76.438310d0, 267.2772d0, &
           24.4272733d0,  4.280910d0,  78.9281d0, &
           24.0127327d0, 13.218362d0, 123.4722d0, &
           21.7150294d0, 17.818769d0, 188.7132d0, &
          -21.5375347d0,  8.359495d0, 180.1364d0, &
           18.1148363d0, 56.792707d0,  49.1382d0, &
          -16.9603104d0,  8.448301d0, 152.5268d0, &
          -16.1765215d0,  1.978796d0,  98.2198d0, &
           15.5567653d0,  8.863925d0,  97.4808d0, &
           15.4846529d0,  0.186365d0, 221.5376d0, &
           15.2150632d0,  8.996212d0, 168.2438d0, &
           14.5047426d0,  6.771027d0, 161.1199d0, &
          -14.3873316d0, 45.815258d0,  55.0196d0, &
           13.1351419d0, 12.002811d0, 262.6495d0, &
           12.8776311d0, 75.278220d0, 200.3284d0, &
           11.9867234d0, 65.241008d0, 201.6651d0, &
           11.9385578d0, 18.870667d0, 294.6547d0, &
           11.7030822d0, 22.009553d0,  99.8233d0, &
           11.6018181d0, 64.604291d0, 213.5577d0, &
          -11.2617293d0, 11.498094d0, 154.1631d0, &
          -10.4664199d0,  0.578834d0, 232.7153d0, &
           10.4333970d0,  9.237738d0, 138.3034d0, &
          -10.2377466d0, 49.747842d0, 204.6609d0, &
           10.1934446d0,  2.147012d0, 106.5938d0, &
          -10.1280191d0,  1.196895d0, 250.4676d0, &
           10.0289441d0,  2.133898d0, 332.3345d0, &
          -10.0034259d0,  0.173168d0,  27.3039d0/
!
    pi = 4.0d0*datan(1.0d0)
    twopi = 2.0d0*pi
    piz180 = twopi/360.0d0
    
    ! NOTE:  change Year CE to Year BP for calculations
    ym1950 = year-1950.0d0

    !  Obliquity from Table 1 (2):
    !    OBLIQ# = 23.320556 (degrees)             Equation 5.5 (15)
    !    OBLIQD = OBLIQ# + sum[A cos(ft+delta)]   Equation 1 (5)
    
    sumc = 0.
    do  i=1,47
        arg = piz180*(ym1950*table1(2,i)/3600.0d0+table1(3,i))
        sumc = sumc + table1(1,i)*dcos(arg)
    end do
    obliqd = 23.320556d0 + sumc/3600.0d0
    obliq = obliqd*piz180

    !  Eccentricity from Table 4 (1):
    !    ECCEN sin(pi) = sum[M sin(gt+beta)]           Equation 4 (1)
    !    ECCEN cos(pi) = sum[M cos(gt+beta)]           Equation 4 (1)
    !    ECCEN = ECCEN sqrt[sin(pi)^2 + cos(pi)^2]
    
    esinpi = 0.d0
    ecospi = 0.d0
    do i=1,19
        arg = piz180*(ym1950*table4(2,i)/3600.0d0+table4(3,i))
        esinpi = esinpi + table4(1,i)*dsin(arg)
        ecospi = ecospi + table4(1,i)*dcos(arg)
    end do
    eccen = sqrt(esinpi*esinpi+ecospi*ecospi)

    !  Perihelion from Equation 4,6,7 (9) and Table 4,5 (1,3):
    !    PSI# = 50.439273 (seconds of degree)         Equation 7.5 (16)
    !    ZETA =  3.392506 (degrees)                   Equation 7.5 (17)
    !    PSI = PSI# t + ZETA + sum[F sin(ft+delta)]   Equation 7 (9)
    !    PIE = atan[ECCEN sin(pi) / ECCEN cos(pi)]
    !    OMEGVP = PIE + PSI + 3.14159                 Equation 6 (4.5)
  
    pie = atan2(esinpi,ecospi)
    fsinfd = 0.d0
    do i=1,78
        arg = piz180*(ym1950*table5(2,i)/3600.+table5(3,i))
        fsinfd = fsinfd + table5(1,i)*dsin(arg)
    end do
    psi = piz180*(3.392506d0+(ym1950*50.439273d0+fsinfd)/3600.d0)
    omegvp = modulo(pie+psi+0.5d0*twopi, twopi)

end subroutine orbpar

end module giss_orbpar_subs