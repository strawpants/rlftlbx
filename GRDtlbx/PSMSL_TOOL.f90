!!Coded by Roelof Rietbroek, Wed Oct 14 14:17:54 2009
!!GeoForschungZentrum Potsdam,
!!Telegrafenberg A17,
!!D-17743, Potsdam, Germany
!!email: roelof@gfz-potsdam.de

!set of tools to read the psmldatabase
!!Updated by Roelof Rietbroek, Thu Jul 29 08:59:41 2010
!!major rewrite associated with the update of the PSMSL databse format of april 2010

!!Updated by Roelof Rietbroek, Tue Oct 26 16:06:58 2010
!! updated psmsl_tindex (round times to months/years instead of truncation)


module psmsl_tool
implicit none
integer::max_stat
parameter(max_stat=3000)
!   This derived type contains all the data for 1 particular station:
type  PSMSL_stat
   !station
   character(5)::sid
   character(3)::scode
   character(40)::name
   character(60)::coast
   character(1)::qcflag ! quality check flag Y means check documentation
   
   !longitude and latitude of the station
   double precision::lon,lat

   !start and end time
   double precision::t0,tn


   !character arrays holding comments
   character(80),pointer,dimension(:)::comm=>null()
   integer::ncom

   !data arrays
   double precision,pointer,dimension(:)::time=>null()
   integer,pointer,dimension(:)::H=>null()
   integer,pointer,dimension(:)::misd=>null()
   integer,pointer,dimension(:)::flag=>null()
   character(1),pointer,dimension(:)::yrflag=>null()
   integer::ndat





     ! character(1)::mrs
!      integer::iyrlr
!      !authority and frequency code
!      character(2)::acode,fcode
!      !country code, station code, gloss code
!      character(3)::ccode,scode,gloss


     !start end end time of the series
!     integer::t0,tn     
     
!      !character arrays holding comments
!      character(80),pointer,dimension(:)::a_d=>null()
!      character(80),pointer,dimension(:)::s_d=>null()
!      character(80),pointer,dimension(:)::c_d=>null()
!      integer::comc,coms,coma
     
!      integer::na ! number of valid years of data
     
!      !monthly data
!      integer,pointer,dimension(:)::mmet=>null()
!      !yearly data
!      integer,pointer,dimension(:)::amet=>null()
!      integer,pointer,dimension(:)::arlr=>null()
     
!      integer,pointer,dimension(:)::yrs
!      !missing days per month
!      integer,pointer,dimension(:)::mmis=>null()
     
!      !missing days per year
!      integer,pointer,dimension(:)::amis=>null()
     
  end type PSMSL_stat

  
  type PSMSL_struc
     character(120)::dir=''
     integer::nstat=0 ! index of last station read in
     type(PSMSL_stat)::stat(max_stat) !default may take up to max_stat stations
     integer::ncount=204
     character(60)::coast(204)
     character(8)::ext! data filename extension
     logical::yearly=.false.
     !timestep
     double precision::step
  end type PSMSL_struc

 contains


   subroutine psmsl_rdstat(psmsl,pos,readdoc)
     use FORTtlbx
     type(PSMSL_struc),target,intent(inout)::psmsl
     integer,intent(in)::pos !specify the position of the station in the list
     logical,intent(in)::readdoc
     type(PSMSL_stat),pointer::stat=>null() ! pointer for readability

     character(200)::file
     integer::unit,last,chunk,mem
     
     !allocate data vectors per chunk
     chunk=100

     !pointer to the relevant stat structure (shortcut for readability)
     stat=>psmsl%stat(pos)
     unit=13
     

     if(readdoc)then! read documentation
        file=trim(psmsl%dir)//'docu/'//trim(adjustl(stat%sid))//'.txt'
        
        !read documentation from file
        stat%ncom=0
        last=0
        open(unit=unit,form='formatted',file=file)
        do while(last .eq. 0)
           stat%ncom=stat%ncom+1
           call realloc_ptr(stat%comm,1)
           read(unit=unit,fmt='(A80)',iostat=last)stat%comm(stat%ncom)
        end do
        close(unit)
        stat%ncom=stat%ncom-1
     else! read data file
        file=trim(psmsl%dir)//'data/'//trim(adjustl(stat%sid))//psmsl%ext
        !read data 
        stat%ndat=0
        mem=chunk
        last=0
        open(unit=unit,form='formatted',file=file)
        if(psmsl%yearly)then ! yearly data
           allocate(stat%H(chunk),stat%yrflag(mem),stat%time(mem),stat%flag(mem))
           do while(last .eq. 0)
              stat%ndat=stat%ndat+1
              
              if(stat%ndat>mem)then ! reallocate if necessary
                 call realloc_ptr(stat%time,chunk)
                 call realloc_ptr(stat%H,chunk)
                 call realloc_ptr(stat%flag,chunk)
                 call realloc_ptr(stat%yrflag,chunk)
                 mem=mem+chunk
              end if
              
              read(unit=unit,fmt='(1x,f4.0,1x,i6,1x,a1,1x,i3)',iostat=last)stat%time(stat%ndat),&
                   stat%H(stat%ndat),stat%yrflag(stat%ndat),stat%flag(stat%ndat)
              
           end do
          ! add half a year to center yearly timtags
           stat%time(1:stat%ndat-1)=stat%time(1:stat%ndat-1)+0.5
        else ! monthly
           allocate(stat%H(chunk),stat%misd(mem),stat%time(mem),stat%flag(mem))
           do while(last .eq. 0)
              stat%ndat=stat%ndat+1
              
              if(stat%ndat>mem)then ! reallocate if necessary
                 call realloc_ptr(stat%time,chunk)
                 call realloc_ptr(stat%H,chunk)
                 call realloc_ptr(stat%flag,chunk)
                 call realloc_ptr(stat%misd,chunk)
                 mem=mem+chunk
           end if
           
           read(unit=unit,fmt='(1x,f10.4,1x,i6,1x,i2,1x,i3)',iostat=last)stat%time(stat%ndat),&
                stat%H(stat%ndat),stat%misd(stat%ndat),stat%flag(stat%ndat)
        end do
     end if

     stat%ndat=stat%ndat-1 !correct the amount of data for the last invalid line
     close(unit)

  end if
     
end subroutine psmsl_rdstat
     
!function to retrieve the index of a time
   !if the year is missing -1 will be returned
    function psmsl_tindex(psmsl,pos,yr)
      implicit none
      type(PSMSL_struc),intent(in)::psmsl
      double precision,intent(in)::yr
      integer,intent(in)::pos
      integer::psmsl_tindex
      
      psmsl_tindex=int((yr-psmsl%stat(pos)%t0)/psmsl%step+0.5)+1

   end function psmsl_tindex

   subroutine psmsl_init(psmsl)
     type(PSMSL_struc),intent(inout)::psmsl
     integer::err,stderr,unit,last,i
     integer::n
     stderr=0
     unit=13

     !coastline codes
     psmsl%coast(1)='     010    ICELAND                                  '
     psmsl%coast(2)='     012    JAN MAYEN                                '
     psmsl%coast(3)='     015    FAEROE ISLANDS                           '
     psmsl%coast(4)='     025    SPITSBERGEN                              '
     psmsl%coast(5)='     030    RUSSIAN FEDERATION (ARCTIC)              '
     psmsl%coast(6)='     040    NORWAY                                   '
     psmsl%coast(7)='     050    SWEDEN                                   '
     psmsl%coast(8)='     060    FINLAND                                  '
     psmsl%coast(9)='     080    BALTIC STATES (FORMER USSR)              '
     psmsl%coast(10)='     110    POLAND                                   '
     psmsl%coast(11)='     120    GERMANY (FORMER DDR) BALTIC              '
     psmsl%coast(12)='     125    GERMANY (FORMER FRG) BALTIC              '
     psmsl%coast(13)='     130    DENMARK                                  '
     psmsl%coast(14)='     140    GERMANY (NORTH SEA)                      '
     psmsl%coast(15)='     150    NETHERLANDS                              '
     psmsl%coast(16)='     160    BELGIUM                                  '
     psmsl%coast(17)='     170    UNITED KINGDOM                           '
     psmsl%coast(18)='     175    IRELAND                                  '
     psmsl%coast(19)='     180    CHANNEL ISLANDS                          '
     psmsl%coast(20)='     190    FRANCE (ATLANTIC)                        '
     psmsl%coast(21)='     200    SPAIN (ATLANTIC)                         '
     psmsl%coast(22)='     210    PORTUGAL                                 '
     psmsl%coast(23)='     215    GIBRALTAR                                '
     psmsl%coast(24)='     220    SPAIN (MEDITERRANEAN)                    '
     psmsl%coast(25)='     225    SPAIN (BALEARIC ISLANDS)                 '
     psmsl%coast(26)='     230    FRANCE (MEDITERRANEAN)                   '
     psmsl%coast(27)='     232    CORSICA                                  '
     psmsl%coast(28)='     233    MONACO                                   '
     psmsl%coast(29)='     240    ITALY (SARDINIA)                         '
     psmsl%coast(30)='     250    ITALY (MEDITERRANEAN)                    '
     psmsl%coast(31)='     260    ITALY (SICILY)                           '
     psmsl%coast(32)='     265    MALTA                                    '
     psmsl%coast(33)='     270    ITALY (ADRIATIC)                         '
     psmsl%coast(34)='     279    SLOVENIA                                 '
     psmsl%coast(35)='     280    CROATIA                                  '
     psmsl%coast(36)='     281    MONTENEGRO                               '
     psmsl%coast(37)='     290    GREECE                                   '
     psmsl%coast(38)='     295    BULGARIA                                 '
     psmsl%coast(39)='     297    ROMANIA                                  '
     psmsl%coast(40)='     298    UKRAINE                                  '
     psmsl%coast(41)='     300    RUSSIAN FEDERATION (BLACK SEA)           '
     psmsl%coast(42)='     305    GEORGIA                                  '
     psmsl%coast(43)='     310    TURKEY                                   '
     psmsl%coast(44)='     315    CYPRUS                                   '
     psmsl%coast(45)='     320    ISRAEL (MEDITERRANEAN)                   '
     psmsl%coast(46)='     330    EGYPT                                    '
     psmsl%coast(47)='     340    SPANISH N. AFRICA                        '
     psmsl%coast(48)='     350    MOROCCO                                  '
     psmsl%coast(49)='     360    PORTUGAL (AZORES)                        '
     psmsl%coast(50)='     365    PORTUGAL (MADEIRA)                       '
     psmsl%coast(51)='     370    SPAIN (CANARY ISLANDS)                   '
     psmsl%coast(52)='     380    CAPE VERDE ISLANDS                       '
     psmsl%coast(53)='     390    SENEGAL                                  '
     psmsl%coast(54)='     396    GUINEA                                   '
     psmsl%coast(55)='     400    SIERRA LEONE                             '
     psmsl%coast(56)='     402    ASCENSION                                '
     psmsl%coast(57)='     405    IVORY COAST                              '
     psmsl%coast(58)='     410    GHANA                                    '
     psmsl%coast(59)='     420    NIGERIA                                  '
     psmsl%coast(60)='     421    FERNANDO POO                             '
     psmsl%coast(61)='     423    SAO TOME AND PRINCIPE                    '
     psmsl%coast(62)='     424    CONGO                                    '
     psmsl%coast(63)='     425    ST. HELENA                               '
     psmsl%coast(64)='     426    ANGOLA                                   '
     psmsl%coast(65)='     427    NAMIBIA                                  '
     psmsl%coast(66)='     430    SOUTH AFRICA                             '
     psmsl%coast(67)='     432    MOZAMBIQUE                               '
     psmsl%coast(68)='     433    CROZET IS.                               '
     psmsl%coast(69)='     434    KERGUELEN ISLAND                         '
     psmsl%coast(70)='     436    SAINT PAUL ISLAND                        '
     psmsl%coast(71)='     438    COMOROS                                  '
     psmsl%coast(72)='     440    MADAGASCAR                               '
     psmsl%coast(73)='     441    ALDABRA                                  '
     psmsl%coast(74)='     442    SEYCHELLES                               '
     psmsl%coast(75)='     450    MAURITIUS                                '
     psmsl%coast(76)='     451    REUNION IS.                              '
     psmsl%coast(77)='     453    CHAGOS ARCHIPELAGO                       '
     psmsl%coast(78)='     454    MALDIVES                                 '
     psmsl%coast(79)='     455    LACCADIVE ISLANDS                        '
     psmsl%coast(80)='     460    TANZANIA                                 '
     psmsl%coast(81)='     470    KENYA                                    '
     psmsl%coast(82)='     475    DJIBOUTI                                 '
     psmsl%coast(83)='     477    SUDAN                                    '
     psmsl%coast(84)='     480    ISRAEL (GULF OF AQABA)                   '
     psmsl%coast(85)='     482    GULF                                     '
     psmsl%coast(86)='     485    YEMEN PDR                                '
     psmsl%coast(87)='     487    MUSCAT & OMAN                            '
     psmsl%coast(88)='     489    IRAN                                     '
     psmsl%coast(89)='     490    PAKISTAN                                 '
     psmsl%coast(90)='     500    INDIA                                    '
     psmsl%coast(91)='     510    BANGLADESH                               '
     psmsl%coast(92)='     520    SRI LANKA                                '
     psmsl%coast(93)='     530    MYANMAR                                  '
     psmsl%coast(94)='     540    ANDAMAN ISLANDS                          '
     psmsl%coast(95)='     545    THAILAND (ANDAMAN SEA)                   '
     psmsl%coast(96)='     550    MALAYSIA                                 '
     psmsl%coast(97)='     555    SINGAPORE                                '
     psmsl%coast(98)='     560    INDONESIA                                '
     psmsl%coast(99)='     563    CHRISTMAS ISLAND                         '
     psmsl%coast(100)='     570    KALIMANTAN                               '
     psmsl%coast(101)='     580    SULAWESI                                 '
     psmsl%coast(102)='     590    MALUKU                                   '
     psmsl%coast(103)='     600    THAILAND (GULF OF THAILAND)              '
     psmsl%coast(104)='     605    VIET NAM                                 '
     psmsl%coast(105)='     609    MACAU                                    '
     psmsl%coast(106)='     610    CHINA                                    '
     psmsl%coast(107)='     611    HONG KONG, CHINA                         '
     psmsl%coast(108)='     612    TAIWAN (FORMOSA)                         '
     psmsl%coast(109)='     620    KOREA (SOUTH)                            '
     psmsl%coast(110)='     625    KOREA (NORTH)-SEA OF JAPAN               '
     psmsl%coast(111)='     630    RUSSIAN FEDERATION (PACIFIC OCEAN)       '
     psmsl%coast(112)='     640    JAPAN (KARAFUTO)                         '
     psmsl%coast(113)='     641    JAPAN (HOKKAIDO)                         '
     psmsl%coast(114)='     642    JAPAN (HONSHU-PACIFIC)                   '
     psmsl%coast(115)='     643    JAPAN (HONSHU-INLAND SEA)                '
     psmsl%coast(116)='     644    JAPAN (SHIKOKU)                          '
     psmsl%coast(117)='     645    JAPAN (KYUSHU)                           '
     psmsl%coast(118)='     646    JAPAN (AMAMI GUNTO)                      '
     psmsl%coast(119)='     647    JAPAN (HONSHU-JAPAN SEA)                 '
     psmsl%coast(120)='     648    JAPAN (OGASAWARA GUNTO)                  '
     psmsl%coast(121)='     649    JAPAN (MINAMI-TORI-SHIMA)                '
     psmsl%coast(122)='     660    PHILIPPINES                              '
     psmsl%coast(123)='     663    SARAWAK                                  '
     psmsl%coast(124)='     665    SABAH                                    '
     psmsl%coast(125)='     670    PAPUA NEW GUINEA                         '
     psmsl%coast(126)='     680    AUSTRALIA                                '
     psmsl%coast(127)='     690    NEW ZEALAND                              '
     psmsl%coast(128)='     700    NORTHERN MARIANA ISLANDS                 '
     psmsl%coast(129)='     701    GUAM                                     '
     psmsl%coast(130)='     710    CAROLINE IS (FED. STATES OF MICRONESIA)  '
     psmsl%coast(131)='     711    PALAU ISLANDS                            '
     psmsl%coast(132)='     715    NAURU                                    '
     psmsl%coast(133)='     720    MARSHALL ISLANDS                         '
     psmsl%coast(134)='     730    KIRIBATI                                 '
     psmsl%coast(135)='     732    TUVALU                                   '
     psmsl%coast(136)='     734    SOLOMON ISLANDS                          '
     psmsl%coast(137)='     740    NEW CALEDONIA                            '
     psmsl%coast(138)='     741    VANUATU                                  '
     psmsl%coast(139)='     742    FIJI                                     '
     psmsl%coast(140)='     744    TONGA                                    '
     psmsl%coast(141)='     745    AMERICAN SAMOA                           '
     psmsl%coast(142)='     746    WESTERN SAMOA                            '
     psmsl%coast(143)='     750    PHOENIX ISLANDS (KIRIBATI)               '
     psmsl%coast(144)='     760    HAWAIIAN ISLANDS                         '
     psmsl%coast(145)='     770    LINE ISLANDS                             '
     psmsl%coast(146)='     775    PENRHYN ISLAND                           '
     psmsl%coast(147)='     780    ILES DE LA SOCIETE                       '
     psmsl%coast(148)='     785    COOK ISLANDS                             '
     psmsl%coast(149)='     790    TUBUAI ISLANDS                           '
     psmsl%coast(150)='     800    TUAMOTU ARCHIPELAGO                      '
     psmsl%coast(151)='     805    MARQUESAS                                '
     psmsl%coast(152)='     808    GAMBIER ISLAND                           '
     psmsl%coast(153)='     810    EASTER ISLAND                            '
     psmsl%coast(154)='     820    USA (ALEUTIAN ISLANDS)                   '
     psmsl%coast(155)='     821    USA (ALASKA)                             '
     psmsl%coast(156)='     822    CANADA (PACIFIC COAST)                   '
     psmsl%coast(157)='     823    USA (PACIFIC COAST)                      '
     psmsl%coast(158)='     830    MEXICO (PACIFIC)                         '
     psmsl%coast(159)='     831    CLIPPERTON ISLAND                        '
     psmsl%coast(160)='     832    GUATEMALA (PACIFIC)                      '
     psmsl%coast(161)='     833    EL SALVADOR                              '
     psmsl%coast(162)='     836    COSTA RICA (PACIFIC)                     '
     psmsl%coast(163)='     840    PANAMA (PACIFIC)                         '
     psmsl%coast(164)='     842    COLOMBIA (PACIFIC)                       '
     psmsl%coast(165)='     845    ECUADOR                                  '
     psmsl%coast(166)='     848    PERU                                     '
     psmsl%coast(167)='     850    CHILE                                    '
     psmsl%coast(168)='     860    ARGENTINA                                '
     psmsl%coast(169)='     863    FALKLAND ISLANDS (MALVINAS)              '
     psmsl%coast(170)='     866    SOUTH GEORGIA                            '
     psmsl%coast(171)='     870    URUGUAY                                  '
     psmsl%coast(172)='     874    BRAZIL                                   '
     psmsl%coast(173)='     876    FRENCH GUIANA                            '
     psmsl%coast(174)='     880    GUYANA                                   '
     psmsl%coast(175)='     890    TRINIDAD & TOBAGO                        '
     psmsl%coast(176)='     893    GRENADA                                  '
     psmsl%coast(177)='     895    ST. VINCENT                              '
     psmsl%coast(178)='     900    VENEZUELA                                '
     psmsl%coast(179)='     902    COLOMBIA (CARIBBEAN)                     '
     psmsl%coast(180)='     904    PANAMA (CARIBBEAN)                       '
     psmsl%coast(181)='     906    COSTA RICA (CARIBBEAN)                   '
     psmsl%coast(182)='     908    HONDURAS                                 '
     psmsl%coast(183)='     910    BARBADOS                                 '
     psmsl%coast(184)='     911    ST. LUCIA                                '
     psmsl%coast(185)='     912    MARTINIQUE                               '
     psmsl%coast(186)='     914    GUADELOUPE                               '
     psmsl%coast(187)='     916    GUATEMALA (CARIBBEAN)                    '
     psmsl%coast(188)='     918    BELIZE                                   '
     psmsl%coast(189)='     920    MEXICO (GULF)                            '
     psmsl%coast(190)='     930    CUBA                                     '
     psmsl%coast(191)='     931    CAYMAN ISLANDS                           '
     psmsl%coast(192)='     932    JAMAICA                                  '
     psmsl%coast(193)='     934    HAITI                                    '
     psmsl%coast(194)='     936    DOMINICAN REPUBLIC                       '
     psmsl%coast(195)='     937    ST. KITTS                                '
     psmsl%coast(196)='     938    PUERTO RICO                              '
     psmsl%coast(197)='     939    VIRGIN ISLANDS                           '
     psmsl%coast(198)='     940    USA (GULF)                               '
     psmsl%coast(199)='     941    BAHAMAS                                  '
     psmsl%coast(200)='     950    BERMUDA                                  '
     psmsl%coast(201)='     960    USA (ATLANTIC)                           '
     psmsl%coast(202)='     970    CANADA (ATLANTIC AND ARCTIC)             '
     psmsl%coast(203)='     980    GREENLAND                                '
     psmsl%coast(204)='     999    ANTARCTICA                               '
     

     !set file extensions for later data files
     if(index(psmsl%dir,'met_monthly') >0 )then
        psmsl%ext='.metdata'
        psmsl%step=1/12.d0
     else if(index(psmsl%dir,'rlr_annual') >0)then
        psmsl%ext='.rlrdata'
        psmsl%yearly=.true.
        psmsl%step=1
     else ! monthly RLR
        psmsl%ext='.rlrdata'
        psmsl%step=1/12.d0
     end if
     
     !open file inventory ( extended version)
     open(unit=unit,file=trim(psmsl%dir)//'filelist.txt_extnd',form='formatted',iostat=err)
     if(err .ne. 0)then
       write(stderr,*)"ERROR:opening file:",trim(psmsl%dir)//'filelist.txt_extnd'
       stop
    end if
    last=0
    n=0
    do while(last .eq. 0)
       n=n+1
       !       read(unit=unit,fmt='(A200)',iostat=last)dum
       read(unit=unit,fmt='(a5,2x,f10.6,2x,f11.6,2x,a40,2x,a3,2x,a3,2x,a1,2x,f9.4,2x,f9.4)',iostat=last)&
            psmsl%stat(n)%sid,psmsl%stat(n)%lat,psmsl%stat(n)%lon,&
            psmsl%stat(n)%name,psmsl%stat(n)%coast(1:3),psmsl%stat(n)%scode&
            ,psmsl%stat(n)%qcflag,psmsl%stat(n)%t0,psmsl%stat(n)%tn
       if(last .ne. 0)exit
       !lookup coastline code
       do,i=1,psmsl%ncount
          if(index(psmsl%coast(i),psmsl%stat(n)%coast(1:3)) > 0)then
             psmsl%stat(n)%coast=psmsl%coast(i)
             exit
          end if
       end do
    end do
    psmsl%nstat=n-1 ! set amount of stations found
    close(unit)
    

  end subroutine psmsl_init


  
 end module psmsl_tool


! c
! c example job to read the file of PSMSL monthly mean values.
! c for a description of each variable, see 'psmsldat.hel' and the
! c 'Data Holdings of the PSMSL' report.
! c
!       CHARACTER FILEIN*80
!       CHARACTER SNAME*40,SLAT*8,SLON*8,ACODE*2,FCODE*2
!       CHARACTER*80 TEXTS(999),TEXTA(999),TEXTC(999)
!       CHARACTER*3 CCODE,SCODE,GLOSS
!       CHARACTER*26 MISSDAYS,AMISSDAYS(200)
!       CHARACTER*1 RLRMET,MR,MRS
!       CHARACTER*2 MDAYS(13)
!       CHARACTER*1 IDOCFLS,IDOCFLY(3000)
!       INTEGER YR,MEANM(13),MEANR(13),RLRFAC,IYRLR
!       INTEGER IYEAR(200),AMEANM(13,200),AMEANR(13,200)
! C
! C OPEN THE FILE - THIS OPEN STATEMENT MAY BE SYSTEM DEPENDENT
! C
!       WRITE(6,*)
!       WRITE(6,*) 'ENTER INPUT FILE NAME'
!       READ(5,906) FILEIN
!   906 FORMAT(A)
!       OPEN(1,FILE=FILEIN,STATUS='OLD',FORM='FORMATTED',IOSTAT=IST)
!       IF(IST.NE.0) WRITE(6,*) ' ERROR FILE1: IST =',IST
!       IF(IST.NE.0) STOP
!       NNYEAR=0
!       NNCOMS=0
!       NNCOMC=0
!       NNCOMA=0
! C
!     1 CONTINUE
!       READ(1,901,END=9) SNAME,CCODE,SCODE,SLAT,SLON,ACODE,FCODE,IYRLR,
!      & GLOSS,MRS
!   901 FORMAT(A40,2A3,2A8,2A2,I4,A3,A1,6X)
!       WRITE(6,*) ' READING STATION ',CCODE,SCODE,' CALLED ',SNAME
! c
! c SNAME is the station name
! c CCODE,SCODE is the country,station code
! c SLAT,SLON are latitude,longitude
! c ACODE,FCODE are the authority,frequency codes
! c IYRLR will be 9999 if the station is not RLR i.e.Metric only
! c GLOSS will be '   ' if the station is not in GLOSS
! c MRS flags that there should be an entry for this station
! c in the documentation which follows.
! c
!       IDOCFLS=MRS
!       DO 1530 I=1,3000
!       IDOCFLY(I) = ' '
!  1530 CONTINUE
! C
! c Subroutine LATLON converts SLAT,SLON to REAL*4 parameters
! c
!       CALL LATLON(SLAT,SLON,ALAT,ALON)
! C
!       READ(1,1901) NYEAR,NCOMS,NCOMC,NCOMA
!  1901 FORMAT(4I3,68X)
!       NNYEAR=NNYEAR+NYEAR
!       NNCOMS=NNCOMS+NCOMS
!       NNCOMC=NNCOMC+NCOMC
!       NNCOMA=NNCOMA+NCOMA
! c
! c NYEAR is the number of years of data
! c NCOMS is the number of lines of station comments
! c NCOMC is the number of lines of country comments
! c NCOMA is the number of lines of authority comments
! c
! c the IYEAR(IY),AMISSDAYS(IY),AMEANM(I,IY),AMEANR(I,IY) arrays
! c store the year, missing days, metric and rlr values for
! c month I =(1,13) and year counter IY =(1,NYEAR)
! c
!       DO 7 IY=1,200
!       IYEAR(IY)=0
!       DO 17 I=1,13
!       AMEANM(I,IY)=99999
!       AMEANR(I,IY)=99999
!    17 CONTINUE
!     7 CONTINUE
! C
!       IF(NYEAR.GT.200) THEN
!           WRITE(6,*) ' ARRAY SIZE EXCEEDED'
!           STOP
!       ENDIF
! C  
!       IF(NYEAR.EQ.0) GOTO 5
!       DO 6 IY=1,NYEAR
!       READ(1,911) YR,MISSDAYS,MR
!   911 FORMAT(I4,6X,A26,4X,A1,39X)
!       READ(MISSDAYS,913) MDAYS
!   913 FORMAT(13A2)
! c
!       IYEAR(IY)=YR
!       AMISSDAYS(IY)=MISSDAYS
! c
! c YR is the year of this station-year
! c MISSDAYS contains 2 bytes of missing days information for each
! c month in the format as described in 'Data Holdings'.
! c MR flags that there should be an entry for this station-year
! c in the documentation which follows.
! C
!       IDOCFLY(YR)=MR
! C
!       READ(1,912) MEANM,RLRFAC
!   912 FORMAT(13I5,I10,5X)
! c
! c MEANM(1-12) are the 12 Metric monthly mean values. MEANM(13) is
! c the annual mean value. A value of any of these of 99999 flags
! c a missing monthly or annual mean. RLRFAC is the RLR factor for
! c the station-year. A value of 99999 flags that this year is
! c not RLR.
! c
! c To convert Metric values to RLR values, add the RLR factor. RLR
! c values should be approximately 7000. Units are millimetres.
! c
!       DO 11 I=1,13
!       MEANR(I)=99999
!       IF(RLRFAC.NE.99999.AND.MEANM(I).NE.99999)
!      &    MEANR(I)=MEANM(I) + RLRFAC
! C
!       AMEANM(I,IY)=MEANM(I)
!       AMEANR(I,IY)=MEANR(I)
!    11 CONTINUE
! C
!     6 CONTINUE
! c
! c now read the station, country and authority comments
! c
!     5 IF(NCOMS.EQ.0) GOTO 2
!       DO 21 I=1,NCOMS
!       READ(1,903) TEXTS(I)
!   903 FORMAT(A80)
!    21 CONTINUE
!     2 IF(NCOMC.EQ.0) GOTO 3
!       DO 31 I=1,NCOMC
!       READ(1,903) TEXTC(I)
!    31 CONTINUE
!     3 IF(NCOMA.EQ.0) GOTO 4
!       DO 41 I=1,NCOMA
!       READ(1,903) TEXTA(I)
!    41 CONTINUE
!     4 CONTINUE
! c
! c the subroutine PRTSTN prints data in the PSMSL 'Publication Format'
! c - see the 'Data Holdings' report for more information. The example
! c below prints RLR values (hence RLRMET ='r') for Walvis Bay. If RLRMET
! c is not 'r' or 'R', then the Metric data will be printed.
! c
!       IF(CCODE.EQ.'427'.AND.SCODE.EQ.'001') THEN
!          RLRMET='R'
!          CALL PRTSTN(RLRMET,SNAME,CCODE,SCODE,SLAT,SLON,ACODE,
!      &   FCODE,IYRLR,GLOSS,NYEAR,IYEAR,AMISSDAYS,AMEANM,AMEANR,
!      &   NCOMS,TEXTS,NCOMC,TEXTC,NCOMA,TEXTA,
!      &   IDOCFLS,IDOCFLY)
!       ENDIF
! C
!       GOTO 1
!     9 WRITE(6,*) 
!      &' END OF FILE: NNYEAR (Total No. of station-years) = ',NNYEAR
!       WRITE(6,*) 'NNCOMS (No.station comments) = ',NNCOMS
!       WRITE(6,*) 'NNCOMC (No.country comments) = ',NNCOMC
!       WRITE(6,*) 'NNCOMA (No.authority comments) = ',NNCOMA
!       CLOSE (1)
!       STOP
!       END
!       SUBROUTINE LATLON(CLAT,CLON,ALAT,ALON)
! C
!       SAVE
! C
! C CLAT AND CLON ARE 8 BYTE CHARACTERS CONTAINING LAT AND LON
! C E.G. ' 51 44 N' OR '123 33 E'
! C ALAT AND ALON ARE LAT (RANGE +90 TO -90 NORTH POSITIVE)
! C               AND LON (RANGE 0 TO 360 EAST)
! C
!       CHARACTER*8 CLAT,CLON
!       CHARACTER*1 CH(8),CC
! C
!       READ(CLAT,901) CH
!   901 FORMAT(8A1)
!       IF(CH(1).NE.' '.OR.CH(4).NE.' '.OR.CH(7).NE.' '.OR.
!      &  (CH(8).NE.'N'.AND.CH(8).NE.'S')) THEN
!        WRITE(6,*) 
!      & ' ILLEGAL FORMAT IN S/R LATLON: CLAT =',CLAT
!        STOP
!       ENDIF
! C
!       READ(CLON,901) CH
!       IF(CH(4).NE.' '.OR.CH(7).NE.' '.OR.
!      &  (CH(8).NE.'W'.AND.CH(8).NE.'E')) THEN
!        WRITE(6,*)
!      & ' ILLEGAL FORMAT IN S/R LATLON: CLON =',CLON
!        STOP
!       ENDIF
! C
!       READ(CLAT,902) II,JJ,CC
!   902 FORMAT(2I3,1X,A1)
!       ALAT = FLOAT(II) + FLOAT(JJ)/60.
!       IF(CC.EQ.'S') ALAT = - ALAT
! C
!       READ(CLON,902) II,JJ,CC
!       ALON = FLOAT(II) + FLOAT(JJ)/60.
!       IF(CC.EQ.'W') ALON = 360. - ALON
! C
!       RETURN
!       END
!       SUBROUTINE PRTSTN(RLRMET,SNAME,CCODE,SCODE,SLAT,SLON,ACODE,
!      &   FCODE,IYRLR,GLOSS,NYEAR,IYEAR,AMISSDAYS,AMEANM,AMEANR,
!      &   NCOMS,TEXTS,NCOMC,TEXTC,NCOMA,TEXTA,
!      &   IDOCFLS,IDOCFLY)
! C
!       SAVE
! C
!       CHARACTER*1 RLRMET,METRLR
!       CHARACTER CCODE*3,SCODE*3,SNAME*40,SLAT*8,SLON*8
!       CHARACTER ACODE*2,FCODE*2,GLOSS*3
!       CHARACTER*80 TEXTS(999),TEXTC(999),TEXTA(999)
!       CHARACTER*1 IDOCFLS,IDOCFLY(3000)
!       CHARACTER*26 MISSD,AMISSDAYS(200)
!       INTEGER YR,IYRLR,IYEAR(200),AMEANM(13,200),AMEANR(13,200)
!       CHARACTER*2 CHMD(13),BL,DASH
!       CHARACTER*4 CHMDB(13),BL4
!       CHARACTER*5 DOT,ACON(13)
!       DIMENSION ICON(13)
!       DATA DOT/'   ..'/,BL/'  '/
!       DATA BL4/'    '/
!       DATA DASH/' -'/
! C
!       METRLR=RLRMET
!       IF(METRLR.EQ.'r') METRLR='R'
! C
!       IF(IYRLR.EQ.9999.AND.METRLR.EQ.'R') GOTO 7002
!       WRITE(6,800)
!   800 FORMAT(1H1)
!       WRITE(6,401) SNAME,SLAT,SLON
!   401 FORMAT(38X,A40,A8,1X,A8)
!       WRITE(6,412) CCODE,SCODE,ACODE,FCODE
!   412 FORMAT(/,26X,'COUNTRY CODE: ',A3,3X,'STATION CODE: ',A3,3X,
!      & 'AUTHORITY CODE: ',A2,3X,'FREQUENCY CODE: ',A2)
!       IF(GLOSS.NE.'   ') WRITE(6,1412) GLOSS
!  1412 FORMAT(/,26X,29X,'GLOSS CODE  : ',A3)
!       IF(METRLR.EQ.'R') WRITE(6,1403) IYRLR
!  1403 FORMAT(//,45X,'VALUES ARE MEASURED TO DATUM OF RLR ',I5)
!       IF(METRLR.NE.'R') THEN
!                WRITE(6,404)
!   404          FORMAT(//,38X,'SUPPLIED DATA VALUES ONLY -',
!      &                   ' NOT MEASURED TO A COMMON DATUM')
!                WRITE(6,4404)
!  4404          FORMAT(38X,17X,
!      &         '(I.E. A "METRIC" RECORD)')
!                WRITE(6,4405)
!  4405          FORMAT(33X,'THESE VALUES SHOULD NOT BE USED FOR',
!      &                    ' MULTI-YEAR TIME SERIES ANALYSIS')
!       ENDIF
! C
!       WRITE(6,405)
!   405 FORMAT(/,1H ,37X,'MONTHLY & ANNUAL MEAN HEIGHTS OF SEA LEVEL',
!      &' IN MILLIMETRES.')
! C
!       WRITE(6,406)
!   406 FORMAT(/,1H ,9X,'I',8X,'II',7X,'III',6X,'IV',7X,'V',
!      &8X,'VI',7X,'VII',5X,'VIII',6X,'IX',7X,'X',8X,'XI',7X,
!      &'XII',6X,'Y')
! C
!       IF(NYEAR.LE.0) GOTO 403
!       DO 402 IY=1,NYEAR
!       YR=IYEAR(IY)
!       MISSD=AMISSDAYS(IY)
! C
! C SET UP THE METRIC OR RLR VALUES.
!       DO 43 I=1,13
!       ICON(I)=AMEANM(I,IY)
!       IF(METRLR.EQ.'R') ICON(I)=AMEANR(I,IY)
!       WRITE(ACON(I),943) ICON(I)
!   943 FORMAT(I5)
!       IF(ICON(I).EQ.99999) ACON(I)=DOT
!    43 CONTINUE
! C
! C READ THE MISSING DAYS INFORMATION
!       DO 246 I=1,13
!       CHMD(I)=BL
!       CHMDB(I)=BL4
!   246 CONTINUE
!       READ(MISSD,902) CHMD
!   902 FORMAT(13A2)
!       DO 146 I=1,13
!       IF(I.EQ.13.AND.CHMD(I).EQ.DASH) CHMD(I)=BL
!       IF(CHMD(I).EQ.BL) GOTO 146
!       WRITE(CHMDB(I),946) CHMD(I)
!   946 FORMAT('(',A2,')')
!   146 CONTINUE
!    46 CONTINUE
! C
!       WRITE(6,409) YR,(ACON(I),CHMDB(I),I=1,13),YR
!   409 FORMAT(2X,I4,2X,13(A5,A4),I4)
! C
! C IF A YEAR OF DATA IS 'METRIC ONLY' AND THIS IS AN RLR PRINTOUT
! C THEN EACH MONTH'S INFORMATION WILL BE A SET OF DOTS.
! C
!   402 CONTINUE
!   403 CONTINUE
! C
! C PRINT EXPLANATION OF ABOVE PRINTOUT
! C
!       WRITE(6,4601)
!  4601 FORMAT(/,' VALUES IN BRACKETS SHOW NUMBER OF MISSING DAYS',
!      &         ' EACH MONTH',
!      &         ' WITH NO INTERPOLATION MADE IN COMPUTING THE MEAN;')
!       WRITE(6,4602)
!  4602 FORMAT(' "XX" SIGNIFIES MISSING OBSERVATIONS WERE INTERPOLATED',
!      &       ' BEFORE COMPUTING THE MONTHLY MEAN;')
!       WRITE(6,4603)
!  4603 FORMAT(' "XX" FOR AN ANNUAL MEAN SIGNIFIES A VALUE LIKELY TO BE',
!      & ' MATERIALLY AFFECTED BY MISSING DATA;')
!       WRITE(6,4604)
!  4604 FORMAT(' YEARS WITH MORE THAN ONE MISSING MONTH HAVE ANNUAL',
!      &       ' MEANS DROPPED.')
! C
!       WRITE(6,6600) ACODE
!  6600 FORMAT(/,' DATA COME FROM AUTHORITY "',A2,
!      & '" - SEE FILE indexa.html ON PSMSL DISK FOR FULL ADDRESS')
! C
! C READ DOCUMENTATION FOR THIS STATION.
!       WRITE(6,411)
!   411 FORMAT(/,1X,'ANY COMMENTS FOR THIS STATION ARE GIVEN BELOW:',/)
!       IF(FCODE.EQ.' C'.OR.FCODE.EQ.'C ') WRITE(6,6611)
!  6611 FORMAT(' FREQUENCY',
!      &       ' CODE "C " IMPLIES DATA OBTAINED FROM INTEGRATION',
!      & ' FROM CONTINUOUS RECORDS')
!       IF(FCODE.EQ.'HL') WRITE(6,6612)
!  6612 FORMAT(' FREQUENCY',
!      &       ' CODE "HL" IMPLIES MEAN TIDE LEVEL (I.E. HIGH AND',
!      & ' LOW WATERS)')
!       IF(FCODE.NE.'C '.AND.FCODE.NE.' C'.AND.FCODE.NE.'HL')
!      & WRITE(6,6613) FCODE,FCODE
!  6613 FORMAT(' FREQUENCY',
!      &       ' CODE "',A2,'" IMPLIES DATA OBTAINED FROM ',A2,
!      & ' READINGS PER DAY')
!       WRITE(6,*)
! C
! C PRINT DOCUMENTATION FLAG WARNINGS
! C
!       ISKIP=0
!       IF(IDOCFLS.NE.' ') THEN
!          ISKIP=1
!          WRITE(6,9502)
!  9502 FORMAT(' WARNING:',
!      &' DOCUMENTATION FLAG SET FOR ENTIRE STATION - SEE COMMENTS BELOW')
!       ENDIF
!       DO 503 I=1,3000
!       IF(IDOCFLY(I).NE.' ') THEN
!           ISKIP=1
!           WRITE(6,9503) I
!  9503 FORMAT(' WARNING: DOCUMENTATION FLAG SET FOR YEAR ',I5,'  - SEE',
!      &    ' COMMENTS BELOW')
!       ENDIF
!   503 CONTINUE
!       IF(ISKIP.EQ.1) WRITE(6,*)
! C
! C STATION COMMENTS
!       IF(NCOMS.LE.0) GOTO 502
!       DO 501 I=1,NCOMS
!       WRITE(6,410) TEXTS(I)
!   410 FORMAT(1X,A80)
!   501 CONTINUE
!   502 CONTINUE
! C
! C COUNTRY COMMENTS
!       IF(NCOMC.LE.0) GOTO 1502
!       WRITE(6,1411)
!  1411 FORMAT(/,1X,'ANY COMMENTS FOR THIS COUNTRY ARE GIVEN BELOW:',/)
!       DO 1501 I=1,NCOMC
!       WRITE(6,410) TEXTC(I)
!  1501 CONTINUE
!  1502 CONTINUE
! C
! C AUTHORITY COMMENTS
!       IF(NCOMA.LE.0) GOTO 2502
!       WRITE(6,2411)
!  2411 FORMAT(/,1X,'ANY COMMENTS FOR THIS AUTHORITY ARE GIVEN BELOW:',/)
!       DO 2501 I=1,NCOMA
!       WRITE(6,410) TEXTA(I)
!  2501 CONTINUE
!  2502 CONTINUE
! C
!   700 CONTINUE
!       RETURN
! C
!  7002 WRITE(6,97002) CCODE,SCODE
! 97002 FORMAT(' STATION ',A3,'/',A3,' IS NOT AN RLR STATION',
!      &       ' NO PRINTOUT PRODUCED')
!       GOTO 700
!       END
