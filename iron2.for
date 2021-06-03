      PROGRAM IRON2
* Copyrights Roberto Moretti, 2021                                        *
* V1.0                                                                    *
* ATTENZIOEN A X10 e X23 SE USI MISFIT
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* GAS SYSTEM IS H-O-S. Carbon is not accounted for.....                   *
*                                                                         *
*_________________________________________________________________________*
* based on the simplified polymeric approach of
* Ottonello et al. (2000) for the evaluation of redox conditions
*                                                                         *
* ancillary informations                                                  *
*                                                                         *
* Z: basicity moderating parameters (eq. 29-II in Ottonello et al. (2000))*
*     exp. data of Young et al. (1992) with modifications given in        *
*     Moretti & Ottonello (2002) (e.g. water set to 2.56)                 *
*                                                                         *
*  Y: oxides gfw                                                          *
*                                                                         *
*  Zcat: n. cations per formula unit                                      *
*                                                                         *
* volx: molar volumes of oxides in the melt phase (cc/mol)                *
*                                                                         *
* eoxm: expansivity of oxides in melt phase (cc/k mol)                    *
*                                                                         *
* coxm: compressibility of oxides in melt phase                           *
*                                                                         *
* Reference standard state is 1 bar @ T of interest                       *
*                                                                         *
* Fugacities are employed. Activity terms in the liquid phase are         *
* corrected for pressure greater than 1 bar.                              *
*                                                                         *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     DECLARATIONS
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      common /setmod/answ
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),  XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,     A,foga(ndi),  XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         regg(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &        volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi), sobs(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ,ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
 
*     UNIT NUMBERS
 
      ISYSRD = 17   ! Read unit
      IVIDEO =  6
      ISYSWR =  8
      ISYSIN =  9
      ISYSCK = 12
      ISYSQZ = 13
      ISYSWZ = 14
      ISYSPR = 15
      ISYSPU = 16
      ISYSMF = 19
      ISYSVO = 18
 
      open(file='INPUT.txt'   ,unit=ISYSRD,status='UNKNOWN')
      open(File='OUTPUT.txt'  ,unit=ISYSWZ,status='UNKNOWN')

      write (ISYSWZ,*) 'TK    Pbars   loFe2/Fe3exp   logFe2/Fe3calc
     $logfO2exp     logfO2calc '
 
      answ=0.
 
      CALL NEWSULF
 
      END
 
      BLOCK DATA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       sumfe(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
 
      DATA Z/2.09,1.54,2.5,1.67,2.09,1.82,1.354,1.40,1.,1.28,.87,
     & .71,2.50/
c N.B gamma per Mn ä 1.69 nel lavoro precedente e per H2O 2.56!
      DATA Y/60.085,79.899,141.94,101.96,159.69,151.99,71.846,
     & 70.937,56.079,40.311,61.979,94.203,18.015/
 
      DATA Zcat/1.,1.,2.,2.,2.,2.,1.,1.,1.,1.,2.,2.,2./
  
**********************************************************************
c    GOOD VALUES ARE BELOW (27 november 2003)
***********************************************************************
c   VERY IMPORTANT: Verify x27 and x227 (i.e. the annealing term -alfaV)
 
c  OXIDE-SULFIDE disproportions: new values 23 March 2004 for H2O
      DATA X/12.752,-4.7934,1.7132,7.70870,15.0328,0.,5.4414,-0.917,
c  Below Papale value without P effect
c    $ -4.0239,-4.1698,14.5300,11.266,-3.341,-36418.,
c  Below Papale value with P effect (S annealing)
     $ -4.0239,-4.1698,14.5300,11.266,-3.341,-36418.,
c  Below Burnham value without P effect
c    $ -4.0239,-4.1698,14.5300,11.266,-4.069,-36418.,
c  Below Burnham value with P effect (S annealing)
*    $ -4.0239,-4.1698,14.5300,11.266,-4.069,-36418.,
     $-16575.9,-24786.,-28926.,-20000.,0.,-15265.,-8391.7,-6804.2,
     $-15181.,-16396.,-16432.,-18764./
 
c   OXIDE_SULFATE disproportions: new values (23 March 2004)
      DATA XX/-26.436,-27.0906,0.,-32.98960,-33.89291,0.,-40.958569,
     $-60.128854,-42.38865,-36.408443,-22.089,13.70000,
c  below Papale model derived value without P effect
c    $-65.95324,77291.1,98230.0,0.,
c  below Papale model derived value with P effect (S annealing)
     $-121.62,77291.1,98230.0,0.,
c  below Burnham model derived value without P effect
c    $-60.8531,77291.100,98230.000,0.,
c  below Burnham model derived value with  P effect (S annealing)
*    $-105.488,77291.100,98230.000,0.,
     $85879.800,85879.800,0.,99540.800,106414.10,107754.894,95831.,
     $100934.,100741.,67292./


c NEW VOLUME REGRESSION (file Volumes_melts.xls) Expansivities for SiO2 and Al2O3
c have been adjusted on the basis of the OptBas correlation. Water Expansivity obtained
c by considering H2O estimated V at 1673K and Richet Volume at room T (12 cc/mol).
c About expansivities:
c SiO2 : use 6.24 (recalculated) or 0.10 (Lange)
c Al2O3: use 8.10 (recalculated) or 2.62 (Lange)
c P2O5:  use 18.3 (recalculated) or 2.62 (P2O5 as Al2O3 of Lange)
      data voxm/26.9,23.16,82.16,37.11,42.13,36.36,13.65,11.62,16.57,
     $          11.45,28.78,45.84,16.79/
      data eoxm/0.10,7.24,2.62,2.62,9.09,8.22,2.92,2.73,2.92,2.62,7.41,
     $          11.91,3.48/
      data coxm/-1.89,-2.31,-8.93,-2.26,-2.53,-1.96,-0.45,-0.37,-1.34,
     $          -0.40,-2.4,-6.75,-1.88/
***********************************************************************
c  segue per calcolo porositÖ ionica
 
      data noxy/2,2,5,3,3,3,1,1,1,1,1,1,1/
      data radi/0.26,0.605,0.17,0.39,0.49,0.615,0.78,0.67,1.,0.72,1.02,
     $          1.37,0./
 
      END
 
      SUBROUTINE NEWSULF
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* this is the main program, disguised as a subroutine for *
* reasons of compatibility between systems.               *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      common/setmod/answ
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         reg3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
 
      Read (isysrd,*) ncomp
	  Read (isysrd,*) ! This is to skip the header row
C     rapmed=0.0
      do 111 i=1,ncomp
C     jflag(i)=0
      Read (isysrd,*) poss(i,1),poss(i,2),poss(i,4),poss(i,5),poss(i,6),
     &poss(i,7),poss(i,8),poss(i,10),
     &poss(i,9),poss(i,11),poss(i,12),poss(i,3), poss(i,13),
     &tcent(i),pbar(i),xossi(i),iaut(I),kflag(i)

c Nota Bene:
c kflag=0. il programma itera per Fe2+/Fe3+
c kflag=2. il programma accetta Fe2+/Fe3+ di input e log fO2
 
 
      write(6,*)'I read ',i
      tkelv(i)=tcent(i)+273.15
      WRITE (*,*) '  '
      WRITE (*,*)'---------------- COMPOSITION ',I,'-------------------'
 
      if (poss(i,5).ne.0.or.poss(i,7).ne.0) then
      write (*,*) 'Composition n. ',i,' holds Iron !!!'
      write (*,*) 'Fe2O3 =',poss(i,5)
      write (*,*) 'FeO  =',poss(i,7)
      endif
      jflag=0.
 
c  calcola le proporzioni molari a 100% catione e ossido
c  poi si assume acqua misurata indipendentemente (tipo FT-IR)
c  per cui non entra nella normalizzazione.....
      summa = 0.
      summb = 0.
      do 961 lk=1,12
      summa=poss(i,lk)+summa
961   continue
      do 871 lk=1,12
      poss(i,lk)=poss(i,lk)*100./summa
871   continue
      wt(i)=summa
      summc=0.
      do 9623 lh = 1,12
      poss(i,lh)= poss(i,lh)*(1.-poss(i,13)/100.)
      write (*,*) poss(i,lh)
      summc=summc+poss(i,lh)
9623  continue
c     write (*,*) poss(i,13),summc+poss(i,13)
c     pause
 
      sumfe(i)=poss(i,5)+poss(i,7)
      NIT = 50
      do 686 kj=1,NIT
 
      write (*,*) 'Poss 5 = ',poss(i,5)
      write (*,*) 'Poss 7 = ',poss(i,7)
 
      if (poss(i,5).ne.0.or.poss(i,7).ne.0) then
      write (*,*) 'Iteration N. ',kj,' @ composition ',i
      endif
 
c  Procedura di normalizzazione su FeOtotale (somma degli ossidi wt%)
      if (sumfe(i).ne.0) then
      summb=poss(i,5)+poss(i,7)
      poss(i,5)=poss(i,5)*sumfe(i)/summb
      poss(i,7)=poss(i,7)*sumfe(i)/summb
      endif
 
      d(i)=0
      dd(i)=0
      do 100 j=1,13
c     write(*,*) j,poss(i,j),zcat(j),y(j)
      xcat(j)=poss(i,j)*zcat(j)/y(j)
c     pause
      d(i) = d(i)+ xcat(j)
      dd(i)=dd(i)+poss(i,j)/y(j)
 100  continue
      som(i)=d(i)
      xcat0(i)=xcat(13)
c considero tutto H2O come OH (o H!)
      dd(i)=dd(i)-poss(i,13)/y(13)
      xw(i)=2.*poss(i,13)/y(13)
      dd(i)=dd(i)+xw(i)
      somm(i)=dd(i)
      do 200 j=1,13
      xcat(j)=xcat(j)/d(i)
 200  continue
      xfetot(i)=xcat(5)+xcat(7)
      if (kj.eq.1) xftot=xfetot(i)
 
CCC   if (xcat(5).ne.0.and.xcat(7).ne.0.and.kj.eq.1) jflag=2
      IF (KFLAG(I).EQ.2) JFLAG=2
 
       if (xcat(5).ne.0.and.xcat(7).ne.0.) then
        write (*,*) 'redoz is = ',dlog10(xcat(7)/xcat(5))
        redoz(i)=dlog10(xcat(7)/xcat(5))
       endif
       if (xcat(5).ne.0.and.xcat(7).eq.0.) then
        write (*,*) 'redoz is UNDEFINED'
        redoz(i)=1.0
       endif
       if (xcat(5).eq.0.and.xcat(7).ne.0.) then
        write (*,*) 'redoz is UNDEFINED'
        redoz(i)=1.00
       endif
       if (xfetot(i).eq.0) then
        write (*,*) 'redoz is UNDEFINED'
        xcat(5)=0.0000000001
        xcat(7)=0.0000000001
       endif
 
       if (kj.eq.1) then
        refirst(i)=redoz(i)
       endif
 
519   continue
c
c  Verifica la consistenza.
      FRIT=0.
      KWFLAG=0
      jwflag=0
      ITER=50
      if (poss(I,13).eq.0.) then
       iter=1
       kwflag=3
      endif
      do 789 jk=1,iter
      jwflag=jwflag+1
 
      CALL CALCOMP
      CALL BASOPT
      CALL TOOPSAMIS
      CALL ACTION
 
      frit=dabs(aossi(i)-aossiz(i))
        if (poss(i,13).ne.0.and.frit.lt.0.00001.and.jk.gt.1) then
         kwflag = 3
        endif
 
      if (kwflag.eq.3) goto 790
      if (jk.eq.iter) then
      pause
      goto 790
 
      endif
 
 
789   continue
 
790   CALL TOOPSAMIS2
      CALL ACTION2
      CALL FEREDOX
 
      if (xfetot(i).eq.0.or.jflag.eq.1.) then
      goto 696
      endif
 
      write (*,*) 'redox last (redoX) = ',redox(i)
      write (*,*) 'redox before (redoZ) = ',redoz(i)
 
      zxf3=(xfetot(i)/(10**redoz(i)+1))
      zxf2=(xfetot(i)/(10**redoz(i)+1))*(10**redoz(i))
 
      xf3=(xfetot(i)/(10**redox(i)+1))
      xf2=(xfetot(i)/(10**redox(i)+1))*(10**redox(i))
 
c      xff3=xf3*d(i)
c      xff2=xf2*d(i)
      xff3=xf3*d(i)
      xff2=xf2*d(i)
      pfe3=xff3*y(5)/zcat(5)
      pfe2=xff2*y(7)
 
      if (poss(i,5).eq.0.) then
       poss(i,5)=0.01
       poss(i,7)=poss(i,7)-0.01
      endif
 
      if (poss(i,7).eq.0.) then
       poss(i,7)=0.01
       poss(i,5)=poss(i,5)-0.01
      endif
      write (*,*) 'Fe2/Fetot calc  = ', xf2/(xf2+xf3)
      write (*,*) 'Fe2/Fetot before= ', zxf2/(zxf2+zxf3)
 
      chi=dabs(bossi(i)-xossi(i))
      chio=dabs(xf2/(xf2+xf3)-zxf2/(zxf2+zxf3))
      write (*,*) 'bossi = ',bossi(i), 'xossi = ',xossi(i)
      write (*,*) 'REMINDER: bossi ä ricalcolato da redox!'
      write (*,*) 'REMINDER: bossi = xossi se INPUT logfO2 = 0'
      write (*,*) ' JFLAG is ',jflag
      if (jflag.eq.2) then
      write (*,*)' !!! I HAVE USED DEFAULT FeO/Fe2O3 AND INPUT fO2 !!!'
      write (*,*)' !!!           I MADE NO ITERATIONS              !!!'
c     pause
      go to 696
      endif
 
c L'ho spostato qua ma stava sopra!!!!!
 
      poss(i,5)=dsqrt(pfe3*poss(i,5))
      poss(i,7)=dsqrt(pfe2*poss(i,7))
 
      if (chi.lt.0.0001.and.chio.lt.0.0001) then
      jflag=1
      write (*,*) 'JFLAG= ',jflag
c     pause
      endif
 
      if (kj.eq.nit) then
      write (*,*) 'This was the last iteration. '
c     pause
      endif
c     if (i.lt.ncomp) pause
686   continue
 
696   continue
 
C     rapmed=dlog10(cstart(i))+rapmed
 
*********************************************************************
 
c  calcola le proporzioni molari a 100% catione e ossido
      summit = 0.
      do 991 lk=1,13
       summit=poss(i,lk)/y(lk)+summit
991   continue
      xh2o(i)=(poss(i,13)/y(13))/summit
***************************************************************************
 
      CALL CALCPROP
      CALL DATAOUT
111   CONTINUE
      END
 
      SUBROUTINE CALCOMP
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     computes stoichiometric amounts                             *
*     charge balance for melt complexes                           *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      common/setmod/answ
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         reg3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
 
c     if (poss(i,13).ne.0) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c NOTA BENE............                                                      c
c INIZIA DA QUI IL RICALCOLO DI H2O SPECIATA SECONDO MORETTI & OTTONELLO     c
c (2003) IN PREP. TRA H+ E OH-. OH- ä ANIONE LIBERO SULLA MATRICE ANIONICA   c
c MENTRE H+ CONTRIBUISCE ALLA MATRICE W CON QUANTO NE CONSEGUE.              c
c QUINDI H+ ENTRA NEL COMPUTO DELLA DEPOLIMERIZZAZIONE SECONDO QUANTO EVINTO c
c DA FRASER.I RUOLI STRUTTURALI NON CI INTERESSANO (CHIUSURA TERMINAZIONI    c
c POLIMERICHE O QUANT'ALTRO....)                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      xcat(1)=poss(i,1)*zcat(1)/y(1)
      xcat(3)=poss(i,3)*zcat(3)/y(3)
      xcat(4)=poss(i,4)*zcat(4)/y(4)
      xcat(5)=poss(i,5)*zcat(5)/y(5)
      xcat(6)=poss(i,6)*zcat(6)/y(6)
      xcat(7)=poss(i,7)*zcat(7)/y(7)
      xcat(2)=poss(i,2)*zcat(2)/y(2)
      xcat(13)=poss(i,13)*zcat(13)/y(13)
      xcat(8)=poss(i,8)*zcat(8)/y(8)
      xcat(9)=poss(i,9)*zcat(9)/y(9)
      xcat(10)=poss(i,10)*zcat(10)/y(10)
      xcat(11)=poss(i,11)*zcat(11)/y(11)
      xcat(12)=poss(i,12)*zcat(12)/y(12)
      xcat(1)=xcat(1)/som(i)
      xcat(3)=xcat(3)/som(i)
      xcat(4)=xcat(4)/som(i)
      xcat(5)=xcat(5)/som(i)
      xcat(6)=xcat(6)/som(i)
      xcat(2)=xcat(2)/som(i)
      xcat(7)=xcat(7)/som(i)
      xcat(8)=xcat(8)/som(i)
      xcat(9)=xcat(9)/som(i)
      xcat(10)=xcat(10)/som(i)
      xcat(11)=xcat(11)/som(i)
      xcat(12)=xcat(12)/som(i)
      xcat(13)=xcat(13)/som(i)
      aossiz(i)=aossi(i)
      if (jk.eq.1.and.iter.gt.1.) then
      aossiz(i)=1.
      xoh(i)=xcat0(i)
      endif
 
      xh(i)=xcat0(i)-xoh(i)
      xohz(i)=xoh(i)
      xcat(13)=xh(i)/som(i)
      xwd(i)=xoh(i)/som(i)
 
c     endif
 
       w=3.0*xcat(6)+2.0*(xcat(7)+xcat(8)+xcat(9)+xcat(10))+xcat(11)+
     & xcat(12)+xcat(13)+xcat(3)+4.0*xcat(2)
 
      IF (Xcat(4).GT.W) GO TO 540
      xfor=xcat(1)+xcat(4)/2.0
      w=w-xcat(4)
      xallu=xcat(4)/2.0
      xcat(4)=0.
      GO TO 590
 540  u=xcat(4)-w
      xfor=xcat(1)+w/2.0
      xcat(4)=u
      XALLU=W/2.0
      w=0.
      GO TO 680
 590  IF (xcat(5).GT.W) GO TO 640
      xfor=xfor+xcat(5)/2.0
      w=w-xcat(5)
      xfe2o3=xcat(5)/2.0
c      feo2m=xcat(5)
      afeo2(i)=xcat(5)
      afeo22(i)=xcat(5)/2.0
c      fe3p=0
      xcat(5)=0.
      GO TO 680
 640  v=xcat(5)-w
      xfor=xfor+w/2.0
      afeo2(i)=w
      afeo22(i)=w/2.0
      xcat(5)=v
      xfe2o3=w/2.0
c      feo2m=w
c      fe3p=v
      w=0.
c     pause
 680  xfor=xfor+xcat(3)/2.0
      afeo2(i)=afeo2(i)/xfor
 1691 CONTINUE
 
      END
 
      SUBROUTINE BASOPT
 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     computes optical basicity of network formers and network          *
*     modifiers;                                                        *
*     defines oxides stoichiometric proportions for oxide-sulfide       *
*     reactions                                                         *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 k,mnom,mgom,nao05m,ko05m
      PARAMETER (ndi=250)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),    XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
C
      basfor=xcat(1)*z(1)/xfor+xcat(3)*z(3)/xfor+xallu*z(4)/xfor+
     & xfe2o3*z(5)/xfor
c
      scat(1)=xcat(1)
      scat(3)=xcat(3)
      scat(2)=0.0
      scat(4)=xallu
      scat(5)=xfe2o3
      xmod=0.0
      do 6969 ij=2,13
      xmod=xmod+xcat(ij)/zcat(ij)
6969  continue
      basmod=0.0
      do 6996 ij=2,13
      xcat(ij)=xcat(ij)/(zcat(ij)*xmod)
      basmod=basmod+xcat(ij)*z(ij)
6996  continue
      afe2(i)=xcat(7)
      afe22(i)=xcat(7)
c     afe2(i)=xcat(7)/xcat(10) ! usato solo per fornire [Fe]/[Mg]
      afe3(i)=xcat(5)
      afe32(i)=xcat(5)
      ana(i)=xcat(11)
      ak(i)=xcat(12)
      amg(i)=xcat(10)
      amn(i)=xcat(8)
      aca(i)=xcat(9)
      acr(i)=xcat(6)
      ati(i)=xcat(2)
      aprot(i)=xcat(13)
      ap(i)=xcat(3)
      k=dexp(((basmod-basfor)/0.2145)-1.1445)
      polcos(i)=k
      basdif(i)=basmod-basfor
      a=1.0-4.0*k
      xcat(1)=xfor/(xfor+xmod)
      acidic(i)=xcat(1)
      e=1.0-xcat(1)
      totcat(i)=e
      afe22(i)=afe22(i)/(xfor+xmod)
      afe32(i)=afe32(i)/(xfor+xmod)
      afe22(i)=afe22(i)/totcat(i)*xmod
      afe32(i)=afe32(i)/totcat(i)*xmod
      afeo2(i)=afeo2(i)/(xfor+xmod)
      afeo22(i)=afeo22(i)/(xfor+xmod)
      END
 
      SUBROUTINE TOOPSAMIS
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  polymeric aproach to melt structure (see Ottonello et al., 2000; *
*  Ottonello, 2001; for details)                                    *
*  computes singly bonded, doubly bonded and free oxygen ion on a   *
*  1-mole of melt basis and furnishes structural details such as    *
*  mean extension of polymeric units, acidity etc. etc.             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
 
      atoop=-a
      btoop=2.0+2.0*xcat(1)
      ctoop=8.0*xcat(1)*(xcat(1)-1.0)
      o(i)=(-btoop+dsqrt(btoop**2.-4.0*atoop*ctoop))/(2.0*atoop)
 
      o1(i)=1.0-xcat(1)-o(i)/2.0
      o3(i)=(4.0*xcat(1)-o(i))/2.0
      o4(i)=o(i)/(o(i)+o3(i)+xcat(1))
 
      s=dexp(-1.7165*dlog(o4(i))+2.8776)
      s1=acidic(i)/s
      o2=o1(i)/(o1(i)+s1+xwd(i))
c     aossi(i)=o2
      totani(i)=o1(i)+s1+xwd(i)
      totali(i)=totani(i)
      ah(i)=totcat(i)/totani(i)
      END
 
      SUBROUTINE ACTION
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     computes the activity of sulphide and of iron ionic species   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),    XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
C
      totold(i)=o1(i)+s1
      struct(i)=s1
      aossi(i)=o2
      aozzi(i)=aossi(i)

c  I volumi sotto sono in joule/bar ...fattore 0.1...
c  Raggi ionici di Shannon
       voh=(4./3.)*3.14159*(1.40+0./10000*(tkelv(i)-298.15))**3
      voh=voh*0.6022045
      vo2=(4./3.)*3.14159*(1.40+0./10000*(Tkelv(i)-298.15))**3
      vo2=vo2*0.6022045
      vh=(4./3.)*3.14159*(0.+0./10000*(tkelv(i)-298.15))**3
      vh=vh*0.6022045
      if (VH.LT.0) vh=0.
c Calcolo delv assumendo espansione termica = 0.
      delv0=vh+vo2-voh+0./1000.*(tkelv(i)-298.15)
      delv(i)=delv0
      delv0=delv0*0.1/(8.3147*2.303)
      delv0=delv0*(pbar(i)-1)/tkelv(i)

c   19 maggio 2004
c        partw=-1.05834853820627+607.514123515489/tkelv(i)
        partw=-1.835+1304.65/tkelv(i)
        partw=10.**partw
      ratius=(aossi(i)/ah(i))/partw
      xoh(i)=xcat0(i)*ratius/(ratius+1)
      xh(I)=XCAT0(I)-XOH(I)
      xoh(i)=dsqrt(xohz(i)*xoh(i))
 
      END
 
      SUBROUTINE TOOPSAMIS2
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  polymeric aproach to melt structure (see Ottonello et al., 2000; *
*  Ottonello, 2001; for details)                                    *
*  computes singly bonded, doubly bonded and free oxygen ion on a   *
*  1-mole of melt basis and furnishes structural details such as    *
*  mean extension of polymeric units, acidity etc. etc.             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi), AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
 
c speculazione su H+ e O-
      voh=(4./3.)*3.14159*((1.40+0./1.d6*(tkelv(i)-298.15))**3)
      voh=voh*0.6022045
c  sarebbe x(7) il coeff lineare di espansione
      vo=(4./3.)*3.14159*((1.40+0./1.d6*(Tkelv(i)-298.15))**3)
      vo=vo*0.6022045
      vh=(4./3.)*3.14159*((0.+0./1.d6*(tkelv(i)-298.15))**3)
      vh=vh*0.6022045
      if (VH.LT.0) vh=0.
      dv1=-vh-vo+voh
      dv(i)=dv1
      dv1=dv1*0.1/(8.3147*2.303)
      dv1=dv1*(pbar(i)-1)/tkelv(i)

      costx=10**(-1.335)
      aoin=o(i)
      ahin=xh(i)/som(i)
      aspec=costx
      bspec=-(costx*ahin+costx*aoin+totcat(i))
      cspec=costx*ahin*aoin
      root1=(-bspec-dsqrt(bspec**2-4*aspec*cspec))/(2*aspec)
c     root2=(-bspec+dsqrt(bspec**2-4*aspec*cspec))/(2*aspec)

c ATTENZIONE EFFETTUO UN CAMBIO IMPORTANTE:
c      root1=0.
      o(i)=o(i)-root1
 


c ATTENZIONE EFFETTUO UN CAMBIO IMPORTANTE:
      root1=0.
      o(i)=o(i)-root1
 
      o1(i)=1.0-xcat(1)-o(i)/2.0
      o3(i)=(4.0*xcat(1)-o(i))/2.0
      o4(i)=o(i)/(o(i)+o3(i)+xcat(1)+root1)
 
      root(i)=root1
 
      s=dexp(-1.7165*dlog(o4(i))+2.8776)
      sii(i)=s
      s1=acidic(i)/s
      o2=o1(i)/(o1(i)+s1+xwd(i))
c     aossi(i)=o2
      totani(i)=o1(i)+s1+xwd(i)
      afeo2(i)=afeo2(I)/s
      aFeO2(i)=afeo2(i)/totani(i)
      afeo22(i)=afeo22(i)/s
      afeo22(i)=afeo22(i)/totani(i)
c     write (*,*) 's = ',s,' afeo22 = ',afeo22(i)
      ah(i)=totcat(i)/totani(i)
      END
 
      SUBROUTINE ACTION2
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*     computes the activity of sulphide and of iron ionic species   *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),    XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
C
      totold(i)=o1(i)+s1
      struct(i)=s1
      aossi(i)=o2
 
      END
 
      SUBROUTINE FEREDOX
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* computes the redox state of iron in terms of oxides proportions *
* and/or the ensuying oxygen partial pressure at equilibrium      *
* (utilized whenever the data are unsufficient or not precise)    *
c the adopted constants are from Ottonello et al. (2000)          *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),  XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       sumfe(ndi),
     &         reg3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
 
c Pressure effect on cost 54 (cost 21 in published work) accounting
c We assume that for the other constants -equilibrium among components
c in the liquid phase- the Dvolume of reaction is 0.
c calcola la fo2 in base al rapporto di ossidazione del ferro
 
c  cost54 => FeO + 0.25O2 <=> FeO1.5
c  cost57 => Fe2O3 + O2- <=> 2FeO2-
c  cost58 => Fe2O3 <=> 2Fe3+ + 3O2-
c  cost59 => FeO <=> Fe2+ + O2-
 
********* OLD VALUES (Chemical Geology) **************
 
      cost54=-2.8792+6364.8/Tkelv(i)
      cost57=5.3392-3357.4/Tkelv(i)
      cost58=1.8285-4100.2/Tkelv(i)
      cost59=1.1529-1622.4/Tkelv(i)
 
ccc   Now we use cost54 at pressures greater than 1 bar. FeO and FeO1.5
ccc   activities will be corrected for volume terms.
ccc   We use the values given by Lange in Carroll and Webster, 1994
ccc   They have also been recomputed with correlations with optical basicity.
ccc   First, we must recalculate for every T at 1 bar pressure
ccc   Values are in cc/mole.
ccc   These are parameters recommended for metaluminous and peralkaline
ccc   liquid !!!!!!  But our extrapolation makes them good for every melt.
 
c   VFeO2 usando i CR (assumo struttura vincolata e non free anions)
c   0.63 ä HS CR di Shannon per Fe3+ (IR sarebbe 0.49) in coord IV
c   0.55 ä LS IR di Shannon per Fe3+ (HS sarebbe 0.645) in coord VI
c   0.61 ä LS IR di Shannon per Fe2+ (HS sarebbe 0.780) in coord VI
c   1.24 ä CR di Shannon per O2- (IR sarebbe 1.38) in coord IV
c   1.26 ä CR di Shannon per O2- (IR ä 1.40...) in coord VI
      vfeo2m=4./3.*3.14159*(0.49+2*1.40+
     $ 0./10000*(Tkelv(i)-298.15))**3
      vo2m=4./3.*3.14159*((1.40+0./10000*(Tkelv(i)-298.15))**3)
      vfe3p=4./3.*3.14159*((0.645+0./10000*
     $ (Tkelv(i)-298.15))**3)
      vfe2p=4./3.*3.14159*((0.78+0./10000*(Tkelv(i)-298.15))**3)
      conv=0.6022045
      vfeo2m=vfeo2m*conv
      vfeo2m=vfeo2m+0.*0.001*(Tkelv(i)-298.15)
      vo2m=vo2m*conv
      vo2m=vo2m+0.*0.001*(Tkelv(i)-298.15)
c Di seguito equivale ad assumere che DV di FeO2- + 2O2- = FeO45- ä 0.
      vfeo2m=vfeo2m-2*vo2m
      vfe3p=vfe3p*conv
      vfe3p=vfe3p+0.*0.001*(Tkelv(i)-298.15)
      vfe2p=vfe2p*conv
      vfe2p=vfe2p+0.*0.001*(Tkelv(i)-298.15)
c sopra ho ripetuto: x(8)=x(9). Cioä le espansioni di Fe3 e Fe2 sono uguali
c espansione libera: sui volumi molari espressione intera
c     vfeO15=21.065+4.545*0.001*(Tkelv(i)-1673)
c     vfeO=13.65+2.92*0.001*(Tkelv(i)-1673)
c espansione fissata: T = 298.15
      vfeO15=21.065+4.545*0.001*(298.15-1673)
      vfeO=13.65+2.92*0.001*(298.15-1673)
 
      dv57=2*(vfeo2m-vfeo15-0.5*vo2m)
      dv58=2*(vfe3p+1.5*vo2m-vfeo15)
      dv59=vfe2p+vo2m-vfeo
c da scommentare per verificare consistenza con gaslast
c     dv57=0
c     dv58=0
c     dv59=0
 
      cost57=cost57-(dv57*0.1)*(pbar(i)-1)/(tkelv(i)*8.3147*2.303)
      cost58=cost58-(dv58*0.1)*(pbar(i)-1)/(tkelv(i)*8.3147*2.303)
      cost59=cost59-(dv59*0.1)*(pbar(i)-1)/(tkelv(i)*8.3147*2.303)
 
c Therefore:
      cost54=10.**cost54
      cost57=10.**cost57
      cost58=10.**cost58
      cost59=10.**cost59
 
*********************************************************
 
      V0II1 =voxm(7)+eoxm(7)*0.001*(tkelv(i)-1673.)
      V0III1=voxm(5)+eoxm(5)*0.001*(tkelv(i)-1673.)
c     v0iip =v0ii1 +(pbar(i)-1.)*0.0001*(-1.8)
c     v0iiip=(v0iii1+(pbar(i)-1.)*0.0001*3.1)*0.5
 
c     terms above are converted in joule/bar
      v0ii1 =0.1*v0ii1
      v0iii1=0.1*v0iii1
 
      cFeii =coxm(7)*0.0001*0.1
      cfeiii=coxm(5)*0.0001*0.1
 
      dvs = (0.5*v0iii1-v0ii1)*(pbar(i)-1.)+(0.5*cfeiii-cfeii)*
     $(pbar(i)**2./2.-pbar(i)+0.5)
 
c   Dvs becomes dvs/RT
      dvs=dvs/(tkelv(i)*8.3147)
c   Correction term dvs entered below...
 
c serve solo quando non sono dati o sono imprecisi
c     rappox(i)=dlog10(xcat(7)/xcat(5))
      raptru=10.0**redoz(i)
      den1=cost54*raptru
      den1=den1/(dexp(dvs))
      den2=cost57**(0.5)*aossi(i)**2*totani(i)+cost58**(0.5)*totcat(i)
      razio=aossi(i)**(0.5)*cost59*totcat(i)/(den1*den2)
      bossi(i)=razio**4.
      bossi(i)=dlog10(bossi(i))
 
***  Vale nel caso che fO2 non sia data in input e si ha Fe2+/Fe3+ ***
      if (xossi(i).eq.0.and.xftot.gt.0.) xossi(i)=bossi(i)
**************************************************************************
      XOSS =10.0**XOSSI(I)
      den1=cost54*xoss**0.25
      den1=den1/(dexp(dvs))
      den2=cost57**(0.5)*aossi(i)**2*totani(i)+cost58**(0.5)*totcat(i)
      ratio=aossi(i)**(0.5)*cost59*totcat(i)/(den1*den2)
      a3a2(i)=den1
      ca3ca2(i)=aossi(i)**(0.5)*cost59*totcat(i)/den2
c  sotto abbiamo il rapporto [FeO2-]/[Fe2+] come log
      afeo2m(i)=((cost57**0.5)/cost59)*(aossi(i)**1.5)*den1
      afeo2m(i)=dlog10(afeo2m(i))
c     afeo2m(i)=afeo2m(i)*afe22(i)
c  sotto abbiamo il rapporto Fe3+/Fe2+
      afe3m(i)=((cost58**0.5)/cost59)*den1/(aossi(i)**1.5)
      afe3m(i)=dlog10(afe3m(i))
c     afe3m(i)=afe3m(i)*afe22(i)
      redox(i)=dlog10(ratio)
c      write (*,*) 'superstronzo!!!!!'
c      pause
 
      END
 
 
      SUBROUTINE CALCPROP
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* assigns network-former and network-modifier fractions for     *
* sulfidation reactions of type                                 *
*   SiO2 + 0.5S2(g)  <=>  SiOS + 0.5O2(g)                       *
*   TiO2 +0.5S2(g)   <=>  TiOS + 0.5O2(g)                       *
*   P2O5 +0.5S2(g)   <=>  P2O4S + 0.5O2(g)                      *
*   Al2O3 +0.5S2(g)  <=>  Al2O2S + 0.5O2(g)                     *
*   Fe2O3 +0.5S2(g)  <=>  Fe2O2S + 0.5O2(g)                     *
*   Cr2O3 +0.5S2(g)  <=>  Cr2O2S + 0.5O2(g)                     *
*   FeO + 0.5S2(g)   <=>  FeS + 0.5O2(g)                        *
*   MnO + 0.5S2(g)   <=>  MnS + 0.5O2(g)                        *
*   CaO + 0.5S2(g)   <=>  CaS + 0.5O2(g)                        *
*   MgO + 0.5S2(g)   <=>  MgS + 0.5O2(g)                        *
*   Na2O +0.5S2(g)   <=>  Na2S + 0.5O2(g)                       *
*   K2O + 0.5S2(g)   <=>  K2S +  0.5O2(g)                       *
*   H2O + 0.5S2(g)   <=>  H2S +  0.5O2(g)                       *
*   values are normalized to the total number of cations in one *
*   mole of melt                                                *
*                                                               *
*   We HAVE ALSO ANALOGOUS SULPHATATION REACTIONS               *
* * * * * * * * ** * * * * * * * * * ** * * * * * * * * * * * * *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi), XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO, ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK,  ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag
 
      write (*,*)'Composition treated is ',i
c     pause
 
      si=2.0*scat(1)
      p=2.5*scat(3)
      al=1.50*scat(4)
      fe=1.50*scat(5)
      wat=0.5*xwd(i)
      totale=si+p+al+fe+wat
      scat(1)=si/totale
      scat(2)=wat/totale
      scat(3)=p/totale
      scat(4)=al/totale
      scat(5)=fe/totale
      sio2f(i)=(scat(1)*(1.0-totcat(i)))
      po25f(i)=(scat(3)*(1.0-totcat(i)))
      alo15f(i)=(scat(4)*(1.0-totcat(i)))
      feo15f(i)=(scat(5)*(1.0-totcat(i)))
      ho05f(i)=(scat(2)*(1.0-totcat(i)))
C     totfor=(scat(1)+scat(3)+scat(4)+scat(5))*(1.0-totcat(i))
      ti=2.0*xcat(2)
      p=2.5*xcat(3)
      al=1.5*xcat(4)
      fe3=1.5*xcat(5)
      cr=1.5*xcat(6)
      fe2=xcat(7)*1.
      xmn=xcat(8)*1.
      ca=xcat(9)*1.
      xmg=xcat(10)*1.
      xna=xcat(11)*0.5
      xk=xcat(12)*0.5
      xhh=xcat(13)*0.5
      totale=ti+p+al+fe3+cr+fe2+xmn+ca+xmg+xna+xk+xhh
      xcat(2)=ti/totale
      xcat(3)=p/totale
      xcat(4)=al/totale
      xcat(5)=fe3/totale
      xcat(6)=cr/totale
      xcat(7)=fe2/totale
      xcat(8)=xmn/totale
      xcat(9)=ca/totale
      xcat(10)=xmg/totale
      xcat(11)=xna/totale
      xcat(12)=xk/totale
      xcat(13)=xhh/totale
      tio2m(i)=(xcat(2)*totcat(i))
      po25m(i)=(xcat(3)*totcat(i))
      alo15m(i)=(xcat(4)*totcat(i))
      feo15m(i)=(xcat(5)*totcat(i))
      cro15m(i)=(xcat(6)*totcat(i))
      feom(i)=xcat(7)*totcat(i)
      mnom(i)=xcat(8)*totcat(i)
      caom(i)=xcat(9)*totcat(i)
      mgom(i)=xcat(10)*totcat(i)
      nao05m(i)=(xcat(11)*totcat(i))
      ko05m(i)=(xcat(12)*totcat(i))
      ho05m(i)=(xcat(13)*totcat(i))
C     totmod=(xcat(2)+xcat(3)+xcat(4)+xcat(5)+xcat(6)+xcat(7)+
C    $ xcat(8)+xcat(9)+xcat(10)+xcat(11)+xcat(12)+xcat(13))*totcat(i)
 
      END
 
      SUBROUTINE DATAOUT
 
* prints informations on melt composition and structure
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8   mnom,mgom,nao05m,ko05m
      parameter (ndi=250)
      COMMON
     1/DATINT/ Z(13), Y(13) , X(26), ZCAT(13),   XX(26),
     &         voxm(13),eoxm(13),coxm(13),radi(13),noxy(13)
     2/DATEXT/ TCENT(ndi), PBAR(ndi),  POSS(ndi,13), IAUT(ndi),
     &         SO2(ndi),   XOSSI(ndi), FS2(ndi),     NCOMP,wt(ndi)
     3/PARINT/ TKELV(ndi), RAPPOX(ndi),REDOX(ndi),   BOSSI(ndi),
     &         SIO2F(ndi), PO25F(ndi), ALO15F(ndi),  FEO15F(ndi),
     &         TIO2M(ndi), PO25M(ndi), ALO15M(ndi),  FEO15M(ndi),
     &         CRO15M(ndi),FEOM(ndi),  MNOM(ndi),    CAOM(ndi),
     &         MGOM(ndi),  NAO05M(ndi),KO05M(ndi),   HO05M(ndi),
     &         REG(ndi),   XSO2(ndi),  KFLAG(ndi),   SPPM(ndi),
     &         D(ndi),     DD(ndi),    BEPPE(ndi),   REG2(ndi),
     &         REDOZ(ndi), XFETOT(ndi),SCAT(5),      XCAT(13),
     &         S1,HO05f(ndi),         O2,         A,foga(ndi),XFOR,
     &         XFE2O3,     XALLU,      RAPMED,       SUMFE(ndi),
     &         REG3(ndi),  xftot,      jflag,        BEPPI(ndi),
     &         REGG(ndi),  SPPB(ndi),  reg4(ndi),    regg3(ndi),
     &         xh2o(ndi),  a3a2(ndi),  ca3ca2(ndi),refirst(ndi),
     &         som(ndi),somm(ndi),xcat0(ndi),xw(ndi),basdif(ndi),
     &         volion(ndi),volm(ndi),por(ndi),xmol(ndi,13),xxcat(ndi,13)
 
     4/DATOUT/ ACIDIC(ndi),POLCOS(ndi),STRUCT(ndi),  TOTOLD(ndi),
     &         TOTANI(ndi),totali(ndi),TOTCAT(ndi),O4(ndi),      O(ndi),
     &         O1(ndi),    O3(ndi),    S,sii(ndi),            STOT(ndi),
     &         AOSSI(ndi),aozzi(ndi),  SOBS(ndi),  AFE2(ndi),AFE3(ndi),
     &         CSCAL,      CSOBS,      CSC,          CSO,
     &         AFEO2(ndi), STOTT(ndi), CSTART(ndi),  CSTAR(ndi),
     &         CSCAL2,     CSOBS2,     CSC2,         CSO2,
     &         S2MENO(ndi),SO4(ndi),   OXSOLF(ndi),  TOTS(ndi),CS(ndi),
     &         s2k(ndi),so4k(ndi),afeo22(ndi),afeo2m(ndi),
     &         afe32(ndi),afe3m(ndi),afe22(ndi),sulf(ndi),
     &         aossiz(ndi),xh(ndi),xoh(ndi),xohz(ndi),xwd(ndi),ah(ndi),
     &         ana(ndi),aprot(ndi),ak(ndi),amg(ndi),amn(ndi),aca(ndi),
     &         acr(ndi),ati(ndi),ap(ndi),dv(ndi),delv(ndi),root(ndi)
     5/UNIT/   ISYSRD,     ISYSWR,     ISYSPU,       IVIDEO,  ISYSQZ,
     &         ISYSPR,     ISYSIN,     ISYSCK, ISYSWZ, ISYSMF,isysvo
     6/COUNT/  I, kj,jk,iter
     7/FLAG/   jwflag

      WRITE(ISYSWZ,123) tkelv(I),pbar(I),refirst(i),redoz(i),xossi(i),
     $bossi(i)
123   format (f8.2,f8.1,x,4(x,f11.5))


      END
 
