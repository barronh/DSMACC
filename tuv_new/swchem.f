* This file contains the following subroutines, related to specifying 
* chemical spectral weighting functions (cross sections x quantum yields)
*     swphys

*=============================================================================*

      SUBROUTINE swchem(nw,wl,nz,tlev,airden,
     $     j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Load various "weighting functions" (products of cross section and        =*
*=  quantum yield at each altitude and each wavelength).  The altitude       =*
*=  dependence is necessary to ensure the consideration of pressure and      =*
*=  temperature dependence of the cross sections or quantum yields.          =*
*=  The actual reading, evaluation and interpolation is done in separate     =*
*=  subroutines for ease of management and manipulation.  Please refer to    =*
*=  the inline documentation of the specific subroutines for detail          =*
*=  information.                                                             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section * quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw, this_j
      REAL wl(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* local:
      REAL wc(kw)
      INTEGER iw
*_______________________________________________________________________

* complete wavelength grid

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

*____________________________________________________________________________

* O2 + hv -> O + O
* reserve first position.  Cross section parameterization in Schumman-Runge and 
* Lyman-alpha regions are zenith-angle dependent, will be written in 
* subroutine seto2.f.
 
      j = 1
      jlabel(j) = 'O2 -> O + O'

* O3 + hv ->  (both channels)
      CALL r01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* HO2 + hv -> OH + O
      CALL r39(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* H2O2 + hv -> 2 OH
      CALL r08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* NO2 + hv -> NO + O(3P)
      CALL r02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* NO3 + hv ->  (both channels)
      CALL r03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
     
* N2O5 + hv -> (both channels)
      CALL r04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* N2O + hv -> N2 + O(1D)
      CALL r44(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* HNO2 + hv -> OH + NO
      CALL r05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* HNO3 + hv -> OH + NO2
      CALL r06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* HNO4 + hv -> HO2 + NO2
      CALL r07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* NO3-(aq) + hv -> NO2 + O-     (for snow)
* NO3-(aq) + hv -> NO2- + O(3P) (for snow)
      CALL r118(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH2O + hv -> (both channels)
      CALL r10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3CHO + hv -> (all three channels)
      CALL r11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* C2H5CHO + hv -> C2H5 + HCO
      CALL r12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH2(OH)CHO + hv -> Products
      CALL r101(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH2=CHCHO + hv -> Products
      CALL r122(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH2=C(CH3)CHO + hv -> Products
      CALL r104(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3COCH3 + hv -> Products
      CALL r15(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3COCHCH2 + hv -> Products
      CALL r103(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3COCH2CH3 -> CH3CO + CH2CH3
            CALL r119(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH2(OH)COCH3 -> CH3CO + CH2(OH)
* CH2(OH)COCH3 -> CH2(OH)CO + CH3
      CALL r112(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3(OOH) + hv -> CH3O + OH
      CALL r16(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* HOCH2OOH -> HOCH2O. + OH
      CALL r121(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3CO(OOH) + hv -> Products
      CALL r123(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3(ONO2) + hv -> CH3O + NO2
      CALL r17(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3CH2(ONO2) -> CH3CH2O + NO2
      CALL r106(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3CH(ONO2)CH3 -> CH3CHOCH3 + NO2
      CALL r107(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2
      CALL r108(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
      CALL r109(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
      CALL r110(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* PAN + hv -> CH3CO(OO) + NO2
* PAN + hv -> CH3CO(O) + NO3
      CALL r18(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3CH2COO2NO2 -> CH3CH2CO(OO) + NO2
* CH3CH2COO2NO2 -> CH2CH2CO(O) + NO3
            CALL r120(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CHOCHO + hv -> Products
      CALL r13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3COCHO + hv -> Products
      CALL r14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3COCOCH3 + hv -> Products
      CALL r102(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3COCO(OH) + hv -> Products
      CALL r105(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* Cl2 + hv -> Cl + Cl
      CALL r47(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* ClOOCl -> Cl + ClOO
      CALL r111(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* ClOO + hv -> Products
      CALL r31(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* ClONO2 + hv -> Products
      CALL r45(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* Br2 -> Br + Br
      CALL r115(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* BrO -> Br + O
      CALL r114(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* HOBr -> OH + Br
      CALL r113(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* BrONO2 + hv -> Products
      CALL r46(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3Cl + hv -> Products
      CALL r30(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CCl2O + hv -> Products
      CALL r19(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CCl4 + hv -> Products
      CALL r20(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CClFO + hv -> Products
      CALL r21(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CCF2O + hv -> Products
      CALL r22(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF2ClCFCl2 (CFC-113) + hv -> Products
      CALL r23(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF2ClCF2Cl (CFC-114) + hv -> Products
      CALL r24(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF3CF2Cl (CFC-115) + hv -> Products
      CALL r25(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CCl3F (CFC-111) + hv -> Products
      CALL r26(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CCl2F2 (CFC-112) + hv -> Products
      CALL r27(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3CCl3 + hv -> Products
      CALL r29(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF3CHCl2 (HCFC-123) + hv -> Products
      CALL r32(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF3CHFCl (HCFC-124) + hv -> Products
      CALL r33(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3CFCl2 (HCFC-141b) + hv -> Products
      CALL r34(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3CF2Cl (HCFC-142b) + hv -> Products
      CALL r35(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF3CF2CHCl2 (HCFC-225ca) + hv -> Products
      CALL r36(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF2ClCF2CHFCl (HCFC-225cb) + hv -> Products
      CALL r37(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CHClF2 (HCFC-22) + hv -> Products
      CALL r38(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3Br + hv -> Products
      CALL r28(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CHBr3 + hv -> Products
      CALL r09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF3Br (Halon-1301) + hv -> Products
      CALL r42(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF2BrCF2Br (Halon-2402) + hv -> Products
      CALL r43(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF2Br2 (Halon-1202) + hv -> Products
      CALL r40(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CF2BrCl (Halon-1211) + hv -> Products
      CALL r41(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

******************************************
* Additions from the DSMACC version of TUV
* nC3H7CHO -> nC3H7 + HCO
* nC3H7CHO -> C2H4 + CH3CHO
      CALL rm116(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* iC3H7CHO -> iC3H7 + HCO
      CALL rm117(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* nC4H9CHO -> nC4H9 + HCO
* nC4H9CHO -> C3H6 + CH3CHO
      CALL rm118(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* iC4H9CHO -> iC4H9 + HCO
* iC4H9CHO -> C3H6 + CH3CHO
      CALL rm119(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* tC4H9CHO -> tC4H9 + HCO
      CALL rm120(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* nC5H11CHO -> nC5H11 + HCO
* nC5H11CHO -> C4H8 + CH3CHO
      CALL rm121(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* nC6H13CHO -> nC6H13 + HCO
* nC6H13CHO -> C5H10 + CH3CHO
      CALL rm122(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3COC2H5 -> C2H5 + CH3CO
      CALL rm123(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* pinonaldehyde -> products
      CALL rm124(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH2CHCHO -> CH2CHCO + H
* CH2CHCHO -> CH2CH + HCO
      CALL rm125(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* nC3H7ONO2 -> nC3H7O + NO2
      CALL rm126(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* nC4H9ONO2 -> nC4H9O + NO2
      CALL rm127(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* iC4H9ONO2 -> iC4H9O + NO2
      CALL rm128(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* 2-C4H9ONO2 -> 2-C4H9O + NO2
      CALL rm129(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* nC5H11ONO2 -> nC5H11O + NO2
      CALL rm130(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* 2-C5H11ONO2 -> 2-C5H11O + NO2
      CALL rm131(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* 3-C5H11ONO2 -> 3-C5H11O + NO2
      CALL rm132(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* i-C5H11ONO2 -> i-C5H11O + NO2
      CALL rm133(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH3CH2C(O)CH2(ONO2) -> CH3CH2C(O)CH2(O.) + NO2
      CALL rm134(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH3CH(ONO2)C(O)CH3 -> CH3CH(O.)C(O)CH3 + NO2
      CALL rm135(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH3CH(ONO2)CH2(ONO2) -> CH3CH(O.)CH2(ONO2) + NO2
* CH3CH(ONO2)CH2(ONO2) -> CH3CH(ONO2)CH2(O.) + NO2
      CALL rm136(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* C2H5CH(ONO2)CH2(ONO2) -> C2H5CH(O.)CH2(ONO2) + NO2
* C2H5CH(ONO2)CH2(ONO2) -> C2H5CH(ONO2)CH2(O.) + NO2
      CALL rm137(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* CH3CH(ONO2)CH(ONO2)CH3 -> CH3CH(O.)CH(ONO2)CH3 + NO2           
      CALL rm138(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH2=CHCH(ONO2)CH2(ONO2) -> C2H3CH(O.)CH2(ONO2) + NO2
* CH2=CHCH(ONO2)CH2(ONO2) -> C2H3CH(ONO2)CH2(O.) + NO2
      CALL rm139(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH2(ONO2)C2H2CH2(ONO2) -> CH2(ONO2)C2H2CH2(O.) + NO2
      CALL rm140(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH3CH=CHCHO -> CH3CH=CHCO + H
* CH3CH=CHCHO -> CH3CH=CH + HCO
* CH3CH=CHCHO -> CH3CH=CH2 + CO
      CALL rm141(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH3COC3H7 -> C3H7 + CH3CO
      CALL rm142(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* C2H5COC2H5 -> C2H5 + C2H5CO
      CALL rm143(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* ((CH3)2CH)2C(O) + hv -> (CH3)2C + (CH3)2CHCO
      CALL rm144(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH3C(O)CH2CH(CH3)2 + hv -> CH3CO + (CH3)2CHCH2
* CH3C(O)CH2CH(CH3)2 + hv -> CH3C(O)CH3 + CH2=CHCH3
      CALL rm145(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH3C(O)C3H5(CH3)2 -> CH3CO + (CH3)2C3H5
* CH3C(O)C3H5(CH3)2 -> CH3C(O)CH3 + CH2=CHCH(CH3)2
      CALL rm146(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      
* CH3O2NO2 -> CH3O2 + NO2
* CH3O2NO2 -> CH3O + NO3
      CALL rm147(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)     
****************************************************************

      IF (j .GT. kj) STOP '1002'
      RETURN
      END
