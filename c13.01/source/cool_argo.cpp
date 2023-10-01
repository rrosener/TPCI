/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolArgo compute argon cooling */
#include "cddefines.h"
#include "coolheavy.h"
#include "phycon.h"
#include "ligbar.h"
#include "taulines.h"
#include "dense.h"
#include "thermal.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"

void CoolArgo(void)
{
	realnum a12, 
	  a13, 
	  a14, 
	  a15, 
	  a23, 
	  a24, 
	  a25, 
	  a34, 
	  a35, 
	  a45, 
	  cs12, 
	  cs13, 
	  cs14, 
	  cs15, 
	  cs23, 
	  cs24, 
	  cs25, 
	  cs34, 
	  cs35, 
	  cs45, 
	  pop2, 
	  popn3;
	double cs, p[5];
	static double gAr4[5]={4.,4.,6.,2.,4.};
	static double exAr4[4]={21090.4,128.9,13636.2,177.1};

	DEBUG_ENTRY( "CoolArgo()" );

	/* Argon II 6.98 micron
	 * trans prob from 
	 * >>refer	Ar2	AS	Nussbaumer, H., & Storey, P. J. 1988, A&A, 200, L25
	 * Y(ik) from
	 * >>refer	Ar2	CS	Pelan, J., & Berrington, K. A. 1995, A&AS, 110, 209 */
	PutCS(3.1,TauLines[ipTAr7]);
	atom_level2(TauLines[ipTAr7]);

	/* A III 7136, 7751, 3109, 5192 CS
	 * >>refer	Ar3	CS	Galavis, M. E., Mendoza, C., & Zeippen, C. J. 1995, A&AS, 111, 347
	 * >>chng 97 mar 19, break up into three level atom */
	cs12 = 4.825f;
	cs13 = 0.841f;
	cs23 = (realnum)MIN2(1.30,3.296/(phycon.te10*phycon.te01/phycon.te001/
	  phycon.te001));
	a12 = 0.3963f;
	a23 = 2.59f;
	a13 = 3.952f;
	/* POP3(G1,G2,G3,O12,O13,O23,A21,A31,A32,E12,E23,P2,ABUND,GAM2) */
	popn3 = (realnum)(atom_pop3(9.,5.,1.,cs12,cs13,cs23,a12,a13,a23,1.955e4,2.770e4,
	  &pop2,dense.xIonDense[ipARGON][2],0.,0.,0.));
	CoolHeavy.c7136 = pop2*a12*2.7e-12;
	CoolHeavy.c5192 = popn3*a23*3.83e-12;
	CoolHeavy.c3109 = popn3*a13*6.40e-12;
	CoolAdd("Ar 3",7136,CoolHeavy.c7136);
	CoolAdd("Ar 3",5192,CoolHeavy.c5192);
	CoolAdd("Ar 3",3109,CoolHeavy.c3109);
	/* add to deriv */
	/* >>chng 01 mar 10, did not have second pair of lines included */
	thermal.dCooldT += CoolHeavy.c7136*1.955e4*thermal.tsq1 + 
		(CoolHeavy.c5192+CoolHeavy.c3109)*4.73e4*thermal.tsq1;

	/* Ar III 21.8(J=0,1), 9.0 (J=1,2) mircon lines,
	 * >>refer	Ar3	CS	Galavis, M. E., Mendoza, C., & Zeippen, C. J. 1995, A&AS, 111, 347 */
	PutCS(3.1,TauLines[ipTAr9]);
	cs = MIN2(1.384,3.110/(phycon.te10/phycon.te001*phycon.te001));
	PutCS(cs,TauLines[ipTAr22]);
	if( phycon.te < 1e4 )
	{
		cs = 0.671;
	}
	else
	{
		cs = MIN2(0.906,0.150*phycon.te20/phycon.te02/phycon.te02*
		  phycon.te003);
	}
	PutCS(cs,*TauDummy);
	atom_level3(TauLines[ipTAr9],TauLines[ipTAr22],*TauDummy);

	/* Argon IV 4711+4740, 7335 lines (O II like)
	 * CS from
	 * >>refer	Ar4	CS	Zeippen, C. J., Le Bourlot, J., & Butler, K. 1987, A&A, 188, 251
	 * >>chng 97 jan 31, increase to full 5 level atom
	 * Ar IV, cs data from 
	 * >>refer	Ar4	CS	Ramsbottom, C. A., Bell., K. L., & Keenan, F. P. 1997, MNRAS 284, 754
	 * differs by 2-3x from older values
	 * temp dependence form 
	 *>>refer	Ar4	CS	Ramsbottom, C. A., & Bell, K. L. 1997, At. Data Nucl. Data Tables, 66, 65 */
	cs12 = (realnum)MAX2(0.761,0.481*phycon.te05);
	cs12 = (realnum)MIN2(0.853,cs12);
	a12 = 2.23e-2f;

	cs13 = (realnum)MAX2(1.14,0.719*phycon.te05);
	cs13 = (realnum)MIN2(1.3,cs13);
	a13 = 1.77e-3f;

	cs14 = (realnum)MAX2(0.39,0.108*phycon.te10*phycon.te02*phycon.te02);
	cs14 = (realnum)MIN2(0.5,cs14);
	a14 = 0.862f;

	cs15 = (realnum)MAX2(0.78,0.216*phycon.te10*phycon.te02*phycon.te02);
	cs15 = (realnum)MIN2(1.0,cs15);
	a15 = 2.11f;

	cs23 = 7.06f;
	a23 = 2.30e-5f;

	cs24 = (realnum)MAX2(1.53,0.346*phycon.te10*phycon.te05);
	cs24 = (realnum)MIN2(1.96,cs24);
	a24 = 0.603f;

	cs25 = (realnum)MAX2(2.18,0.664*phycon.te10*phycon.te02);
	cs25 = (realnum)MIN2(2.65,cs25);
	a25 = 0.789f;

	cs34 = (realnum)MAX2(1.56,0.475*phycon.te10*phycon.te02);
	cs34 = (realnum)MIN2(1.89,cs34);
	a34 = 0.119f;

	cs35 = (realnum)MAX2(4.01,1.00*phycon.te10*phycon.te02*phycon.te02);
	cs35 = (realnum)MIN2(5.03,cs35);
	a35 = 0.598f;

	cs45 = (realnum)(0.0359*phycon.te20*phycon.te20*phycon.te02*phycon.te02);
	a45 = 4.94e-5f;

	/* FIVEL( G(1-5) , ex(wn,1-5), cs12,cs13,14,15,23,24,25,34,35,45,
	 *  A21,31,41,51,32,42,52,43,53,54, pop(1-5), abund) */
	double Cooling;
	double CoolingDeriv;
	atom_pop5(gAr4,exAr4,cs12,cs13,cs14,cs15,cs23,cs24,cs25,cs34,cs35,
	  cs45,a12,a13,a14,a15,a23,a24,a25,a34,a35,a45,p,
	  dense.xIonDense[ipARGON][3],&Cooling , &CoolingDeriv, 0.,0.,0.,0.);
	CoolHeavy.Ar4740 = p[1]*a12*4.20e-12;
	CoolHeavy.Ar4711 = p[2]*a13*4.20e-12;
	CoolHeavy.Ar2868 = p[3]*a14*6.94e-12;
	CoolHeavy.Ar2854 = p[4]*a15*6.94e-12;
	CoolHeavy.Ar7263 = p[3]*a24*2.74e-12;
	CoolHeavy.Ar7171 = p[4]*a25*2.74e-12;
	CoolHeavy.Ar7331 = p[3]*a34*2.74e-12;
	CoolHeavy.Ar7237 = p[4]*a35*2.74e-12;

	// total cooling from 5-level atom
	CoolAdd("Ar 4",4740,Cooling);
	thermal.dCooldT +=  CoolingDeriv;

	/* Argon V 6435+7007,  
	 * >>refer	Ar5	AS	Mendoza, C., & Zeippen, C. J. 1982, MNRAS, 199, 1025
	 * >>refer	Ar5	CS	Galavis, M. E., Mendoza, C., & Zeippen, C. J. 1995, A&AS, 111, 347
	 * POPEXC( O12,g1,g2,A21,excit,abund); result already*a21, excit in Kelvin */
	if( phycon.te < 1e4 )
	{
		cs12 = 3.09f;
	}
	else
	{
		cs12 = (realnum)MIN2(4.454,0.634*phycon.te20/phycon.te03*phycon.te001*
		  phycon.te001);
	}
	cs13 = 0.56f;
	cs23 = 1.65f;
	a12 = 0.68f;
	a13 = 6.55f;
	a23 = 3.35f;

	/* >>chng 01 mar 10, convert from 2 to 3 level atom */
	popn3 = (realnum)(atom_pop3(9.,5.,1.,cs12,cs13,cs23,a12,a13,a23,2.055e4,3.110e4,
	  &pop2,dense.xIonDense[ipARGON][4],0.,0.,0.));
	CoolHeavy.c7007 = pop2*a12*2.84e-12;
	CoolHeavy.c4626 = popn3*a23*4.30e-12;
	CoolHeavy.c2691 = popn3*a13*7.39e-12;

	CoolAdd("Ar 5",7007,CoolHeavy.c7007);
	CoolAdd("Ar 5",4626,CoolHeavy.c4626);
	CoolAdd("Ar 5",2691,CoolHeavy.c2691);

	/* add to deriv */
	/* >>chng 01 mar 10, did not have second pair of lines included */
	thermal.dCooldT += CoolHeavy.c7007*2.055e4*thermal.tsq1 + 
		(CoolHeavy.c4626+CoolHeavy.c2691)*5.17e4*thermal.tsq1;

	/* Ar V 3P fine structure lines , A from
	 * >>refer	Ar5	AS	Mendoza, C. 1982, in IAU Symp. 103, Planetary 
	 * >>refercon	Nebulae, ed. D. R. Flower, (Dordrecht, Holland: D. Reidel Publishing Co.), 143
	 * >>refer	Ar5	CS	Galavis, M. E., Mendoza, C., & Zeippen, C. J. 1995, A&AS, 111, 347 */
	cs = MIN2(3.26,26.27/(phycon.te20*phycon.te03*phycon.te01*
	  phycon.te003));
	PutCS(cs,TauLines[ipTAr13]);

	cs = MIN2(8.47,44.31/(phycon.te20/phycon.te003/phycon.te003));
	PutCS(cs,TauLines[ipTAr8]);

	cs = MIN2(1.95,7.280/(phycon.te10*phycon.te05*phycon.te005));
	PutCS(cs,*TauDummy);

	atom_level3(TauLines[ipTAr13],TauLines[ipTAr8],*TauDummy);

	/* [Ar VI] 4.53 micron, cs from 
	 * >>refer	Ar6	CS	Saraph, H. E., & Storey, P. J. 1996, A&AS, 115, 151
	 * >>chng 96 dec 11, cs should have been 6.33, caught by Simon Casassus */
	PutCS(6.33,TauLines[ipAr06453]);

	atom_level2(TauLines[ipAr06453]);

	return;
}
