/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolCalc compute calcium cooling */
#include "cddefines.h"
#include "taulines.h"
#include "doppvel.h"
#include "phycon.h"
#include "ca.h"
#include "dense.h"
#include "thermal.h"
#include "opacity.h"
#include "rfield.h"
#include "ligbar.h"
#include "lines_service.h"
#include "atoms.h"
#include "cooling.h"
#include "iso.h"

void CoolCalc(void)
{
	double  a21, 
	  a31, 
	  a41, 
	  a42, 
	  a51, 
	  a52, 
	  a53, 
	  c21, 
	  Ca2pop[5] ,
	  cs,
	  cs14, 
	  cs15, 
	  d41, 
	  d42, 
	  d51, 
	  d52, 
	  d53, 
	  hlgam, 
	  op41, 
	  op51, 
	  opckh, 
	  opcxyz, 
	  PhotoRate2,
	  PhotoRate3, 
	  PhotoRate4, 
	  PhotoRate5, 
	  r21, 
	  r31, 
	  r41, 
	  r42, 
	  r51, 
	  r52, 
	  r53;
	static double gCa2[5]={2.,4.,6.,2.,4.};
	static double exCa2[4]={13650.2,60.7,11480.6,222.9};
	static realnum opCax = 0.f;
	static realnum opCay = 0.f;
	static realnum opCaz = 0.f;

	DEBUG_ENTRY( "CoolCalc()" );

	/* Ca I 4228 */
	MakeCS(TauLines[ipCaI4228]);
	atom_level2(TauLines[ipCaI4228]);

	/* photoionization of evcited levels by Ly-alpha */
	hlgam = rfield.otslin[ iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).ipCont() -1];
	PhotoRate5 = 1.7e-18*hlgam;
	PhotoRate4 = 8.4e-19*hlgam;
	PhotoRate3 = 7.0e-18*hlgam;
	PhotoRate2 = 4.8e-18*hlgam;

	/* spontaneous decays
	 * frist two trans prob from 
	 * >>refer	Ca2	AS	Zeippen, C.J. 1990, A&A, 229, 248 */
	a21 = 1.02*TauLines[ipT7324].Emis().Pesc();
	a31 = 1.05*TauLines[ipT7291].Emis().Pesc();
	a41 = 1.4e8*TauLines[ipT3969].Emis().Pesc();
	a51 = 1.4e8*TauLines[ipT3934].Emis().Pesc();
	a42 = 7.9e6*TauLines[ipT8662].Emis().Pesc();
	a52 = 8.2e5*TauLines[ipT8498].Emis().Pesc();
	a53 = 7.48e6*TauLines[ipT8542].Emis().Pesc();

	/* destruction of IR triplet by continuous opacities */
	opcxyz = opac.opacity_abs[ TauLines[ipT7324].ipCont() -1];

	/* opcxyz = opac(icaxyz) */
	if( opcxyz > 0. )
	{
		d52 = 5.6*opcxyz/(opcxyz + opCax)*(1. - TauLines[ipT8498].Emis().Pesc());
		d53 = 5.6*opcxyz/(opcxyz + opCay)*(1. - TauLines[ipT8542].Emis().Pesc());
		d42 = 5.6*opcxyz/(opcxyz + opCaz)*(1. - TauLines[ipT8662].Emis().Pesc());
	}
	else
	{
		d52 = 0.;
		d53 = 0.;
		d42 = 0.;
	}

	/* near UV dest of KH by background continuum */
	opckh = opac.opacity_abs[ TauLines[ipT3969].ipCont() -1];

	/* opckh = opac(icakh) */
	if( opckh > 0. )
	{
		op51 = dense.xIonDense[ipCALCIUM][1]*3.89e-7/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]);
		d51 = 5.6*opckh/(opckh + op51);
		op41 = dense.xIonDense[ipCALCIUM][1]*1.96e-7/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]);
		d41 = 5.6*opckh/(opckh + op41);
	}
	else
	{
		op51 = 0.;
		d51 = 0.;
		op41 = 0.;
		d41 = 0.;
	}
	/* net rates */
	r21 = PhotoRate2 + a21;
	r31 = PhotoRate3 + a31;
	r41 = a41 + PhotoRate4 + d41;
	r51 = a51 + PhotoRate5 + d51;
	r42 = a42 + d42;
	r52 = a52 + d52;
	r53 = a53 + d53;
	cs14 = 0.923*phycon.te10*phycon.te10;
	cs15 = cs14*2.;
	TauLines[ipT3969].Coll().col_str() = (realnum)cs14;
	TauLines[ipT3934].Coll().col_str() = (realnum)cs15;

	/* following used to correct rec contribution
	 * fcakh = a51 / ( a51 + eden*1.5e-5 / sqrte ) 
	 * cs 1-2 from 
	 * >>refer	Ca2	CS	Saraph, H.E. 1970, J. Phys. B, 3, 952
	 * other 
	 * >>refer	Ca2	CS	Chidichimo, M.C. 1981, J. Phys. B, 14, 4149 */
	double Cooling , CoolingDeriv;
	atom_pop5(gCa2,exCa2,5.8,8.6,cs14,cs15,20.6,22.9,9.8,3.4,44.4,1.0,
		r21,r31,r41,r51,0.,r42,r52,0.,r53,0.,Ca2pop,
		dense.xIonDense[ipCALCIUM][1],&Cooling , &CoolingDeriv, 0.,0.,0.,0.);

	/* CDSQTE = 8.629E-6*EDEN/SQRTE */
	c21 = 5.8/4.*dense.cdsqte;

	/* remember largest ratio of Ly-al removal to total */
	if( dense.xIonDense[ipCALCIUM][1] > 0. )
		ca.Ca2RmLya = MAX2(ca.Ca2RmLya,(realnum)(PhotoRate2/(PhotoRate2+a21+c21)));

	ca.Cak = (realnum)(Ca2pop[4]*a51*5.06e-12);
	ca.Cah = (realnum)(Ca2pop[3]*a41*5.01e-12);
	ca.Cax = (realnum)(Ca2pop[4]*a52*2.34e-12);
	ca.Cay = (realnum)(Ca2pop[4]*a53*2.33e-12);
	ca.Caz = (realnum)(Ca2pop[3]*a42*2.30e-12);
	ca.Caf1 = (realnum)(Ca2pop[2]*a31*2.73e-12);
	ca.Caf2 = (realnum)(Ca2pop[1]*a21*2.72e-12);
	ca.popca2ex = (realnum)(Ca2pop[1] + Ca2pop[2] + Ca2pop[3] + Ca2pop[4]);

	/* this is the total cooling due to the model atom */
	ca.Cair = ca.Cax + ca.Cay + ca.Caz;
	ca.c7306 = ca.Caf1 + ca.Caf2;
	ca.Cakh = ca.Cak + ca.Cah;

	// total cooling from 5-level atom
	CoolAdd("Ca 2",7306,Cooling);
	thermal.dCooldT += CoolingDeriv;

	/*fprintf(ioQQQ,"DEBUG ca2\t%.2f\t%.5e\t%.4e\t%.4e\n",
		fnzone, phycon.te,ca.Cakh,dense.xIonDense[ipCALCIUM][1]);*/

	/* level populations that will be used for excited state photoionization */
	ca.dstCala = (realnum)(Ca2pop[4]*PhotoRate5 + Ca2pop[3]*PhotoRate4);
	ca.dCakh = (realnum)(ca.dstCala*5.03e-12);
	ca.dCaf12 = (realnum)((Ca2pop[2]*PhotoRate3 + Ca2pop[1]*PhotoRate2)*2.31e-12);
	opCax = (realnum)(Ca2pop[1]*1.13e-8/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]));
	opCay = (realnum)(Ca2pop[2]*6.87e-8/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]));
	opCaz = (realnum)(Ca2pop[1]*5.74e-8/GetDopplerWidth(dense.AtomicWeight[ipCALCIUM]));

	/* total rate Lalpha destroys CaII,
	 * this is only used in ioncali to increase ionization rate by
	 * adding it to the ct vector */
	if( dense.xIonDense[ipCALCIUM][1] > 0. )
	{
		ca.dstCala = (realnum)(
			(ca.dstCala + ca.dCaf12/2.31e-12)/dense.xIonDense[ipCALCIUM][1]);
		{
			/*@-redef@*/
			enum {DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ," hlgam is %e\n", hlgam);
			}
		}
	}
	else
	{
		ca.dstCala = 0.;
	}
	ca.Ca3d = (realnum)(Ca2pop[1] + Ca2pop[2]);
	ca.Ca4p = (realnum)(Ca2pop[3] + Ca2pop[4]);

	/* incl stimulated emission for Calcium II 5-level atom */
	TauLines[ipT3934].Emis().PopOpc() = (Ca2pop[0] - Ca2pop[4]/2.);
	(*TauLines[ipT3934].Hi()).Pop() = Ca2pop[4];
	(*TauLines[ipT3934].Lo()).Pop() = Ca2pop[0];
	TauLines[ipT3969].Emis().PopOpc() = (Ca2pop[0] - Ca2pop[3]);
	(*TauLines[ipT3969].Hi()).Pop() = Ca2pop[3];
	(*TauLines[ipT3969].Lo()).Pop() = Ca2pop[0];

	TauLines[ipT8498].Emis().PopOpc() = (Ca2pop[1] - Ca2pop[4]);
	(*TauLines[ipT8498].Hi()).Pop() = Ca2pop[4];
	(*TauLines[ipT8498].Lo()).Pop() = Ca2pop[1];
	TauLines[ipT8542].Emis().PopOpc() = (Ca2pop[2] - Ca2pop[4]*1.5);
	(*TauLines[ipT8542].Hi()).Pop() = Ca2pop[4];
	(*TauLines[ipT8542].Lo()).Pop() = Ca2pop[2];
	TauLines[ipT8662].Emis().PopOpc() = (Ca2pop[1] - Ca2pop[3]*2.);
	(*TauLines[ipT8662].Hi()).Pop() = Ca2pop[3];
	(*TauLines[ipT8662].Lo()).Pop() = Ca2pop[1];
	TauLines[ipT7291].Emis().PopOpc() = dense.xIonDense[ipCALCIUM][1];
	(*TauLines[ipT7291].Hi()).Pop() = 0.;
	(*TauLines[ipT7291].Lo()).Pop() = dense.xIonDense[ipCALCIUM][1];
	TauLines[ipT7324].Emis().PopOpc() = dense.xIonDense[ipCALCIUM][1];
	(*TauLines[ipT7324].Hi()).Pop() = 0.;
	(*TauLines[ipT7324].Lo()).Pop() = dense.xIonDense[ipCALCIUM][1];

	/* Ca IV 3.2 micron; data from
	 * >>refer	Ca4	AS	Mendoza, C. 1982, in IAU Symp. 103, Planetary
	 * >>refercon	Nebulae, ed. D.R. Flower, (Dordrecht, Holland: D. Reidel Publishing Co.), 143
	 * Y(ik) from 
	 * >>refer	Ca4	CS	Pelan, J., & Berrington, K. A. 1995, A&AS, 110, 209 */
	if( phycon.te <= 1e5 )
	{
		cs = MAX2(1.0,8.854e-3*phycon.sqrte);
	}
	else if( phycon.te < 2.512e5 )
	{
		cs = 2.8;
	}
	else
	{
		cs = 641.1/(phycon.te30*phycon.te10*phycon.te02*phycon.te02/
		  phycon.te003);
	}
	PutCS(cs,TauLines[ipTCa3]);
	atom_level2(TauLines[ipTCa3]);

	return;
}
