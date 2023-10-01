/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_lv1_k_zn place lines of elements potassium and heavier into lines storage stack */
#include "cddefines.h"
#include "cddrive.h"
#include "coolheavy.h"
#include "ca.h"
#include "fe.h"
#include "rfield.h"
#include "dense.h"
#include "phycon.h"
#include "radius.h"
#include "taulines.h"
#include "trace.h"
#include "lines_service.h"
#include "rt.h"
#include "atomfeii.h"
#include "lines.h"
#include "iso.h"

void lines_lv1_k_zn(void)
{
	long int i, 
	  ipnt,
	  ilo,
	  ihi;

	double eff, fela;

	DEBUG_ENTRY( "lines_lv1_k_zn()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_lv1_k_zn called\n" );
	}

	PutLine(TauLines[ipKI7745],
		"  potassium K I 7745 ");

	PutLine(TauLines[ipxK03462],
		"  [K III] 4.62 microns ");

	PutLine(TauLines[ipxK04598],
		"  [KIV] 5.983 min  ");

	PutLine(TauLines[ipxK04154],
		"  [KIV] 15.39 micron ");

	PutLine(TauLines[ipxK07319],
		"  [K VII] 3.189 microns ");

	PutLine(TauLines[ipCaI4228],
		"  calcium Ca I 4228 ");

	linadd(ca.Cakh,3933,"Ca 2",'i',
		" coll excited calcium k+h " );

	linadd(ca.Cair,8579,"Ca 2",'i' ,
		" infrared triplet ");

	linadd(ca.c7306,7306,"Ca 2",'i',
		" forbidden lines, 7291+7324 together " );

	linadd(ca.dCakh,3933,"Phot",'i' ,
		" fraction H Ly-alpha destruction of excited levels ");

	linadd(ca.dCaf12,7306,"Phot",'i' ,
		" fraction H Ly-alpha destruction of excited levels ");

	PntForLine(3934.,"Ca2K",&ipnt);
	lindst(ca.Cak,3934,"Ca2K",ipnt,'c',true,
		" individual lines from five level atom");


	PntForLine(3969.,"Ca2H",&ipnt);
	lindst(ca.Cah,3969,"Ca2H",ipnt,'c',true,
		" individual lines from five level atom" );


	PntForLine(8498.,"Ca2X",&ipnt);
	lindst(ca.Cax,8498,"Ca2X",ipnt,'c',true,
		" individual lines from five level atom " );


	PntForLine(8542.,"Ca2Y",&ipnt);
	lindst(ca.Cay,8542,"Ca2Y",ipnt,'c',true,
		"  individual lines from five level atom" );


	PntForLine(8662.,"Ca2Z",&ipnt);
	lindst(ca.Caz,8662,"Ca2Z",ipnt,'c',true,
		" individual lines from five level atom" );


	PntForLine(7291.,"CaF1",&ipnt);
	lindst(ca.Caf1,7291,"CaF1",ipnt,'c',true,
		" individual lines from five level atom" );


	PntForLine(7324.,"CaF2",&ipnt);
	lindst(ca.Caf2,7324,"CaF2",ipnt,'c',true,
		" individual lines from five level atom" );

	eff = dense.eden*dense.xIonDense[ipCALCIUM][2]*5.4e-21/(phycon.te/
	  phycon.te10/phycon.te10);
	linadd(eff,3933,"Rec ",'i',
		" recombination contribution to CaII emission" );

	PutLine(TauLines[ipTCa3],
		"  Ca IV 3.2 micron ");


	PntForLine(22.08e4,"Sc 2",&ipnt);
	lindst(CoolHeavy.Sc22p08m,22.08e4,"Sc 2",ipnt,'c',true,
		" Sc II 2.08 (1-3) " );


	PntForLine(24.1e4,"Sc 2",&ipnt);
	lindst(CoolHeavy.Sc24p1m,24.1e4,"Sc 2",ipnt,'c',true,
		" Sc II 4.1 micron (1-2)" );


	PntForLine(24.2e4,"Sc 2",&ipnt);
	lindst(CoolHeavy.Sc24p2m,24.2e4,"Sc 2",ipnt,'c',true,
		"  Sc II 4.22 (2-3)" );


	PntForLine(3933.,"Sc 3",&ipnt);
	lindst(CoolHeavy.Sc33936,3933,"Sc 3",ipnt,'c',true,
		" Sc III 3936" );

	PutLine(TauLines[ipSc05231],
		"  [Sc V] 1.46 microns ");


	PntForLine(5054.,"Sc 6",&ipnt);
	lindst(CoolHeavy.Sc45058,5054,"Sc 6",ipnt,'c',true ,
		" Sc VI 5054 (1-2)");


	PntForLine(3592.,"Sc 6",&ipnt);
	lindst(CoolHeavy.Sc43595,3592,"Sc 6",ipnt,'c',true,
		" Sc VI 3595 (2-3)" );


	PntForLine(2100.,"Sc 6",&ipnt);
	lindst(CoolHeavy.Sc42100,2100,"Sc 6",ipnt,'c',true,
		"  Sc VI 2100 (1-3)" );

	PutLine(TauLines[ipSc13264],
		"  [Sc 13] 2637.97A");

	PntForLine(8823.,"V  3",&ipnt);
	lindst(CoolHeavy.V38830,8823,"V  3",ipnt,'c',true ,
		"  V III 8823 ");


	PntForLine(8507.,"V  3",&ipnt);
	lindst(CoolHeavy.V38507,8507,"V  3",ipnt,'c',true,
		"  V III 8507" );

	PntForLine(5828.,"Cr 3",&ipnt);
	lindst(CoolHeavy.Cr3l21,5828,"Cr 3",ipnt,'c',true,
		" [CrIII] multiplet blend at 5828A" );

	PntForLine(7267.,"Cr 4",&ipnt);
	lindst(CoolHeavy.Cr4l21,7267,"Cr 4",ipnt,'c',true,
		" [CrIV] 2 - 1 multiplet blend at 7272" );


	PntForLine(6801.,"Cr 4",&ipnt);
	lindst(CoolHeavy.Cr4l31,6801,"Cr 4",ipnt,'c',true,
		" [CrIV] 3 - 1 multiplet blend at 6806" );


	PntForLine(7979.,"Cr 5",&ipnt);
	lindst(CoolHeavy.Cr5l21,7979,"Cr 5",ipnt,'c',true,
		" [CrV] 2 - 1 multiplet blend at 7985" );

	PntForLine(6577.,"Cr 5",&ipnt);
	lindst(CoolHeavy.Cr5l31,6577,"Cr 5",ipnt,'c',true,
		"  [CrV] 3 - 1 multiplet blend at 6582" );


	PntForLine(3.75e4,"Cr 5",&ipnt);
	lindst(CoolHeavy.Cr5l32,3.75e4,"Cr 5",ipnt,'c',true,
		" [CrV] 3 - 2 multiplet blend at 3.75 microns " );

	/* bob Rubin's UV line
	 * f2 = dense.xIonDense(26,4)*sexp(50 764./te)*0.45*cdsqte/6.*7.01e-12
	 * call linadd( f2 , 2837 , 'BobR' , 'i')
	 * f2 = dense.xIonDense(26,4)*sexp(55 989./te)*0.384*cdsqte/6.*7.74e-12
	 * call linadd( f2 , 2568 , 'BobR' , 'i') */

	/* iron */

	PutLine(TauLines[ipFe1_24m],
		"  Fe 1 24m ");

	PutLine(TauLines[ipFe1_35m],
		"  Fe 1 35m ");

	PutLine(TauLines[ipFe1_54m],
		"  Fe 1 54m ");

	PutLine(TauLines[ipFe1_111m],
		"  Fe 1 111m ");

	PutLine(TauLines[ipFeI3884],
		"  Fe 1 3884 ");

	PutLine(TauLines[ipFeI3729],
		"  Fe 1 3729 ");

	PutLine(TauLines[ipFeI3457],
		"  Fe 1 3457 ");

	PutLine(TauLines[ipFeI3021],
		"  Fe 1 3021 ");

	PutLine(TauLines[ipFeI2966],
		"  Fe 1 2966 ");

	linadd(MAX2(0.,FeII.Fe2_large_cool+FeII.Fe2_UVsimp_cool),0,"Fe2c",'c' ,
		"total of all Fe 2 cooling, both simple UV and large atom together ");

	linadd(MAX2(0.,-FeII.Fe2_large_cool-FeII.Fe2_UVsimp_cool),0,"Fe2h",'h' ,
		"total of all Fe 2 heating, both simple UV and large atom together ");

	linadd(FeII.for7,4300,"Fe 2",'i' ,
		" Fe 2 forbidden 2-1 transition from Netzer's atom ");

	PutLine(TauLines[ipTuv3],
		" 2400 in simple Wills, Netzer, Wills FeII");
	PutLine(TauLines[ipTr48],
		" 6200 in simple Wills, Netzer, Wills FeII");
	PutLine(TauLines[ipTFe16],
		" 1080 in simple Wills, Netzer, Wills FeII");
	PutLine(TauLines[ipTFe26],
		" 1500  in simple Wills, Netzer, Wills FeII");
	PutLine(TauLines[ipTFe34],
		" 11500 in simple Wills, Netzer, Wills FeII");
	PutLine(TauLines[ipTFe35],
		" 2500 in simple Wills, Netzer, Wills FeII");
	PutLine(TauLines[ipTFe46],
		" 2300 in simple Wills, Netzer, Wills FeII");
	PutLine(TauLines[ipTFe56],
		" 8900 in simple Wills, Netzer, Wills FeII");

	/* option to save all intensities predicted by large FeII atom,
	 * code is in FeIILevelPops */
	FeIIAddLines();
	/* we were called by lines, and we want to zero out Fe2SavN */
	for( long ipLo=0; ipLo < (FeII.nFeIILevel_malloc - 1); ipLo++ )
	{
		for( long ipHi=ipLo + 1; ipHi < FeII.nFeIILevel_malloc; ipHi++ )
		{
			/* only evaluate real transitions */
			TransitionProxy tr = Fe2LevN[ipFe2LevN[ipHi][ipLo]];
			if( tr.ipCont() > 0 ) 
				PutLine( tr ," Fe II emission" );
		}
	}

	/* emission bands from the model Fe II atom */
	for( i=0; i < nFeIIBands; i++ )
	{
		double SumBandInward;
		/* [i][0] is center wavelength, [i][1] and [i][2] are upper and
		 * lower bounds in Angstroms.  These are set in FeIIZero 
		 * units are erg s-1 cm-3 */
		eff = FeIISumBand(FeII_Bands[i][1],FeII_Bands[i][2],
			&SumBandInward);

		linadd(eff,FeII_Bands[i][0],"Fe2b",'i' ,
			" total Fe II emission in Fe II bands, as defined in FeII_bands.ini ");
		linadd(SumBandInward,FeII_Bands[i][0],"Inwd",'i' ,
			" inward Fe II emission in Fe II bands, as defined in FeII_bands.ini ");
	}

	// integrate the pseudo continuum of FeII emission
	if( LineSave.ipass > 0 )
	{
		// initialize
		if( nzone == 1 )
		{
			for( i=0; i < nFeIIConBins; i++ )
			{
				/* initialize arrays */
				FeII_Cont[i][1] = 0.;
				FeII_Cont[i][2] = 0.;
			}
		}

		// integrate
		for( i=0; i < nFeIIConBins; i++ )
		{
			double SumBandInward;
			/* [i][0] is total intensity in cell, [i][1] and [i][2] are lower and
			 * upper bounds in Angstroms.  these are set in FeIIZero *
			 * find total emission from large FeII atom, integrated over band */
			double TotalFeII = FeIISumBand(FeII_Cont[i][0],FeII_Cont[i+1][0],
				&SumBandInward);
			FeII_Cont[i][1] += (realnum)(SumBandInward*radius.dVeffAper);
			FeII_Cont[i][2] += (realnum)(MAX2(0.,TotalFeII-SumBandInward)*radius.dVeffAper);
			/*fprintf(ioQQQ,"DEBUG feii\t%li\t%.2e\n", i, FeII_Cont[i][0]);*/
		}
	}

	PutLine(TauLines[ipT191],
		"  anomalous Fe 2 transition at 1787, RMT 191");

	linadd(fe.Fe3CoolTot,0,"Fe3c",'c' ,
		" chng 05 dec 16, FeIII code created by Kevin Blagrave  Fe3c 0 - total cooling due to 14-level Fe 3 atom ");

	/* Fe 3 14-level atom 
	 * following from print statements within loop */
	/* Fe 3 22.92m from Blagrave 14-level atom */
 	/* Fe 3 13.53m from Blagrave 14-level atom */
 	/* Fe 3 33.03m from Blagrave 14-level atom */
 	/* Fe 3 10.72m from Blagrave 14-level atom */
 	/* Fe 3 20.15m from Blagrave 14-level atom */
 	/* Fe 3 51.67m from Blagrave 14-level atom */
 	/* Fe 3 9.732m from Blagrave 14-level atom */
 	/* Fe 3 16.91m from Blagrave 14-level atom */
 	/* Fe 3 34.66m from Blagrave 14-level atom */
 	/* Fe 3 105.3m from Blagrave 14-level atom */
 	/* Fe 3  5152A from Blagrave 14-level atom */
 	/* Fe 3  5271A from Blagrave 14-level atom */
 	/* Fe 3  5356A from Blagrave 14-level atom */
 	/* Fe 3  5412A from Blagrave 14-level atom */
 	/* Fe 3  5440A from Blagrave 14-level atom */
 	/* Fe 3  4986A from Blagrave 14-level atom */
 	/* Fe 3  5097A from Blagrave 14-level atom */
 	/* Fe 3  5177A from Blagrave 14-level atom */
 	/* Fe 3  5230A from Blagrave 14-level atom */
 	/* Fe 3  5256A from Blagrave 14-level atom */
 	/* Fe 3 15.47m from Blagrave 14-level atom */
 	/* Fe 3  4925A from Blagrave 14-level atom */
 	/* Fe 3  5033A from Blagrave 14-level atom */
 	/* Fe 3  5111A from Blagrave 14-level atom */
 	/* Fe 3  5162A from Blagrave 14-level atom */
 	/* Fe 3  5188A from Blagrave 14-level atom */
 	/* Fe 3 11.16m from Blagrave 14-level atom */
 	/* Fe 3 40.04m from Blagrave 14-level atom */
 	/* Fe 3  4881A from Blagrave 14-level atom */
 	/* Fe 3  4988A from Blagrave 14-level atom */
 	/* Fe 3  5064A from Blagrave 14-level atom */
 	/* Fe 3  5114A from Blagrave 14-level atom */
 	/* Fe 3  5139A from Blagrave 14-level atom */
 	/* Fe 3 9.282m from Blagrave 14-level atom */
 	/* Fe 3 23.21m from Blagrave 14-level atom */
 	/* Fe 3 55.20m from Blagrave 14-level atom */
 	/* Fe 3  4833A from Blagrave 14-level atom */
 	/* Fe 3  4937A from Blagrave 14-level atom */
 	/* Fe 3  5012A from Blagrave 14-level atom */
 	/* Fe 3  5061A from Blagrave 14-level atom */
 	/* Fe 3  5085A from Blagrave 14-level atom */
 	/* Fe 3 7.789m from Blagrave 14-level atom */
 	/* Fe 3 15.69m from Blagrave 14-level atom */
 	/* Fe 3 25.79m from Blagrave 14-level atom */
 	/* Fe 3 48.41m from Blagrave 14-level atom */
 	/* Fe 3  4714A from Blagrave 14-level atom */
 	/* Fe 3  4813A from Blagrave 14-level atom */
 	/* Fe 3  4884A from Blagrave 14-level atom */
 	/* Fe 3  4931A from Blagrave 14-level atom */
 	/* Fe 3  4954A from Blagrave 14-level atom */
 	/* Fe 3 5.543m from Blagrave 14-level atom */
 	/* Fe 3 8.638m from Blagrave 14-level atom */
 	/* Fe 3 11.01m from Blagrave 14-level atom */
 	/* Fe 3 13.76m from Blagrave 14-level atom */
 	/* Fe 3 19.22m from Blagrave 14-level atom */
 	/* Fe 3  4659A from Blagrave 14-level atom */
 	/* Fe 3  4755A from Blagrave 14-level atom */
 	/* Fe 3  4825A from Blagrave 14-level atom */
 	/* Fe 3  4870A from Blagrave 14-level atom */
 	/* Fe 3  4893A from Blagrave 14-level atom */
 	/* Fe 3 4.859m from Blagrave 14-level atom */
 	/* Fe 3 7.085m from Blagrave 14-level atom */
 	/* Fe 3 8.608m from Blagrave 14-level atom */
 	/* Fe 3 10.20m from Blagrave 14-level atom */
 	/* Fe 3 12.92m from Blagrave 14-level atom */
 	/* Fe 3 39.41m from Blagrave 14-level atom */
 	/* Fe 3  4608A from Blagrave 14-level atom */
 	/* Fe 3  4702A from Blagrave 14-level atom */
 	/* Fe 3  4770A from Blagrave 14-level atom */
 	/* Fe 3  4814A from Blagrave 14-level atom */
 	/* Fe 3  4836A from Blagrave 14-level atom */
 	/* Fe 3 4.356m from Blagrave 14-level atom */
 	/* Fe 3 6.063m from Blagrave 14-level atom */
 	/* Fe 3 7.146m from Blagrave 14-level atom */
 	/* Fe 3 8.208m from Blagrave 14-level atom */
 	/* Fe 3 9.884m from Blagrave 14-level atom */
 	/* Fe 3 20.34m from Blagrave 14-level atom */
 	/* Fe 3 42.06m from Blagrave 14-level atom */
 	/* Fe 3  4574A from Blagrave 14-level atom */
 	/* Fe 3  4668A from Blagrave 14-level atom */
 	/* Fe 3  4734A from Blagrave 14-level atom */
 	/* Fe 3  4778A from Blagrave 14-level atom */
 	/* Fe 3  4800A from Blagrave 14-level atom */
 	/* Fe 3 4.077m from Blagrave 14-level atom */
 	/* Fe 3 5.535m from Blagrave 14-level atom */
 	/* Fe 3 6.423m from Blagrave 14-level atom */
 	/* Fe 3 7.269m from Blagrave 14-level atom */
 	/* Fe 3 8.554m from Blagrave 14-level atom */
 	/* Fe 3 15.41m from Blagrave 14-level atom */
 	/* Fe 3 25.31m from Blagrave 14-level atom */
 	/* Fe 3 63.56m from Blagrave 14-level atom */
	for( ihi=1; ihi<NLFE3; ++ihi )
	{
		for( ilo=0; ilo<ihi; ++ilo )
		{
			/* emission in these lines */
			PntForLine(fe.Fe3_wl[ihi][ilo],"Fe 3",&ipnt);
#			if 0
			fprintf( ioQQQ,"\t/* FeIII ");
			prt_wl( ioQQQ , (realnum)(fe.Fe3_wl[ihi][ilo]+0.5) );
			fprintf( ioQQQ," from Blagrave 14-level atom */\n" );
#			endif
			lindst( fe.Fe3_emiss[ihi][ilo] , (realnum)(fe.Fe3_wl[ihi][ilo]+0.5) , "Fe 3",ipnt,'c',true,
				" " );
		}
	}

	/*>>chng 05 dec 18, following are now in the above */
	/* sum of 3p and 3g states together */
	/*	linadd(CoolHeavy.c5270,0,"Fe 3",'c' ); */

	/* Fe 3 5270, predictions from Garstang et al 78
	PntForLine(5270.,"Fe 3",&ipnt);
	lindst(CoolHeavy.c5270*0.2090,5270,"Fe 3",ipnt,'c',true );*/

	/* Fe 3 5270, predictions from Garstang et al 78 
	PntForLine(4658.,"Fe 3",&ipnt);
	lindst(CoolHeavy.c5270*0.3667,4658,"Fe 3",ipnt,'c',true ); */

	PutLine(TauLines[ipT1122]," Fe 3 1122 entire multiplet");

	linadd(fe.Fe4CoolTot,0,"Fe4c",'i',
		" Fe4c 0 - total cooling due to 12-level Fe 4 atom " );


	PntForLine(3096.,"Fe 4",&ipnt);
	lindst(fe.fe40401,3096,"Fe 4",ipnt,'c',true,
		" Fe 4 3096.A, 4-1 and 5-1 transitions together"  );


	PntForLine(2836.,"Fe 4",&ipnt);
	lindst(fe.fe42836,2836,"Fe 4",ipnt,'c',true,
		" Fe 4 2835.7A, 6-1 transition, 4P5/2 - 6S5/2 "  );


	PntForLine(2829.,"Fe 4",&ipnt);
	lindst(fe.fe42829,2829,"Fe 4",ipnt,'c',true,
		"   Fe 4 2829.4A, 7-1 transition, 4P3/2 - 6S5/2"  );


	PntForLine(2567.,"Fe 4",&ipnt);
	lindst(fe.fe42567,2567,"Fe 4",ipnt,'c',true,
		"  Fe 4 2567.6+ 2567.4. 11-1 and 12-1 transitions"  );


	PntForLine(2.774e4,"Fe 4",&ipnt);
	lindst(fe.fe41207,2.774e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.774 microns 12-7 transition "  );


	PntForLine(2.714e4,"Fe 4",&ipnt);
	lindst(fe.fe41206,2.714e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.714 microns 12-6 transition "  );


	PntForLine(2.716e4,"Fe 4",&ipnt);
	lindst(fe.fe41106,2.716e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.716 microns 11-6 transition"  );


	PntForLine(2.806e4,"Fe 4",&ipnt);
	lindst(fe.fe41007,2.806e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.806 microns 10-7 transition "  );


	PntForLine(2.865e4,"Fe 4",&ipnt);
	lindst(fe.fe41008,2.865e4,"Fe 4",ipnt,'c',true ,
		"  Fe 4 2.865 microns 10-8 transition");


	PntForLine(2.836e4,"Fe 4",&ipnt);
	lindst(fe.fe40906,2.836e4,"Fe 4",ipnt,'c',true,
		" Fe 4 2.836 microns 9-6 transition" );


	PntForLine(3892.,"Fe 5",&ipnt);
	lindst(CoolHeavy.c3892,3892,"Fe 5",ipnt,'c',true,
		" Fe 5  3892+3839" );
	/* recombination Ka */
	if( dense.lgElmtOn[ipIRON] )
	{
		/* these lines added to outlin in metdif - following must be false
		 * fela = xLyaHeavy(nelem,nelem)*dense.xIonDense(nelem,nelem+1) */
		fela = iso_sp[ipH_LIKE][ipIRON].trans(ipH2p,ipH1s).Emis().xIntensity();
	}
	else
	{
		fela = 0.;
	}

	/* >>chng 02 jan 14, add grain fe to this sum */
	/* total intensity of K-alpha line */
	/*linadd((fe.fekcld+fe.fegrain)*1.03e-8+(fe.fekhot+fela)*1.11e-8,2,"FeKa",'i' );*/
	if( dense.lgElmtOn[ipIRON] )
	{
		lindst((fe.fekcld+fe.fegrain)*1.03e-8+(fe.fekhot+fela)*1.11e-8,1.78f,"FeKa",
			iso_sp[ipH_LIKE][ipIRON].trans(ipH2p,ipH1s).ipCont(),'i',false,
			   "total intensity of K-alpha line" );
	}

	linadd(fela*1.11e-8,2,"FeLr",'i' ,
		" recombination from fully stripped ion ");

	/* >>chng 03 aug 14, label changed from TotH to AugH to be like rest total hot iron Ka; */
	linadd((fe.fekhot+fela)*1.11e-8,2,"AugH",'i' ,
		"  Auger hot iron, assumes case b for H and He-like ");

	linadd(fe.fekcld*1.03e-8,2,"AugC",'i',
		" Auger production of cold iron, less than or 17 times ionized " );

	linadd(fe.fegrain*1.03e-8,2,"AugG",'i' ,
		" grain production of cold iron ");

	PutLine(TauLines[ipCo11527],
		"  [Co XI] 5168. A ");

	PutLine(TauLines[ipNi1_7m],
		"  nickel  [Ni I] 7m ");

	/* nickel*/


	PutLine(TauLines[ipNi1_11m],
		"  [Ni I] 11m ");

	/* copper */

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_lv1_k_zn returns\n" );
	}
	return;
}
