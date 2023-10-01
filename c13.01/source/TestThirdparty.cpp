/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include <UnitTest++.h>
#include "cddefines.h"
#include "thirdparty.h"

namespace {

	// constants below were calculated with xmaxima 5.17.0 with fpprec set to 30

	TEST(TestFactorial)
	{
		CHECK( fp_equal( factorial(0), 1. ) );
		CHECK( fp_equal( factorial(1), 1. ) );
		CHECK( fp_equal( factorial(2), 2. ) );
		CHECK( fp_equal( factorial(10), 3628800. ) );
		CHECK( fp_equal( factorial(170), 7.257415615307998967396728211129316e306 ) );
	}

	TEST(TestLogFactorial)
	{
		CHECK( fp_equal( lfactorial(0), 0. ) );
		CHECK( fp_equal( lfactorial(1), 0. ) );
		CHECK( fp_equal( lfactorial(2), 3.01029995663981195213738894725e-1, 10 ) );
		CHECK( fp_equal( lfactorial(10), 6.55976303287679375117476123996e0, 10 ) );
		CHECK( fp_equal( lfactorial(170), 3.06860781994828320192380273111e2, 10 ) );
		CHECK( fp_equal( lfactorial(1700), 4.75547688279135071178908638597e3, 10 ) );
	}

	TEST(TestCDGamma)
	{
		complex<double> y = cdgamma( complex<double>(1.,0.) );
		CHECK( fp_equal( y.real(), 1., 10 ) && fp_equal( y.imag(), 0. ) );
		y = cdgamma( complex<double>(2.,0.) );
		CHECK( fp_equal( y.real(), 1., 10 ) && fp_equal( y.imag(), 0. ) );
		y = cdgamma( complex<double>(11.,0.) );
		CHECK( fp_equal( y.real(), 3628800., 50 ) && fp_equal( y.imag(), 0. ) );	
		y = cdgamma( complex<double>(-0.5,0.) );
		CHECK( fp_equal( y.real(), -3.544907701811032054596334966682277e0, 10 ) &&
		       fp_equal( y.imag(), 0. ) );
		y = cdgamma( complex<double>(0.,1.) );
		CHECK( fp_equal( y.real(), -1.549498283018106851249551304838863e-1, 10 ) &&
		       fp_equal( y.imag(), -4.980156681183560427136911174621973e-1, 10 ) );
		y = cdgamma( complex<double>(-1.,-2.) );
		CHECK( fp_equal( y.real(), -3.23612885501927256868232032760969e-2, 30 ) &&
		       fp_equal( y.imag(), -1.122942423463261735043406872030743e-2, 30 ) );		       
	}

	// constants below were taken from Abramowitz & Stegun, Handbook of Mathematical Functions

	TEST(TestBesselJ0)
	{
		CHECK( fp_equal( bessel_j0(0.), 1., 10 ) );
		CHECK( fp_equal_tol( bessel_j0(1.), 0.765197686557967, 1.e-15 ) );
		CHECK( fp_equal_tol( bessel_j0(2.9), -0.224311545791968, 1.e-15 ) );
		CHECK( fp_equal_tol( bessel_j0(4.8), -0.240425327291183, 1.e-15 ) );
		CHECK( fp_equal_tol( bessel_j0(8.3), 0.096006100895010, 1.e-15 ) );
		CHECK( fp_equal_tol( bessel_j0(17.4), -0.118955856336348, 1.e-15 ) );
	}

	TEST(TestBesselJ1)
	{
		CHECK( fp_equal( bessel_j1(0.), 0., 10 ) );
		CHECK( fp_equal( bessel_j1(1.e-30), 5.e-31, 10 ) );
		CHECK( fp_equal_tol( bessel_j1(1.e-15), 5.e-16, 1.e-27 ) );
		CHECK( fp_equal_tol( bessel_j1(1.), 0.4400505857, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_j1(2.9), 0.3754274818, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_j1(4.8), -0.2984998581, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_j1(8.3), 0.2657393020, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_j1(17.4), -0.1532161760, 1.e-10 ) );
	}

	TEST(TestBesselJn)
	{
		CHECK( fp_equal( bessel_jn(2,0.), 0., 10 ) );
		CHECK( fp_equal( bessel_jn(2,1.e-30), 1.25e-61, 10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,2.e-15), 5.e-31, 1.e-42 ) );
		CHECK( fp_equal_tol( bessel_jn(2,2.e-10), 5.e-21, 1.e-32 ) );
		CHECK( fp_equal_tol( bessel_jn(2,0.099999999999), 0.0012489587, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,0.100000000001), 0.0012489587, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,1.), 0.1149034849, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,2.9), 0.4832270505, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,4.8), 0.1160503864, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,8.3), -0.0319725341, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,17.4), 0.1013448016, 1.e-10 ) );

		CHECK( fp_equal( bessel_jn(8,0.), 0., 10 ) );
		CHECK( fp_equal( bessel_jn(8,2.e-20), 1.e-160/40320., 10 ) );
		CHECK( fp_equal_tol( bessel_jn(8,2.e-15), 1.e-120/40320., 1.e-136 ) );
		CHECK( fp_equal_tol( bessel_jn(8,1.), 9.4223e-8, 1.e-12 ) );
		CHECK( fp_equal_tol( bessel_jn(8,2.8), 2.9367e-4, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_jn(8,4.8), 1.4079e-2, 1.e-6 ) );
		CHECK( fp_equal_tol( bessel_jn(8,8.2), 0.24257, 1.e-5 ) );

		CHECK( fp_equal_tol( bessel_jn(20,1.), 3.873503009e-25, 1.e-34 ) );
		CHECK( fp_equal_tol( bessel_jn(50,1.), 2.906004948e-80, 1.e-89 ) );
		CHECK( fp_equal_tol( bessel_jn(100,1.), 8.431828790e-189, 1.e-198 ) );
		CHECK( fp_equal_tol( bessel_jn(100,100.), 9.636667330e-2, 1.e-11 ) );
	}

	TEST(TestBesselY0)
	{
		CHECK( fp_equal_tol( bessel_y0(1.e-30), -44.0499402279, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(1.), 0.0882569642, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(2.9), 0.4079117692, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(4.8), -0.2723037945, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(8.3), 0.2595149638, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(17.4), -0.1497391883, 1.e-10 ) );
	}

	TEST(TestBesselY1)
	{
		CHECK( fp_equal_tol( bessel_y1(1.e-30), -6.36619772368e29, 1.e18 ) );
		CHECK( fp_equal_tol( bessel_y1(1.), -0.7812128213, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y1(2.9), 0.2959400546, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y1(4.8), 0.2135651673, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y1(8.3), -0.0805975035, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y1(17.4), 0.1147053859, 1.e-10 ) );
	}

	TEST(TestBesselYn)
	{
		CHECK( fp_equal_tol( bessel_yn(2,1.e-30), -1.27323954474e60, 1.e49 ) );
		CHECK( fp_equal_tol( bessel_yn(2,1.), -1.65068261, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_yn(2,2.9), -0.20381518, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_yn(2,4.8), 0.36128928, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_yn(2,8.3), -0.27893605, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_yn(2,17.4), 0.16292372, 1.e-8 ) );

		CHECK( fp_equal_tol( bessel_yn(8,2.e-20), -1.60428182637e163, 1.e152 ) );
		CHECK( fp_equal_tol( bessel_yn(8,1.), -4.2567e5, 1.e1 ) );
		CHECK( fp_equal_tol( bessel_yn(8,2.8), -1.4486e2, 1.e-2 ) );
		CHECK( fp_equal_tol( bessel_yn(8,4.8), -3.5855, 1.e-4 ) );
		CHECK( fp_equal_tol( bessel_yn(8,8.2), -0.35049, 1.e-5 ) );

		CHECK( fp_equal_tol( bessel_yn(20,1.), -4.113970315e22, 1.e13 ) );
		CHECK( fp_equal_tol( bessel_yn(50,1.), -2.191142813e77, 1.e68 ) );
		CHECK( fp_equal_tol( bessel_yn(100,1.), -3.775287810e185, 1.e176 ) );
		CHECK( fp_equal_tol( bessel_yn(100,100.), -1.669214114e-1, 1.e-10 ) );
	}

	TEST(TestBesselI0)
	{
		CHECK( fp_equal( bessel_i0_scaled(0.), 1., 10 ) );
		CHECK( fp_equal_tol( bessel_i0_scaled(1.), 0.4657596077, 2.e-10 ) );
		CHECK( fp_equal_tol( bessel_i0_scaled(2.9), 0.2477557304, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_i0_scaled(4.8), 0.1875862042, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_i0_scaled(8.3), 0.1407239098, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_i0_scaled(17.4), 0.0963498277, 1.e-10 ) );

		CHECK( fp_equal( bessel_i0_scaled(1.), exp(-1.)*bessel_i0(1.) ) );
		CHECK( fp_equal( bessel_i0_scaled(17.4), exp(-17.4)*bessel_i0(17.4) ) );
	}

	TEST(TestBesselI1)
	{
		CHECK( fp_equal( bessel_i1_scaled(0.), 0., 10 ) );
		CHECK( fp_equal_tol( bessel_i1_scaled(1.), 0.2079104154, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_i1_scaled(2.9), 0.1987772816, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_i1_scaled(4.8), 0.1666757058, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_i1_scaled(8.3), 0.1319524362, 2.e-10 ) );
		CHECK( fp_equal_tol( bessel_i1_scaled(17.4), 0.0935388542, 1.e-10 ) );

		CHECK( fp_equal( bessel_i1_scaled(1.), exp(-1.)*bessel_i1(1.) ) );
		CHECK( fp_equal( bessel_i1_scaled(17.4), exp(-17.4)*bessel_i1(17.4) ) );
	}

	TEST(TestBesselK0)
	{
		CHECK( fp_equal_tol( bessel_k0(2.e-30), 68.5003371249, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_k0_scaled(1.), 1.1444630797, 2.e-10 ) );
		CHECK( fp_equal_tol( bessel_k0_scaled(2.9), 0.7089049774, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_k0_scaled(4.8), 0.5586133194, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_k0_scaled(8.3), 0.4288766329, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_k0_scaled(17.4), 0.2983665276, 1.e-10 ) );

		CHECK( fp_equal( bessel_k0_scaled(1.), exp(1.)*bessel_k0(1.) ) );
		CHECK( fp_equal( bessel_k0_scaled(17.4), exp(17.4)*bessel_k0(17.4) ) );
	}

	TEST(TestBesselK1)
	{
		CHECK( fp_equal_tol( bessel_k1(2.e-30), 5e29, 1.e18 ) );
		CHECK( fp_equal_tol( bessel_k1_scaled(1.), 1.6361534863, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_k1_scaled(2.9), 0.8230420403, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_k1_scaled(4.8), 0.6142566003, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_k1_scaled(8.3), 0.4540139001, 2.e-10 ) );
		CHECK( fp_equal_tol( bessel_k1_scaled(17.4), 0.3068236027, 1.e-10 ) );

		CHECK( fp_equal( bessel_k1_scaled(1.), exp(1.)*bessel_k1(1.) ) );
		CHECK( fp_equal( bessel_k1_scaled(17.4), exp(17.4)*bessel_k1(17.4) ) );
	}

	TEST(TestEllpk)
	{
		CHECK( fp_equal( ellpk(1.), PI/2., 10 ) );
		CHECK( fp_equal_tol( ellpk(0.86), 1.630575548881754, 1.e-15 ) );
		CHECK( fp_equal_tol( ellpk(0.56), 1.806327559107699, 1.e-15 ) );
		CHECK( fp_equal_tol( ellpk(0.36), 1.995302777664729, 1.e-15 ) );
		CHECK( fp_equal_tol( ellpk(0.13), 2.455338028321380, 1.e-15 ) );
		CHECK( fp_equal_tol( ellpk(0.01), 3.695637362989875, 1.e-15 ) );
	}

	TEST(TestExpn)
	{
		CHECK( fp_equal_tol( expn(1,0.25), 0.25*0.9408157528 - log(0.25) - EULER, 1.e-10 ) );
		CHECK( fp_equal_tol( expn(1,0.79), 0.316277004, 1.e-9 ) );
		CHECK( fp_equal_tol( expn(1,1.97), 0.050976988, 1.e-9 ) );

		CHECK( fp_equal_tol( expn(2,0.25), 0.8643037 + log(0.25)/4., 1.e-7 ) );
		CHECK( fp_equal_tol( expn(2,0.79), 0.2039860, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(2,1.97), 0.0390322, 1.e-7 ) );

		CHECK( fp_equal( expn(3,0.), 0.5, 10 ) );
		CHECK( fp_equal_tol( expn(3,0.15), 0.3822761, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(3,0.79), 0.1463479, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(3,1.97), 0.0312817, 1.e-7 ) );

		CHECK( fp_equal( expn(4,0.), 1./3., 10 ) );
		CHECK( fp_equal_tol( expn(4,0.15), 0.2677889, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(4,0.79), 0.1127433, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(4,1.97), 0.0259440, 1.e-7 ) );

		CHECK( fp_equal( expn(10,0.), 1./9., 10 ) );
		CHECK( fp_equal_tol( expn(10,0.15), 0.0938786, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(10,0.79), 0.0459453, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(10,1.97), 0.0124964, 1.e-7 ) );
	}

	TEST(TestErf)
	{
		/* constants calculated with xmaxima 5.22.1 */
		CHECK( fp_equal_tol( erf(0.), 0., 1.e-22 ) );
		// erf(x) loses some precision around 1.e-10, but should still be plenty good...
		CHECK( fp_equal_tol( erf(1.e-10), 1.1283791671081724525e-10, 3.e-21 ) );
		CHECK( fp_equal_tol( erf(1.e-5), 1.1283791670579000307e-5, 1.e-17 ) );
		CHECK( fp_equal_tol( erf(0.1), 1.1246291601828489259e-1, 1.e-13 ) );
		CHECK( fp_equal_tol( erf(0.5), 5.2049987781304653768e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erf(1.), 8.4270079294971486934e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erf(2.), 9.9532226501895273416e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erf(10.), 1.0, 1.e-12 ) );
		CHECK( fp_equal( erf(-1.), -erf(1.) ) );
	}

	TEST(TestErfc)
	{
		/* constants calculated with xmaxima 5.22.1 */
		CHECK( fp_equal_tol( erfc(-1.), 1.8427007929497148693e0, 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(0.), 1., 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(1.e-5), 9.99988716208329421e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(0.1), 8.8753708398171510741e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(0.5), 4.7950012218695346232e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(1.), 1.5729920705028513066e-1, 1.e-13 ) );
		CHECK( fp_equal_tol( erfc(2.), 4.6777349810472658396e-3, 1.e-15 ) );
		CHECK( fp_equal_tol( erfc(10.), 2.088487583762544757e-45, 1.e-57 ) );
	}

	TEST(TestErfce)
	{
		/* constants taken from Finn G.D. & Mugglestone D., 1965, MNRAS 129, 221 */
		/* this is the voigt function H(a,0) normalized according to 9-44 of Mihalas */
		CHECK( fp_equal_tol( erfce(0.), 1., 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.01), 9.88815e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.02), 9.77826e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.05), 9.45990e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.10), 8.96456e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.20), 8.09019e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.55), 5.90927e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(1.00), 4.27583e-1, 1.e-6 ) );
	}

	TEST(TestVoigtH0)
	{
		/* constants taken from Finn G.D. & Mugglestone D., 1965, MNRAS 129, 221 */
		/* this is the voigt function H(a,0) normalized according to 9-44 of Mihalas */
		CHECK( fp_equal_tol( VoigtH0(0.), 1., 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.01), 9.88815e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.02), 9.77826e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.05), 9.45990e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.10), 8.96456e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.20), 8.09019e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.55), 5.90927e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(1.00), 4.27583e-1, 1.e-6 ) );
	}

	TEST(TestVoigtU0)
	{
		/* constants taken from Finn G.D. & Mugglestone D., 1965, MNRAS 129, 221 */
		/* this is the voigt function U(a,0) normalized according to 9-45 of Mihalas */
		CHECK( fp_equal_tol( VoigtU0(0.), 1./SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.01), 9.88815e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.02), 9.77826e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.05), 9.45990e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.10), 8.96456e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.20), 8.09019e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.55), 5.90927e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(1.00), 4.27583e-1/SQRTPI, 1.e-6 ) );
	}

	TEST(TestVoigtH)
	{
		// check that the Voigt profile returned by VoigtH() is properly normalized
		const int NP = 200;
		realnum v[NP], a, y[NP];
		// for a > 0.1, VoigtH() calls humlik(), for smaller values it
		// calls FastVoigtH().
		//
		// humlik() is set up for a relative precision of 1e-4 over its
		// entire range, but looks to be more precise in practice (at
		// least at v=0)
		// FastVoigtH() is set up for a rel. precision of 2.5e-3 over its
		// entire range, but should be better than 1e-4 for a < 0.0235
		a = realnum(0.0002);
		for( int i=0; i < 9; ++i )
		{
			// test both humlik() and FastVoigtH()
			for( int i=0; i < NP; ++i )
				v[i] = realnum(i)*max(a,1.f)/realnum(5.);
			VoigtH( a, v, y, NP );
			realnum integral = realnum(0.);
			// We need the integral from -infinity to +infinity, which is simply
			// two times the integral from 0 to +infinity. Hence we omit the
			// division by 2 in the trapezoidal rule
			for( int i=1; i < NP; ++i )
				integral += (y[i]+y[i-1])*(v[i]-v[i-1]);
			// add in the integral over the Lorentz wings assuming H(a,v) = c/v^2
			integral += realnum(2.)*v[NP-1]*y[NP-1];
			// VoigtH() calculates H(a,v), so integral should be sqrt(pi)
			CHECK( fp_equal_tol( integral, realnum(SQRTPI), realnum(1.e-4) ) );
			// also check the central value...
			realnum h0 = realnum(VoigtH0(a));
			CHECK( fp_equal_tol( y[0], h0, realnum(1.e-4)*h0 ) );
			a *= realnum(10.);
		}
		// now do some spot checks
		// constants taken from
		// >>refer	Zaghloul M.R. & Ali A.N. 2011, ACM Transactions on Mathematical Software, 38, 15
		// available from http://arxiv.org/abs/1106.0151
		v[0] = realnum(5.76);
		VoigtH( realnum(1.e-20), v, y, 1 );
		// this constant comes from page 21 ("Present algorithm").
		CHECK( fp_equal_tol( y[0], realnum(3.900779639194698e-015), realnum(4.e-19) ) );
		v[0] = realnum(6.3e-2);
		v[1] = realnum(6.3e+0);
		v[2] = realnum(6.3e+2);
		VoigtH( realnum(1.e-20), v, y, 3 );
		// these constants come from Table 2 of the same paper ("Present algorithm").
		CHECK( fp_equal_tol( y[0], realnum(9.960388660702479e-001), realnum(1.e-4) ) );
		CHECK( fp_equal_tol( y[1], realnum(5.792460778844116e-018), realnum(6.e-22) ) );
		CHECK( fp_equal_tol( y[2], realnum(1.421495882582395e-026), realnum(1.4e-30) ) );
		VoigtH( realnum(1.e-14), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(9.960388660702366e-001), realnum(1.e-4) ) );
		CHECK( fp_equal_tol( y[1], realnum(1.536857621303163e-016), realnum(1.5e-20) ) );
		CHECK( fp_equal_tol( y[2], realnum(1.421495882582395e-020), realnum(1.4e-24) ) );
		VoigtH( realnum(1.e-12), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(9.960388660691284e-001), realnum(1.e-4) ) );
		CHECK( fp_equal_tol( y[1], realnum(1.479513723737753e-014), realnum(1.5e-18) ) );
		CHECK( fp_equal_tol( y[2], realnum(1.421495882582395e-018), realnum(1.4e-22) ) );
		VoigtH( realnum(1.e-10), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(9.960388659583033e-001), realnum(1.e-4) ) );
		CHECK( fp_equal_tol( y[1], realnum(1.478940284762099e-012), realnum(1.5e-16) ) );
		CHECK( fp_equal_tol( y[2], realnum(1.421495882582395e-016), realnum(1.4e-20) ) );
		VoigtH( realnum(1.e-6), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(9.960377466254801e-001), realnum(1.e-4) ) );
		CHECK( fp_equal_tol( y[1], realnum(1.478934493028404e-008), realnum(1.5e-12) ) );
		CHECK( fp_equal_tol( y[2], realnum(1.421495882582395e-012), realnum(1.4e-16) ) );
		VoigtH( realnum(1.e-2), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(9.849424862549039e-001), realnum(1.e-4) ) );
		CHECK( fp_equal_tol( y[1], realnum(1.478930389133934e-004), realnum(1.5e-8) ) );
		CHECK( fp_equal_tol( y[2], realnum(1.421495882224242e-008), realnum(1.4e-12) ) );
		VoigtH( realnum(1.e+1), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(5.613881832823886e-002), realnum(6.e-6) ) );
		CHECK( fp_equal_tol( y[1], realnum(4.040671157393835e-002), realnum(4.e-6) ) );
		CHECK( fp_equal_tol( y[2], realnum(1.421137820009847e-005), realnum(1.4e-9) ) );
		VoigtH( realnum(1.2e+1), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(4.685295149211637e-002), realnum(5.e-6) ) );
		CHECK( fp_equal_tol( y[1], realnum(3.684277239564798e-002), realnum(4.e-6) ) );
		CHECK( fp_equal_tol( y[2], realnum(1.705176395541707e-005), realnum(1.7e-9) ) );
		VoigtH( realnum(1.5e+1), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(3.752895161491574e-002), realnum(4.e-6) ) );
		CHECK( fp_equal_tol( y[1], realnum(3.194834330452605e-002), realnum(3.e-6) ) );
		CHECK( fp_equal_tol( y[2], realnum(2.131035743074598e-005), realnum(2.e-9) ) );
		VoigtH( realnum(2.e+2), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(2.820912377324508e-003), realnum(3.e-7) ) );
		CHECK( fp_equal_tol( y[1], realnum(2.818116555672206e-003), realnum(3.e-7) ) );
		CHECK( fp_equal_tol( y[2], realnum(2.582702147491469e-004), realnum(3.e-8) ) );
		VoigtH( realnum(1.e+5), v, y, 3 );
		CHECK( fp_equal_tol( y[0], realnum(5.641895835193228e-006), realnum(6.e-10) ) );
		CHECK( fp_equal_tol( y[1], realnum(5.641895812802746e-006), realnum(6.e-10) ) );
		CHECK( fp_equal_tol( y[2], realnum(5.641671917237128e-006), realnum(6.e-10) ) );
		v[0] = realnum(1.e+0);
		VoigtH( realnum(1.e-20), v, y, 1 );
		CHECK( fp_equal_tol( y[0], realnum(3.678794411714423e-001), realnum(4.e-5) ) );
		v[0] = realnum(5.5e+0);
		VoigtH( realnum(1.e-14), v, y, 1 );
		CHECK( fp_equal_tol( y[0], realnum(7.307386729528773e-014), realnum(7.e-18) ) );
		v[0] = realnum(3.9e+4);
		VoigtH( realnum(1.e+0), v, y, 1 );
		CHECK( fp_equal_tol( y[0], realnum(3.709333226385423e-010), realnum(4.e-14) ) );
		v[0] = realnum(1.e+0);
		VoigtH( realnum(2.8e+4), v, y, 1 );
		CHECK( fp_equal_tol( y[0], realnum(2.014962794529686e-005), realnum(2.e-9) ) );
	}

	TEST(TestVoigtU)
	{
		// check that the Voigt profile returned by VoigtU() is properly normalized
		const int NP = 200;
		realnum v[NP], a, y[NP];
		// for a > 0.1, VoigtU() calls humlik(), for smaller values it
		// calls FastVoigtH() and divides by sqrt(pi).
		//
		// humlik() is set up for a relative precision of 1e-4 over its
		// entire range, but looks to be more precise in practice (at
		// least at v=0)
		// FastVoigtH() is set up for a rel. precision of 2.5e-3 over its
		// entire range, but should be better than 1e-4 for a < 0.0235
		a = realnum(0.0002);
		for( int i=0; i < 9; ++i )
		{
			// test both humlik() and FastVoigtH()
			for( int i=0; i < NP; ++i )
				v[i] = realnum(i)*max(a,1.f)/realnum(5.);
			VoigtU( a, v, y, NP );
			realnum integral = realnum(0.);
			// We need the integral from -infinity to +infinity, which is simply
			// two times the integral from 0 to +infinity. Hence we omit the
			// division by 2 in the trapezoidal rule
			for( int i=1; i < NP; ++i )
				integral += (y[i]+y[i-1])*(v[i]-v[i-1]);
			// add in the integral over the Lorentz wings assuming U(a,v) = c/v^2
			integral += realnum(2.)*v[NP-1]*y[NP-1];
			// VoigtU() calculates U(a,v), so integral should be 1
			CHECK( fp_equal_tol( integral, realnum(1.), realnum(1.e-4) ) );
			// also check the central value...
			CHECK( fp_equal_tol( y[0], realnum(VoigtU0(a)), realnum(1.e-4) ) );
			a *= realnum(10.);
		}
	}

	TEST(TestMD5string)
	{
		string test;
		// md5sum of an empty file...
		CHECK( MD5string( test ) == "d41d8cd98f00b204e9800998ecf8427e" );
		CHECK( MD5string( test ).length() == NMD5 );
		// check if padding is done correctly
		// an extra block of padding needs to be added when length%64 == 56
		test = "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 55 );
		CHECK( MD5string( test ) == "426ec4ac35ad38d125f6efb39da03098" );
		test = "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 56 );
		CHECK( MD5string( test ) == "d03607b2c89adc0c4abf5a0f1d9e40c9" );
		test = "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 57 );
		CHECK( MD5string( test ) == "bac1b47748411cb6eee0cae3befb8377" );
		string test64 = "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		test = test64 + "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 64+55 );
		CHECK( MD5string( test ) == "10d49aad1fc69976376fbe7c8c5ed118" );
		test = test64 + "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 64+56 );
		CHECK( MD5string( test ) == "61ec7da14576f3b585038c6d72cd5bd5" );
		test = test64 + "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 64+57 );
		CHECK( MD5string( test ) == "f17a0475a26d0930e2a35bb320c10e0d" );
		// check that leading zeros are printed correctly
		test = "ghijklmn";
		CHECK( MD5string( test ) == "0256b9cea63bc1f97b8c5aea92c24a98" );
	}

}
