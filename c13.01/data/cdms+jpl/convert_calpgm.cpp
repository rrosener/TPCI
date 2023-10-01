#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <cmath>

using namespace std;

struct state;

void process_file(const string& filename);

void InsertState( state** s, vector<state*>& v );

template<class T>
inline T pow2(T a) { return a*a; }
template<class T>
inline T pow3(T a) { return a*a*a; }

template<class T>
inline void read_data(const string& str,
		      T& x);

void read_Jval(const string& str, int& x);

inline string TrimSpaces(const string& str);
inline bool FindAndErase(string& str,
			 const string& substr);
inline bool FindAndReplace(string& str,
			   const string& substr,
			   const string& newstr);

template<class T>
class EnergyCmp
{
public:
	inline bool operator() ( const T* x, const T* y ) const
	{
		return ( x->energy < y->energy );
	}
};

struct state
{
	int index;
	int nQuant;
	int quant[6];
	bool lgHalfInt[6];
	int glev;
	double energy;
	double unc;
	// flag that experimental frequenecy was used to calculate level energy
	// this implies that the energy may not match the energy calculated by the code
	bool lgExperimental;
	bool operator< (const state& s) const
	{
		for( int i=0; i < nQuant; ++i )
		{
			if( quant[i] < s.quant[i] )
				return true;
			if( quant[i] > s.quant[i] )
				return false;
		}
		return false;
	}
	bool operator== (const state& s) const
	{
		for( int i=0; i < nQuant; ++i )
		{
			if( quant[i] != s.quant[i] )
				return false;
		}
		return true;
	}
	state() : index(-1), nQuant(0), glev(0), energy(0.), unc(0.), lgExperimental(false)
	{
		for( int i=0; i < 6; ++i )
		{
			quant[i] = 0;
			lgHalfInt[i] = false;
		}
	}
	state( int n, const string q[], bool lgHI[], int g, double E, double u, bool lgExp )
	{
		index = -1;
		nQuant = n;
		for( int i=0; i < nQuant; ++i )
		{
			read_Jval( q[i], quant[i] );
			lgHalfInt[i] = lgHI[i];
		}
		for( int i=nQuant; i < 6; ++i )
		{
			quant[i] = 0;
			lgHalfInt[i] = false;
		}
		glev = g;
		energy = E;
		unc = u;
		lgExperimental = lgExp;
	}
};

struct line
{
	const state *lo;
	const state *hi;
	double freq;
	double Aki;
	line() : lo(0), hi(0), freq(0.), Aki(0.) {}
	line(const state* slo, const state* shi, double f, double A) : lo(slo), hi(shi), freq(f), Aki(A) {}
};

int main(int argc, char* argv[])
{
	if( argc != 2 )
	{
		cerr << "usage: " << argv[0] << " <filename>" << endl;
		return 1;
	}

	try
	{
		process_file(argv[1]);
		return 0;
	}
	catch( exception& e )
	{
		cerr << "An exception was caught, what: " << e.what() << endl;
		return 1;
	}
}

void process_file(const string& filename)
{
	string textline;
	map<string,double> Q300;
	map<string,string> names;

	bool lgJPL = ( filename[4] == '0' );

	ifstream pft;
	if( lgJPL )
		pft.open( "partition_function_jpl.html" );
	else
		pft.open( "partition_function_cdms.html" );

	if( pft == 0 )
		throw runtime_error( "could not open partition function file" );
	    
	if( !lgJPL )
	{
		// skip html header
		while( getline(pft,textline) && textline.find( "<pre>" ) == string::npos ) {}
		getline(pft,textline);
		getline(pft,textline);
	}
	while( getline(pft,textline) && textline.find( "</pre>" ) == string::npos )
	{
		string tag, name, Qrs_str;
		tag = TrimSpaces(textline.substr(0,6));
		if( lgJPL )
		{
			name = TrimSpaces(textline.substr(7,13));
			Qrs_str = textline.substr(26,7);
		}
		else
		{
			name = TrimSpaces(textline.substr(7,25));
			Qrs_str = textline.substr(64,13);
		}
		double Qrs;
		if( Qrs_str.find("---") == string::npos )
		{
			read_data( Qrs_str, Qrs );
			Qrs = pow( 10., Qrs );
		}
		else
			Qrs = -1.;
		Q300[tag] = Qrs;
		names[tag] = name;
	}

	vector<state*> states;
	vector<state*>::iterator s;

	vector<line> lines;
	vector<line>::iterator l;

	int molweight;

	int nQCode = -1;
	string tag = "";

	ifstream ldt( filename.c_str() );
	if( ldt == 0 )
		throw runtime_error( "could not open data file" );

	while( getline(ldt,textline) )
	{
		// pad with spaces if necessary
		if( textline.length() < 80 )
			for( string::size_type i=textline.length(); i < 80; ++i )
				textline += ' ';

		string field[20];
		// frequency of the line in MHz or cm^-1
		field[0] = textline.substr(0,13);
		// uncertainty in the frequency
		// a negative value indicates units cm^-1
		// a positive value indicates units MHz
		field[1] = textline.substr(13,8);
		// base-10 log of the line intensity at 300 K
		// can also be S*mu^2 for selected molecules !!!
		// in those cases the partition function will not
		// be set in partition_function_cdms.html
		field[2] = textline.substr(21,8);
		// degree of freedom, not used
		field[3] = textline.substr(29,2);
		// energy of the lower level in cm^-1
		field[4] = textline.substr(31,10);
		// upper state degeneracy g_up
		field[5] = textline.substr(41,3);
		// molecule tag, a negative value indicates that
		// the line frequency and uncertainty are experimental
		field[6] = textline.substr(44,7);
		// coding of the quantum numbers
		field[7] = textline.substr(51,4);
		// up to 6 quantum numbers for the upper state
		field[8] = textline.substr(55,2);
		field[9] = textline.substr(57,2);
		field[10] = textline.substr(59,2);
		field[11] = textline.substr(61,2);
		field[12] = textline.substr(63,2);
		field[13] = textline.substr(65,2);
		// up to 6 quantum numbers for the lower state
		field[14] = textline.substr(67,2);
		field[15] = textline.substr(69,2);
		field[16] = textline.substr(71,2);
		field[17] = textline.substr(73,2);
		field[18] = textline.substr(75,2);
		field[19] = textline.substr(77,2);

		bool lgFreqInMHz = ( field[1].find( '-' ) == string::npos );
		bool lgExperimentalFreq = ( field[6].find( '-' ) != string::npos );
		string ltag = TrimSpaces(field[6]);
		FindAndErase( ltag, "-" );
		if( tag.length() == 0 )
			tag = ltag;
		if( ltag != tag )
			throw runtime_error( "tag does not match" );

		// read the number of quantum numbers for each state
		int lnQCode;
		read_data( field[7], lnQCode );
		if( nQCode == -1 )
			nQCode = lnQCode;
		if( nQCode != lnQCode )
			throw runtime_error( "quantum number code does not match" );
		int nQuant = nQCode%10;

		// read coding for half-integer quantum numbers
		bool lgHalfInt[6];
		for( int i=0; i < 6; ++i )
			lgHalfInt[i] = false;

		int CodeHI = (nQCode/10)%10;
		if( CodeHI > 7 )
			throw runtime_error( "the half-integer code must be <= 7" );
		for( int i=0; i < 3; ++i )
		{
			if( (CodeHI&1) == 1 )
				lgHalfInt[nQuant-i-1] = true;
			CodeHI /= 2;
		}

		read_data( field[6].substr(1,3), molweight );

		// The data file only contains g_up, so we need to calculate g_lo.
		// The statistical weight is defined as g = g_l * g_F = g_l*(2*F+1).
		// When spin splitting is present, F is the last quantum number...
		// When there is no spin splitting, N = J = F, and F is omitted. But
		// then N is always the first quantum number, so use that instead...
		// There is no general formula for g_l. However, we can use the
		// fact that g_l never changes in an allowed transition, so we can
		// compute g_l by comparing g_up in the file with 2*F_up+1.

		int gup2;
		read_data( field[5], gup2 );

		int glo, gup;
		// most likely more codes need to be listed here...!
		if( nQCode == 202 || nQCode == 303 || nQCode == 1202 || nQCode == 1303 )
		{
			// no spin splitting, so F == N and we should use the first quantum number
			read_Jval( field[14], glo );
			read_Jval( field[8], gup );
		}
		else
		{
			// spin splitting included, so F is the last quantum number
			read_Jval( field[14+nQuant-1], glo );
			read_Jval( field[8+nQuant-1], gup );
		}
		// convert to statistical weight
		if( lgHalfInt[nQuant-1] )
		{
			// quantum numbers were already rounded up to the next integer
			glo = 2*glo;
			gup = 2*gup;
		}
		else
		{
			glo = 2*glo+1;
			gup = 2*gup+1;
		}

		// compute g_l by comapring the two values of g_up
		int gl = gup2/gup;

		glo *= gl;
		gup *= gl;

		if( gup != gup2 )
			throw runtime_error( "g_up insanity" );

		double Elo, Ehi, freq, freqInvCM, freqMHz;
		// read line frequency
		read_data( field[0], freq );
		// and convert
		if( lgFreqInMHz )
		{
			freqMHz = freq;
			freqInvCM = freq/2.99792458e4;
		}
		else
		{
			freqMHz = freq*2.99792458e4;
			freqInvCM = freq;
		}
		// read lower level energy, already in cm^-1
		read_data( field[4], Elo );
		Ehi = Elo + freqInvCM;

		double unc;
		read_data( field[1], unc );
		if( lgFreqInMHz )
			unc /= 2.99792458e4;
		unc = sqrt( pow2(unc) + 1.e-8 );

		state* slo = new state( nQuant, &field[14], lgHalfInt, glo, Elo, 1.e-4, false );
		InsertState( &slo, states );

		state* shi = new state( nQuant, &field[8], lgHalfInt, gup, Ehi, unc, lgExperimentalFreq );
		InsertState( &shi, states );

		// kT at 300K in cm^-1
		double kT300 = 300.*0.6950356;
		double intensity;
		read_data( field[2], intensity );
		intensity = pow( 10., intensity );
		if( Q300.find(tag) == Q300.end() )
			throw runtime_error( "tag not found" );
		double Qrs = Q300[tag];
		double EinsteinA;
		if( Qrs > 0. )
			EinsteinA = -intensity*pow2(freqMHz)*Qrs/gup/
				(exp(-Elo/kT300)*expm1(-freqInvCM/kT300))*2.7964e-16;
		else
			EinsteinA = 1.16395e-20*pow3(freqMHz)*intensity/double(gup);
		lines.push_back( line( slo, shi, freqMHz*1.e-3, EinsteinA ) );
	}

	sort( states.begin(), states.end(), EnergyCmp<state>() );

	int n = 0;
	// give each state an index number...
	for( s=states.begin(); s != states.end(); ++s )
		(*s)->index = ++n;

	// create LAMDA compatible output
	cout << "!MOLECULE" << endl;
	if( lgJPL )
	{
		cout << names[tag] << " (Levels and lines obtained from the ";
		cout << "JPL Molecular Spectroscopy Database)" << endl;
	}
	else
	{
		cout << names[tag] << " (Levels and lines obtained from the ";
		cout << "Cologne Database for Molecular Spectroscopy)" << endl;
	}
	cout << "!MOLECULAR WEIGHT" << endl;
	// this is rounded to the nearest integer, but doesn't appear to be used by Cloudy...
	cout << molweight << ".0" << endl;
	cout << "!NUMBER OF ENERGY LEVELS" << endl;
	cout << states.size() << endl;
	cout << "!LEVEL + ENERGIES(cm^-1) + WEIGHT + (quantum numbers)" << endl;
	for( s=states.begin(); s != states.end(); ++s )
	{
		cout << setw(5) << (*s)->index;
		cout << setw(16) << setprecision(9) << fixed << abs((*s)->energy);
		cout << setw(4) << (*s)->glev << ".0    ";
		for( int i=0; i < (*s)->nQuant; ++i )
		{
			if( (*s)->lgHalfInt[i] )
				cout << (*s)->quant[i]-1 << ".5";
			else
				cout << (*s)->quant[i];
			if( i != (*s)->nQuant-1 )
				cout << "_";
		}
		cout << endl;
	}
	cout << "!NUMBER OF RADIATIVE TRANSITIONS" << endl;
	cout << lines.size() << endl;
	cout << "!TRANS + UP + LOW + EINSTEIN_A[s-1] + FREQ[GHz] + E_u[K]" << endl;
	n = 0;
	for( l=lines.begin(); l != lines.end(); ++l )
	{
		cout << setw(5) << ++n;
		cout << setw(5) << l->hi->index;
		cout << setw(5) << l->lo->index;
		cout << setw(12) << setprecision(3) << scientific << l->Aki;
		cout << setw(17) << setprecision(7) << fixed << l->freq;
		cout << setw(11) << setprecision(2) << fixed << l->hi->energy/0.6950356;
		cout << endl;
	}
	cout << "!NUMBER OF COLL PARTNERS" << endl;
	// -1 means: do the molecule in LTE
	// 0 means: do the molecule in NLTE with guesses for the collision strengths
	cout << "-1" << endl;

	// free the memory
	for( s=states.begin(); s != states.end(); ++s )
		delete *s;
}

void InsertState( state** s, vector<state*>& v )
{
	vector<state*>::iterator p;

	for( p = v.begin(); p != v.end(); ++p )
		if( **p == **s )
			break;

	if( p == v.end() )
	{
		v.push_back( *s );
	}
	else
	{
		double cunc = sqrt( pow2((*s)->unc) + pow2((*p)->unc) );
		if( !(*s)->lgExperimental && !(*p)->lgExperimental &&
		    abs((*s)->energy - (*p)->energy) > cunc )
			throw runtime_error( "found 2 states with duplicate quantum numbers" );
		if( (*p)->glev != (*s)->glev )
			throw runtime_error( "statistical weight mismatch" );
		// copy over level energy if it comes from the code...
		if( (*p)->lgExperimental && !(*s)->lgExperimental )
		{
			(*p)->lgExperimental = (*s)->lgExperimental;
			(*p)->energy = (*s)->energy;
			(*p)->unc = (*s)->unc;
		}
		delete *s;
		*s = *p;
	}
}

template<class T>
inline void read_data(const string& str,
		      T& x)
{
	istringstream iss( str );
	iss >> x;
	if( iss.fail() || !iss.eof() )
		throw invalid_argument( "read_data() failed" );
}

void read_Jval(const string& str, int& x)
{
	if( str.length() != 2 )
		throw invalid_argument( "a quantum number must be exactly 2 characters long!" );

	if( !isdigit(str[1]) )
		throw invalid_argument( "illegal 2nd character found in quantum number" );

	if( str[0] == ' ' || str[0] == '-' || isdigit(str[0]) )
		read_data( str, x );
	else if( islower(str[0]) )
	{
		// special coding for negative quantum numbers
		int tens, units;
		tens = str[0] - 'a' + 1;
		units = str[1] - '0';
		x = -(tens*10 + units);
	}
	else if( isupper(str[0]) )
	{
		// special coding for positive quantum numbers
		int tens, units;
		tens = str[0] - 'A' + 10;
		units = str[1] - '0';
		x = tens*10 + units;
	}
	else
		throw invalid_argument( "illegal 1st character found in quantum number" );
}

inline string TrimSpaces(const string& str)
{
	string::size_type s1,s2;
	for( s1=0; s1 < str.length(); s1++ ) {
		if( str[s1] != ' ' )
			break;
	}
	for( s2=str.length()-1; s2 >= s1; s2-- ) {
		if( str[s2] != ' ' )
			break;
	}
	return ( s2 < s1 ) ? "" : str.substr( s1, s2-s1+1 );
}

inline bool FindAndErase(string& str,
			 const string& substr)
{
	return FindAndReplace( str, substr, "" );
}

inline bool FindAndReplace(string& str,
			   const string& substr,
			   const string& newstr)
{
	string::size_type ptr = str.find( substr );
	if( ptr != string::npos )
		str.replace( ptr, substr.length(), newstr );
	return ptr != string::npos;
}
