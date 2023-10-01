/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CO_Init called from cdInit to initialize co routines */
/*CO_update_chem_rates update rate coefficients, only temp part - in mole_co_etc.c */
#include "cdstd.h"
#include <cctype>
#include <string.h>
#include <algorithm>
#include <stdlib.h>
#include "cddefines.h"
#include "co.h"
#include "colden.h"
#include "conv.h"
#include "deuterium.h"
#include "doppvel.h"
#include "elementnames.h"
#include "h2.h"
#include "iso.h"
#include "phycon.h"
#include "physconst.h"
#include "mole.h"
#include "mole_priv.h"
#include "hmi.h"
#include "radius.h"
#include "rfield.h"
#include "rt.h"
#include "secondaries.h"
#include "dense.h"
#include "ionbal.h"
#include "grainvar.h"
#include "timesc.h"
#include "taulines.h"
#include "trace.h"
/*lint -e778 constant expression evaluates to 0 in operation '-' */

/* CO_update_chem_rates update rate coefficients, only temp part - in
 * mole_co_etc.c called in conv_base before any chemistry or
 * ionization is done */

enum spectype {MOLECULE, OTHER};

STATIC void read_species_file( string filename, bool lgCreateIsotopologues = true );
STATIC void newelement(const char label[], int ipion);
STATIC void newisotope( const count_ptr<chem_element> &el, int massNumberA,
	realnum mass_amu, double frac );
STATIC realnum MeanMassOfElement( const count_ptr<chem_element> &el );
// newspecies is overloaded.  The first one just calls the second with the additional lgCreateIsotopologues set true
STATIC molecule *newspecies(const char label[], enum spectype type,
	enum mole_state state, realnum form_enthalpy);
STATIC molecule *newspecies(const char label[], enum spectype type,
	enum mole_state state, realnum form_enthalpy, bool lgCreateIsotopologues);
STATIC count_ptr<chem_atom> findatom(const char buf[]);
STATIC bool isactive(const molecule &mol);
STATIC bool ispassive(const molecule &mol);
STATIC void ReadIsotopeFractions( const vector<bool>& lgResolveNelem );
//STATIC void create_isotopologues(count_ptr<molecule> mol, ChemAtomList& atoms, vector<int>& numAtoms);

namespace mole_priv {
	map <string,count_ptr<molecule> > spectab;
	map <string,count_ptr<mole_reaction> > reactab;
	map <string,count_ptr<chem_element> > elemtab;
	map <string,count_ptr<chem_atom> > atomtab;
	map <string,count_ptr<mole_reaction> > functab;
};
molecule *null_mole;
molezone *null_molezone;
chem_element *null_element;
chem_atom *null_atom;
vector< count_ptr <chem_element> > element_list;
ChemAtomList unresolved_atom_list;
ChemAtomList atom_list;
molecule **groupspecies;

#include <functional>

namespace
{
	class MoleCmp : public binary_function<const count_ptr<molecule>,
														const count_ptr<molecule>,bool>
	{
	public:
		bool operator()(const count_ptr<molecule> &mol1, 
							 const count_ptr<molecule> &mol2) const
			{
				return mol1->compare(*mol2) < 0;
			}
		bool operator()(const molecule *mol1, const molecule *mol2) const
			{
				return mol1->compare(*mol2) < 0;
			}
	};
}

void t_mole_global::sort(t_mole_global::MoleculeList::iterator start, 
					 t_mole_global::MoleculeList::iterator end)
{
	std::sort(start,end,MoleCmp());	
}
void t_mole_global::sort(molecule **start, molecule **end)
{
	std::sort(start,end,MoleCmp());	
}

STATIC void ReadIsotopeFractions( const vector<bool>& lgResolveNelem )
{
	DEBUG_ENTRY( "ReadIsotopeFractions()" );

	char chLine[INPUT_LINE_LENGTH];
	long i;
	bool lgEOL;

	static const char chFile[] = "isotope_fractions.dat";
	FILE *ioDATA = open_data( chFile, "r" );
	ASSERT( ioDATA != NULL );

	while( read_whole_line( chLine, (int)sizeof(chLine), ioDATA ) != NULL )
	{
		if( chLine[0] == '#' )
			continue;

		i = 1;
		long Z = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
		long A = (long)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);
		// file has fractionation as percent; the 0.01 converts to fraction
		double frac = 0.01 * (double)FFmtRead(chLine,&i,INPUT_LINE_LENGTH,&lgEOL);

		fixit(); // need real masses here instead of just A (protons + neutrons)
		if( (unsigned)Z <= lgResolveNelem.size() && lgResolveNelem[Z-1] )
			newisotope( element_list[Z-1], A, (realnum)A, frac );
		// always do this to continue history of predicting 13CO	
		else if( Z-1==ipCARBON )
			newisotope( element_list[Z-1], A, (realnum)A, frac );
	}

	fclose( ioDATA );

	return;
}

void t_mole_global::make_species(void)
{
	DEBUG_ENTRY( "mole_global::make_species()" );

	long int i;
	molecule *sp;

	null_element = new chem_element(0,"") ;
	null_atom = new( chem_atom );
	null_molezone = new( molezone );
	null_molezone->den = 0.;

	/* set up concordance of elemental species to external Cloudy indices */
	for (i=0;i<LIMELM;i++) 
	{
		newelement(elementnames.chElementSym[i], i);
	}

	// flip this to treat deuterium
	if( deut.lgElmtOn )
		lgTreatIsotopes[ipHYDROGEN] = true;

	// read and define isotopes, set default fractionations
	ReadIsotopeFractions( lgTreatIsotopes );

	// special code to maintain effect of "set 12C13C" command
	{
		count_ptr<chem_atom> atom12C = findatom("^12C");
		count_ptr<chem_atom> atom13C = findatom("^13C");
		
		if( co.C12_C13_isotope_ratio >= 0. )
		{
			double ratio = co.C12_C13_isotope_ratio;
			atom12C->frac = ratio / (ratio + 1.);
			atom13C->frac = 1. / (ratio + 1.);
		}
		else
		{
			double ratio = atom12C->frac / SDIV( atom13C->frac );
			co.C12_C13_isotope_ratio = ratio;
		}
		// Make sure only two isotopes are defined: 12C, 13C.
		// The mean-abundance isotope is defined below.
		ASSERT( element_list[ipCARBON]->isotopes.size() == 2 );
	}


	if( lgTreatIsotopes[ipHYDROGEN] )
	{
		SetDeuteriumFractionation( element_list[ipHYDROGEN]->isotopes[2]->frac );
		SetGasPhaseDeuterium( dense.gas_phase[ipHYDROGEN] );
		InitDeuteriumIonization();
	}

	for( long nelem=0; nelem<LIMELM; ++nelem )
	{
		realnum mass_amu = MeanMassOfElement( element_list[nelem] );
		//define generic mean-abundance isotopes
		newisotope( element_list[nelem], -1, mass_amu, 1.0 );
	}
	
	ASSERT( (long) unresolved_atom_list.size() == LIMELM );
	
	/* set up properties of molecular species -- chemical formulae,
		 array indices, elementary components (parsed from formula), 
		 status within CO network, location of stored value external 
		 to CO network (must be floating point). */

	/* Sizes of different parts of network are calculated by increments in newspecies */
	mole_global.num_total = mole_global.num_calc = 0;
	/* Enthalpies of formation taken from
	 * >>refer    mol     Le Teuff, Y. E., Millar, T. J., & Markwick, A. J.,2000, A&AS, 146, 157
	 */
	
	/* Zero density pseudo-species to return when molecule is switched off */
	null_mole = newspecies("      ",OTHER,MOLE_NULL, 0.);
	null_mole->index = -1;

	read_species_file( "chem_species.dat" );

	if(gv.lgDustOn() && mole_global.lgGrain_mole_deplete )
	{		
		/* What are formation enthalpies of COgrn, H2Ogrn, OHgrn?  For 
			 present, take grn as standard state, and neglect adsorption enthalpy.

			 -- check e.g. http://www.arxiv.org/abs/astro-ph/0702322 for CO adsorption energy.
		*/
		if (0)
		{
			read_species_file( "chem_species_grn.dat" );
		}
		else
		{
			sp = newspecies("COgrn ",MOLECULE,MOLE_ACTIVE, -113.8f);
			sp = newspecies("H2Ogrn",MOLECULE,MOLE_ACTIVE, -238.9f);
			sp = newspecies("OHgrn ",MOLECULE,MOLE_ACTIVE, 38.4f);
			//sp = newspecies("Hgrn ",MOLECULE,MOLE_ACTIVE, 216.f);
		}
	}

	/* Add passive species to complete network */
	sp = newspecies("e-    ",OTHER,MOLE_PASSIVE, 0.0f);
	sp->charge = -1;	sp->mole_mass = (realnum)ELECTRON_MASS; /* Augment properties for this non-molecular species */
	sp = newspecies("grn   ",OTHER,MOLE_PASSIVE, 0.0f);
	sp = newspecies("PHOTON",OTHER,MOLE_PASSIVE, 0.0f);
	sp = newspecies("CRPHOT",OTHER,MOLE_PASSIVE, 0.0f);
	sp = newspecies("CRP   ",OTHER,MOLE_PASSIVE, 0.0f);

	if (!mole_global.lgLeidenHack)
		sp = newspecies("H-    ",MOLECULE,MOLE_ACTIVE, 143.2f);
	if (hmi.lgLeiden_Keep_ipMH2s) 
	{
		sp = newspecies("H2*   ",MOLECULE,MOLE_ACTIVE, 
							 h2.ENERGY_H2_STAR * KJMOL1CM);  
	}

	if( deut.lgElmtOn )
	{
		read_species_file( "chem_species_deuterium.dat", false );
	}

	/* Add species for all other elements and their first ions -- couple to network at least via H- */
	for( ChemAtomList::iterator atom = unresolved_atom_list.begin();
		  atom != unresolved_atom_list.end(); ++atom) 
	{
		long int nelem = (*atom)->el->Z-1;

		for( long ion=0; ion<=nelem+1; ion++ )
		{
			char str[CHARS_ISOTOPE_SYM+CHARS_ION_STAGE+1];
			char temp[CHARS_ION_STAGE+1];
			if( ion==0 )
				temp[0] = '\0';
			else if( ion==1 )
				sprintf( temp, "+" );
			else
				sprintf( temp, "+%ld", ion );
			sprintf( str, "%s%s", (*atom)->label().c_str(), temp );
			if (findspecies(str) == null_mole)
			{
				sp = newspecies(str,MOLECULE,MOLE_ACTIVE, 0.f);
				fixit(); // populate these in a local update
#if 0
				if( sp != NULL )
				{
					sp->levels = NULL;
					sp->numLevels = 0;
				}
#endif
			}
		}
	}

	return;
}

void mole_make_list()
{
	DEBUG_ENTRY( "mole_make_list()" );

	/* Create linear list of species and populate it... */
	mole_global.list.resize(mole_global.num_total);

	/* ...first active species */
	long int i = 0;
	for (molecule_i p = mole_priv::spectab.begin(); p != mole_priv::spectab.end(); ++p) 
	{
		if (isactive(*(p->second)))
			mole_global.list[i++] = p->second;
	}
	ASSERT (i == mole_global.num_calc); 

	/* Sort list into a standard ordering */
	t_mole_global::sort(mole_global.list.begin(),mole_global.list.begin()+mole_global.num_calc);

	for (molecule_i p = mole_priv::spectab.begin(); p != mole_priv::spectab.end(); ++p) 
	{
		if (ispassive(*(p->second)))
			mole_global.list[i++] = p->second;
	}
	ASSERT (i == mole_global.num_total); 

	/* Set molecule indices to order of list just created */
	for(i=0;i<mole_global.num_total;i++) 
	{
		mole_global.list[i]->index = i;
	}

	for(i=0;i<mole_global.num_total;i++)
	{
		if( !mole_global.list[i]->parentLabel.empty() )
		{
			long parentIndex = findspecies( mole_global.list[i]->parentLabel.c_str() )->index;
			mole_global.list[i]->parentIndex = parentIndex;
		}
		else
			mole_global.list[i]->parentIndex = -1;
	}

	/* Register the atomic ladders */
	for(i=0;i<mole_global.num_total;i++) 
	{
		molecule *sp = &(*mole_global.list[i]);
		if (sp->isMonatomic())
		{
			ASSERT( (int)sp->nAtom.size() == 1 );
			count_ptr<chem_atom> atom = sp->nAtom.begin()->first;
			ASSERT(sp->charge <= atom->el->Z);
			if(sp->charge >= 0 && sp->lgGas_Phase) 
			{
				atom->ipMl[sp->charge] = i;
			}
		}
	}

	return;

}

STATIC void read_species_file( string filename, bool lgCreateIsotopologues )
{
	DEBUG_ENTRY( "read_sepcies_file()" );
	
	fstream ioDATA;
	open_data( ioDATA, filename.c_str(), mode_r );
	string line;

	while( getline( ioDATA,line ) )
	{
		if( line.empty() )
			break;
		if( line[0] == '#' )
			continue;
		istringstream iss( line );
		string species;
		double formation_enthalpy;
		iss >> species;
		iss >> formation_enthalpy;
		ASSERT( iss.eof() );
		newspecies( species.c_str(), MOLECULE,MOLE_ACTIVE, formation_enthalpy, lgCreateIsotopologues );
		//fprintf( ioQQQ, "DEBUGGG read_species_file %s\t%f\n", species.c_str(), formation_enthalpy );
	}

	return;
}
/*lint +e778 constant expression evaluates to 0 in operation '-' */

void create_isotopologues_one(
	ChemAtomList& atoms,
	vector< int >& numAtoms,
	string atom_old,
	string atom_new,
	string embellishments,
	vector<string>& newLabels )
{
	fixit(); // make sure atom_new and atom_old are isotopes
	fixit(); // make sure atom_new is not already present

	//for( ChemAtomList::iterator it = atoms.begin(); it != atoms.end(); ++it )
	for( unsigned position = 0; position < atoms.size(); ++position )
	{
		stringstream str;
		if( atoms[position]->label() == atom_old )
		{
			for( unsigned i=0; i<position; ++i )
			{
				str << atoms[i]->label();
				if( numAtoms[i]>1 )
					str << numAtoms[i];
			}
			
			if( numAtoms[position] > 1 )
			{
				str << atom_old;
				if( numAtoms[position] > 2 )
					str << numAtoms[position]-1;
			}
			

			if( position+1 == atoms.size() )
				str << atom_new;

			for( unsigned i=position+1; i<atoms.size(); ++i )
			{
				if( i==position+1 )
				{
					// remove adjacent duplicates
					if( atom_new == atoms[i]->label() )
					{
						str << atoms[i]->label();
						ASSERT( numAtoms[i] + 1 > 1 );
						str << numAtoms[i] + 1;
					}
					else
					{
						str << atom_new;
						str << atoms[i]->label();
						if( numAtoms[i] > 1 )
							str << numAtoms[i];
					}
				}
				else
				{
					str << atoms[i]->label();
					if( numAtoms[i] > 1 )
						str << numAtoms[i];
				}
			}

			// add on charge, grn, and excitation embellishments
			str << embellishments;

			newLabels.push_back( str.str() );
			//fprintf( ioQQQ, "DEBUGGG create_isotopologues_one %s\n", newLabels.back().c_str() );
		}
	}

	return;
}

/* Fill element linking structure */
STATIC void newelement(const char label[], int ipion)
{
	char *s;

	DEBUG_ENTRY("newelement()");

	/* Create private workspace for label; copy and strip trailing whitespace */
	int len = strlen(label)+1;
	auto_vec<char> mylab_v(new char[len]);
	char *mylab = mylab_v.data();
	strncpy(mylab,label,len);
	s = strchr(mylab,' ');
	if (s)
		*s = '\0';

	int exists = (mole_priv::elemtab.find(mylab) != mole_priv::elemtab.end());
	if (!exists)
	{
		count_ptr<chem_element> element(new chem_element(ipion+1,mylab));
		mole_priv::elemtab[element->label] = element;
		element_list.push_back(element);
	}
	return;
}

/* Fill isotope lists */
STATIC void newisotope( const count_ptr<chem_element>& el, int massNumberA, 
	realnum mass_amu, double frac )
{

	DEBUG_ENTRY("newisotope()");

	ASSERT( massNumberA < 3 * LIMELM && ( massNumberA > 0 || massNumberA == -1 ) );
	ASSERT( mass_amu < 3. * LIMELM && mass_amu > 0. );
	ASSERT( frac <= 1. + FLT_EPSILON && frac >= 0. );

	count_ptr<chem_atom> isotope(new chem_atom);
	isotope->A = massNumberA;
	isotope->mass_amu = mass_amu;
	isotope->frac = frac;
	isotope->ipMl.resize(el->Z+1);
	isotope->el = el.get_ptr();
	for (long int ion = 0; ion < el->Z+1; ion++)
		isotope->ipMl[ion] = -1; 	/* Chemical network species indices not yet defined */

	//int exists = (mole_priv::atomtab.find( isotope->label() ) != mole_priv::atomtab.end());
	mole_priv::atomtab[ isotope->label() ] = isotope;
	atom_list.push_back(isotope);
	if( isotope->lgMeanAbundance() )
		unresolved_atom_list.push_back(isotope);
	// register with 'parent' element
	el->isotopes[massNumberA] = isotope;
}

STATIC realnum MeanMassOfElement( const count_ptr<chem_element>& el )
{
	DEBUG_ENTRY("MeanMassOfElement()");

	// if no isotopes have been defined, just use the mean mass stored elsewhere
	if( el->isotopes.size()==0 )
		return dense.AtomicWeight[el->Z-1];

	realnum averageMass = 0., fracsum = 0.;
	for( isotopes_i it = el->isotopes.begin(); it != el->isotopes.end(); ++it )
	{
		fracsum += it->second->frac;
		averageMass += it->second->mass_amu * it->second->frac;
	}
	ASSERT( fp_equal( fracsum, 1.f ) );

	return averageMass;
}

STATIC molecule *newspecies(
	const char label[],
	enum spectype type,
	enum mole_state state,
	realnum form_enthalpy)
{
	return newspecies( label, type, state, form_enthalpy, true);
}

/* Parse species string to find constituent atoms, charge etc. */
STATIC molecule *newspecies(
	const char label[], 
	enum spectype type, 
	enum mole_state state, 
	realnum form_enthalpy,
	bool lgCreateIsotopologues )/* formation enthalpy at 0K */
{
	DEBUG_ENTRY("newspecies()");

	int exists;
	ChemAtomList atomsLeftToRight;
	vector< int > numAtoms;
	string embellishments;
	char *s;
	count_ptr<molecule> mol(new molecule);
	
	mol->parentLabel.clear();
	mol->isEnabled = true;
	mol->charge = 0;
	mol->lgExcit = false;
	mol->mole_mass = 0.0;
	mol->state = state;
	mol->lgGas_Phase = true;
	mol->form_enthalpy = form_enthalpy;
	mol->groupnum = -1;

	/* Create private workspace for label; copy and strip trailing whitespace */
	int len = strlen(label)+1;
	auto_vec<char> mylab_v(new char[len]);
	char *mylab = mylab_v.data();
	strncpy(mylab,label,len);
	s = strchr(mylab,' ');
	if (s)
		*s = '\0';
	mol->label = mylab;

	if(type == MOLECULE)
	{
		if( parse_species_label( mylab, atomsLeftToRight, numAtoms, embellishments, mol->lgExcit, mol->charge, mol->lgGas_Phase ) == false )
			return NULL;		
		
		for( unsigned i = 0; i < atomsLeftToRight.size(); ++i )
		{
			mol->nAtom[ atomsLeftToRight[i] ] += numAtoms[i];
			mol->mole_mass += numAtoms[i] * atomsLeftToRight[i]->mass_amu * (realnum)ATOMIC_MASS_UNIT;
		}
	}
	
	// we also kill H- if molecules are disabled.  This is less than ideal,
	// but physically motivated by the fact that one of the strongest H- sinks
	// involves formation of H2 (disabled by "no molecules"), while one the 
	// fastest sources is e- recombination (which would still be allowed).
	if ( (mol->n_nuclei() > 1 || (mol->isMonatomic() && mol->charge==-1) ) && mole_global.lgNoMole) 
	{
		if( trace.lgTraceMole )
			fprintf(ioQQQ,"No species %s as molecules off\n",label);
		return NULL;
	}

	if (mol->n_nuclei() > 1 && mole_global.lgNoHeavyMole)
	{
		for( nAtoms_ri it=mol->nAtom.rbegin(); it != mol->nAtom.rend(); --it )
		{
			if( it->first->el->Z-1 != ipHYDROGEN )
			{
				ASSERT( it->second > 0 );
				if( trace.lgTraceMole )
					fprintf(ioQQQ,"No species %s as heavy molecules off\n",label);
				return NULL;
			}
		}
	}

	/* Insert species into hash table */
	exists = (mole_priv::spectab.find(mol->label) != mole_priv::spectab.end());
	if( exists )
	{
		fprintf( ioQQQ,"Warning: duplicate species %s - using first one\n", 
					mol->label.c_str() );
		return NULL;
	}
	else
		mole_priv::spectab[mol->label] = mol;

	// all map entries should have strictly positive number of nuclei
	for( nAtoms_i it=mol->nAtom.begin(); it != mol->nAtom.end(); ++it )
		ASSERT( it->second > 0 );

	if (state != MOLE_NULL)
		mole_global.num_total++;
	if (state == MOLE_ACTIVE)
		mole_global.num_calc++;

	// this is a special case to always treat 13CO (as has long been done) 
	if( lgCreateIsotopologues && type==MOLECULE && mol->label.compare("CO")==0 ) 
	{
		molecule *sp = newspecies( "^13CO", MOLECULE, mol->state, mol->form_enthalpy, false );
		sp->parentLabel = mol->label;
	}

	// create singly-substituted isotopologues	
	if( lgCreateIsotopologues && type==MOLECULE && !mol->isMonatomic() )
	{
		for( nAtoms_i it1 = mol->nAtom.begin(); it1 != mol->nAtom.end(); ++it1 )
		{
			for( map<int, count_ptr<chem_atom> >::iterator it2 = it1->first->el->isotopes.begin(); 
				it2 != it1->first->el->isotopes.end(); ++it2 )
			{
				// we don't want to create ^1H isotopologues (only substitute D for H)
				if( it1->first->el->Z-1 == ipHYDROGEN && it2->second->A != 2 )
					continue;
				if( !mole_global.lgTreatIsotopes[it1->first->el->Z-1] )
					continue;
				if( it2->second->lgMeanAbundance() )
					continue;
				vector<string> isoLabs;
				create_isotopologues_one( atomsLeftToRight, numAtoms, it1->first->label(), it2->second->label(), embellishments, isoLabs );
				//fprintf( ioQQQ, " DEBUGGG %10s isotopologues of %10s:", it1->first->label().c_str(), mol->label.c_str() );
				//for( vector<string>::iterator lab = isoLabs.begin(); lab != isoLabs.end(); ++ lab )
				//	fprintf( ioQQQ, " %10s", lab->c_str() );
				//fprintf( ioQQQ, "\n" );
				for( vector<string>::iterator newLabel = isoLabs.begin(); newLabel != isoLabs.end(); ++newLabel )
				{
					molecule *sp = newspecies( newLabel->c_str(), MOLECULE, mol->state, mol->form_enthalpy, false );
					// D is special -- don't set parentLabel
					if( sp!=NULL && it2->second->A != 2 )
						sp->parentLabel = mol->label;
				}
			}
		}
	}
		
	return &(*mol);
}
bool parse_species_label( const char label[], ChemAtomList &atomsLeftToRight, vector<int> &numAtoms, string &embellishments ) 
{
	bool lgExcit, lgGas_Phase; 
	int charge;
	bool lgOK = parse_species_label( label, atomsLeftToRight, numAtoms, embellishments, lgExcit, charge, lgGas_Phase );
	return lgOK; 
}
bool parse_species_label( const char label[], ChemAtomList &atomsLeftToRight, vector<int> &numAtoms, string &embellishments, 
	bool &lgExcit, int &charge, bool &lgGas_Phase )
{
	long int i, n, ipAtom;
	char thisAtom[CHARS_ISOTOPE_SYM];
	count_ptr<chem_atom> atom;
	char mylab[CHARS_SPECIES];
	char *s;

	strncpy( mylab, label, CHARS_SPECIES );
	
	/* Excitation... */
	s = strpbrk(mylab,"*");
	if(s)
	{
		lgExcit = true;
		embellishments = s;
		*s = '\0';
	} 

	/* ...Charge */
	s = strpbrk(mylab,"+-");
	if(s)
	{
		if(isdigit(*(s+1))) 
			n = atoi(s+1);
		else
			n = 1;
		if(*s == '+')
			charge = n;
		else
			charge = -n;
		embellishments = s + embellishments;
		*s = '\0';
	}
	/* ...Grain */
	s = strstr(mylab,"grn");
	if(s) 
	{
		lgGas_Phase = false;
		embellishments = s + embellishments;
		*s = '\0';
	} 
	else 
	{
		lgGas_Phase = true;
	}
	//fprintf( ioQQQ, "DEBUGGG parse_species_label %s\t%s\t%s\n", label, mylab, embellishments.c_str() );

	/* Now analyse chemical formula */
	i = 0;
	while (mylab[i] != '\0' && mylab[i] != ' ' && mylab[i] != '*') 
	{
		/* Select next element in species, matches regexp [A-Z][a-z]? */
		ipAtom = 0;
		/* look for isotope prefix */
		if(mylab[i]=='^')
		{
			thisAtom[ipAtom++] = mylab[i++];
			ASSERT( isdigit(mylab[i]) );
			thisAtom[ipAtom++] = mylab[i++];
			if(isdigit(mylab[i]))
			{
				thisAtom[ipAtom++] = mylab[i++];
			}
		}
		// should be first character of an element symbol
		thisAtom[ipAtom++] = mylab[i++];
		if(islower(mylab[i])) 
		{
			thisAtom[ipAtom++] = mylab[i++];
		}
		thisAtom[ipAtom] = '\0';
		ASSERT(ipAtom <= CHARS_ISOTOPE_SYM);

		atom = findatom(thisAtom);
		if(atom.get_ptr() == NULL) 
		{
			fprintf(stderr,"Did not recognize atom at %s in \"%s \"[%ld]\n",thisAtom,mylab,i);
			exit(-1);
		}
		if(!dense.lgElmtOn[atom->el->Z-1])
		{
			if( trace.lgTraceMole )
				fprintf(ioQQQ,"No species %s as element %s off\n",mylab,atom->el->label.c_str() );
			return false;
		}
			
		if(isdigit(mylab[i])) /* If there is >1 of this atom */
		{
			n = 0;
			do {
				n = 10*n+(long int)(mylab[i]-'0');
				i++;
			} while (isdigit(mylab[i]));
		}
		else
		{
			n = 1;
		}
		atomsLeftToRight.push_back( atom );
		numAtoms.push_back( n );
	}

	return true;
}
STATIC bool isactive(const molecule &mol)
{
	DEBUG_ENTRY("isactive()");
	return mol.state == MOLE_ACTIVE;
}
STATIC bool ispassive(const molecule &mol)
{

	DEBUG_ENTRY("ispassive()");
	return mol.state == MOLE_PASSIVE;
}

bool lgDifferByExcitation( const molecule &mol1, const molecule &mol2 )
{
	if( mol1.label == mol2.label + "*" )
		return true;
	else if( mol2.label == mol1.label + "*" )
		return true;
	else 
		return false;
}

molecule *findspecies(const char buf[])
{
	DEBUG_ENTRY("findspecies()");

	// strip string of the first space and anything after it
	string s;
	for (const char *pb = buf; *pb && *pb != ' '; ++pb)
	{
		s += *pb;
	}

	const molecule_i p = mole_priv::spectab.find(s);

	if(p != mole_priv::spectab.end()) 
		return &(*p->second);
	else 
		return null_mole;
}

molezone *findspecieslocal(const char buf[])
{
	DEBUG_ENTRY("findspecieslocal()");

	// strip string of the first space and anything after it
	string s;
	for (const char *pb = buf; *pb && *pb != ' '; ++pb)
	{
		s += *pb;
	}

	const molecule_i p = mole_priv::spectab.find(s);

	if(p != mole_priv::spectab.end()) 
		return &mole.species[ p->second->index ];
	else 
		return null_molezone;
}

STATIC count_ptr<chem_atom> findatom(const char buf[])
{
	chem_atom_i p;

	DEBUG_ENTRY("findatom()");

	p = mole_priv::atomtab.find(buf);

	if(p != mole_priv::atomtab.end())  
		return p->second;
	else 
		return count_ptr<chem_atom>(NULL);
}

void mole_update_species_cache(void)
{
	int i;
	double den_times_area, den_grains, adsorbed_density;

	DEBUG_ENTRY("mole_update_species_cache()");

	enum { DEBUG_MOLE = false };

	for (i=0;i<mole_global.num_total;i++) 
	{
		if(mole.species[i].location != NULL) 
		{
			ASSERT( mole_global.list[i]->parentLabel.empty() );
			mole.species[i].den = *(mole.species[i].location);
			if (DEBUG_MOLE)
				fprintf(ioQQQ,"%s: %f\n",mole_global.list[i]->label.c_str(),mole.species[i].den);
		}
	}
	
	mole.set_isotope_abundances();
 
	/* For rates that are dependent on grain physics.  This includes grain density, 
	cross sectional area, and dust temperature of each constituent.  Note that 
	
	gv.bin[nd]->IntArea/4.*gv.bin[nd]->cnv_H_pCM3
	
	is the integrated projected grain surface area per cm^3 of gas for each grain size bin */

	/* >>chng 06 feb 28, turn off this rate when no grain molecules */
	/* >>chng 06 dec 05 rjrw: do this in newreact rather than rate */
	if( gv.lgDustOn() )
	{
		den_times_area = den_grains = 0.0;
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			/* >>chng 06 mar 04, update expression for projected grain surface area, PvH */
			den_times_area += gv.bin[nd]->IntArea/4.*gv.bin[nd]->cnv_H_pCM3;
			den_grains += gv.bin[nd]->cnv_GR_pCM3;
		}
		
		adsorbed_density = 0.0;
		for (i=0;i<mole_global.num_total;i++) 
		{
			if( !mole_global.list[i]->lgGas_Phase && mole_global.list[i]->parentLabel.empty() ) 
				adsorbed_density += mole.species[i].den;
		}

		mole.grain_area = (realnum) den_times_area;
		mole.grain_density = (realnum) den_grains;

		double mole_cs = 1e-15;
		if (4*den_times_area <= mole_cs*adsorbed_density)
			mole.grain_saturation = 1.0;
		else
			mole.grain_saturation = (realnum)((mole_cs*adsorbed_density)/(4.*den_times_area));
	}
	else
	{
		mole.grain_area = 0.0;
		mole.grain_density = 0.0;
		mole.grain_saturation = 1.0;
	}
	if (DEBUG_MOLE)
		fprintf(ioQQQ,"Dust: %f %f %f\n",
						mole.grain_area,mole.grain_density,mole.grain_saturation);

}

realnum mole_return_cached_species(const GroupMap &) 
// Finding the total atom density from MoleMap.molElems w
{
	int i;

	/* These two updates should together maintain the abundance invariant */

	// Assert invariant
	ASSERT(lgElemsConserved());

	// Update total of non-ladder species
	dense.updateXMolecules();
	total_molecule_deut( deut.xMolecules );

	/* charge on molecules */
	mole.elec = 0.;
	for(i=0;i<mole_global.num_calc;i++)
	{
		if (mole.species[i].location == NULL && mole_global.list[i]->parentLabel.empty())
			mole.elec += mole.species[i].den*mole_global.list[i]->charge;
	}

	// Update ionization ladder species
	realnum delta = 0.0;
	long ncpt = 0;

	for (i=0;i<mole_global.num_total;i++) 
	{
		if(mole.species[i].location && mole_global.list[i]->state == MOLE_ACTIVE) 
		{
			realnum new_pop = mole.species[i].den;

			if( mole_global.list[i]->isMonatomic() )
			{
				realnum old_pop = *(mole.species[i].location);
				long nelem = mole_global.list[i]->nAtom.begin()->first->el->Z-1;
				realnum frac = (new_pop-old_pop)/SDIV(new_pop+old_pop+1e-8*
																  dense.gas_phase[nelem]);
				delta += frac*frac;
				++ncpt;
			}

			//if( iteration >= 3 && nzone >= 100  )
			//	fprintf( ioQQQ, "DEBUGGG mole_return_ %i\t%s\t%.12e\t%.12e\t%.12e\t%.12e\t%li\n", 
			//	i, mole_global.list[i]->label.c_str(), new_pop, old_pop, frac, delta, ncpt );
			*(mole.species[i].location) = new_pop;
		}
	}

	// Assert invariant
	ASSERT(lgElemsConserved());
	return ncpt>0 ? sqrt(delta/ncpt) : 0.f;
}

void t_mole_local::set_isotope_abundances( void )
{
	DEBUG_ENTRY("t_mole_local::set_isotope_abundances()");

	// loop over unresolved elements
	for(ChemAtomList::iterator atom = unresolved_atom_list.begin(); atom != unresolved_atom_list.end(); ++atom)
	{
		// loop over all isotopes of each element
		for( isotopes_i it = (*atom)->el->isotopes.begin(); it != (*atom)->el->isotopes.end(); ++it )
		{
			// skip mean-abundance "isotopes"
			if( it->second->lgMeanAbundance() )
				continue;
			
			// loop over all ions of the isotope
			for( unsigned long ion=0; ion<it->second->ipMl.size(); ++ion )
			{
				if( it->second->ipMl[ion] != -1 && 
					(species[ it->second->ipMl[ion] ].location == NULL ) && (*atom)->ipMl[ion] != -1 )
				{
					species[ it->second->ipMl[ion] ].den = species[ (*atom)->ipMl[ion] ].den * it->second->frac;
				}
			}
		}
	}
	
	return;
}

void t_mole_local::set_location( long nelem, long ion, double *density )
{
	DEBUG_ENTRY( "t_mole_local::set_location()" );

	ASSERT( nelem < LIMELM );
	ASSERT( ion < nelem + 2 );
	long mole_index = unresolved_atom_list[nelem]->ipMl[ion];
	// element not enabled if index is -1
	if( mole_index == -1 )
		return;
	ASSERT( mole_index < mole_global.num_total );
	species[mole_index].location = density;

	return;
}

void total_molecule_deut( realnum &total_f )
{
	DEBUG_ENTRY( "total_molecule_deut()" );
	
	double total = 0.;

	if( !deut.lgElmtOn )
		return;
	
	for (long int i=0;i<mole_global.num_calc;++i) 
	{
		if (mole.species[i].location == NULL && mole_global.list[i]->parentLabel.empty() )
		{
			for( molecule::nAtomsMap::iterator atom=mole_global.list[i]->nAtom.begin();
				  atom != mole_global.list[i]->nAtom.end(); ++atom)
			{
				long int nelem = atom->first->el->Z-1;
				if( nelem==0 && atom->first->A==2 )
				{
					total += mole.species[i].den*atom->second;
				}
			}
		}
	}

	total_f = (realnum)total;

	return;
}
void total_molecule_elems(realnum total[LIMELM])
{
	DEBUG_ENTRY( "total_molecule_elems()" );

	/* now set total density of each element locked in molecular species */
	for ( long int nelem=ipHYDROGEN;nelem<LIMELM; ++nelem )
	{
		total[nelem] = 0.;
	}
	for (long int i=0;i<mole_global.num_calc;++i) 
	{
		if (mole.species[i].location == NULL && mole_global.list[i]->parentLabel.empty() )
		{
			for( molecule::nAtomsMap::iterator atom=mole_global.list[i]->nAtom.begin();
				  atom != mole_global.list[i]->nAtom.end(); ++atom)
			{
				ASSERT( atom->second > 0 );
				long int nelem = atom->first->el->Z-1;
				if( atom->first->lgMeanAbundance() )
					total[nelem] += (realnum) mole.species[i].den*atom->second;
			}
		}
	}
}
void total_network_elems(double total[LIMELM])
{
	DEBUG_ENTRY( "total_network_elems()" );

	/* now set total density of each element locked in molecular species */
	for ( long int nelem=ipHYDROGEN;nelem<LIMELM; ++nelem )
	{
		total[nelem] = 0.;
	}
	for (long int i=0;i<mole_global.num_calc;++i) 
	{
		if (mole_global.list[i]->parentLabel.empty())
		{
			for( molecule::nAtomsMap::iterator atom=mole_global.list[i]->nAtom.begin();
				  atom != mole_global.list[i]->nAtom.end(); ++atom)
			{
				long int nelem = atom->first->el->Z-1;
				total[nelem] += (realnum) mole.species[i].den*atom->second;
			}
		}
	}
}
realnum total_molecules(void)
{
	long int i;
	realnum total;

	DEBUG_ENTRY( "total_molecules()" );

	total = 0.;
	for (i=0;i<mole_global.num_calc;++i) 
	{
		if (mole.species[i].location == NULL && mole_global.list[i]->parentLabel.empty())
			total += (realnum) mole.species[i].den;
	}
	return total;
}
realnum total_molecules_gasphase(void)
{
	long int i;
	realnum total;

	DEBUG_ENTRY( "total_molecules_gasphase()" );

	total = 0.;
	for (i=0;i<mole_global.num_calc;++i) 
	{
		if (mole_global.list[i]->lgGas_Phase && mole.species[i].location == NULL && mole_global.list[i]->parentLabel.empty())
			total += (realnum) mole.species[i].den;
	}
	return total;
}
/*lint +e778 const express eval to 0 */
/*lint +e725 expect positive indentation */

void mole_make_groups(void)
{
	long int i, j;
	/* Neutrals and positive ions will be treated as single species inside 
		 molecular equilibrium solver, to facilitate coupling with ionization
		 solvers */
	DEBUG_ENTRY ("mole_make_groups()");
	if (mole_global.num_calc == 0)
	{
		groupspecies = NULL;
		mole_global.num_compacted = 0;
		return;
	}
	groupspecies = (molecule **) MALLOC(mole_global.num_calc*sizeof(molecule *));
	for (i=0,j=0;i<mole_global.num_calc;i++) 
	{
		if( mole_global.list[i]->parentLabel.empty() && ( !mole_global.list[i]->isMonatomic() || mole_global.list[i]->charge <= 0 || ! mole_global.list[i]->lgGas_Phase ) )
		{
			/* Compound molecules and negative ions are represented individually */
			mole_global.list[i]->groupnum = j;
			groupspecies[j++] = &(*mole_global.list[i]);
		}
		else
		{
			/* All positive ions are collapsed into single macrostate (i.e. H+ -> H) */
			/* Need to increase constant if higher ions are included */
			ASSERT (mole_global.list[i]->charge < LIMELM+1);
			ASSERT (mole_global.list[i]->groupnum == -1);
		}
	}
	mole_global.num_compacted = j;
	groupspecies = (molecule **) REALLOC((void *)groupspecies,
		mole_global.num_compacted*sizeof(molecule *));

	for (i=0;i<mole_global.num_calc;i++) 
	{
		if (mole_global.list[i]->groupnum == -1)
		{
			if( mole_global.list[i]->isMonatomic() && mole_global.list[i]->parentLabel.empty() )
			{
				for( nAtoms_i it = mole_global.list[i]->nAtom.begin(); it != mole_global.list[i]->nAtom.end(); ++it )
				{
					ASSERT( it->second > 0 );
					mole_global.list[i]->groupnum = mole_global.list[it->first->ipMl[0]]->groupnum;
					break;
				}
			}
			else
			{
				ASSERT( !mole_global.list[i]->parentLabel.empty() );
				mole_global.list[i]->groupnum = mole_global.list[ mole_global.list[i]->parentIndex ]->groupnum;
			}
		}
	
		ASSERT( mole_global.list[i]->groupnum != -1);
	}

	return;
}

