/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef QUANTUMSTATE_H_
#define QUANTUMSTATE_H_
#include "energy.h"
#include "proxy_iterator.h"

class quantumStateLabels
{
	char m_chLabel[5];
	char m_chConfig[11];
public:
	char *chLabel()
	{
		return m_chLabel;
	}
	const char *chLabel() const
	{
		return m_chLabel;
	}
	char *chConfig()
	{
		return m_chConfig;
	}
	const char *chConfig() const
	{
		return m_chConfig;
	}
};

// Quantum state properties as structure of lists -- packing
// like physical properties in vector storage is more cache efficient
// than a list-of-structures layout.
class qStateProxy;
class qStateConstProxy;
void Junk(qStateProxy);
void Zero(qStateProxy);

class qList
{
	vector<quantumStateLabels> m_labels;
	vector<double> m_ConBoltz;
	vector<double> m_Boltzmann;
	vector<Energy> m_energy;
	vector<realnum> m_g;
	vector<long> m_j;
	vector<long> m_J;
	vector<int> m_IonStg;
	vector<int> m_nelem;
	vector<long> m_l;
	vector<double> m_lifetime;
	vector<long> m_n;
	vector<double> m_ColDen;
	vector<double> m_Pop;
	vector<double> m_DestCollBB;
	vector<double> m_DestPhotoBB;
	vector<double> m_CreatCollBB;
	vector<long> m_S;
	vector<long> m_v;
	friend class qStateProxy;
	friend class qStateConstProxy;
public:
	typedef ProxyIterator<qStateProxy,qStateConstProxy> iterator;
	typedef ProxyIterator<qStateConstProxy,qStateConstProxy> const_iterator;
	typedef qStateProxy reference;
	typedef qStateConstProxy const_reference;
	explicit qList()
	{
		resize(0);
	}
	explicit qList(size_t i)
	{
		resize(i);
	}
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	reference operator[](int i);
	const_reference operator[](int i) const;
	// Must resize *all* the list members
	void resize(size_t i)
	{
		size_t old_size = size();
		m_labels.resize(i);
		m_ConBoltz.resize(i);
		m_Boltzmann.resize(i);
		m_energy.resize(i);
		m_g.resize(i);
		m_IonStg.resize(i);
		m_j.resize(i);
		m_J.resize(i);
		m_lifetime.resize(i);
		m_l.resize(i);
		m_n.resize(i);
		m_nelem.resize(i);
		m_ColDen.resize(i);
		m_Pop.resize(i);
# ifdef USE_NLTE7
		m_DestCollBB.resize(i);
		m_DestPhotoBB.resize(i);
		m_CreatCollBB.resize(i);
# endif
		m_S.resize(i);
		m_v.resize(i);
		for (size_t n=old_size; n<i; ++n)
		{
			reset(n);
		}
	}
	void reset(int n);
	// The size of any of the lists will do, as they should all be
	// the same.  No point really in asserting this, as the assert
	// will be at least as fragile as the original resize.
	size_t size() const
	{
		return m_labels.size();
	}
};

// Quantum state proxy object.  This is used to give access to
// structure-of-lists class qList in an 'object-like' manner, e.g.
//
// qs = List[element]; qs.IonStg() = ???; qs.S() = ???; update(qs);
//
// The member functions of the proxy class must be the same as the 
// lists which are included in the qList class.
class qStateConstProxy;
class qStateProxy
{
public:
	typedef qList list_type;
	typedef ProxyIterator<qStateProxy,qStateConstProxy> iterator;
private:
	friend class  ProxyIterator<qStateProxy,qStateConstProxy>;
	list_type *m_list;
	int m_index;
public:
   explicit qStateProxy(list_type* list, int index) : 
	m_list(list), m_index(index) {}
   explicit qStateProxy(void) : m_list(NULL), m_index(0) {}
// Proxy functions for members of qList below
	char *chLabel() const
	{
		return m_list->m_labels[m_index].chLabel();
	}
	char *chConfig() const
	{
		return m_list->m_labels[m_index].chConfig();
	}
	/** energy of the state */
	Energy &energy() const
	{
		return m_list->m_energy[m_index];
	}
	/** statistical weight [dimensionless] */
	realnum  &g() const
	{
		return m_list->m_g[m_index];
	}
	/** population of state [cm-3] */
	double  &Pop() const
	{
		return m_list->m_Pop[m_index];
	}
	/** population of state [cm-3] */
	double  &ColDen() const
	{
		return m_list->m_ColDen[m_index];
	}
# ifdef USE_NLTE7
	/** Destruction Rate due to collisional processes bound-bound */
	double &DestCollBB() const
	{
		return m_list->m_DestCollBB[m_index];
	}
	/** Destruction Rate due to photo processes bound-bound */
	double &DestPhotoBB() const
	{
		return m_list->m_DestPhotoBB[m_index];
	}
	/** Creation Rate due to collisional processes bound-bound */
	double &CreatCollBB() const
	{
		return m_list->m_CreatCollBB[m_index];
	}
# endif
	/** ion stage of element, 1 for atom, 2 ion, etc */
	int &IonStg() const
	{
		return m_list->m_IonStg[m_index];
	}
	/** atomic number of element, 1 for H, 2 for He, etc */
	int &nelem() const
	{
		return m_list->m_nelem[m_index];
	}
	/** ConBoltz excit to continuum */
	double &ConBoltz() const
	{
		return m_list->m_ConBoltz[m_index];
	}
	/** Boltzmann to ground state */
	double &Boltzmann() const
	{
		return m_list->m_Boltzmann[m_index];
	}
	/** Lifetime of the state */
	double &lifetime() const
	{
		return m_list->m_lifetime[m_index];
	}
	long &n() const
	{
		return m_list->m_n[m_index];
	}
	long &l() const
	{
		return m_list->m_l[m_index];
	}
	long &S() const
	{
		return m_list->m_S[m_index];
	}
	long &v() const
	{
		return m_list->m_v[m_index];
	}
	long &j() const
	{
		return m_list->m_j[m_index];
	}
	long &J() const
	{
		return m_list->m_J[m_index];
	}
};
class qStateConstProxy
{
public:
	typedef const qList list_type;
	typedef ProxyIterator<qStateConstProxy,qStateConstProxy> iterator;
private:
	friend class ProxyIterator<qStateConstProxy,qStateConstProxy>;
	const list_type *m_list;
	int m_index;
public:
	explicit qStateConstProxy(const list_type* list, int index) : 
	m_list(list), m_index(index) {}
	explicit qStateConstProxy(void) : m_list(NULL), m_index(0) {}
// Proxy functions for members of qList below
	const char *chLabel() const
	{
		return m_list->m_labels[m_index].chLabel();
	}
	const char *chConfig() const
	{
		return m_list->m_labels[m_index].chConfig();
	}
	Energy energy() const
	{
		return m_list->m_energy[m_index];
	}
	realnum g() const
	{
		return m_list->m_g[m_index];
	}
	double ColDen() const
	{
		return m_list->m_ColDen[m_index];
	}
	double Pop() const
	{
		return m_list->m_Pop[m_index];
	}
# ifdef USE_NLTE7
	double DestCollBB() const
	{
		return m_list->m_DestCollBB[m_index];
	}
	double DestPhotoBB() const
	{
		return m_list->m_DestPhotoBB[m_index];
	}
	double CreatCollBB() const
	{
		return m_list->m_CreatCollBB[m_index];
	}
# endif
	int IonStg() const
	{
		return m_list->m_IonStg[m_index];
	}
	int nelem() const
	{
		return m_list->m_nelem[m_index];
	}
	double ConBoltz() const
	{
		return m_list->m_ConBoltz[m_index];
	}
	double Boltzmann() const
	{
		return m_list->m_Boltzmann[m_index];
	}
	double lifetime() const
	{
		return m_list->m_lifetime[m_index];
	}
	long n() const
	{
		return m_list->m_n[m_index];
	}
	long l() const
	{
		return m_list->m_l[m_index];
	}
	long S() const
	{
		return m_list->m_S[m_index];
	}
	long v() const
	{
		return m_list->m_v[m_index];
	}
	long j() const
	{
		return m_list->m_j[m_index];
	}
	long J() const
	{
		return m_list->m_J[m_index];
	}
};

inline qList::iterator qList::begin()
{
	return iterator(this,0);
}
inline qList::const_iterator qList::begin() const
{
	return const_iterator(this,0);
}
inline qList::iterator qList::end()
{
	return iterator(this,m_labels.size());
}
inline qList::const_iterator qList::end() const
{
	return const_iterator(this,m_labels.size());
}
inline qList::reference qList::operator[](int i)
{
	return begin()[i];
}
inline qList::const_reference qList::operator[](int i) const
{
	return begin()[i];
}
inline void qList::reset(int n)
{
	Junk((*this)[n]);
	Zero((*this)[n]);
}
	

#endif
