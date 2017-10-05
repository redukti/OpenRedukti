/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */
/**
 * Portions Copyright (C) 1994-2016 Lua.org, PUC-Rio.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ******************************************************************************/

#ifndef _REDUKTI_HASHTABLE_H
#define _REDUKTI_HASHTABLE_H

#include <cassert>
#include <cinttypes>
#include <cstdio>

#include <functional>
#include <string>

namespace redukti
{

// The hashtable implementation below is based upon
// the Lua table implementation. Its main features are
// a) it stores the data in an array
// b) the key and value types must be concrete types
// c) you cannot delete keys
// d) both key and value must have a special nil value

template <typename T> struct ValueConverter;
template <typename T> struct NilHelper;

// Wraps the key type
template <typename K, typename KNH> class TKey
{
	int32_t next_; // index of next node in the hash chain, -1 == null
	K tvk_;

      public:
	TKey() : next_(-1) { setnil(); }
	bool isnil() const { return KNH::is_nil(&tvk_); }
	void setnil() { KNH::set_nil(&tvk_); }
	void setkey(const K &k) { tvk_ = k; }
	const K &key() const { return tvk_; }
	int32_t next() const { return next_; }
	void setnext(int32_t n)
	{
		assert(n >= -1);
		next_ = n;
	}
};

template <typename V, typename VNH> class TValue
{
	V value_;

      public:
	TValue() { setnil(); }
	bool isnil() const { return VNH::is_nil(&value_); }
	void setnil() { VNH::set_nil(&value_); }
	void setval(const V &v) { value_ = v; }
	const V &val() const { return value_; }
};

template <typename K, typename V, typename KNH, typename VNH> struct Node {
	typedef TKey<K, KNH> key_type;
	typedef TValue<V, VNH> value_type;

	key_type i_key;
	value_type i_val;

	const K &key() const { return i_key.key(); }
	const V &value() const { return i_val.val(); }
};

template <typename node_type> class NodeIterator : public std::iterator<std::input_iterator_tag, node_type>
{
	node_type *p_;
	node_type *last_;

      public:
	NodeIterator(node_type *p, node_type *last) : p_(p), last_(last) {}
	NodeIterator(const NodeIterator &other) : p_(other.p_), last_(other.last_) {}
	NodeIterator &operator++()
	{
		if (p_ < last_)
			++p_;
		adjust();
		return *this;
	}
	NodeIterator operator++(int)
	{
		NodeIterator tmp(*this);
		if (p_ < last_)
			operator++();
		adjust();
		return tmp;
	}
	bool operator==(const NodeIterator &rhs) const { return p_ == rhs.p_; }
	bool operator!=(const NodeIterator &rhs) const { return p_ != rhs.p_; }
	node_type &operator*() { return *p_; }
	node_type *operator->() { return p_; }
	// moves the current node to next non-nil node
	NodeIterator &adjust()
	{
		while (p_ < last_ && (p_->i_key.isnil() || p_->i_val.isnil()))
			p_++;
		return *this;
	}
};

template <typename K, typename V, typename KNH = NilHelper<K>, typename VNH = NilHelper<V>,
	  typename HashFunc = std::hash<K>>
struct Table {
	typedef Node<K, V, KNH, VNH> node_type;
	typedef typename node_type::key_type key_type;
	typedef typename node_type::value_type value_type;
	typedef NodeIterator<node_type> iterator;

      private:
	node_type *node;    // start of the node vector
	int32_t lastfree;   // points to one past last available free node
	uint32_t lsizenode; // size of the node vector stored as power of 2

	node_type dummy_node;
	HashFunc hash;

	// modulus - when size is power of 2
	int32_t lmod(int32_t s, int32_t size) const { return s & (size - 1); }

	// determine the size of the node vector
	// following is equivalent to 2^lsizenode
	int32_t sizenode() const
	{
		int32_t n = 1 << lsizenode;
		assert(n >= 0);
		return n;
	}

	// get node at position i
	node_type *gnode(int32_t i) const
	{
		assert(i >= 0);
		return &node[i];
	}

	/*
	** returns the `main' position of an element in a table (that is, the
	*index
	** of its hash value)
	*/
	int32_t mainposition(const K &key) const { return lmod(hash(key), sizenode()); }

	// Gets an unused node
	int32_t getfreepos()
	{
		// Note lastfree points to one past the last free position
		while (lastfree > 0) {
			lastfree--;
			if (gnode(lastfree)->i_key.isnil())
				return lastfree;
		}
		return -1;
	}

	bool isdummy(node_type *n) const { return n == &dummy_node; }

	enum { MAXBITS = 30 };

	// Counts the number of used elements in the
	// table
	int32_t numusehash() const
	{
		int32_t totaluse = 0;
		int32_t i = sizenode();
		while (i-- > 0) {
			node_type *n = gnode(i);
			if (!n->i_val.isnil()) {
				totaluse++;
			}
		}
		return totaluse;
	}

	int32_t ceillog2(uint32_t x) const
	{
		static const unsigned char log_2[256] = {
		    0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
		int32_t l = 0;
		x--;
		while (x >= 256) {
			l += 8;
			x >>= 8;
		}
		return l + log_2[x];
	}

	// Allocates the node vector to required size
	// Note that the elements in node vector are
	// initialised to default values - i.e. nil
	void setnodevector(int32_t size)
	{
		int32_t lsize;
		if (size == 0) { /* no elements to hash part? */
			node = &dummy_node;
			lsize = 0;
		} else {
			if (size < 4)
				size = 4;
			lsize = ceillog2(size);
			assert(lsize <= MAXBITS);
			size = 1 << lsize;
			node = new node_type[size];
		}
		lsizenode = (unsigned char)lsize;
		lastfree = size;
	}

	// resize the hash table and rehash the elements
	void resize(int32_t nhsize)
	{
		int32_t oldhsize = lsizenode;
		node_type *nold = node; /* save old hash */
		/* create new hash part with appropriate size */
		setnodevector(nhsize);
		/* re-insert elements from hash part */
		int32_t n = (1 << oldhsize) - 1; // index of last node in old vector
		for (int32_t i = n; i >= 0; i--) {
			node_type *old = nold + i;
			if (!old->i_val.isnil()) {
				// hash to new position and copy value
				*set(old->i_key.key()) = old->i_val;
			}
		}
		if (!isdummy(nold))
			delete[] nold;
	}

	value_type *set(const K &key)
	{
		value_type *p = get(key);
		if (p != nullptr)
			return p;
		else
			return newkey(key);
	}

	void rehash(K ek)
	{
		int32_t totaluse = numusehash() + 1; // count in-use keys in the table, add 1
		// resize the table to new computed sizes
		resize(totaluse);
	}

      public:
	/*
	** main search function - if key is not found
	** nullptr is returned
	*/
	value_type *get(const K &key)
	{
		int32_t pos = mainposition(key);
		do {
			auto n = gnode(pos);
			if (!n->i_key.isnil() && n->i_key.key() == key)
				return &n->i_val;
			else {
				pos = n->i_key.next();
			}
		} while (pos != -1);
		return nullptr;
	}

	// Insert a new key - if key already exists then
	// the value node associated with the key is returned, else
	// the key is inserted and a new value node is returned
	// If necessary the table is grown.
	value_type *newkey(const K &key)
	{
		int32_t mainpos = mainposition(key);
		node_type *mp = gnode(mainpos);
		if (!mp->i_val.isnil() || isdummy(mp)) {
			node_type *othern;
			int32_t freepos = getfreepos(); // get a free place
			if (freepos == -1) {		// cannot find a free place?
				rehash(key);		// grow table
				return set(key);	// insert key into grown table
			}
			node_type *n = gnode(freepos);
			// mp is the colliding node
			assert(!isdummy(n));
			// get the main node at the key associated with mp
			int32_t otherpos = mainposition(mp->i_key.key());
			othern = gnode(otherpos);
			if (otherpos != mainpos) { // is colliding node out of
						   // its main position?
				// yes; move colliding node into free position
				// as 'othern' points to the main position
				// and 'mp' is in the chain starting at 'othern'
				// so search the chain until we get to the node
				// before 'mp' objective is to find previous
				// node to 'mp'
				while (othern->i_key.next() != mainpos) {
					otherpos = othern->i_key.next();
					othern = gnode(otherpos);
				}
				// Now we need to free 'mp' so link othern->next
				// = n
				othern->i_key.setnext(freepos); // redo the
								// chain with
								// 'n' in place
								// of mp
				*n = *mp;			// copy colliding node into free pos.
								// (mp->next also goes)
				mp->i_val.setnil();		// now 'mp' is free
			} else {				// colliding node is in its own main position
				// new node will go into free position
				// insert n into the chain
				n->i_key.setnext(mp->i_key.next());
				mp->i_key.setnext(freepos);
				mp = n;
			}
		}
		mp->i_key.setkey(key);
		assert(mp->i_val.isnil());
		return &mp->i_val;
	}

	Table(int32_t initsize = 0)
	{
		assert(initsize >= 0);
		setnodevector(initsize);
	}

	~Table()
	{
		if (!isdummy(node))
			delete[] node;
	}

	iterator begin()
	{
		// return iterator adjusted to first non-nil node
		return iterator(gnode(0), gnode(sizenode())).adjust();
	}

	iterator end() { return iterator(gnode(sizenode()), gnode(sizenode())); }

	void dump()
	{
		typedef ValueConverter<K> kc;
		typedef ValueConverter<V> vc;
		printf("dump of table\n");
		char buf[80];
		for (int32_t i = 0; i < sizenode(); i++) {
			node_type *n = gnode(i);
			if (n->i_key.isnil()) {
				printf("[%d] key is nil\n", i);
			} else {
				printf("[%d] key = %s", i, kc::as_string(buf, sizeof buf, n->i_key.key()));
				printf(" mainpos = %d", mainposition(n->i_key.key()));
				printf(" next = %d", n->i_key.next());
				printf(" value = %s\n", vc::as_string(buf, sizeof buf, n->i_val.val()));
			}
		}
	}
};

template <> struct ValueConverter<int> {
	static const char *as_string(char *buf, size_t n, int i)
	{
		snprintf(buf, n, "%d", i);
		return buf;
	}
};

template <> struct ValueConverter<double> {
	static const char *as_string(char *buf, size_t n, double i)
	{
		snprintf(buf, n, "%0.6f", i);
		return buf;
	}
};

template <> struct ValueConverter<std::string> {
	static const char *as_string(char *buf, size_t n, const std::string &i)
	{
		snprintf(buf, n, "%.*s", (int)(n - 1), i.c_str());
		return buf;
	}
};

extern int test_hashtable();

} // namespace redukti

#endif
