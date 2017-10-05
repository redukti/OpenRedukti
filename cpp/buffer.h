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

#ifndef _REDUKTI_BUFFER_H_
#define _REDUKTI_BUFFER_H_

#include <algorithm>

#include <cassert>
#include <cstddef>
#include <cstdint>

namespace redukti
{

// Similar to vector<T> but
// do not grow dynamically
// and has a reference count
template <typename _T> struct buffer {
	typedef _T value_type;
	typedef _T *pointer;
	typedef const _T *const_pointer;
	typedef _T &reference;
	typedef const _T &const_reference;
	typedef pointer iterator;
	typedef const_pointer const_iterator;
	typedef pointer reverse_iterator;
	typedef const_pointer const_reverse_iterator;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

	explicit buffer(size_type n) noexcept
	{
		_buf = ::new _T[n];
		std::fill(_buf, _buf + n, _T());
		_size = n;
		_count = 1;
	}

	~buffer() noexcept { ::delete[] _buf; }

	reference operator[](size_type n) noexcept
	{
		assert(n < _size && n >= 0);
		return *(_buf + n);
	}

	const_reference operator[](size_type n) const noexcept
	{
		assert(n < _size && n >= 0);
		return *(_buf + n);
	}

	iterator begin() noexcept { return _buf; }

	const_iterator begin() const noexcept { return _buf; }

	iterator end() noexcept { return _buf + _size; }

	const_iterator end() const noexcept { return _buf + _size; }

	size_type size() const noexcept { return _size; }

	int _count;
	pointer _buf;
	size_type _size;

      private:
	buffer(const buffer &) = delete;
	buffer &operator=(const buffer &) = delete;
};

// reference counted buffer
// does not grow dynamically
template <typename _T> class buffer_ptr
{
      public:
	typedef _T value_type;
	typedef _T *pointer;
	typedef const _T *const_pointer;
	typedef _T &reference;
	typedef const _T &const_reference;
	typedef pointer iterator;
	typedef const_pointer const_iterator;
	typedef pointer reverse_iterator;
	typedef const_pointer const_reverse_iterator;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

	explicit buffer_ptr(size_type n) noexcept { _ref = new buffer<_T>(n); }
	~buffer_ptr() noexcept
	{
		if (--_ref->_count == 0)
			delete _ref;
	}
	buffer_ptr(const buffer_ptr &p) noexcept : _ref(p._ref) { ++_ref->_count; }
	buffer_ptr &operator=(const buffer_ptr &p) noexcept
	{
		if (this == &p)
			return *this;
		buffer<_T> *old = _ref;
		_ref = p._ref;
		++_ref->_count;
		if (--old->_count == 0)
			delete old;
		return *this;
	}
	reference operator[](size_type n) noexcept
	{
		buffer<_T> &b = *_ref;
		return b[n];
	}
	const_reference operator[](size_type n) const noexcept
	{
		const buffer<_T> &b = *_ref;
		return b[n];
	}
	iterator begin() noexcept { return _ref->begin(); }
	const_iterator begin() const noexcept { return _ref->begin(); }
	iterator end() noexcept { return _ref->end(); }
	const_iterator end() const noexcept { return _ref->end(); }
	size_type size() const noexcept { return _ref->size(); }
	int ref_count() const noexcept { return _ref->_count; }

      private:
	buffer<_T> *_ref;
};

} // namespace redukti

#endif
