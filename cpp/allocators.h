/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017-2019 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */
#ifndef _REDUKTI_ALLOCATORS_H
#define _REDUKTI_ALLOCATORS_H

#include <logger.h>

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <iterator>
#include <type_traits>

// IMPORTANT
//
// The allocators defined below are NOT thread safe
// You must ensure that an allocator (other than the
// MallocAllocator to be accurate) is never shared across
// threads
//
// Secondly these allocators are fine tuned to requirements
// in this project and are not general purpose.

namespace redukti
{

// Generic allocator interface
class Allocator
{
	public:
	virtual ~Allocator() = default;

	// Allocate at least size bytes
	// A size of 0 will result in nullptr being returned
	virtual void *allocate(size_t size) noexcept = 0;

	// Allocate at least size bytes
	// A size of 0 will result in nullptr being returned
	// Memory allocation failure will be treated as fatal error and the
	// program will be terminated
	// This is based on the philosophy that memory allocation
	// failures are unrecoverable and it is better to fail fast.
	void *safe_allocate(size_t size) noexcept
	{
		if (!size)
			return nullptr;
		void *p = allocate(size);
		if (p == nullptr) {
			die("Out of memory, exiting\n");
		}
		return p;
	}

	// Depending upon the type of allocator a deallocate may
	// not do anything
	virtual void deallocate(void *address) noexcept = 0;
};

// Utility for associating a deleter with a
// unique_ptr when memory was allocated using an allocator.
// When freeing up the object its destructor will be called.
// If the allocator is not null, then after the object is
// destroyed it will be passed to the deallocate() method
// of the allocator.
//
// Example:
//  Allocator *A;
//  std::unique_ptr<YieldCurve, Deleter<YieldCurve>>(
//	     new (*A) YieldCurve(), Deleter<YieldCurve>(A));
//
template <typename T> class Deleter
{
	public:
	explicit Deleter(Allocator *A = nullptr) : A_(A) {}

	void operator()(T *p)
	{
		if (p != nullptr) {
			// Call destructor
			p->~T();
			// Free memory
			assert(A_ != nullptr);
			if (A_ != nullptr)
				A_->deallocate(p);
		}
	}

	private:
	Allocator *A_;
};

// This allocator uses calloc/free
// This is the only allocator that is thread
// safe!
// Memory allocation falures are treated as fatal
// errors as they are generally unrecoverable anyway
// and it is better to fail fast.
class MallocAllocator : public Allocator
{
	public:
	MallocAllocator() = default;

	void *allocate(size_t size) noexcept final
	{
		if (size == 0)
			return nullptr;
		void *p = ::calloc(1, size);
		if (p == nullptr) {
			die("Out of memory, exiting\n");
		}
		return p;
	}

	void deallocate(void *address) noexcept final
	{
		if (address)
			::free(address);
	}

	public:
	MallocAllocator(const MallocAllocator &) = delete;

	MallocAllocator &operator=(const MallocAllocator &) = delete;
};

// Returns the default allocator which should normally be the
// malloc allocator
extern Allocator* get_default_allocator();

// Allocator interface where it is not necessary
// to destroy or free individual objects
//
// IMPORTANT
//
// Do not use for objects requiring destruction
//
class RegionAllocator : public Allocator
{
	public:
	// When a RegionAllocator is destroyed all memory allocated
	// may be released depending upon how the allocator
	// acquired that memory. User does not need to call
	// deallocate() explicitly on objects.
	// Note therefore that this allocator is unsuitable for
	// objects with destructors!
	~RegionAllocator() override = default;

	void *allocate(size_t size) noexcept override = 0;

	// Deallocate does nothing
	void deallocate(void *address) noexcept final {}

	// Resets the allocator so that all memory
	// is either freed and available for reuse
	virtual void release() noexcept = 0;
};

union alignment_type {
	char c;
	double d;
	uint64_t i;
	void *p;
};

// This is an allocator that returns memory from a fixed
// sized memory buffer. The buffer may be externally provided or
// owned. When the buffer is exhausted any allocation requests
// will fail and allocate() will return nullptr.
//
// As it is a RegionAllocator, deallocate() is a no-op
struct FixedRegionAllocator : public RegionAllocator {
	// memory externally supplied
	FixedRegionAllocator(char *start, char *end) noexcept
	    : data_(start), offset_(0), size_(end - start), delete_(false)
	{
		std::fill(data_, data_ + size_, 0);
	}

	// memory externally supplied
	FixedRegionAllocator(void *start, size_t n) noexcept
	    : data_(static_cast<char *>(start)), offset_(0), size_(n), delete_(false)
	{
		std::fill(data_, data_ + size_, 0);
	}

	// Acquire memory
	// Memory will be owned by this instance
	// It is a fatal error if memory cannot be allocated
	explicit FixedRegionAllocator(size_t n) noexcept
	{
		assert(n != 0);
		data_ = (char *)calloc(n, sizeof(char));
		if (!data_) {
			die("Out of memory, exiting\n");
		}
		size_ = n;
		offset_ = 0;
		delete_ = true;
	}

	// Base pointer
	const char *ptr() const noexcept { return data_; }

	// Current position
	size_t pos() const noexcept { return offset_; }

	// Sets current position
	// This is useful for scenarios where the user
	// wants to use the allocator in a stack like fashion
	// This is used by FixedRegionAllocatorGuard to
	// undo allocation upon destruction
	void pos(size_t i) noexcept
	{
		assert(i < size_);
		if (i < size_)
			offset_ = i;
	}

	size_t size() const noexcept { return size_; }

	// Memory is allocated from the buffer we
	// were given - we simply bump the current
	// position
	void *allocate(size_t size) noexcept final
	{
		if (!size)
			return nullptr;
		// round to alignment multiple
		constexpr size_t alignment = std::alignment_of<alignment_type>::value;
		auto allocsize = (size + alignment - 1u) & ~(alignment - 1u);
		assert(allocsize % alignment == 0);
		assert(allocsize >= size);
		assert(allocsize < size + alignment);
		// Is there enough space?
		if (remaining() < allocsize) {
			return nullptr;
		}
		void *result = static_cast<void *>(&data_[offset_]);
		offset_ += allocsize;
		return result;
	}

	// Returns remaining size of internal array
	// As memory may be allocated elsewhere this is
	// only useful for testing
	size_t remaining() const noexcept { return size() - pos(); }

	// Resets the memory used to 0 and makes
	// all memory available for reuse
	void release() noexcept final { offset_ = 0; }

	~FixedRegionAllocator() noexcept override
	{
		if (delete_ && data_) {
			::free(data_);
		}
	}

	public:
	FixedRegionAllocator(const FixedRegionAllocator &) = delete;

	FixedRegionAllocator &operator=(const FixedRegionAllocator &) = delete;

	private:
	char *data_;    // buffer to use for memory allocations
	size_t offset_; // Current position up to which memory is allocated,
	// when offset_ == size_ memory is exhausted
	size_t size_; // buffer size
	bool delete_; // if set the memory will be freed
};

// This allocator uses a fixed size array internal array
// if the requested memory allocation request can fit. Else
// memory is allocated via a fallback allocator. Note
// that the fallback allocator is a RegionAllocator,
// i.e. there is no deallocation or destruction of individual
// objects
//
// It is intended that this allocator is created on the stack
// which ensures that the internal array is stack allocated.
//
// Note that this allocator does not track memory that was
// allocated via the fallback allocator. Hence if the fallback
// allocator is used then any objects allocated with that
// will only be freed when the fallback allocator deallocates them.
template <int N> struct StackRegionWithFallbackAllocator : public RegionAllocator {
	explicit StackRegionWithFallbackAllocator(RegionAllocator *backup_allocator) noexcept
	    : my_allocator_(buf_, N), fallback_allocator_(backup_allocator), overflowed_(false)
	{
	}

	// Allocate at least size bytes
	// A size of 0 will result in nullptr
	// If allocation fails, it is treated as fatal error
	void *allocate(size_t size) noexcept final
	{
		if (!size)
			return nullptr;
		void *p = my_allocator_.allocate(size);
		if (!p && fallback_allocator_) {
			p = fallback_allocator_->allocate(size);
			overflowed_ = overflowed_ || p != nullptr;
			if (!p) {
				die("Out of memory, exiting\n");
			}
		}
		return p;
	}

	// Resets the internal array pointer
	void release() noexcept final
	{
		if (overflowed_) {
			warn("WARNING: possible memory leak\n");
		}
		my_allocator_.release();
	}

	public:
	StackRegionWithFallbackAllocator(const StackRegionWithFallbackAllocator &) = delete;

	StackRegionWithFallbackAllocator &operator=(const StackRegionWithFallbackAllocator &) = delete;

	private:
	FixedRegionAllocator my_allocator_;
	RegionAllocator *fallback_allocator_;
	char buf_[N]; // static memory
	bool overflowed_;
};

// Sequential buffer allocator for use cases where a bunch of small objects
// of varying sizes are needed, but all of the objects can then be thrown
// away once finished. This allocator uses a finite amount of memory -
// 32 * buffer size. By default buffer size is 5MB. If it runs out of memory
// the program will be aborted. Maximum individual allocation is
// restricted to buffer size. This allocator is designed for use by a single
// thread - i.e. it is not safe to share across threads.
class DynamicRegionAllocator : public RegionAllocator
{
	private:
	enum { N = 32 };
	FixedRegionAllocator *buffers_[N];
	int current_;
	size_t buffer_size_;
	size_t used_;
	size_t allocated_;

	public:
	explicit DynamicRegionAllocator(size_t buffer_size = 5 * 1024 * 1024) noexcept
	    : current_(0), buffer_size_(buffer_size), used_(0), allocated_(0)
	{
		std::fill(std::begin(buffers_), std::end(buffers_), nullptr);
	}

	void *allocate(size_t size) noexcept final
	{
		if (!size)
			return nullptr;
		if (size > buffer_size_) {
			die("attempt to allocate more than allowed size: requested %lld, max size %lld\n",
			    (long long)size, (long long)buffer_size_);
		}
	try_again:
		if (current_ >= N) {
			die("Out of memory: %lld requested\n", (long long)size);
		}
		if (!buffers_[current_]) {
			buffers_[current_] = ::new FixedRegionAllocator(buffer_size_);
			allocated_ += buffer_size_;
		}
		FixedRegionAllocator *buf_ = buffers_[current_];
		void *result = buf_->allocate(size);
		if (result == nullptr) {
			current_++;
			goto try_again;
		}
		used_ += size;
		return result;
	}

	void release() noexcept final
	{
		for (size_t i = 0; i < N; i++) {
			FixedRegionAllocator *buf_ = buffers_[i];
			if (buf_)
				buf_->release();
		}
		current_ = used_ = 0;
	}

	~DynamicRegionAllocator() noexcept override
	{
		for (size_t i = 0; i < N; i++) {
			FixedRegionAllocator *buf_ = buffers_[i];
			::delete buf_;
		}
	}

	public:
	DynamicRegionAllocator(const DynamicRegionAllocator &) = delete;

	DynamicRegionAllocator &operator=(const DynamicRegionAllocator &) = delete;
};

// This guard can be used to restore a FixedRegionAllocator to
// its previous allocation state. It relies on the fact that
// a FixedRegionAllocator is a bump the pointer allocator, and
// can be restored by simply resetting the pointer to the previous
// position
class FixedRegionAllocatorGuard
{
	public:
	explicit FixedRegionAllocatorGuard(FixedRegionAllocator *A) : A_(A)
	{
		// The the current position
		saved_pos_ = A_->pos();
	}

	~FixedRegionAllocatorGuard()
	{
		// Restore original position
		A_->pos(saved_pos_);
	}

	private:
	FixedRegionAllocator *A_;
	size_t saved_pos_;
};

// Each thread is given a set of allocators to use
// To obtain the thread specific allocator set call
// get_threadspecific_allocators().
struct AllocatorSet {
	RegionAllocator *cashflow_allocator;
	RegionAllocator *sensitivities_allocator;
	FixedRegionAllocator *tempspace_allocator;

	void reset()
	{
		cashflow_allocator->release();
		sensitivities_allocator->release();
		tempspace_allocator->release();
	}
};

// Retrieves the thread specific allocator set.
extern AllocatorSet *get_threadspecific_allocators();

extern int test_allocators();

// TODO move to utils
inline void string_copy(char *buf, const char *src, size_t buflen)
{
	if (buflen != 0) {
		strncpy(buf, src, buflen);
		buf[buflen - 1] = 0;
	}
}

} // namespace redukti

inline void *operator new(size_t n, redukti::Allocator &alloc) { return alloc.allocate(n); }

void *operator new(size_t n, redukti::Allocator *alloc) = delete;

inline void operator delete(void *ptr, redukti::Allocator &alloc) { alloc.deallocate(ptr); }

void operator delete(void *ptr, redukti::Allocator *alloc) = delete;

inline void *operator new[](size_t n, redukti::Allocator &alloc) { return alloc.allocate(n); }

inline void operator delete[](void *ptr, redukti::Allocator &alloc) { alloc.deallocate(ptr); }

#endif
