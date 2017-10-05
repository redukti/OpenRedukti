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

#include <allocators.h>
#include <logger.h>

#include <cstring>
#include <mutex>
#include <thread>

namespace redukti
{

MallocAllocator GlobalAllocator;

static thread_local AllocatorSet *allocSet = nullptr;

// Returns the thread specific allocators.
// If they do not exist they are created
// As these objects are specific to the thread there is no need to
// synchronize
AllocatorSet *get_threadspecific_allocators()
{
	if (allocSet == nullptr) {
		inform("Initializing thread specific allocators\n");
		allocSet = new AllocatorSet();
		allocSet->cashflow_allocator = new DynamicRegionAllocator(64 * 1024);
		allocSet->sensitivities_allocator = new DynamicRegionAllocator(256 * 1024);
		allocSet->tempspace_allocator = new FixedRegionAllocator(2 * 1024 * 1024);
	}
	return allocSet;
}

static void thread_func(AllocatorSet *mainAlloc, volatile int *result)
{
	*result = 0;
	auto allocSet = get_threadspecific_allocators();
	if (!allocSet) {
		printf("Failed to get allocSet\n");
		*result = 1;
		return;
	}

	if (get_threadspecific_allocators() != allocSet) {
		printf("AllocSet mismatched\n");
		*result = 1;
		return;
	}

	if (allocSet == mainAlloc) {
		printf("AllocSet not thread specific\n");
		*result = 1;
		return;
	}
}

// Test memory allocators
int test_allocators()
{
	char stackmem[1024];
	FixedRegionAllocator balloc(stackmem, sizeof stackmem);
	if (balloc.allocate(0) != nullptr)
		return 1;
	if (balloc.remaining() != sizeof stackmem)
		return 1;
	void *p2 = balloc.allocate(63); // actual allocation should be 64
	if (balloc.remaining() != 1024 - 64)
		return 1;
	if (!p2)
		return 1;
	{
		FixedRegionAllocatorGuard guard(&balloc);
		if (!balloc.allocate(512 - 64))
			return 1;
		if (!balloc.allocate(512))
			return 1;
		if (balloc.remaining() != 0)
			return 1;
	}

	if (balloc.remaining() != 1024 - 64)
		return 1;

	auto allocSet = get_threadspecific_allocators();
	if (!allocSet)
		return 1;

	if (get_threadspecific_allocators() != allocSet)
		return 1;

	int result = 2;
	std::thread t = std::thread(thread_func, allocSet, &result);
	t.join();

	if (result != 0)
		return 1;

	printf("Memory Allocations Tests OK\n");
	return 0;
}

} // namespace redukti
