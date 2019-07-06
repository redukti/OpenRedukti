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

#include <allocators.h>
#include <logger.h>

#include <cstring>
#include <mutex>
#include <thread>

namespace redukti
{

static MallocAllocator g_GlobalAllocator;
static Allocator *g_DefaultAllocator = &g_GlobalAllocator;

Allocator* get_default_allocator()
{
    return g_DefaultAllocator;
}

static thread_local AllocatorSet *g_PAllocatorSet = nullptr;

// Returns the thread specific allocators.
// If they do not exist they are created
// As these objects are specific to the thread there is no need to
// synchronize
AllocatorSet *get_threadspecific_allocators()
{
	if (g_PAllocatorSet == nullptr) {
		debug("Initializing thread specific allocators\n");
		g_PAllocatorSet = new AllocatorSet();
		g_PAllocatorSet->cashflow_allocator = new DynamicRegionAllocator(64 * 1024);
		g_PAllocatorSet->sensitivities_allocator = new DynamicRegionAllocator(256 * 1024);
		g_PAllocatorSet->tempspace_allocator = new FixedRegionAllocator(2 * 1024 * 1024);
	}
	return g_PAllocatorSet;
}

static void thread_func(AllocatorSet *main_thread_allocator_set, volatile int *result)
{
	*result = 0;
	auto this_thread_allocator_set = get_threadspecific_allocators();
	if (!this_thread_allocator_set) {
		printf("Failed to get AllocatorSet\n");
		*result = 1;
		return;
	}

	if (get_threadspecific_allocators() != this_thread_allocator_set) {
		printf("AllocatorSet mismatched\n");
		*result = 1;
		return;
	}

	if (this_thread_allocator_set == main_thread_allocator_set) {
		printf("AllocatorSet not thread specific\n");
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

	auto main_thread_allocator_set = get_threadspecific_allocators();
	if (!main_thread_allocator_set)
		return 1;

	if (get_threadspecific_allocators() != main_thread_allocator_set)
		return 1;

	int result = 2;
	std::thread t = std::thread(thread_func, main_thread_allocator_set, &result);
	t.join();

	if (result != 0)
		return 1;

	printf("Memory Allocations Tests OK\n");
	return 0;
}

} // namespace redukti
