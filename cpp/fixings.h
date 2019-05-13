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

#ifndef _REDUKTI_FIXINGS_H_
#define _REDUKTI_FIXINGS_H_

#include <index.h>
#include <timeseries.h>

#include <map>
#include <memory>

namespace redukti
{

class FixingDataService
{
	private:
	std::map<IndexId, std::unique_ptr<TimeSeries>> fixings_;

	public:
	FixingDataService() {}
	~FixingDataService() {}

	void set_fixings(IndexId id, std::unique_ptr<TimeSeries> &&ts) { fixings_.emplace(id, std::move(ts)); }
	TimeSeries *get_fixings(IndexId id)
	{
		auto const &fixing = fixings_.find(id);
		if (fixing == fixings_.cend())
			return nullptr;
		return fixing->second.get();
	}

	void reset() { fixings_.clear(); }
};

int test_fixings();

} // namespace redukti

#endif
