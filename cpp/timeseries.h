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

#ifndef _REDUKTI_TIMESERIES_H_
#define _REDUKTI_TIMESERIES_H_

#include <date.h>

#include <algorithm>
#include <memory>

namespace redukti
{

class Value
{
	private:
	Date date_;
	double value_;

	public:
	Value() : value_(0.0) {}
	Value(Date d, double v = 0.0) : date_(d), value_(v) {}
	Value(const Value &other) : date_(other.date_), value_(other.value_) {}
	Value &operator=(const Value &other)
	{
		if (this != &other) {
			date_ = other.date_;
			value_ = other.value_;
		}
		return *this;
	}
	Date date() const { return date_; }
	double value() const { return value_; }
};

struct ValueComparator {
	bool operator()(const Value &a, const Value &b) { return a.date() < b.date(); }
};

class TimeSeries
{
	private:
	std::vector<Value> values_;

	public:
	TimeSeries() {}
	TimeSeries(size_t num_values, Value *data) : values_(num_values) { add(num_values, data); }
	void add(size_t num_values, Value *data)
	{
		for (auto i = 0; i < num_values; i++)
			values_.push_back(data[i]);
		std::sort(std::begin(values_), std::end(values_), ValueComparator());
	}
	~TimeSeries() {}
	Value *begin() { return values_.data(); }
	const Value *cbegin() const { return values_.data(); }
	Value *end() { return values_.data() + values_.size(); }
	const Value *cend() const { return values_.data() + values_.size(); }
	bool find(Date d, double &value)
	{
		auto iter = std::lower_bound(std::begin(values_), std::end(values_), Value(d), ValueComparator());
		if (iter != std::end(values_) && iter->date() == d) {
			value = iter->value();
			return true;
		}
		value = 0.0; // avoid undefined value
		return false;
	}

	private:
	TimeSeries(const TimeSeries &) = delete;
	TimeSeries &operator=(const TimeSeries &) = delete;
};

} // namespace redukti
#endif
