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

#include <index.h>
#include <index.pb.h>

#include <dayfractions.h>

#include <map>
#include <mutex>

namespace redukti
{

struct InterestRateIndexTemplate {
	IsdaIndex isda_index;
	IndexFamily index_family;
	Currency currency;
	Tenor tenor;
	Tenor tenor_threshold;
	int fixing_lag;
	BusinessDayConvention short_tenor_bdc;
	BusinessDayConvention long_tenor_bdc;
	bool eom;
	const Calendar *fixing_calendar;
	const Calendar *value_date_calendar;
	const Calendar *index_calendar;
	const DayFraction *day_fraction;

	InterestRateIndexTemplate(IsdaIndex isda_index, IndexFamily index_family, Currency currency, Tenor tenor,
				  Tenor short_tenor_threshold, int fixing_lag, BusinessDayConvention short_tenor_bdc,
				  BusinessDayConvention long_tenor_bdc, bool eom, const Calendar *fixing_calendar,
				  const Calendar *value_date_calendar, const Calendar *index_calendar,
				  const DayFraction *day_fraction)
	    : isda_index(isda_index), index_family(index_family), currency(currency), tenor(tenor),
	      tenor_threshold(short_tenor_threshold), fixing_lag(fixing_lag), short_tenor_bdc(short_tenor_bdc),
	      long_tenor_bdc(long_tenor_bdc), eom(eom), fixing_calendar(fixing_calendar),
	      value_date_calendar(value_date_calendar), index_calendar(index_calendar), day_fraction(day_fraction)
	{
	}
};

class InterestRateIndexImpl : public InterestRateIndex
{
	private:
	IndexId id_;
	IsdaIndex isda_index_;
	IndexFamily index_family_;
	Currency currency_;
	Tenor tenor_;
	int fixing_lag_;
	BusinessDayConvention bdc_;
	bool eom_;
	const Calendar *fixing_calendar_;
	const Calendar *value_date_calendar_;
	const Calendar *index_calendar_;
	const DayFraction *day_fraction_;

	public:
	InterestRateIndexImpl(IndexId id, IsdaIndex isda_index, IndexFamily index_family, Currency currency,
			      Tenor tenor, int fixing_lag, BusinessDayConvention bdc, bool eom,
			      const Calendar *fixing_calendar, const Calendar *value_date_calendar,
			      const Calendar *index_calendar, const DayFraction *day_fraction)
	    : id_(id), isda_index_(isda_index), index_family_(index_family), currency_(currency), tenor_(tenor),
	      fixing_lag_(fixing_lag), bdc_(bdc), eom_(eom), fixing_calendar_(fixing_calendar),
	      value_date_calendar_(value_date_calendar), index_calendar_(index_calendar), day_fraction_(day_fraction)
	{
	}
	Currency currency() const override final { return currency_; };
	IndexFamily family() const override final { return index_family_; }
	Tenor tenor() const override final { return tenor_; }
	IsdaIndex isda_index() const override final { return isda_index_; }
	virtual IndexId id() const override final { return id_; }

	Date fixing_date(Date accrual_start_date) const override final
	{
		// Following moves actual business days so dayconvention is not
		// relevant
		return fixing_calendar_->advance(accrual_start_date, -fixing_lag_, PeriodUnit::DAYS);
	}

	Date value_date(Date fixing_date) const override final
	{
		Date d = fixing_calendar_->advance(fixing_date, fixing_lag_, PeriodUnit::DAYS);
		// If the value date is not a good business day as per
		// value date calendar then adjust
		if (value_date_calendar_->is_holiday(d)) {
			d = value_date_calendar_->adjust(d, bdc_);
		}
		return d;
	}

	Date maturity_date(Date value_date) const override final
	{
		return index_calendar_->advance(value_date, Period::tenor_to_period(tenor_), bdc_, eom_);
	}
	bool is_valid_fixing_date(Date date) const override final { return fixing_calendar_->is_businessday(date); }
	const Calendar *fixing_calendar() const override final { return fixing_calendar_; }
	const DayFraction *day_fraction() const override final { return day_fraction_; }
	BusinessDayConvention day_convention() const override final { return bdc_; }
};

class IndexServiceImpl : public IndexService
{
	private:
	std::mutex lock_;
	std::map<IndexId, std::unique_ptr<InterestRateIndexTemplate>> cached_templates_;
	std::map<IndexId, std::unique_ptr<InterestRateIndexImpl>> cached_indices_;
	std::map<int, IsdaIndex> index_family_to_isdaindex_;

	private:
	public:
	bool register_index(const IndexDefinition &definition) override final
	{
		auto calendar_service = get_calendar_factory();
		auto fixing_calendar = build_calendar(calendar_service, definition.fixing_calendars(),
						      definition.fixing_calendars_join_rule());
		auto value_date_calendar = build_calendar(calendar_service, definition.value_date_calendars(),
							  definition.value_date_calendars_join_rule());
		auto index_calendar = build_calendar(calendar_service, definition.index_calendars(),
						     definition.index_calendars_join_rule());
		if (fixing_calendar == nullptr) {
			fprintf(stderr, "Error: fixing calendar must be specified\n");
			return false;
		}
		if (value_date_calendar == nullptr)
			value_date_calendar = fixing_calendar;
		if (index_calendar == nullptr)
			index_calendar = fixing_calendar;
		const DayFraction *dfc = get_day_fraction(definition.day_count_fraction());
		if (dfc == nullptr) {
			fprintf(stderr, "Error: Day count fraction must be specified\n");
			return false;
		}
		std::unique_ptr<InterestRateIndexTemplate> index =
		    std::unique_ptr<InterestRateIndexTemplate>(new InterestRateIndexTemplate(
			definition.isda_index(), definition.index_family(), definition.currency(), definition.tenor(),
			definition.short_tenor_threshold(), definition.fixing_lag(),
			definition.short_tenor_convention(), definition.long_tenor_convention(), definition.eom(),
			fixing_calendar, value_date_calendar, index_calendar, dfc));
		std::lock_guard<std::mutex> guard(lock_);
		cached_templates_.emplace(make_index_id(definition.isda_index(), definition.tenor()), std::move(index));
		if (definition.default_for_index_family()) {
			int family_id = definition.currency() | definition.index_family() << 8;
			if (index_family_to_isdaindex_.find(family_id) == index_family_to_isdaindex_.end())
				index_family_to_isdaindex_[family_id] = definition.isda_index();
		}
		return true;
	}

	InterestRateIndex *get_index(IsdaIndex isda_index, Tenor tenor) override final
	{
		std::lock_guard<std::mutex> guard(lock_);
		IndexId id = make_index_id(isda_index, tenor);
		auto const &iter1 = cached_indices_.find(id);
		if (iter1 != cached_indices_.end())
			return iter1->second.get();
		InterestRateIndexTemplate *index_template = nullptr;
		auto const &iter2 = cached_templates_.find(id);
		if (iter2 != cached_templates_.end())
			index_template = iter2->second.get();
		if (index_template == nullptr) {
			auto const &iter3 = cached_templates_.find(make_index_id(isda_index, TENOR_UNSPECIFIED));
			if (iter3 != cached_templates_.end())
				index_template = iter3->second.get();
			else
				return nullptr;
		}
		BusinessDayConvention bdc = index_template->long_tenor_bdc;
		if (tenor <= index_template->tenor_threshold)
			bdc = index_template->short_tenor_bdc;
		std::unique_ptr<InterestRateIndexImpl> index = std::unique_ptr<InterestRateIndexImpl>(
		    new InterestRateIndexImpl(id, isda_index, index_template->index_family, index_template->currency,
					      tenor, tenor == TENOR_1D ? 0 : index_template->fixing_lag, bdc,
					      tenor >= TENOR_1M ? index_template->eom : false,
					      index_template->fixing_calendar, index_template->value_date_calendar,
					      index_template->index_calendar, index_template->day_fraction));
		InterestRateIndex *ptr = index.get();
		cached_indices_.emplace(id, std::move(index));
		return ptr;
	}

	InterestRateIndex *get_index(Currency currency, IndexFamily index_family, Tenor tenor) override final
	{
		int family_id = currency | index_family << 8;
		IsdaIndex index = ISDA_INDEX_UNSPECIFIED;
		{
			std::lock_guard<std::mutex> guard(lock_);
			auto iter = index_family_to_isdaindex_.find(family_id);
			if (iter != index_family_to_isdaindex_.end())
				index = iter->second;
		}
		if (index == ISDA_INDEX_UNSPECIFIED)
			return nullptr;
		return get_index(index, tenor);
	}
};

std::unique_ptr<IndexService> default_index_service;
std::mutex default_index_service_lock;

IndexService *get_default_index_service()
{
	std::lock_guard<std::mutex> guard(default_index_service_lock);
	if (!default_index_service) {
		default_index_service = std::make_unique<IndexServiceImpl>();
		IndexDefinition definition;

		definition.set_isda_index(IsdaIndex::USD_Federal_Funds_H_15_OIS_COMPOUND);
		definition.set_index_family(IndexFamily::FEDFUND);
		definition.set_currency(Currency::USD);
		definition.set_tenor(TENOR_1D);
		definition.set_fixing_lag(-1);
		definition.set_short_tenor_threshold(TENOR_2W);
		definition.set_short_tenor_convention(BusinessDayConvention::FOLLOWING);
		definition.set_long_tenor_convention(BusinessDayConvention::MODIFIED_FOLLOWING);
		definition.add_fixing_calendars(USNY);
		definition.set_fixing_calendars_join_rule(JOIN_HOLIDAYS);
		definition.add_value_date_calendars(USNY);
		definition.set_value_date_calendars_join_rule(JOIN_HOLIDAYS);
		definition.add_index_calendars(USNY);
		definition.set_index_calendars_join_rule(JOIN_HOLIDAYS);
		definition.set_day_count_fraction(DayCountFraction::ACT_360);
		definition.set_eom(false);
		definition.set_default_for_index_family(true);
		default_index_service->register_index(definition);

		definition.set_tenor(TENOR_UNSPECIFIED);
		default_index_service->register_index(definition);

		definition.set_isda_index(IsdaIndex::USD_LIBOR_BBA);
		definition.set_index_family(IndexFamily::LIBOR);
		definition.set_currency(Currency::USD);
		definition.set_tenor(TENOR_1D);
		definition.set_fixing_lag(0);
		definition.set_short_tenor_threshold(TENOR_2W);
		definition.set_short_tenor_convention(BusinessDayConvention::FOLLOWING);
		definition.set_long_tenor_convention(BusinessDayConvention::MODIFIED_FOLLOWING);
		definition.clear_fixing_calendars();
		definition.add_fixing_calendars(GBLO);
		definition.clear_value_date_calendars();
		definition.add_value_date_calendars(GBLO);
		definition.add_value_date_calendars(USNY);
		definition.clear_index_calendars();
		definition.add_index_calendars(GBLO);
		definition.add_index_calendars(USNY);
		definition.set_default_for_index_family(true);
		default_index_service->register_index(definition);

		definition.set_tenor(TENOR_UNSPECIFIED);
		definition.set_fixing_lag(2);
		definition.set_eom(true);
		default_index_service->register_index(definition);

		definition.set_isda_index(IsdaIndex::GBP_LIBOR_BBA);
		definition.set_index_family(IndexFamily::LIBOR);
		definition.set_currency(Currency::GBP);
		definition.set_tenor(TENOR_UNSPECIFIED);
		definition.set_fixing_lag(0);
		definition.clear_fixing_calendars();
		definition.add_fixing_calendars(GBLO);
		definition.clear_value_date_calendars();
		definition.add_value_date_calendars(GBLO);
		definition.clear_index_calendars();
		definition.add_index_calendars(GBLO);
		definition.set_day_count_fraction(DayCountFraction::ACT_365_FIXED);
		definition.set_eom(true);
		definition.set_default_for_index_family(true);
		default_index_service->register_index(definition);

		definition.set_isda_index(IsdaIndex::GBP_WMBA_SONIA_COMPOUND);
		definition.set_index_family(IndexFamily::SONIA);
		definition.set_tenor(TENOR_1D);
		definition.set_eom(false);
		definition.set_default_for_index_family(true);
		default_index_service->register_index(definition);

		definition.set_tenor(TENOR_UNSPECIFIED);
		default_index_service->register_index(definition);

		definition.set_isda_index(IsdaIndex::EUR_LIBOR_BBA);
		definition.set_index_family(IndexFamily::LIBOR);
		definition.set_currency(Currency::EUR);
		definition.set_tenor(TENOR_1D);
		definition.set_fixing_lag(0);
		definition.clear_fixing_calendars();
		definition.add_fixing_calendars(EUTA);
		definition.clear_value_date_calendars();
		definition.add_value_date_calendars(EUTA);
		definition.add_value_date_calendars(GBLO);
		definition.set_value_date_calendars_join_rule(JOIN_BUSINESS_DAYS);
		definition.clear_index_calendars();
		definition.add_index_calendars(EUTA);
		definition.add_index_calendars(GBLO);
		definition.set_index_calendars_join_rule(JOIN_BUSINESS_DAYS);
		definition.set_day_count_fraction(DayCountFraction::ACT_360);
		definition.set_eom(false);
		definition.set_default_for_index_family(true);
		default_index_service->register_index(definition);

		definition.set_tenor(TENOR_UNSPECIFIED);
		definition.set_fixing_lag(2);
		definition.set_eom(true);
		default_index_service->register_index(definition);

		definition.set_isda_index(IsdaIndex::EUR_EURIBOR_Reuters);
		definition.set_index_family(IndexFamily::EURIBOR);
		definition.set_currency(Currency::EUR);
		definition.set_tenor(TENOR_UNSPECIFIED);
		definition.set_fixing_lag(2);
		definition.clear_fixing_calendars();
		definition.add_fixing_calendars(EUTA);
		definition.clear_value_date_calendars();
		definition.add_value_date_calendars(EUTA);
		definition.set_value_date_calendars_join_rule(JOIN_HOLIDAYS);
		definition.clear_index_calendars();
		definition.add_index_calendars(EUTA);
		definition.set_index_calendars_join_rule(JOIN_HOLIDAYS);
		definition.set_day_count_fraction(DayCountFraction::ACT_360);
		definition.set_eom(true);
		definition.set_default_for_index_family(true);
		default_index_service->register_index(definition);

		definition.set_isda_index(IsdaIndex::EUR_EONIA_OIS_COMPOUND);
		definition.set_index_family(IndexFamily::EONIA);
		definition.set_tenor(TENOR_1D);
		definition.set_fixing_lag(0);
		definition.set_eom(false);
		definition.set_default_for_index_family(true);
		default_index_service->register_index(definition);

		definition.set_tenor(TENOR_UNSPECIFIED);
		default_index_service->register_index(definition);
	}
	return default_index_service.get();
}

////////////////////////////////////////////// TESTS
#if 0
struct IndexData {
	Currency currency;
	const char *fixing_date;
	const char *on_value_date;
	const char *other_value_date;
};

static IndexData data[] = {
    {Currency::USD, "2009-01-01", "No Fix", "No Fix"},
    {Currency::USD, "2009-01-02", "2009-01-02", "2009-01-06"},
    {Currency::USD, "2009-01-05", "2009-01-05", "2009-01-07"},
    {Currency::USD, "2009-01-06", "2009-01-06", "2009-01-08"},
    {Currency::USD, "2009-01-07", "2009-01-07", "2009-01-09"},
    {Currency::USD, "2009-01-13", "2009-01-13", "2009-01-15"},
    {Currency::USD, "2009-01-14", "2009-01-14", "2009-01-16"},
    {Currency::USD, "2009-01-15", "2009-01-15", "2009-01-20"},
    {Currency::USD, "2009-01-16", "2009-01-16", "2009-01-20"},
    {Currency::USD, "2009-01-19", "No Fix", "2009-01-21"},
    {Currency::USD, "2009-01-20", "2009-01-20", "2009-01-22"},
    {Currency::USD, "2009-01-21", "2009-01-21", "2009-01-23"},
    {Currency::USD, "2009-01-22", "2009-01-22", "2009-01-26"},
    {Currency::USD, "2009-01-23", "2009-01-23", "2009-01-27"},
    {Currency::USD, "2009-02-09", "2009-02-09", "2009-02-11"},
    {Currency::USD, "2009-02-10", "2009-02-10", "2009-02-12"},
    {Currency::USD, "2009-02-11", "2009-02-11", "2009-02-13"},
    {Currency::USD, "2009-02-12", "2009-02-12", "2009-02-17"},
    {Currency::USD, "2009-02-13", "2009-02-13", "2009-02-17"},
    {Currency::USD, "2009-02-16", "No Fix", "2009-02-18"},
    {Currency::USD, "2009-02-17", "2009-02-17", "2009-02-19"},
    {Currency::USD, "2009-02-18", "2009-02-18", "2009-02-20"},
    {Currency::USD, "2009-02-19", "2009-02-19", "2009-02-23"},
    {Currency::USD, "2009-02-20", "2009-02-20", "2009-02-24"},
    {Currency::USD, "2009-02-23", "2009-02-23", "2009-02-25"},
    {Currency::USD, "2009-04-03", "2009-04-03", "2009-04-07"},
    {Currency::USD, "2009-04-06", "2009-04-06", "2009-04-08"},
    {Currency::USD, "2009-04-07", "2009-04-07", "2009-04-09"},
    {Currency::USD, "2009-04-08", "2009-04-08", "2009-04-14"},
    {Currency::USD, "2009-04-09", "2009-04-09", "2009-04-15"},
    {Currency::USD, "2009-04-10", "No Fix", "No Fix"},
    {Currency::USD, "2009-04-13", "No Fix", "No Fix"},
    {Currency::USD, "2009-04-14", "2009-04-14", "2009-04-16"},
    {Currency::USD, "2009-04-15", "2009-04-15", "2009-04-17"},
    {Currency::USD, "2009-04-16", "2009-04-16", "2009-04-20"},
    {Currency::USD, "2009-04-28", "2009-04-28", "2009-04-30"},
    {Currency::USD, "2009-04-29", "2009-04-29", "2009-05-01"},
    {Currency::USD, "2009-04-30", "2009-04-30", "2009-05-05"},
    {Currency::USD, "2009-05-01", "2009-05-01", "2009-05-06"},
    {Currency::USD, "2009-05-04", "No Fix", "No Fix"},
    {Currency::USD, "2009-05-05", "2009-05-05", "2009-05-07"},
    {Currency::USD, "2009-05-06", "2009-05-06", "2009-05-08"},
    {Currency::USD, "2009-05-07", "2009-05-07", "2009-05-11"},
    {Currency::USD, "2009-05-08", "2009-05-08", "2009-05-12"},
    {Currency::USD, "2009-05-19", "2009-05-19", "2009-05-21"},
    {Currency::USD, "2009-05-20", "2009-05-20", "2009-05-22"},
    {Currency::USD, "2009-05-21", "2009-05-21", "2009-05-26"},
    {Currency::USD, "2009-05-22", "2009-05-22", "2009-05-27"},
    {Currency::USD, "2009-05-25", "No Fix", "No Fix"},
    {Currency::USD, "2009-05-26", "2009-05-26", "2009-05-28"},
    {Currency::USD, "2009-05-27", "2009-05-27", "2009-05-29"},
    {Currency::USD, "2009-05-28", "2009-05-28", "2009-06-01"},
    {Currency::USD, "2009-05-29", "2009-05-29", "2009-06-02"},
    {Currency::USD, "2009-06-01", "2009-06-01", "2009-06-03"},
    {Currency::USD, "2009-08-24", "2009-08-24", "2009-08-26"},
    {Currency::USD, "2009-08-25", "2009-08-25", "2009-08-27"},
    {Currency::USD, "2009-08-26", "2009-08-26", "2009-08-28"},
    {Currency::USD, "2009-08-27", "2009-08-27", "2009-09-01"},
    {Currency::USD, "2009-08-28", "2009-08-28", "2009-09-02"},
    {Currency::USD, "2009-08-31", "No Fix", "No Fix"},
    {Currency::USD, "2009-09-01", "2009-09-01", "2009-09-03"},
    {Currency::USD, "2009-09-02", "2009-09-02", "2009-09-04"},
    {Currency::USD, "2009-09-03", "2009-09-03", "2009-09-08"},
    {Currency::USD, "2009-09-04", "2009-09-04", "2009-09-08"},
    {Currency::USD, "2009-09-07", "No Fix", "2009-09-09"},
    {Currency::USD, "2009-09-08", "2009-09-08", "2009-09-10"},
    {Currency::USD, "2009-09-09", "2009-09-09", "2009-09-11"},
    {Currency::USD, "2009-09-10", "2009-09-10", "2009-09-14"},
    {Currency::USD, "2009-09-11", "2009-09-11", "2009-09-15"},
    {Currency::USD, "2009-09-14", "2009-09-14", "2009-09-16"},
    {Currency::USD, "2009-10-05", "2009-10-05", "2009-10-07"},
    {Currency::USD, "2009-10-06", "2009-10-06", "2009-10-08"},
    {Currency::USD, "2009-10-07", "2009-10-07", "2009-10-09"},
    {Currency::USD, "2009-10-08", "2009-10-08", "2009-10-13"},
    {Currency::USD, "2009-10-09", "2009-10-09", "2009-10-13"},
    {Currency::USD, "2009-10-12", "No Fix", "2009-10-14"},
    {Currency::USD, "2009-10-13", "2009-10-13", "2009-10-15"},
    {Currency::USD, "2009-10-14", "2009-10-14", "2009-10-16"},
    {Currency::USD, "2009-11-03", "2009-11-03", "2009-11-05"},
    {Currency::USD, "2009-11-04", "2009-11-04", "2009-11-06"},
    {Currency::USD, "2009-11-05", "2009-11-05", "2009-11-09"},
    {Currency::USD, "2009-11-06", "2009-11-06", "2009-11-10"},
    {Currency::USD, "2009-11-09", "2009-11-09", "2009-11-12"},
    {Currency::USD, "2009-11-10", "2009-11-10", "2009-11-12"},
    {Currency::USD, "2009-11-11", "No Fix", "2009-11-13"},
    {Currency::USD, "2009-11-12", "2009-11-12", "2009-11-16"},
    {Currency::USD, "2009-11-13", "2009-11-13", "2009-11-17"},
    {Currency::USD, "2009-11-16", "2009-11-16", "2009-11-18"},
    {Currency::USD, "2009-11-17", "2009-11-17", "2009-11-19"},
    {Currency::USD, "2009-11-20", "2009-11-20", "2009-11-24"},
    {Currency::USD, "2009-11-23", "2009-11-23", "2009-11-25"},
    {Currency::USD, "2009-11-24", "2009-11-24", "2009-11-27"},
    {Currency::USD, "2009-11-25", "2009-11-25", "2009-11-27"},
    {Currency::USD, "2009-11-26", "No Fix", "2009-11-30"},
    {Currency::USD, "2009-11-27", "2009-11-27", "2009-12-01"},
    {Currency::USD, "2009-11-30", "2009-11-30", "2009-12-02"},
    {Currency::USD, "2009-12-01", "2009-12-01", "2009-12-03"},
    {Currency::USD, "2009-12-02", "2009-12-02", "2009-12-04"},
    {Currency::USD, "2009-12-03", "2009-12-03", "2009-12-07"},
    {Currency::USD, "2009-12-18", "2009-12-18", "2009-12-22"},
    {Currency::USD, "2009-12-21", "2009-12-21", "2009-12-23"},
    {Currency::USD, "2009-12-22", "2009-12-22", "2009-12-24"},
    {Currency::USD, "2009-12-23", "2009-12-23", "2009-12-29"},
    {Currency::USD, "2009-12-24", "2009-12-24", "2009-12-30"},
    {Currency::USD, "2009-12-25", "No Fix", "No Fix"},
    {Currency::USD, "2009-12-28", "No Fix", "No Fix"},
    {Currency::USD, "2009-12-29", "2009-12-29", "2009-12-31"},
    {Currency::USD, "2009-12-30", "2009-12-30", "2010-01-04"},
    {Currency::USD, "2009-12-31", "2009-12-31", "2010-01-05"},
    {Currency::EUR, "2009-01-01", "No Fix", "No Fix"},
    {Currency::EUR, "2009-01-02", "2009-01-02", "2009-01-06"},
    {Currency::EUR, "2009-01-05", "2009-01-05", "2009-01-07"},
    {Currency::EUR, "2009-01-06", "2009-01-06", "2009-01-08"},
    {Currency::EUR, "2009-01-07", "2009-01-07", "2009-01-09"},
    {Currency::EUR, "2009-04-03", "2009-04-03", "2009-04-07"},
    {Currency::EUR, "2009-04-06", "2009-04-06", "2009-04-08"},
    {Currency::EUR, "2009-04-07", "2009-04-07", "2009-04-09"},
    {Currency::EUR, "2009-04-08", "2009-04-08", "2009-04-14"},
    {Currency::EUR, "2009-04-09", "2009-04-09", "2009-04-15"},
    {Currency::EUR, "2009-04-10", "No Fix", "No Fix"},
    {Currency::EUR, "2009-04-13", "No Fix", "No Fix"},
    {Currency::EUR, "2009-04-14", "2009-04-14", "2009-04-16"},
    {Currency::EUR, "2009-04-15", "2009-04-15", "2009-04-17"},
    {Currency::EUR, "2009-04-16", "2009-04-16", "2009-04-20"},
    {Currency::EUR, "2009-04-28", "2009-04-28", "2009-04-30"},
    {Currency::EUR, "2009-04-29", "2009-04-29", "2009-05-04"},
    {Currency::EUR, "2009-04-30", "2009-04-30", "2009-05-05"},
    {Currency::EUR, "2009-05-01", "No Fix", "2009-05-05"},
    {Currency::EUR, "2009-05-04", "2009-05-04", "2009-05-06"},
    {Currency::EUR, "2009-05-05", "2009-05-05", "2009-05-07"},
    {Currency::EUR, "2009-05-06", "2009-05-06", "2009-05-08"},
    {Currency::EUR, "2009-12-21", "2009-12-21", "2009-12-23"},
    {Currency::EUR, "2009-12-22", "2009-12-22", "2009-12-24"},
    {Currency::EUR, "2009-12-23", "2009-12-23", "2009-12-28"},
    {Currency::EUR, "2009-12-24", "2009-12-24", "2009-12-29"},
    {Currency::EUR, "2009-12-25", "No Fix", "No Fix"},
    {Currency::EUR, "2009-12-28", "2009-12-28", "2009-12-30"},
    {Currency::EUR, "2009-12-29", "2009-12-29", "2009-12-31"},
    {Currency::EUR, "2009-12-30", "2009-12-30", "2010-01-04"},
    {Currency::EUR, "2009-12-31", "2009-12-31", "2010-01-05"},
    {Currency::GBP, "2009-01-01", "No Fix", "No Fix"},
    {Currency::GBP, "2009-01-02", "2009-01-02", "2009-01-02"},
    {Currency::GBP, "2009-01-05", "2009-01-05", "2009-01-05"},
    {Currency::GBP, "2009-01-06", "2009-01-06", "2009-01-06"},
    {Currency::GBP, "2009-04-06", "2009-04-06", "2009-04-06"},
    {Currency::GBP, "2009-04-07", "2009-04-07", "2009-04-07"},
    {Currency::GBP, "2009-04-08", "2009-04-08", "2009-04-08"},
    {Currency::GBP, "2009-04-09", "2009-04-09", "2009-04-09"},
    {Currency::GBP, "2009-04-10", "No Fix", "No Fix"},
    {Currency::GBP, "2009-04-13", "No Fix", "No Fix"},
    {Currency::GBP, "2009-04-14", "2009-04-14", "2009-04-14"},
    {Currency::GBP, "2009-04-15", "2009-04-15", "2009-04-15"},
    {Currency::GBP, "2009-04-16", "2009-04-16", "2009-04-16"},
    {Currency::GBP, "2009-04-29", "2009-04-29", "2009-04-29"},
    {Currency::GBP, "2009-04-30", "2009-04-30", "2009-04-30"},
    {Currency::GBP, "2009-05-01", "2009-05-01", "2009-05-01"},
    {Currency::GBP, "2009-05-04", "No Fix", "No Fix"},
    {Currency::GBP, "2009-05-05", "2009-05-05", "2009-05-05"},
    {Currency::GBP, "2009-05-06", "2009-05-06", "2009-05-06"},
    {Currency::GBP, "2009-05-07", "2009-05-07", "2009-05-07"},
    {Currency::GBP, "2009-05-19", "2009-05-19", "2009-05-19"},
    {Currency::GBP, "2009-05-20", "2009-05-20", "2009-05-20"},
    {Currency::GBP, "2009-05-21", "2009-05-21", "2009-05-21"},
    {Currency::GBP, "2009-05-22", "2009-05-22", "2009-05-22"},
    {Currency::GBP, "2009-05-25", "No Fix", "No Fix"},
    {Currency::GBP, "2009-05-26", "2009-05-26", "2009-05-26"},
    {Currency::GBP, "2009-05-27", "2009-05-27", "2009-05-27"},
    {Currency::GBP, "2009-05-28", "2009-05-28", "2009-05-28"},
    {Currency::GBP, "2009-05-29", "2009-05-29", "2009-05-29"},
    {Currency::GBP, "2009-08-25", "2009-08-25", "2009-08-25"},
    {Currency::GBP, "2009-08-26", "2009-08-26", "2009-08-26"},
    {Currency::GBP, "2009-08-27", "2009-08-27", "2009-08-27"},
    {Currency::GBP, "2009-08-28", "2009-08-28", "2009-08-28"},
    {Currency::GBP, "2009-08-31", "No Fix", "No Fix"},
    {Currency::GBP, "2009-09-01", "2009-09-01", "2009-09-01"},
    {Currency::GBP, "2009-09-02", "2009-09-02", "2009-09-02"},
    {Currency::GBP, "2009-09-03", "2009-09-03", "2009-09-03"},
    {Currency::GBP, "2009-12-21", "2009-12-21", "2009-12-21"},
    {Currency::GBP, "2009-12-22", "2009-12-22", "2009-12-22"},
    {Currency::GBP, "2009-12-23", "2009-12-23", "2009-12-23"},
    {Currency::GBP, "2009-12-24", "2009-12-24", "2009-12-24"},
    {Currency::GBP, "2009-12-25", "No Fix", "No Fix"},
    {Currency::GBP, "2009-12-28", "No Fix", "No Fix"},
    {Currency::GBP, "2009-12-29", "2009-12-29", "2009-12-29"},
    {Currency::GBP, "2009-12-30", "2009-12-30", "2009-12-30"},
    {Currency::GBP, "2009-12-31", "2009-12-31", "2009-12-31"},
    {Currency::USD, "2010-01-01", "No Fix", "No Fix"},
    {Currency::USD, "2010-01-04", "2010-01-04", "2010-01-06"},
    {Currency::USD, "2010-01-05", "2010-01-05", "2010-01-07"},
    {Currency::USD, "2010-01-06", "2010-01-06", "2010-01-08"},
    {Currency::USD, "2010-01-13", "2010-01-13", "2010-01-15"},
    {Currency::USD, "2010-01-14", "2010-01-14", "2010-01-19"},
    {Currency::USD, "2010-01-15", "2010-01-15", "2010-01-19"},
    {Currency::USD, "2010-01-18", "No Fix", "2010-01-20"},
    {Currency::USD, "2010-01-19", "2010-01-19", "2010-01-21"},
    {Currency::USD, "2010-01-20", "2010-01-20", "2010-01-22"},
    {Currency::USD, "2010-01-21", "2010-01-21", "2010-01-25"},
    {Currency::USD, "2010-02-10", "2010-02-10", "2010-02-12"},
    {Currency::USD, "2010-02-11", "2010-02-11", "2010-02-16"},
    {Currency::USD, "2010-02-12", "2010-02-12", "2010-02-16"},
    {Currency::USD, "2010-02-15", "No Fix", "2010-02-17"},
    {Currency::USD, "2010-02-16", "2010-02-16", "2010-02-18"},
    {Currency::USD, "2010-02-17", "2010-02-17", "2010-02-19"},
    {Currency::USD, "2010-02-18", "2010-02-18", "2010-02-22"},
    {Currency::USD, "2010-02-19", "2010-02-19", "2010-02-23"},
    {Currency::USD, "2010-03-29", "2010-03-29", "2010-03-31"},
    {Currency::USD, "2010-03-30", "2010-03-30", "2010-04-01"},
    {Currency::USD, "2010-03-31", "2010-03-31", "2010-04-06"},
    {Currency::USD, "2010-04-01", "2010-04-01", "2010-04-07"},
    {Currency::USD, "2010-04-02", "No Fix", "No Fix"},
    {Currency::USD, "2010-04-05", "No Fix", "No Fix"},
    {Currency::USD, "2010-04-06", "2010-04-06", "2010-04-08"},
    {Currency::USD, "2010-04-07", "2010-04-07", "2010-04-09"},
    {Currency::USD, "2010-04-08", "2010-04-08", "2010-04-12"},
    {Currency::USD, "2010-04-28", "2010-04-28", "2010-04-30"},
    {Currency::USD, "2010-04-29", "2010-04-29", "2010-05-04"},
    {Currency::USD, "2010-04-30", "2010-04-30", "2010-05-05"},
    {Currency::USD, "2010-05-03", "No Fix", "No Fix"},
    {Currency::USD, "2010-05-04", "2010-05-04", "2010-05-06"},
    {Currency::USD, "2010-05-05", "2010-05-05", "2010-05-07"},
    {Currency::USD, "2010-05-06", "2010-05-06", "2010-05-10"},
    {Currency::USD, "2010-05-26", "2010-05-26", "2010-05-28"},
    {Currency::USD, "2010-05-27", "2010-05-27", "2010-06-01"},
    {Currency::USD, "2010-05-28", "2010-05-28", "2010-06-02"},
    {Currency::USD, "2010-05-31", "No Fix", "No Fix"},
    {Currency::USD, "2010-06-01", "2010-06-01", "2010-06-03"},
    {Currency::USD, "2010-06-02", "2010-06-02", "2010-06-04"},
    {Currency::USD, "2010-06-03", "2010-06-03", "2010-06-07"},
    {Currency::USD, "2010-08-24", "2010-08-24", "2010-08-26"},
    {Currency::USD, "2010-08-25", "2010-08-25", "2010-08-27"},
    {Currency::USD, "2010-08-26", "2010-08-26", "2010-08-31"},
    {Currency::USD, "2010-08-27", "2010-08-27", "2010-09-01"},
    {Currency::USD, "2010-08-30", "No Fix", "No Fix"},
    {Currency::USD, "2010-08-31", "2010-08-31", "2010-09-02"},
    {Currency::USD, "2010-09-01", "2010-09-01", "2010-09-03"},
    {Currency::USD, "2010-09-02", "2010-09-02", "2010-09-07"},
    {Currency::USD, "2010-09-03", "2010-09-03", "2010-09-07"},
    {Currency::USD, "2010-09-06", "No Fix", "2010-09-08"},
    {Currency::USD, "2010-09-07", "2010-09-07", "2010-09-09"},
    {Currency::USD, "2010-09-08", "2010-09-08", "2010-09-10"},
    {Currency::USD, "2010-09-09", "2010-09-09", "2010-09-13"},
    {Currency::USD, "2010-10-06", "2010-10-06", "2010-10-08"},
    {Currency::USD, "2010-10-07", "2010-10-07", "2010-10-12"},
    {Currency::USD, "2010-10-08", "2010-10-08", "2010-10-12"},
    {Currency::USD, "2010-10-11", "No Fix", "2010-10-13"},
    {Currency::USD, "2010-10-12", "2010-10-12", "2010-10-14"},
    {Currency::USD, "2010-10-13", "2010-10-13", "2010-10-15"},
    {Currency::USD, "2010-10-14", "2010-10-14", "2010-10-18"},
    {Currency::USD, "2010-10-15", "2010-10-15", "2010-10-19"},
    {Currency::USD, "2010-11-08", "2010-11-08", "2010-11-10"},
    {Currency::USD, "2010-11-09", "2010-11-09", "2010-11-12"},
    {Currency::USD, "2010-11-10", "2010-11-10", "2010-11-12"},
    {Currency::USD, "2010-11-11", "No Fix", "2010-11-15"},
    {Currency::USD, "2010-11-12", "2010-11-12", "2010-11-16"},
    {Currency::USD, "2010-11-15", "2010-11-15", "2010-11-17"},
    {Currency::USD, "2010-11-16", "2010-11-16", "2010-11-18"},
    {Currency::USD, "2010-11-22", "2010-11-22", "2010-11-24"},
    {Currency::USD, "2010-11-23", "2010-11-23", "2010-11-26"},
    {Currency::USD, "2010-11-24", "2010-11-24", "2010-11-26"},
    {Currency::USD, "2010-11-25", "No Fix", "2010-11-29"},
    {Currency::USD, "2010-11-26", "2010-11-26", "2010-11-30"},
    {Currency::USD, "2010-11-29", "2010-11-29", "2010-12-01"},
    {Currency::USD, "2010-11-30", "2010-11-30", "2010-12-02"},
    {Currency::USD, "2010-12-01", "2010-12-01", "2010-12-03"},
    {Currency::USD, "2010-12-22", "2010-12-22", "2010-12-24"},
    {Currency::USD, "2010-12-23", "2010-12-23", "2010-12-29"},
    {Currency::USD, "2010-12-24", "2010-12-24", "2010-12-30"},
    {Currency::USD, "2010-12-27", "No Fix", "No Fix"},
    {Currency::USD, "2010-12-28", "No Fix", "No Fix"},
    {Currency::USD, "2010-12-29", "2010-12-29", "2010-12-31"},
    {Currency::USD, "2010-12-30", "2010-12-30", "2011-01-04"},
    {Currency::USD, "2010-12-31", "2010-12-31", "2011-01-05"},
    {Currency::EUR, "2010-01-01", "No Fix", "No Fix"},
    {Currency::EUR, "2010-01-04", "2010-01-04", "2010-01-06"},
    {Currency::EUR, "2010-01-05", "2010-01-05", "2010-01-07"},
    {Currency::EUR, "2010-01-06", "2010-01-06", "2010-01-08"},
    {Currency::EUR, "2010-03-29", "2010-03-29", "2010-03-31"},
    {Currency::EUR, "2010-03-30", "2010-03-30", "2010-04-01"},
    {Currency::EUR, "2010-03-31", "2010-03-31", "2010-04-06"},
    {Currency::EUR, "2010-04-01", "2010-04-01", "2010-04-07"},
    {Currency::EUR, "2010-04-02", "No Fix", "No Fix"},
    {Currency::EUR, "2010-04-05", "No Fix", "No Fix"},
    {Currency::EUR, "2010-04-06", "2010-04-06", "2010-04-08"},
    {Currency::EUR, "2010-04-07", "2010-04-07", "2010-04-09"},
    {Currency::EUR, "2010-04-08", "2010-04-08", "2010-04-12"},
    {Currency::EUR, "2010-12-27", "2010-12-27", "2010-12-29"},
    {Currency::EUR, "2010-12-28", "2010-12-28", "2010-12-30"},
    {Currency::EUR, "2010-12-29", "2010-12-29", "2010-12-31"},
    {Currency::EUR, "2010-12-30", "2010-12-30", "2011-01-03"},
    {Currency::EUR, "2010-12-31", "2010-12-31", "2011-01-04"},
    {Currency::GBP, "2010-01-01", "No Fix", "No Fix"},
    {Currency::GBP, "2010-01-04", "2010-01-04", "2010-01-04"},
    {Currency::GBP, "2010-01-05", "2010-01-05", "2010-01-05"},
    {Currency::GBP, "2010-01-06", "2010-01-06", "2010-01-06"},
    {Currency::GBP, "2010-03-29", "2010-03-29", "2010-03-29"},
    {Currency::GBP, "2010-03-30", "2010-03-30", "2010-03-30"},
    {Currency::GBP, "2010-03-31", "2010-03-31", "2010-03-31"},
    {Currency::GBP, "2010-04-01", "2010-04-01", "2010-04-01"},
    {Currency::GBP, "2010-04-02", "No Fix", "No Fix"},
    {Currency::GBP, "2010-04-05", "No Fix", "No Fix"},
    {Currency::GBP, "2010-04-06", "2010-04-06", "2010-04-06"},
    {Currency::GBP, "2010-04-07", "2010-04-07", "2010-04-07"},
    {Currency::GBP, "2010-04-08", "2010-04-08", "2010-04-08"},
    {Currency::GBP, "2010-04-28", "2010-04-28", "2010-04-28"},
    {Currency::GBP, "2010-04-29", "2010-04-29", "2010-04-29"},
    {Currency::GBP, "2010-04-30", "2010-04-30", "2010-04-30"},
    {Currency::GBP, "2010-05-03", "No Fix", "No Fix"},
    {Currency::GBP, "2010-05-04", "2010-05-04", "2010-05-04"},
    {Currency::GBP, "2010-05-05", "2010-05-05", "2010-05-05"},
    {Currency::GBP, "2010-05-06", "2010-05-06", "2010-05-06"},
    {Currency::GBP, "2010-05-26", "2010-05-26", "2010-05-26"},
    {Currency::GBP, "2010-05-27", "2010-05-27", "2010-05-27"},
    {Currency::GBP, "2010-05-28", "2010-05-28", "2010-05-28"},
    {Currency::GBP, "2010-05-31", "No Fix", "No Fix"},
    {Currency::GBP, "2010-06-01", "2010-06-01", "2010-06-01"},
    {Currency::GBP, "2010-06-02", "2010-06-02", "2010-06-02"},
    {Currency::GBP, "2010-06-03", "2010-06-03", "2010-06-03"},
    {Currency::GBP, "2010-08-25", "2010-08-25", "2010-08-25"},
    {Currency::GBP, "2010-08-26", "2010-08-26", "2010-08-26"},
    {Currency::GBP, "2010-08-27", "2010-08-27", "2010-08-27"},
    {Currency::GBP, "2010-08-30", "No Fix", "No Fix"},
    {Currency::GBP, "2010-08-31", "2010-08-31", "2010-08-31"},
    {Currency::GBP, "2010-09-01", "2010-09-01", "2010-09-01"},
    {Currency::GBP, "2010-09-02", "2010-09-02", "2010-09-02"},
    {Currency::GBP, "2010-12-21", "2010-12-21", "2010-12-21"},
    {Currency::GBP, "2010-12-22", "2010-12-22", "2010-12-22"},
    {Currency::GBP, "2010-12-23", "2010-12-23", "2010-12-23"},
    {Currency::GBP, "2010-12-24", "2010-12-24", "2010-12-24"},
    {Currency::GBP, "2010-12-27", "No Fix", "No Fix"},
    {Currency::GBP, "2010-12-28", "No Fix", "No Fix"},
    {Currency::GBP, "2010-12-29", "2010-12-29", "2010-12-29"},
    {Currency::GBP, "2010-12-30", "2010-12-30", "2010-12-30"},
    {Currency::GBP, "2010-12-31", "2010-12-31", "2010-12-31"},
    {Currency::USD, "2011-01-03", "No Fix", "No Fix"},
    {Currency::USD, "2011-01-04", "2011-01-04", "2011-01-06"},
    {Currency::USD, "2011-01-05", "2011-01-05", "2011-01-07"},
    {Currency::USD, "2011-01-12", "2011-01-12", "2011-01-14"},
    {Currency::USD, "2011-01-13", "2011-01-13", "2011-01-18"},
    {Currency::USD, "2011-01-14", "2011-01-14", "2011-01-18"},
    {Currency::USD, "2011-01-17", "No Fix", "2011-01-19"},
    {Currency::USD, "2011-01-18", "2011-01-18", "2011-01-20"},
    {Currency::USD, "2011-01-19", "2011-01-19", "2011-01-21"},
    {Currency::USD, "2011-01-20", "2011-01-20", "2011-01-24"},
    {Currency::USD, "2011-02-15", "2011-02-15", "2011-02-17"},
    {Currency::USD, "2011-02-16", "2011-02-16", "2011-02-18"},
    {Currency::USD, "2011-02-17", "2011-02-17", "2011-02-22"},
    {Currency::USD, "2011-02-18", "2011-02-18", "2011-02-22"},
    {Currency::USD, "2011-02-21", "No Fix", "2011-02-23"},
    {Currency::USD, "2011-02-22", "2011-02-22", "2011-02-24"},
    {Currency::USD, "2011-02-23", "2011-02-23", "2011-02-25"},
    {Currency::USD, "2011-02-24", "2011-02-24", "2011-02-28"},
    {Currency::USD, "2011-04-19", "2011-04-19", "2011-04-21"},
    {Currency::USD, "2011-04-20", "2011-04-20", "2011-04-26"},
    {Currency::USD, "2011-04-21", "2011-04-21", "2011-04-27"},
    {Currency::USD, "2011-04-22", "No Fix", "No Fix"},
    {Currency::USD, "2011-04-25", "No Fix", "No Fix"},
    {Currency::USD, "2011-04-26", "2011-04-26", "2011-04-28"},
    {Currency::USD, "2011-04-27", "2011-04-27", "2011-05-03"},
    {Currency::USD, "2011-04-28", "2011-04-28", "2011-05-04"},
    {Currency::USD, "2011-04-29", "No Fix", "No Fix"},
    {Currency::USD, "2011-05-02", "No Fix", "No Fix"},
    {Currency::USD, "2011-05-03", "2011-05-03", "2011-05-05"},
    {Currency::USD, "2011-05-04", "2011-05-04", "2011-05-06"},
    {Currency::USD, "2011-05-05", "2011-05-05", "2011-05-09"},
    {Currency::USD, "2011-05-25", "2011-05-25", "2011-05-27"},
    {Currency::USD, "2011-05-26", "2011-05-26", "2011-05-31"},
    {Currency::USD, "2011-05-27", "2011-05-27", "2011-06-01"},
    {Currency::USD, "2011-05-30", "No Fix", "No Fix"},
    {Currency::USD, "2011-05-31", "2011-05-31", "2011-06-02"},
    {Currency::USD, "2011-06-01", "2011-06-01", "2011-06-03"},
    {Currency::USD, "2011-06-02", "2011-06-02", "2011-06-06"},
    {Currency::USD, "2011-06-29", "2011-06-29", "2011-07-01"},
    {Currency::USD, "2011-06-30", "2011-06-30", "2011-07-05"},
    {Currency::USD, "2011-07-01", "2011-07-01", "2011-07-05"},
    {Currency::USD, "2011-07-04", "No Fix", "2011-07-06"},
    {Currency::USD, "2011-07-05", "2011-07-05", "2011-07-07"},
    {Currency::USD, "2011-07-06", "2011-07-06", "2011-07-08"},
    {Currency::USD, "2011-07-07", "2011-07-07", "2011-07-11"},
    {Currency::USD, "2011-08-23", "2011-08-23", "2011-08-25"},
    {Currency::USD, "2011-08-24", "2011-08-24", "2011-08-26"},
    {Currency::USD, "2011-08-25", "2011-08-25", "2011-08-30"},
    {Currency::USD, "2011-08-26", "2011-08-26", "2011-08-31"},
    {Currency::USD, "2011-08-29", "No Fix", "No Fix"},
    {Currency::USD, "2011-08-30", "2011-08-30", "2011-09-01"},
    {Currency::USD, "2011-08-31", "2011-08-31", "2011-09-02"},
    {Currency::USD, "2011-09-01", "2011-09-01", "2011-09-06"},
    {Currency::USD, "2011-09-02", "2011-09-02", "2011-09-06"},
    {Currency::USD, "2011-09-05", "No Fix", "2011-09-07"},
    {Currency::USD, "2011-09-06", "2011-09-06", "2011-09-08"},
    {Currency::USD, "2011-09-07", "2011-09-07", "2011-09-09"},
    {Currency::USD, "2011-09-08", "2011-09-08", "2011-09-12"},
    {Currency::USD, "2011-10-05", "2011-10-05", "2011-10-07"},
    {Currency::USD, "2011-10-06", "2011-10-06", "2011-10-11"},
    {Currency::USD, "2011-10-07", "2011-10-07", "2011-10-11"},
    {Currency::USD, "2011-10-10", "No Fix", "2011-10-12"},
    {Currency::USD, "2011-10-11", "2011-10-11", "2011-10-13"},
    {Currency::USD, "2011-10-12", "2011-10-12", "2011-10-14"},
    {Currency::USD, "2011-10-13", "2011-10-13", "2011-10-17"},
    {Currency::USD, "2011-11-08", "2011-11-08", "2011-11-10"},
    {Currency::USD, "2011-11-09", "2011-11-09", "2011-11-14"},
    {Currency::USD, "2011-11-10", "2011-11-10", "2011-11-14"},
    {Currency::USD, "2011-11-11", "No Fix", "2011-11-15"},
    {Currency::USD, "2011-11-14", "2011-11-14", "2011-11-16"},
    {Currency::USD, "2011-11-15", "2011-11-15", "2011-11-17"},
    {Currency::USD, "2011-11-16", "2011-11-16", "2011-11-18"},
    {Currency::USD, "2011-11-21", "2011-11-21", "2011-11-23"},
    {Currency::USD, "2011-11-22", "2011-11-22", "2011-11-25"},
    {Currency::USD, "2011-11-23", "2011-11-23", "2011-11-25"},
    {Currency::USD, "2011-11-24", "No Fix", "2011-11-28"},
    {Currency::USD, "2011-11-25", "2011-11-25", "2011-11-29"},
    {Currency::USD, "2011-11-28", "2011-11-28", "2011-11-30"},
    {Currency::USD, "2011-11-29", "2011-11-29", "2011-12-01"},
    {Currency::USD, "2011-12-21", "2011-12-21", "2011-12-23"},
    {Currency::USD, "2011-12-22", "2011-12-22", "2011-12-28"},
    {Currency::USD, "2011-12-23", "2011-12-23", "2011-12-29"},
    {Currency::USD, "2011-12-26", "No Fix", "No Fix"},
    {Currency::USD, "2011-12-27", "No Fix", "No Fix"},
    {Currency::USD, "2011-12-28", "2011-12-28", "2011-12-30"},
    {Currency::USD, "2011-12-29", "2011-12-29", "2012-01-03"},
    {Currency::USD, "2011-12-30", "2011-12-30", "2012-01-04"},
    {Currency::USD, "2012-01-02", "No Fix", "No Fix"},
    {Currency::EUR, "2011-01-03", "2011-01-03", "2011-01-05"},
    {Currency::EUR, "2011-01-04", "2011-01-04", "2011-01-06"},
    {Currency::EUR, "2011-01-05", "2011-01-05", "2011-01-07"},
    {Currency::EUR, "2011-04-18", "2011-04-18", "2011-04-20"},
    {Currency::EUR, "2011-04-19", "2011-04-19", "2011-04-21"},
    {Currency::EUR, "2011-04-20", "2011-04-20", "2011-04-26"},
    {Currency::EUR, "2011-04-21", "2011-04-21", "2011-04-27"},
    {Currency::EUR, "2011-04-22", "No Fix", "No Fix"},
    {Currency::EUR, "2011-04-25", "No Fix", "No Fix"},
    {Currency::EUR, "2011-04-26", "2011-04-26", "2011-04-28"},
    {Currency::EUR, "2011-04-27", "2011-04-27", "2011-04-29"},
    {Currency::EUR, "2011-04-28", "2011-04-28", "2011-05-02"},
    {Currency::EUR, "2011-12-20", "2011-12-20", "2011-12-22"},
    {Currency::EUR, "2011-12-21", "2011-12-21", "2011-12-23"},
    {Currency::EUR, "2011-12-22", "2011-12-22", "2011-12-28"},
    {Currency::EUR, "2011-12-23", "2011-12-23", "2011-12-29"},
    {Currency::EUR, "2011-12-26", "No Fix", "No Fix"},
    {Currency::EUR, "2011-12-27", "2011-12-27", "2011-12-29"},
    {Currency::EUR, "2011-12-28", "2011-12-28", "2011-12-30"},
    {Currency::EUR, "2011-12-29", "2011-12-29", "2012-01-03"},
    {Currency::EUR, "2011-12-30", "2011-12-30", "2012-01-04"},
    {Currency::EUR, "2012-01-02", "No Fix", "No Fix"},
    {Currency::GBP, "2011-01-03", "No Fix", "No Fix"},
    {Currency::GBP, "2011-01-04", "2011-01-04", "2011-01-04"},
    {Currency::GBP, "2011-01-05", "2011-01-05", "2011-01-05"},
    {Currency::GBP, "2011-01-06", "2011-01-06", "2011-01-06"},
    {Currency::GBP, "2011-01-07", "2011-01-07", "2011-01-07"},
    {Currency::GBP, "2011-04-18", "2011-04-18", "2011-04-18"},
    {Currency::GBP, "2011-04-19", "2011-04-19", "2011-04-19"},
    {Currency::GBP, "2011-04-20", "2011-04-20", "2011-04-20"},
    {Currency::GBP, "2011-04-21", "2011-04-21", "2011-04-21"},
    {Currency::GBP, "2011-04-22", "No Fix", "No Fix"},
    {Currency::GBP, "2011-04-25", "No Fix", "No Fix"},
    {Currency::GBP, "2011-04-26", "2011-04-26", "2011-04-26"},
    {Currency::GBP, "2011-04-27", "2011-04-27", "2011-04-27"},
    {Currency::GBP, "2011-04-28", "2011-04-28", "2011-04-28"},
    {Currency::GBP, "2011-04-29", "No Fix", "No Fix"},
    {Currency::GBP, "2011-05-02", "No Fix", "No Fix"},
    {Currency::GBP, "2011-05-03", "2011-05-03", "2011-05-03"},
    {Currency::GBP, "2011-05-04", "2011-05-04", "2011-05-04"},
    {Currency::GBP, "2011-05-05", "2011-05-05", "2011-05-05"},
    {Currency::GBP, "2011-05-24", "2011-05-24", "2011-05-24"},
    {Currency::GBP, "2011-05-25", "2011-05-25", "2011-05-25"},
    {Currency::GBP, "2011-05-26", "2011-05-26", "2011-05-26"},
    {Currency::GBP, "2011-05-27", "2011-05-27", "2011-05-27"},
    {Currency::GBP, "2011-05-30", "No Fix", "No Fix"},
    {Currency::GBP, "2011-05-31", "2011-05-31", "2011-05-31"},
    {Currency::GBP, "2011-06-01", "2011-06-01", "2011-06-01"},
    {Currency::GBP, "2011-06-02", "2011-06-02", "2011-06-02"},
    {Currency::GBP, "2011-06-03", "2011-06-03", "2011-06-03"},
    {Currency::GBP, "2011-08-23", "2011-08-23", "2011-08-23"},
    {Currency::GBP, "2011-08-24", "2011-08-24", "2011-08-24"},
    {Currency::GBP, "2011-08-25", "2011-08-25", "2011-08-25"},
    {Currency::GBP, "2011-08-26", "2011-08-26", "2011-08-26"},
    {Currency::GBP, "2011-08-29", "No Fix", "No Fix"},
    {Currency::GBP, "2011-08-30", "2011-08-30", "2011-08-30"},
    {Currency::GBP, "2011-08-31", "2011-08-31", "2011-08-31"},
    {Currency::GBP, "2011-09-01", "2011-09-01", "2011-09-01"},
    {Currency::GBP, "2011-09-02", "2011-09-02", "2011-09-02"},
    {Currency::GBP, "2011-12-20", "2011-12-20", "2011-12-20"},
    {Currency::GBP, "2011-12-21", "2011-12-21", "2011-12-21"},
    {Currency::GBP, "2011-12-22", "2011-12-22", "2011-12-22"},
    {Currency::GBP, "2011-12-23", "2011-12-23", "2011-12-23"},
    {Currency::GBP, "2011-12-26", "No Fix", "No Fix"},
    {Currency::GBP, "2011-12-27", "No Fix", "No Fix"},
    {Currency::GBP, "2011-12-28", "2011-12-28", "2011-12-28"},
    {Currency::GBP, "2011-12-29", "2011-12-29", "2011-12-29"},
    {Currency::GBP, "2011-12-30", "2011-12-30", "2011-12-30"},
    {Currency::GBP, "2012-01-02", "No Fix", "No Fix"},
    {Currency::USD, "2012-01-02", "No Fix", "No Fix"},
    {Currency::USD, "2012-01-03", "2012-01-03", "2012-01-05"},
    {Currency::USD, "2012-01-04", "2012-01-04", "2012-01-06"},
    {Currency::USD, "2012-01-05", "2012-01-05", "2012-01-09"},
    {Currency::USD, "2012-01-11", "2012-01-11", "2012-01-13"},
    {Currency::USD, "2012-01-12", "2012-01-12", "2012-01-17"},
    {Currency::USD, "2012-01-13", "2012-01-13", "2012-01-17"},
    {Currency::USD, "2012-01-16", "NO FIX", "2012-01-18"},
    {Currency::USD, "2012-01-17", "2012-01-17", "2012-01-19"},
    {Currency::USD, "2012-01-18", "2012-01-18", "2012-01-20"},
    {Currency::USD, "2012-01-19", "2012-01-19", "2012-01-23"},
    {Currency::USD, "2012-01-20", "2012-01-20", "2012-01-24"},
    {Currency::USD, "2012-02-15", "2012-02-15", "2012-02-17"},
    {Currency::USD, "2012-02-16", "2012-02-16", "2012-02-21"},
    {Currency::USD, "2012-02-17", "2012-02-17", "2012-02-21"},
    {Currency::USD, "2012-02-20", "NO FIX", "2012-02-22"},
    {Currency::USD, "2012-02-21", "2012-02-21", "2012-02-23"},
    {Currency::USD, "2012-02-22", "2012-02-22", "2012-02-24"},
    {Currency::USD, "2012-02-23", "2012-02-23", "2012-02-27"},
    {Currency::USD, "2012-04-02", "2012-04-02", "2012-04-04"},
    {Currency::USD, "2012-04-03", "2012-04-03", "2012-04-05"},
    {Currency::USD, "2012-04-04", "2012-04-04", "2012-04-10"},
    {Currency::USD, "2012-04-05", "2012-04-05", "2012-04-11"},
    {Currency::USD, "2012-04-06", "No Fix", "No Fix"},
    {Currency::USD, "2012-04-09", "No Fix", "No Fix"},
    {Currency::USD, "2012-04-10", "2012-04-10", "2012-04-12"},
    {Currency::USD, "2012-04-11", "2012-04-11", "2012-04-13"},
    {Currency::USD, "2012-04-12", "2012-04-12", "2012-04-16"},
    {Currency::USD, "2012-05-02", "2012-05-02", "2012-05-04"},
    {Currency::USD, "2012-05-03", "2012-05-03", "2012-05-08"},
    {Currency::USD, "2012-05-04", "2012-05-04", "2012-05-09"},
    {Currency::USD, "2012-05-07", "No Fix", "No Fix"},
    {Currency::USD, "2012-05-08", "2012-05-08", "2012-05-10"},
    {Currency::USD, "2012-05-09", "2012-05-09", "2012-05-11"},
    {Currency::USD, "2012-05-10", "2012-05-10", "2012-05-14"},
    {Currency::USD, "2012-05-23", "2012-05-23", "2012-05-25"},
    {Currency::USD, "2012-05-24", "2012-05-24", "2012-05-29"},
    {Currency::USD, "2012-05-25", "2012-05-25", "2012-05-29"},
    {Currency::USD, "2012-05-28", "NO FIX", "2012-05-30"},
    {Currency::USD, "2012-05-29", "2012-05-29", "2012-05-31"},
    {Currency::USD, "2012-05-30", "2012-05-30", "2012-06-01"},
    {Currency::USD, "2012-05-31", "2012-05-31", "2012-06-06"},
    {Currency::USD, "2012-06-01", "2012-06-01", "2012-06-07"},
    {Currency::USD, "2012-06-04", "No Fix", "No Fix"},
    {Currency::USD, "2012-06-05", "No Fix", "No Fix"},
    {Currency::USD, "2012-06-06", "2012-06-06", "2012-06-08"},
    {Currency::USD, "2012-06-07", "2012-06-07", "2012-06-11"},
    {Currency::USD, "2012-06-08", "2012-06-08", "2012-06-12"},
    {Currency::USD, "2012-06-29", "2012-06-29", "2012-07-03"},
    {Currency::USD, "2012-07-02", "2012-07-02", "2012-06-05"},
    {Currency::USD, "2012-07-03", "2012-07-03", "2012-07-05"},
    {Currency::USD, "2012-07-04", "NO FIX", "2012-07-06"},
    {Currency::USD, "2012-07-05", "2012-07-05", "2012-07-09"},
    {Currency::USD, "2012-07-06", "2012-07-06", "2012-07-10"},
    {Currency::USD, "2012-07-09", "2012-07-09", "2012-07-11"},
    {Currency::USD, "2012-08-21", "2012-08-21", "2012-08-23"},
    {Currency::USD, "2012-08-22", "2012-08-22", "2012-08-24"},
    {Currency::USD, "2012-08-23", "2012-08-23", "2012-08-28"},
    {Currency::USD, "2012-08-24", "2012-08-24", "2012-08-29"},
    {Currency::USD, "2012-08-27", "No Fix", "No Fix"},
    {Currency::USD, "2012-08-28", "2012-08-28", "2012-08-30"},
    {Currency::USD, "2012-08-29", "2012-08-29", "2012-08-31"},
    {Currency::USD, "2012-08-30", "2012-08-30", "2012-09-04"},
    {Currency::USD, "2012-08-31", "2012-08-31", "2012-09-04"},
    {Currency::USD, "2012-09-03", "NO FIX", "2012-09-05"},
    {Currency::USD, "2012-09-04", "2012-09-04", "2012-09-06"},
    {Currency::USD, "2012-09-05", "2012-09-05", "2012-09-07"},
    {Currency::USD, "2012-09-06", "2012-09-06", "2012-09-10"},
    {Currency::USD, "2012-10-02", "2012-10-02", "2012-10-04"},
    {Currency::USD, "2012-10-03", "2012-10-03", "2012-10-05"},
    {Currency::USD, "2012-10-04", "2012-10-04", "2012-10-09"},
    {Currency::USD, "2012-10-05", "2012-10-05", "2012-10-09"},
    {Currency::USD, "2012-10-08", "NO FIX", "2012-10-10"},
    {Currency::USD, "2012-10-09", "2012-10-09", "2012-10-11"},
    {Currency::USD, "2012-10-10", "2012-10-10", "2012-10-12"},
    {Currency::USD, "2012-10-11", "2012-10-11", "2012-10-15"},
    {Currency::USD, "2012-11-06", "2012-11-06", "2012-11-08"},
    {Currency::USD, "2012-11-07", "2012-11-07", "2012-11-09"},
    {Currency::USD, "2012-11-08", "2012-11-08", "2012-11-13"},
    {Currency::USD, "2012-11-09", "2012-11-09", "2012-11-13"},
    {Currency::USD, "2012-11-12", "NO FIX", "2012-11-14"},
    {Currency::USD, "2012-11-13", "2012-11-13", "2012-11-15"},
    {Currency::USD, "2012-11-14", "2012-11-14", "2012-11-16"},
    {Currency::USD, "2012-11-15", "2012-11-15", "2012-11-19"},
    {Currency::USD, "2012-11-16", "2012-11-16", "2012-11-20"},
    {Currency::USD, "2012-11-19", "2012-11-19", "2012-11-21"},
    {Currency::USD, "2012-11-20", "2012-11-20", "2012-11-23"},
    {Currency::USD, "2012-11-21", "2012-11-21", "2012-11-23"},
    {Currency::USD, "2012-11-22", "NO FIX", "2012-11-26"},
    {Currency::USD, "2012-11-23", "2012-11-23", "2012-11-27"},
    {Currency::USD, "2012-11-26", "2012-11-26", "2012-11-28"},
    {Currency::USD, "2012-11-27", "2012-11-27", "2012-11-29"},
    {Currency::USD, "2012-11-28", "2012-11-28", "2012-11-30"},
    {Currency::USD, "2012-11-29", "2012-11-29", "2012-12-03"},
    {Currency::USD, "2012-11-30", "2012-11-30", "2012-12-04"},
    {Currency::USD, "2012-12-03", "2012-12-03", "2012-12-05"},
    {Currency::USD, "2012-12-04", "2012-12-04", "2012-12-06"},
    {Currency::USD, "2012-12-05", "2012-12-05", "2012-12-07"},
    {Currency::USD, "2012-12-06", "2012-12-06", "2012-12-10"},
    {Currency::USD, "2012-12-07", "2012-12-07", "2012-12-11"},
    {Currency::USD, "2012-12-19", "2012-12-19", "2012-12-21"},
    {Currency::USD, "2012-12-20", "2012-12-20", "2012-12-24"},
    {Currency::USD, "2012-12-21", "2012-12-21", "2012-12-27"},
    {Currency::USD, "2012-12-24", "2012-12-24", "2012-12-28"},
    {Currency::USD, "2012-12-25", "No Fix", "No Fix"},
    {Currency::USD, "2012-12-26", "No Fix", "No Fix"},
    {Currency::USD, "2012-12-27", "2012-12-27", "2012-12-31"},
    {Currency::USD, "2012-12-28", "2012-12-28", "2013-01-01"},
    {Currency::USD, "2012-12-31", "2012-12-31", "2013-01-02"},
    {Currency::USD, "2013-01-01", "No Fix", "No Fix"},
    {Currency::EUR, "2012-01-02", "2012-01-02", "2012-01-04"},
    {Currency::EUR, "2012-01-03", "2012-01-03", "2012-01-05"},
    {Currency::EUR, "2012-01-26", "2012-01-26", "2012-01-30"},
    {Currency::EUR, "2012-01-27", "2012-01-27", "2012-01-31"},
    {Currency::EUR, "2012-01-30", "2012-01-30", "2012-02-01"},
    {Currency::EUR, "2012-01-31", "2012-01-31", "2012-02-02"},
    {Currency::EUR, "2012-02-01", "2012-02-01", "2012-02-03"},
    {Currency::EUR, "2012-02-02", "2012-02-02", "2012-02-06"},
    {Currency::EUR, "2012-02-03", "2012-02-03", "2012-02-07"},
    {Currency::EUR, "2012-02-27", "2012-02-27", "2012-02-29"},
    {Currency::EUR, "2012-02-28", "2012-02-28", "2012-03-01"},
    {Currency::EUR, "2012-02-29", "2012-02-29", "2012-03-02"},
    {Currency::EUR, "2012-03-01", "2012-03-01", "2012-03-05"},
    {Currency::EUR, "2012-03-02", "2012-03-02", "2012-03-06"},
    {Currency::EUR, "2012-03-28", "2012-03-28", "2012-03-30"},
    {Currency::EUR, "2012-03-29", "2012-03-29", "2012-04-02"},
    {Currency::EUR, "2012-03-30", "2012-03-30", "2012-04-03"},
    {Currency::EUR, "2012-04-02", "2012-04-02", "2012-04-04"},
    {Currency::EUR, "2012-04-03", "2012-04-03", "2012-04-05"},
    {Currency::EUR, "2012-04-04", "2012-04-04", "2012-04-10"},
    {Currency::EUR, "2012-04-05", "2012-04-05", "2012-04-11"},
    {Currency::EUR, "2012-04-06", "No Fix", "No Fix"},
    {Currency::EUR, "2012-04-09", "No Fix", "No Fix"},
    {Currency::EUR, "2012-04-10", "2012-04-10", "2012-04-12"},
    {Currency::EUR, "2012-04-11", "2012-04-11", "2012-04-13"},
    {Currency::EUR, "2012-04-12", "2012-04-12", "2012-04-16"},
    {Currency::EUR, "2012-04-13", "2012-04-13", "2012-04-17"},
    {Currency::EUR, "2012-04-26", "2012-04-26", "2012-04-30"},
    {Currency::EUR, "2012-04-27", "2012-04-27", "2012-05-02"},
    {Currency::EUR, "2012-04-30", "2012-04-30", "2012-05-02"},
    {Currency::EUR, "2012-05-01", "NO FIX", "2012-05-03"},
    {Currency::EUR, "2012-05-02", "2012-05-02", "2012-05-04"},
    {Currency::EUR, "2012-05-03", "2012-05-03", "2012-05-07"},
    {Currency::EUR, "2012-05-04", "2012-05-04", "2012-05-08"},
    {Currency::EUR, "2012-05-07", "2012-05-07", "2012-05-09"},
    {Currency::EUR, "2012-12-20", "2012-12-20", "2012-12-24"},
    {Currency::EUR, "2012-12-21", "2012-12-21", "2012-12-27"},
    {Currency::EUR, "2012-12-24", "2012-12-24", "2012-12-28"},
    {Currency::EUR, "2012-12-25", "No Fix", "No Fix"},
    {Currency::EUR, "2012-12-26", "No Fix", "No Fix"},
    {Currency::EUR, "2012-12-27", "2012-12-27", "2012-12-31"},
    {Currency::EUR, "2012-12-28", "2012-12-28", "2013-01-01"},
    {Currency::EUR, "2012-12-31", "2012-12-31", "2013-01-02"},
    {Currency::EUR, "2013-01-01", "No Fix", "No Fix"},
    {Currency::GBP, "2012-01-02", "No Fix", "No Fix"},
    {Currency::GBP, "2012-01-03", "2012-01-03", "2012-01-03"},
    {Currency::GBP, "2012-01-04", "2012-01-04", "2012-01-04"},
    {Currency::GBP, "2012-01-05", "2012-01-05", "2012-01-05"},
    {Currency::GBP, "2012-01-06", "2012-01-06", "2012-01-06"},
    {Currency::GBP, "2012-01-27", "2012-01-27", "2012-01-27"},
    {Currency::GBP, "2012-01-30", "2012-01-30", "2012-01-30"},
    {Currency::GBP, "2012-01-31", "2012-01-31", "2012-01-31"},
    {Currency::GBP, "2012-02-01", "2012-02-01", "2012-02-01"},
    {Currency::GBP, "2012-02-02", "2012-02-02", "2012-02-02"},
    {Currency::GBP, "2012-02-24", "2012-02-24", "2012-02-24"},
    {Currency::GBP, "2012-02-27", "2012-02-27", "2012-02-27"},
    {Currency::GBP, "2012-02-28", "2012-02-28", "2012-02-28"},
    {Currency::GBP, "2012-02-29", "2012-02-29", "2012-02-29"},
    {Currency::GBP, "2012-03-01", "2012-03-01", "2012-03-01"},
    {Currency::GBP, "2012-03-02", "2012-03-02", "2012-03-02"},
    {Currency::GBP, "2012-03-05", "2012-03-05", "2012-03-05"},
    {Currency::GBP, "2012-03-28", "2012-03-28", "2012-03-28"},
    {Currency::GBP, "2012-03-29", "2012-03-29", "2012-03-29"},
    {Currency::GBP, "2012-03-30", "2012-03-30", "2012-03-30"},
    {Currency::GBP, "2012-04-02", "2012-04-02", "2012-04-02"},
    {Currency::GBP, "2012-04-03", "2012-04-03", "2012-04-03"},
    {Currency::GBP, "2012-04-04", "2012-04-04", "2012-04-04"},
    {Currency::GBP, "2012-04-05", "2012-04-05", "2012-04-05"},
    {Currency::GBP, "2012-04-06", "No Fix", "No Fix"},
    {Currency::GBP, "2012-04-09", "No Fix", "No Fix"},
    {Currency::GBP, "2012-04-10", "2012-04-10", "2012-04-10"},
    {Currency::GBP, "2012-04-11", "2012-04-11", "2012-04-11"},
    {Currency::GBP, "2012-04-12", "2012-04-12", "2012-04-12"},
    {Currency::GBP, "2012-05-02", "2012-05-02", "2012-05-02"},
    {Currency::GBP, "2012-05-03", "2012-05-03", "2012-05-03"},
    {Currency::GBP, "2012-05-04", "2012-05-04", "2012-05-04"},
    {Currency::GBP, "2012-05-07", "No Fix", "No Fix"},
    {Currency::GBP, "2012-05-08", "2012-05-08", "2012-05-08"},
    {Currency::GBP, "2012-05-09", "2012-05-09", "2012-05-09"},
    {Currency::GBP, "2012-05-10", "2012-05-10", "2012-05-10"},
    {Currency::GBP, "2012-05-30", "2012-05-30", "2012-05-30"},
    {Currency::GBP, "2012-05-31", "2012-05-31", "2012-05-31"},
    {Currency::GBP, "2012-06-01", "2012-06-01", "2012-06-01"},
    {Currency::GBP, "2012-06-04", "No Fix", "No Fix"},
    {Currency::GBP, "2012-06-05", "No Fix", "No Fix"},
    {Currency::GBP, "2012-06-06", "2012-06-06", "2012-06-06"},
    {Currency::GBP, "2012-06-07", "2012-06-07", "2012-06-07"},
    {Currency::GBP, "2012-06-08", "2012-06-08", "2012-06-08"},
    {Currency::GBP, "2012-08-22", "2012-08-22", "2012-08-22"},
    {Currency::GBP, "2012-08-23", "2012-08-23", "2012-08-23"},
    {Currency::GBP, "2012-08-24", "2012-08-24", "2012-08-24"},
    {Currency::GBP, "2012-08-27", "No Fix", "No Fix"},
    {Currency::GBP, "2012-08-28", "2012-08-28", "2012-08-28"},
    {Currency::GBP, "2012-08-29", "2012-08-29", "2012-08-29"},
    {Currency::GBP, "2012-08-30", "2012-08-30", "2012-08-30"},
    {Currency::GBP, "2012-08-31", "2012-08-31", "2012-08-31"},
    {Currency::GBP, "2012-09-03", "2012-09-03", "2012-09-03"},
    {Currency::GBP, "2012-09-04", "2012-09-04", "2012-09-04"},
    {Currency::GBP, "2012-09-05", "2012-09-05", "2012-09-05"},
    {Currency::GBP, "2012-12-19", "2012-12-19", "2012-12-19"},
    {Currency::GBP, "2012-12-20", "2012-12-20", "2012-12-20"},
    {Currency::GBP, "2012-12-21", "2012-12-21", "2012-12-21"},
    {Currency::GBP, "2012-12-24", "2012-12-24", "2012-12-24"},
    {Currency::GBP, "2012-12-25", "No Fix", "No Fix"},
    {Currency::GBP, "2012-12-26", "No Fix", "No Fix"},
    {Currency::GBP, "2012-12-27", "2012-12-27", "2012-12-27"},
    {Currency::GBP, "2012-12-28", "2012-12-28", "2012-12-28"},
    {Currency::GBP, "2012-12-31", "2012-12-31", "2012-12-31"},
    {Currency::GBP, "2013-01-01", "No Fix", "No Fix"},
    {Currency::USD, "2013-01-01", "No Fix", "No Fix"},
    {Currency::USD, "2013-01-02", "2013-01-02", "2013-01-04"},
    {Currency::USD, "2013-01-03", "2013-01-03", "2013-01-07"},
    {Currency::USD, "2013-01-04", "2013-01-04", "2013-01-08"},
    {Currency::USD, "2013-01-07", "2013-01-07", "2013-01-09"},
    {Currency::USD, "2013-01-08", "2013-01-08", "2013-01-10"},
    {Currency::USD, "2013-01-09", "2013-01-09", "2013-01-11"},
    {Currency::USD, "2013-01-17", "2013-01-17", "2013-01-22"},
    {Currency::USD, "2013-01-18", "2013-01-18", "2013-01-22"},
    {Currency::USD, "2013-01-21", "No Fix", "2013-01-23"},
    {Currency::USD, "2013-01-22", "2013-01-22", "2013-01-24"},
    {Currency::USD, "2013-01-23", "2013-01-23", "2013-01-25"},
    {Currency::USD, "2013-01-24", "2013-01-24", "2013-01-28"},
    {Currency::USD, "2013-01-25", "2013-01-25", "2013-01-29"},
    {Currency::USD, "2013-01-28", "2013-01-28", "2013-01-30"},
    {Currency::USD, "2013-01-29", "2013-01-29", "2013-01-31"},
    {Currency::USD, "2013-01-30", "2013-01-30", "2013-02-01"},
    {Currency::USD, "2013-01-31", "2013-01-31", "2013-02-04"},
    {Currency::USD, "2013-02-01", "2013-02-01", "2013-02-05"},
    {Currency::USD, "2013-02-04", "2013-02-04", "2013-02-06"},
    {Currency::USD, "2013-02-05", "2013-02-05", "2013-02-07"},
    {Currency::USD, "2013-02-12", "2013-02-12", "2013-02-14"},
    {Currency::USD, "2013-02-13", "2013-02-13", "2013-02-15"},
    {Currency::USD, "2013-02-14", "2013-02-14", "2013-02-19"},
    {Currency::USD, "2013-02-15", "2013-02-15", "2013-02-19"},
    {Currency::USD, "2013-02-18", "No Fix", "2013-02-20"},
    {Currency::USD, "2013-02-19", "2013-02-19", "2013-02-21"},
    {Currency::USD, "2013-02-20", "2013-02-20", "2013-02-22"},
    {Currency::USD, "2013-02-21", "2013-02-21", "2013-02-25"},
    {Currency::USD, "2013-02-22", "2013-02-22", "2013-02-26"},
    {Currency::USD, "2013-02-25", "2013-02-25", "2013-02-27"},
    {Currency::USD, "2013-02-26", "2013-02-26", "2013-02-28"},
    {Currency::USD, "2013-02-27", "2013-02-27", "2013-03-01"},
    {Currency::USD, "2013-02-28", "2013-02-28", "2013-03-04"},
    {Currency::USD, "2013-03-01", "2013-03-01", "2013-03-05"},
    {Currency::USD, "2013-03-04", "2013-03-04", "2013-03-06"},
    {Currency::USD, "2013-03-05", "2013-03-05", "2013-03-07"},
    {Currency::USD, "2013-03-26", "2013-03-26", "2013-03-28"},
    {Currency::USD, "2013-03-27", "2013-03-27", "2013-04-02"},
    {Currency::USD, "2013-03-28", "2013-03-28", "2013-04-03"},
    {Currency::USD, "2013-03-29", "No Fix", "No Fix"},
    {Currency::USD, "2013-04-01", "No Fix", "No Fix"},
    {Currency::USD, "2013-04-02", "2013-04-02", "2013-04-04"},
    {Currency::USD, "2013-04-03", "2013-04-03", "2013-04-05"},
    {Currency::USD, "2013-04-04", "2013-04-04", "2013-04-08"},
    {Currency::USD, "2013-05-01", "2013-05-01", "2013-05-03"},
    {Currency::USD, "2013-05-02", "2013-05-02", "2013-05-07"},
    {Currency::USD, "2013-05-03", "2013-05-03", "2013-05-08"},
    {Currency::USD, "2013-05-06", "No Fix", "No Fix"},
    {Currency::USD, "2013-05-07", "2013-05-07", "2013-05-09"},
    {Currency::USD, "2013-05-08", "2013-05-08", "2013-05-10"},
    {Currency::USD, "2013-05-09", "2013-05-09", "2013-05-13"},
    {Currency::USD, "2013-05-22", "2013-05-22", "2013-05-24"},
    {Currency::USD, "2013-05-23", "2013-05-23", "2013-05-28"},
    {Currency::USD, "2013-05-24", "2013-05-24", "2013-05-29"},
    {Currency::USD, "2013-05-27", "No Fix", "No Fix"},
    {Currency::USD, "2013-05-28", "2013-05-28", "2013-05-30"},
    {Currency::USD, "2013-05-29", "2013-05-29", "2013-05-31"},
    {Currency::USD, "2013-05-30", "2013-05-30", "2013-06-03"},
    {Currency::USD, "2013-05-31", "2013-05-31", "2013-06-04"},
    {Currency::USD, "2013-07-01", "2013-07-01", "2013-07-03"},
    {Currency::USD, "2013-07-02", "2013-07-02", "2013-07-05"},
    {Currency::USD, "2013-07-03", "2013-07-03", "2013-07-05"},
    {Currency::USD, "2013-07-04", "No Fix", "2013-07-08"},
    {Currency::USD, "2013-07-05", "2013-07-05", "2013-07-09"},
    {Currency::USD, "2013-07-08", "2013-07-08", "2013-07-10"},
    {Currency::USD, "2013-07-09", "2013-07-09", "2013-07-11"},
    {Currency::USD, "2013-08-21", "2013-08-21", "2013-08-23"},
    {Currency::USD, "2013-08-22", "2013-08-22", "2013-08-27"},
    {Currency::USD, "2013-08-23", "2013-08-23", "2013-08-28"},
    {Currency::USD, "2013-08-26", "No Fix", "No Fix"},
    {Currency::USD, "2013-08-27", "2013-08-27", "2013-08-29"},
    {Currency::USD, "2013-08-28", "2013-08-28", "2013-08-30"},
    {Currency::USD, "2013-08-29", "2013-08-29", "2013-09-03"},
    {Currency::USD, "2013-08-30", "2013-08-30", "2013-09-03"},
    {Currency::USD, "2013-09-02", "No Fix", "2013-09-04"},
    {Currency::USD, "2013-09-03", "2013-09-03", "2013-09-05"},
    {Currency::USD, "2013-09-04", "2013-09-04", "2013-09-06"},
    {Currency::USD, "2013-09-05", "2013-09-05", "2013-09-09"},
    {Currency::USD, "2013-09-06", "2013-09-06", "2013-09-10"},
    {Currency::USD, "2013-10-08", "2013-10-08", "2013-10-10"},
    {Currency::USD, "2013-10-09", "2013-10-09", "2013-10-11"},
    {Currency::USD, "2013-10-10", "2013-10-10", "2013-10-15"},
    {Currency::USD, "2013-10-11", "2013-10-11", "2013-10-15"},
    {Currency::USD, "2013-10-14", "No Fix", "2013-10-16"},
    {Currency::USD, "2013-10-15", "2013-10-15", "2013-10-17"},
    {Currency::USD, "2013-10-16", "2013-10-16", "2013-10-18"},
    {Currency::USD, "2013-10-17", "2013-10-17", "2013-10-21"},
    {Currency::USD, "2013-11-06", "2013-11-06", "2013-11-08"},
    {Currency::USD, "2013-11-07", "2013-11-07", "2012-11-12"},
    {Currency::USD, "2013-11-08", "2013-11-08", "2013-11-12"},
    {Currency::USD, "2013-11-11", "No Fix", "2013-11-13"},
    {Currency::USD, "2013-11-12", "2013-11-12", "2013-11-14"},
    {Currency::USD, "2013-11-13", "2013-11-13", "2013-11-15"},
    {Currency::USD, "2013-11-14", "2013-11-14", "2013-11-18"},
    {Currency::USD, "2013-11-25", "2013-11-25", "2013-11-27"},
    {Currency::USD, "2013-11-26", "2013-11-26", "2012-11-29"},
    {Currency::USD, "2013-11-27", "2013-11-27", "2013-11-29"},
    {Currency::USD, "2013-11-28", "No Fix", "2013-12-02"},
    {Currency::USD, "2013-11-29", "2013-11-29", "2013-12-03"},
    {Currency::USD, "2013-12-02", "2013-12-02", "2013-12-04"},
    {Currency::USD, "2013-12-03", "2013-12-03", "2013-12-05"},
    {Currency::USD, "2013-12-20", "2013-12-20", "2013-12-24"},
    {Currency::USD, "2013-12-23", "2013-12-23", "2013-12-27"},
    {Currency::USD, "2013-12-24", "2013-12-24", "2013-12-30"},
    {Currency::USD, "2013-12-25", "No Fix", "No Fix"},
    {Currency::USD, "2013-12-26", "No Fix", "No Fix"},
    {Currency::USD, "2013-12-27", "2013-12-27", "2013-12-31"},
    {Currency::USD, "2013-12-30", "2013-12-30", "2014-01-02"},
    {Currency::USD, "2013-12-31", "2013-12-31", "2014-01-03"},
    {Currency::EUR, "2013-01-01", "No Fix", "No Fix"},
    {Currency::EUR, "2013-01-02", "2013-01-02", "2013-01-04"},
    {Currency::EUR, "2013-01-03", "2013-01-03", "2013-01-07"},
    {Currency::EUR, "2013-01-04", "2013-01-04", "2013-01-08"},
    {Currency::EUR, "2013-01-07", "2013-01-07", "2013-01-09"},
    {Currency::EUR, "2013-03-26", "2013-03-26", "2013-03-28"},
    {Currency::EUR, "2013-03-27", "2013-03-27", "2013-04-02"},
    {Currency::EUR, "2013-03-28", "2013-03-28", "2013-04-03"},
    {Currency::EUR, "2013-03-29", "No Fix", "No Fix"},
    {Currency::EUR, "2013-04-01", "No Fix", "No Fix"},
    {Currency::EUR, "2013-04-02", "2013-04-02", "2013-04-04"},
    {Currency::EUR, "2013-04-03", "2013-04-03", "2013-04-05"},
    {Currency::EUR, "2013-04-04", "2013-04-04", "2013-04-08"},
    {Currency::EUR, "2013-04-26", "2013-04-26", "2013-04-30"},
    {Currency::EUR, "2013-04-29", "2013-04-29", "2013-05-01"},
    {Currency::EUR, "2013-04-30", "2013-04-30", "2013-05-02"},
    {Currency::EUR, "2013-05-01", "No Fix", "2013-05-03"},
    {Currency::EUR, "2013-05-02", "2013-05-02", "2013-05-07"},
    {Currency::EUR, "2013-05-03", "2013-05-03", "2013-05-08"},
    {Currency::EUR, "2013-05-06", "No Fix", "No Fix"},
    {Currency::EUR, "2013-05-07", "2013-05-07", "2013-05-09"},
    {Currency::EUR, "2013-05-08", "2013-05-08", "2013-05-10"},
    {Currency::EUR, "2013-05-09", "2013-05-09", "2013-05-13"},
    {Currency::EUR, "2013-05-22", "2013-05-22", "2013-05-24"},
    {Currency::EUR, "2013-05-23", "2013-05-23", "2013-05-28"},
    {Currency::EUR, "2013-05-24", "2013-05-24", "2013-05-29"},
    {Currency::EUR, "2013-05-27", "No Fix", "No Fix"},
    {Currency::EUR, "2013-05-28", "2013-05-28", "2013-05-30"},
    {Currency::EUR, "2013-05-29", "2013-05-29", "2013-05-31"},
    {Currency::EUR, "2013-05-30", "2013-05-30", "2013-06-03"},
    {Currency::EUR, "2013-08-21", "2013-08-21", "2013-08-23"},
    {Currency::EUR, "2013-08-22", "2013-08-22", "2013-08-27"},
    {Currency::EUR, "2013-08-23", "2013-08-23", "2013-08-28"},
    {Currency::EUR, "2013-08-26", "No Fix", "No Fix"},
    {Currency::EUR, "2013-08-27", "2013-08-27", "2013-08-29"},
    {Currency::EUR, "2013-08-28", "2013-08-28", "2013-08-30"},
    {Currency::EUR, "2013-08-29", "2013-08-29", "2013-09-02"},
    {Currency::EUR, "2013-12-20", "2013-12-20", "2013-12-24"},
    {Currency::EUR, "2013-12-23", "2013-12-23", "2013-12-27"},
    {Currency::EUR, "2013-12-24", "2013-12-24", "2013-12-28"},
    {Currency::EUR, "2013-12-25", "No Fix", "No Fix"},
    {Currency::EUR, "2013-12-26", "No Fix", "No Fix"},
    {Currency::EUR, "2013-12-27", "2013-12-27", "2013-12-31"},
    {Currency::EUR, "2013-12-30", "2013-12-30", "2014-01-02"},
    {Currency::EUR, "2013-12-31", "2013-12-31", "2014-01-03"},
    {Currency::GBP, "2013-01-01", "No Fix", "No Fix"},
    {Currency::GBP, "2013-01-02", "2013-01-02", "2013-01-02"},
    {Currency::GBP, "2013-01-03", "2013-01-03", "2013-01-03"},
    {Currency::GBP, "2013-01-04", "2013-01-04", "2013-01-04"},
    {Currency::GBP, "2013-03-25", "2013-03-25", "2013-03-25"},
    {Currency::GBP, "2013-03-26", "2013-03-26", "2013-03-26"},
    {Currency::GBP, "2013-03-27", "2013-03-27", "2013-03-27"},
    {Currency::GBP, "2013-03-28", "2013-03-28", "2013-03-28"},
    {Currency::GBP, "2013-03-29", "No Fix", "No Fix"},
    {Currency::GBP, "2013-04-01", "No Fix", "No Fix"},
    {Currency::GBP, "2013-04-02", "2013-04-02", "2013-04-02"},
    {Currency::GBP, "2013-04-03", "2013-04-03", "2013-04-03"},
    {Currency::GBP, "2013-04-04", "2013-04-04", "2013-04-04"},
    {Currency::GBP, "2013-05-01", "2013-05-01", "2013-05-01"},
    {Currency::GBP, "2013-05-02", "2013-05-02", "2013-05-02"},
    {Currency::GBP, "2013-05-03", "2013-05-03", "2013-05-03"},
    {Currency::GBP, "2013-05-06", "No Fix", "No Fix"},
    {Currency::GBP, "2013-05-07", "2013-05-07", "2013-05-07"},
    {Currency::GBP, "2013-05-08", "2013-05-08", "2013-05-08"},
    {Currency::GBP, "2013-05-09", "2013-05-09", "2013-05-09"},
    {Currency::GBP, "2013-05-22", "2013-05-22", "2013-05-22"},
    {Currency::GBP, "2013-05-23", "2013-05-23", "2013-05-23"},
    {Currency::GBP, "2013-05-24", "2013-05-24", "2013-05-24"},
    {Currency::GBP, "2013-05-27", "No Fix", "No Fix"},
    {Currency::GBP, "2013-05-28", "2013-05-28", "2013-05-28"},
    {Currency::GBP, "2013-05-29", "2013-05-29", "2013-05-29"},
    {Currency::GBP, "2013-05-30", "2013-05-30", "2013-05-30"},
    {Currency::GBP, "2013-08-21", "2013-08-21", "2013-08-21"},
    {Currency::GBP, "2013-08-22", "2013-08-22", "2013-08-22"},
    {Currency::GBP, "2013-08-23", "2013-08-23", "2013-08-23"},
    {Currency::GBP, "2013-08-26", "No Fix", "No Fix"},
    {Currency::GBP, "2013-08-27", "2013-08-27", "2013-08-27"},
    {Currency::GBP, "2013-08-28", "2013-08-28", "2013-08-28"},
    {Currency::GBP, "2013-08-29", "2013-08-29", "2013-08-29"},
    {Currency::GBP, "2013-12-20", "2013-12-20", "2013-12-20"},
    {Currency::GBP, "2013-12-23", "2013-12-23", "2013-12-23"},
    {Currency::GBP, "2013-12-24", "2013-12-24", "2013-12-24"},
    {Currency::GBP, "2013-12-25", "No Fix", "No Fix"},
    {Currency::GBP, "2013-12-26", "No Fix", "No Fix"},
    {Currency::GBP, "2013-12-27", "2013-12-27", "2013-12-27"},
    {Currency::GBP, "2013-12-30", "2013-12-30", "2013-12-30"},
    {Currency::GBP, "2013-12-31", "2013-12-31", "2013-12-31"},
    {Currency::CURRENCY_UNSPECIFIED, nullptr, nullptr, nullptr}};
#endif

int test_index1()
{
	IndexService *service = get_default_index_service();

	auto index_1d = service->get_index(IsdaIndex::USD_LIBOR_BBA, TENOR_1D);
	if (index_1d == nullptr)
		return 1;
	auto index_1d_2 = service->get_index(Currency::USD, IndexFamily::LIBOR, TENOR_1D);
	if (index_1d_2 != index_1d)
		return 1;
	auto index_1m = service->get_index(IsdaIndex::USD_LIBOR_BBA, TENOR_1M);
	if (index_1m == nullptr)
		return 1;
	auto index_1m_2 = service->get_index(Currency::USD, IndexFamily::LIBOR, TENOR_1M);
	if (index_1m_2 != index_1m)
		return 1;
	Date fixingDate = make_date(23, 12, 2009);
	auto valueDate = index_1d->value_date(fixingDate);
	if (valueDate != fixingDate)
		return 1;
	Date expected_value_date = make_date(29, 12, 2009);
	if (index_1m->value_date(fixingDate) != expected_value_date)
		return 1;
	auto fedfund_1y = service->get_index(Currency::USD, IndexFamily::FEDFUND, TENOR_12M);
	if (fedfund_1y == nullptr)
		return 1;
	return 0;
}

int test_index()
{
	int failure_count = 0;
	failure_count += test_index1();
	if (failure_count != 0)
		printf("Index Tests FAILED\n");
	else
		printf("Index Tests OK\n");
	return failure_count;
}

} // namespace redukti
