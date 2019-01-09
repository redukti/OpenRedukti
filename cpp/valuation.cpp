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
#include <cashflow.h>
#include <cashflow_pricing.h>
#include <converters.h>
#include <curve.h>
#include <date.h>
#include <dayfractions.h>
#include <fixings.h>
#include <logger.h>
#include <matrix.h>
#include <status.h>
#include <valuation.h>

#include <internal/cashflow_internal.h>

#include <cmath>
#include <condition_variable>
#include <map>
#include <mutex>
#include <shared_mutex>

using namespace google::protobuf;

namespace redukti
{

struct Matrix {
	int32_t m; /* rows */
	int32_t n; /* columns */
	double data[];

	Matrix(int32_t rows, int32_t cols)
	{
		m = rows;
		n = cols;
		std::fill_n(data, m * n, 0.0);
		// debug("Matrix created\n");
	}
	~Matrix() { /* debug("Matrix destroyed\n");*/}
};

class PricingEnv
{
	Date business_date_;
	Date payment_cutoff_date_;
	int derivative_order_;
	bool is_todays_fixings_included_;
	int cycle_;
	CurveGroup curve_group_;
	MarketDataQualifier qualifier_;
	int scenario_;
	std::map<CurveId, std::unique_ptr<YieldCurve, Deleter<YieldCurve>>> curves_;
};

class YieldCurveReference : public CurveReference
{
      private:
	std::unique_ptr<YieldCurve, Deleter<YieldCurve>> yield_curve_;

      public:
	YieldCurveReference(std::unique_ptr<YieldCurve, Deleter<YieldCurve>> yield_curve)
	    : yield_curve_(std::move(yield_curve))
	{
	}
	YieldCurve *get() const noexcept override final { return yield_curve_.get(); }
};

class CurvesByGroup
{
      private:
	CurveGroup curve_group_;
	std::map<CurveId, std::unique_ptr<YieldCurveReference>> curves_;
	// A mapping from curve ID to the original definition id,
	// this can be used to return values back to the client
	// referencing the original id
	std::map<CurveId, int> curve_to_definition_id_;
	// reverse map
	std::map<int, CurveId> definition_id_to_curve_;
	std::map<CurveId, std::unique_ptr<Matrix, Deleter<Matrix>>> par_sensitivities_;

      public:
	CurvesByGroup(CurveGroup group) : curve_group_(group) {}
	// Add a curve - this replaces any existing curve of the same id.
	// Also any par sensititivies associated with the curve will be erased.
	void add(std::unique_ptr<YieldCurve, Deleter<YieldCurve>> curve, int curve_definition_id)
	{
		CurveId id = curve->id();
		auto name = curve_id_to_string(id);
		// fprintf(stdout, "INFO Adding curve %s\n", name.c_str());
		auto &&iter = curves_.find(id);
		if (iter != curves_.end()) {
			// fprintf(stderr, "Replace curve %lld\n", (long long)id);
			curves_.erase(id);
		}
		// Also erase any par sensitivities
		auto &&pariter = par_sensitivities_.find(id);
		if (pariter != par_sensitivities_.end()) {
			par_sensitivities_.erase(id);
		}
		curves_.insert(std::pair<CurveId, std::unique_ptr<YieldCurveReference>>(
		    id, std::unique_ptr<YieldCurveReference>(new YieldCurveReference(std::move(curve)))));
		auto iter2 = curve_to_definition_id_.find(id);
		if (iter2 != curve_to_definition_id_.end()) {
			curve_to_definition_id_.erase(id);
		}
		curve_to_definition_id_.insert({id, curve_definition_id});
		auto iter3 = definition_id_to_curve_.find(curve_definition_id);
		if (iter3 != definition_id_to_curve_.end()) {
			definition_id_to_curve_.erase(curve_definition_id);
		}
		definition_id_to_curve_.insert({curve_definition_id, id});
	}
	// Add the provided par sensitivities for the curve
	// Previous par sensitivities will be replaced by the new ones
	// Note that par sensitivities are really in column major matrix
	// format with m instrument rows and n curve pillars.
	void add(CurveId id, std::unique_ptr<Matrix, Deleter<Matrix>> par_sensitivities)
	{
		auto &&iter = par_sensitivities_.find(id);
		if (iter != par_sensitivities_.end()) {
			par_sensitivities_.erase(id);
		}
		par_sensitivities_.insert(
		    std::pair<CurveId, std::unique_ptr<Matrix, Deleter<Matrix>>>(id, std::move(par_sensitivities)));
	}
	const CurveReference *get_curve(CurveId id)
	{
		auto &&iter = curves_.find(id);
		if (iter != curves_.end()) {
			return iter->second.get();
		}
		return nullptr;
	}
	Matrix *get_sensitivities(CurveId id)
	{
		auto &&iter = par_sensitivities_.find(id);
		if (iter != par_sensitivities_.end()) {
			return iter->second.get();
		}
		return nullptr;
	}
	int get_curve_definition_id(CurveId id)
	{
		auto iter = curve_to_definition_id_.find(id);
		if (iter != curve_to_definition_id_.end()) {
			return iter->second;
		}
		return -1;
	}
	CurveId get_curve_id(int curve_definition_id)
	{
		auto iter = definition_id_to_curve_.find(curve_definition_id);
		if (iter != definition_id_to_curve_.end()) {
			return iter->second;
		}
		return 0;
	}
};

class CurveProviderImpl : public CurveProvider
{
      private:
	Date as_of_date_;
	MarketDataQualifier qual_;
	short int cycle_;
	CurvesByGroup *curves_by_group_;
	short int scenario_;

      public:
	CurveProviderImpl(Date as_of_date, MarketDataQualifier qual, short int cycle, CurvesByGroup *curves_by_group,
			  short int scenario = 0)
	    : as_of_date_(as_of_date), qual_(qual), cycle_(cycle), curves_by_group_(curves_by_group),
	      scenario_(scenario)
	{
	}
	const CurveReference *get_curve(PricingCurve curve) const override final
	{
		auto curve_type = curve.curve_type();
		auto currency = curve.currency();
		auto index_family = curve.index_family();
		auto tenor = curve.tenor();
		CurveId curve_id =
		    make_curve_id(curve_type, currency, index_family, tenor, as_of_date_, cycle_, qual_, scenario_);
		CurveId curve_id2 = curve_id;
		CurveId curve_id3 = curve_id;
		const CurveReference *ref = curves_by_group_->get_curve(curve_id);
		if (ref == nullptr && tenor != TENOR_UNSPECIFIED) {
			curve_id3 = curve_id2 = make_curve_id(curve_type, currency, index_family, TENOR_UNSPECIFIED,
							      as_of_date_, cycle_, qual_, scenario_);
			ref = curves_by_group_->get_curve(curve_id2);
		}
		if (ref == nullptr && index_family != INDEX_FAMILY_UNSPECIFIED) {
			curve_id3 = make_curve_id(curve_type, currency, INDEX_FAMILY_UNSPECIFIED, TENOR_UNSPECIFIED,
						  as_of_date_, cycle_, qual_, scenario_);
			ref = curves_by_group_->get_curve(curve_id3);
		}
		if (ref == nullptr) {
			error("ERROR Requested curve not found, we tried following mappings: [%s], "
			      "[%s], [%s]\n",
			      curve_id_to_string(curve_id).c_str(), curve_id_to_string(curve_id2).c_str(),
			      curve_id_to_string(curve_id3).c_str());
		}
		return ref;
	}
};

class ValuationServiceImpl : public ValuationService
{
      private:
	std::map<CurveGroup, std::unique_ptr<SimpleCurveMapper>> curve_mappers_;
	std::map<int, std::unique_ptr<IRCurveDefinition>> curve_definitions_;
	std::map<CurveGroup, std::unique_ptr<CurvesByGroup>> curves_by_group_;
	std::unique_ptr<FixingDataService> fixings_data_service_;
	IndexService *index_service_;
	CalendarService *calendar_service_;
	// Global read/write lock for valuation service
	std::shared_mutex lock_;

      private:
	const IRCurveDefinition *find_curve_definition(int id)
	{
		auto &&iter = curve_definitions_.find(id);
		if (iter == curve_definitions_.end())
			return nullptr;
		return iter->second.get();
	}
	StatusCode validate(int as_of_date, const ZeroCurve &curve)
	{
		if (curve.maturities_size() <= 4) {
			return StatusCode::kVAL_InsufficientMaturities;
		}
		if (curve.maturities_size() != curve.values_size()) {
			return StatusCode::kVAL_MismatchedMaturitiesAndValues;
		}
		if (curve.maturities(curve.maturities_size() - 1) > as_of_date + (61 * 365)) {
			return StatusCode::kVAL_BadMaturityDate;
		}
		for (int j = 0; j < curve.maturities_size(); j++) {
			if ((j > 0 && curve.maturities(j) <= curve.maturities(j - 1)) ||
			    (j == 0 && curve.maturities(j) <= as_of_date)) {
				return StatusCode::kVAL_BadMaturityDate;
			}
			if (curve.values(j) < -0.1 || curve.values(j) > 0.2) {
				error("suspicious curve (definition id %d) value [%d] = %f\n",
				      curve.curve_definition_id(), j, curve.values(j));
				return StatusCode::kVAL_BadZeroRate;
			}
		}
		return StatusCode::kOk;
	}
	StatusCode validate(const ZeroCurve &curve, const ZeroCurveParSensitivities &sens)
	{
		if (sens.values_size() == 0) {
			return StatusCode::kVAL_MissingParSensitivitiesForCurve;
		}
		if (sens.num_maturities() != curve.maturities_size()) {
			return StatusCode::kVAL_ParSensitivitiesForCurveMismatchedMaturities;
		}
		for (auto &iter : sens.values()) {
			uint32_t n = iter.first;
			uint32_t col = n >> 16;
			uint32_t row = n & 0xFFFF;
			if (col >= sens.num_maturities()) {
				return StatusCode::kVAL_ParSensitivitiesForCurveMismatchedMaturities;
			}
			if (row >= sens.num_instruments()) {
				return StatusCode::kVAL_ParSensitivitiesForCurveInstrumentsOutOfRange;
			}
			// TODO check that row is not too much out of line compared to pillars of
			// the curve
		}
		return StatusCode::kOk;
	}

	// std::unique_ptr<YieldCurve, Deleter<YieldCurve>>
	// make_curve(Date as_of_date, const IRCurveDefinition *defn, const ZeroCurve &curve, int deriv_order,
	//		   PricingCurveType type = PRICING_CURVE_TYPE_UNSPECIFIED, MarketDataQualifier mdq = MDQ_NORMAL,
	//		   short int cycle = 0, short int scenario = 0)
	//{
	//	auto n = curve.maturities_size();
	//	int reqsize = n * sizeof(double) * 2;
	//	FixedRegionAllocator *tempalloc = get_threadspecific_allocators()->tempspace_allocator;
	//	FixedRegionAllocatorGuard guard(tempalloc);
	//	StackRegionWithFallbackAllocator<1024> buf(tempalloc);
	//	Date *maturities = (Date *)buf.safe_allocate(sizeof(Date) * n);
	//	double *values = (double *)buf.safe_allocate(sizeof(double) * n);
	//	for (int i = 0; i < n; i++) {
	//		maturities[i] = curve.maturities(i);
	//		values[i] = curve.values(i);
	//	}
	//	if (defn->interpolated_on() == IRRateType::DISCOUNT_FACTOR) {
	//		// Convert the values to discount factors as we need to
	//		// interpolate on discount factors
	//		// FIXME the day count fraction ought to be a parameter in curve
	//		// definition
	//		auto fraction = get_day_fraction(DayCountFraction::ACT_365_FIXED);
	//		for (int i = 0; i < n; i++) {
	//			double t = fraction->year_fraction(as_of_date, maturities[i]);
	//			values[i] = std::exp(-values[i] * t);
	//		}
	//	}
	//	CurveId curveId = make_curve_id(type, defn->currency(), defn->index_family(), defn->tenor(), as_of_date,
	//					cycle, mdq, scenario);
	//	return make_curve(&GlobalAllocator, curveId, as_of_date, maturities, values, n,
	//			  defn->interpolator_type(), defn->interpolated_on(), deriv_order);
	//}
	CurvesByGroup *get_group(CurveGroup group_id)
	{
		auto &&iter = curves_by_group_.find(group_id);
		CurvesByGroup *group = nullptr;
		if (iter != curves_by_group_.end()) {
			group = iter->second.get();
		}
		return group;
	}
	CurveMapper *get_mapper(CurveGroup group_id)
	{
		auto &&iter = curve_mappers_.find(group_id);
		CurveMapper *mapper = nullptr;
		if (iter != curve_mappers_.end()) {
			mapper = iter->second.get();
		}
		return mapper;
	}

	// In this method we are supposed to build the message contents
	// for Zero sensitivities but this is also a convenient place to compute
	// the Par sensitivities.
	// TODO consider refactoring
	StatusCode add_irsensitivities(Arena *arena, CurvesByGroup *curves_by_group, ValuationResult *result,
				       Sensitivities &sens)
	{
		FixedRegionAllocator *tempalloc = get_threadspecific_allocators()->tempspace_allocator;
		FixedRegionAllocatorGuard guard(tempalloc);
		StackRegionWithFallbackAllocator<1024> buffer(tempalloc);
		for (auto deltasens = sens.first_1d_sensitivities(); deltasens != nullptr;
		     deltasens = deltasens->next()) {
			buffer.release();
			auto curve = deltasens->curve();
			auto curveid = curve->id();
			auto defn_id = curves_by_group->get_curve_definition_id(curveid);
			auto curve_par_sensitivities = curves_by_group->get_sensitivities(curveid);
			if (defn_id == -1)
				return StatusCode::kInternalError;
			if (curve_par_sensitivities) {
				if (curve_par_sensitivities->n != curve->last_pillar()) {
					error("Par sensitivities matrix column does not equal number of "
					      "pillars in curve [%s, %d, %d]\n",
					      curve_id_to_string(curveid).c_str(), (int)curve_par_sensitivities->n,
					      (int)curve->last_pillar());
					return StatusCode::kInternalError;
				}
				debug("Par sensitivities available for curve [%s]\n",
				      curve_id_to_string(curveid).c_str());
			}
			PricingCurveType curve_type;
			Currency ccy;
			IndexFamily index_family;
			Tenor tenor;
			Date as_of_date;
			short int scenario;
			short int cycle;
			MarketDataQualifier qual;
			curve_id_components(curveid, curve_type, ccy, index_family, tenor, as_of_date, cycle, qual,
					    scenario);
			if (scenario != 0)
				return StatusCode::kInternalError;
			SensitivityRiskCode risktype =
			    curve_type == PRICING_CURVE_TYPE_DISCOUNT ? SRC_DISCOUNT : SRC_FORWARD;
			auto irsens = result->add_sensitivities();
			irsens->set_curve_definition_id_1(defn_id);
			irsens->set_curve_definition_id_2(0);
			irsens->set_order(SensitivityOrderCode::SOC_DELTA);
			irsens->set_risk_type_1(risktype);
			irsens->set_risk_type_2(SensitivityRiskCode::SRC_UNSPECIFIED);
			irsens->set_sensitivity_type(SensitivityTypeCode::STC_ZERO);
			double *zero_sens = nullptr; /* we need the zero deltas as a vector for
							the matrix multiplication below. */
			if (curve_par_sensitivities) {
				// The zero delta vector has a size of 'n' or pillars on the curve
				// excluding the reference date
				zero_sens = (double *)buffer.safe_allocate(sizeof(double) * curve_par_sensitivities->n);
				std::fill_n(zero_sens, curve_par_sensitivities->n, 0.0);
			}
			// We start at 1 below as 0 is the reference date and not
			// a pillar
			for (int row = 1; row < deltasens->count(); row++) {
				auto v = deltasens->at(row);
				if (v != 0.0) {
					// printf("Returning delta at pillar %d %d/%d/%d\n", row,
					// day_of_month(deltasens->curve()->maturity_date(row)),
					//  month(deltasens->curve()->maturity_date(row)),
					//  year(deltasens->curve()->maturity_date(row)));
					// The delta is scaled by 1e-4 to make it currency units
					// As the curve has pillars from as of date - the maturity offset
					// needs to be adjusted by 1
					// TODO can row=0 have sensitivity?
					uint32_t return_row = row - 1;
					irsens->mutable_values()->insert(
					    google::protobuf::MapPair<google::protobuf::uint32, double>(return_row,
													v * 1e-4));
					if (zero_sens) {
						// Save for par delta calculation
						zero_sens[return_row] = v;
					}
				}
			}
			if (curve_par_sensitivities) {
				// We need to muliply the par matrix with the zero
				// deltas to obtain par deltas
				redukti_matrix_t par_matrix{curve_par_sensitivities->m, curve_par_sensitivities->n,
							    curve_par_sensitivities->data};
				// Zero deltas are presented as a column vector
				redukti_matrix_t zero_delta_vector{curve_par_sensitivities->n, 1, zero_sens};
				double *par_sens =
				    (double *)buffer.safe_allocate(sizeof(double) * curve_par_sensitivities->m);
				std::fill_n(par_sens, curve_par_sensitivities->m, 0.0);
				// PAR deltas will be a column matrix too but different size
				redukti_matrix_t par_delta_vector{curve_par_sensitivities->m, 1, par_sens};
				debug("Multiplying par matrix with zero vector\n");
				redukti_matrix_multiply(&par_matrix, &zero_delta_vector, &par_delta_vector, false,
							false, 1.0, 1.0);
				debug("Obtained Par Sensitivities for curve %s as follows\n",
				      curve_id_to_string(curveid).c_str());
				risktype = curve_type == PRICING_CURVE_TYPE_DISCOUNT ? SRC_DISCOUNT : SRC_FORWARD;
				irsens = result->add_sensitivities();
				irsens->set_curve_definition_id_1(defn_id);
				irsens->set_curve_definition_id_2(0);
				irsens->set_order(SensitivityOrderCode::SOC_DELTA);
				irsens->set_risk_type_1(risktype);
				irsens->set_risk_type_2(SensitivityRiskCode::SRC_UNSPECIFIED);
				irsens->set_sensitivity_type(SensitivityTypeCode::STC_PAR);
				for (int i = 0; i < par_delta_vector.m; i++) {
					if (par_delta_vector.data[i] != 0.0) {
						debug("[%d] = zero %f, par %f\n", i, zero_delta_vector.data[i],
						      par_delta_vector.data[i]);
						irsens->mutable_values()->insert(
						    google::protobuf::MapPair<google::protobuf::uint32, double>(
							i, par_delta_vector.data[i] * 1e-4));
					}
				}
			}
		}
		for (auto gammasens = sens.first_2d_sensitivities(); gammasens != nullptr;
		     gammasens = gammasens->next()) {
			auto curve1 = gammasens->curve1();
			auto curveid1 = curve1->id();
			auto defn_id1 = curves_by_group->get_curve_definition_id(curveid1);
			if (defn_id1 == -1)
				return StatusCode::kInternalError;
			PricingCurveType curve_type;
			Currency ccy;
			IndexFamily index_family;
			Tenor tenor;
			Date as_of_date;
			short int scenario;
			short int cycle;
			MarketDataQualifier qual;
			curve_id_components(curveid1, curve_type, ccy, index_family, tenor, as_of_date, cycle, qual,
					    scenario);
			if (scenario != 0)
				return StatusCode::kInternalError;
			SensitivityRiskCode risktype1 =
			    curve_type == PRICING_CURVE_TYPE_DISCOUNT ? SRC_DISCOUNT : SRC_FORWARD;
			auto curve2 = gammasens->curve2();
			auto curveid2 = curve2->id();
			auto defn_id2 = curves_by_group->get_curve_definition_id(curveid2);
			if (defn_id2 == -1)
				return StatusCode::kInternalError;
			curve_id_components(curveid2, curve_type, ccy, index_family, tenor, as_of_date, cycle, qual,
					    scenario);
			if (scenario != 0)
				return StatusCode::kInternalError;
			SensitivityRiskCode risktype2 =
			    curve_type == PRICING_CURVE_TYPE_DISCOUNT ? SRC_DISCOUNT : SRC_FORWARD;
			auto irsens = result->add_sensitivities();
			irsens->set_curve_definition_id_1(defn_id1);
			irsens->set_curve_definition_id_2(defn_id2);
			irsens->set_order(SensitivityOrderCode::SOC_GAMMA);
			irsens->set_risk_type_1(risktype1);
			irsens->set_risk_type_2(risktype2);
			irsens->set_sensitivity_type(SensitivityTypeCode::STC_ZERO);
			// We start at 1 below as 0 is the as of date and not
			// a pillar
			for (int row = 1; row < gammasens->count1(); row++) {
				for (int col = 1; col < gammasens->count2(); col++) {
					auto v = gammasens->at(row, col);
					if (v != 0.0) {
						// We return offsets that relate to the original curve
						// i.e. pillars start at 0 instead of 1
						// TODO can row=0, col=0 have sensitivity?
						uint32_t return_row = row - 1;
						uint32_t return_col = col - 1;
						uint32_t return_index = (uint32_t)((return_col << 16) | return_row);
						// The gamma values are scaled by 1e-8
						irsens->mutable_values()->insert(
						    google::protobuf::MapPair<google::protobuf::uint32, double>(
							return_index, v * 1e-8));
					}
				}
			}
		}
		return StatusCode::kOk;
	}

      public:
	ValuationServiceImpl(IndexService *index_service, CalendarService *calendar_service)
	    : fixings_data_service_(new FixingDataService()), index_service_(index_service),
	      calendar_service_(calendar_service)
	{
	}
	~ValuationServiceImpl() {}
	CurveInterpolationReply *handle_curve_interpolation_request(Arena *arena,
								    const CurveInterpolationRequest *request)
	{
		auto reply = Arena::CreateMessage<CurveInterpolationReply>(arena);
		ReplyHeader *header = Arena::CreateMessage<ReplyHeader>(arena);
		reply->set_allocated_header(header);
		header->set_response_code(StandardResponseCode::SRC_ERROR);
		header->set_response_message(error_message(StatusCode::kInternalError));

		Date as_of_date = request->business_date();
		if (!is_valid_date(as_of_date)) {
			header->set_response_message(error_message(StatusCode::kVAL_BadBusinessDate));
			return reply;
		}
		StatusCode status = validate(as_of_date, request->curve());
		if (status != StatusCode::kOk) {
			header->set_response_message(error_message(status));
			return reply;
		}
		auto curve = make_curve(as_of_date, &request->definition(), request->curve(), 0);
		PeriodUnit tenor_unit;
		int tenor_mult;
		Period tenor_period;
		Tenor tenor = request->forward_tenor();
		if (tenor == TENOR_UNSPECIFIED)
			tenor = request->definition().tenor();
		if (tenor != TENOR_UNSPECIFIED) {
			get_default_converter()->tenor_to_period_unit_and_multiplier(tenor, &tenor_unit, &tenor_mult);
			tenor_period = Period(tenor_mult, tenor_unit);
		}
		Date last_maturity = curve->last_maturity();
		// TODO we return 0 values to indicate non existant data
		// perhaps we should return NaNs instead?
		for (int i = 0; i < request->dates_size(); i++) {
			Date dt = request->dates(i);
			if (is_valid_date(dt) && dt > as_of_date && dt <= last_maturity) {
				switch (request->rate_type()) {
				case IRRateType::ZERO_RATE:
					reply->add_values(curve->zero_rate(dt));
					break;
				case IRRateType::DISCOUNT_FACTOR:
					reply->add_values(curve->discount(dt));
					break;
				case IRRateType::FORWARD_RATE:
					if (tenor != TENOR_UNSPECIFIED) {
						Date mt = add(dt, tenor_period);
						if (mt <= last_maturity) {
							reply->add_values(curve->forward_rate(dt, mt));
						} else {
							reply->add_values(0);
						}
					} else {
						reply->add_values(0);
					}
					break;
				default:
					reply->add_values(0);
				}
			} else {
				reply->add_values(0);
			}
		}
		header->set_response_code(StandardResponseCode::SRC_OK);
		header->set_response_message("");
		return reply;
	}

	SetCurveMappingsReply *handle_set_curve_mappings_request(Arena *arena,
								 const SetCurveMappingsRequest *request) override final
	{
		SetCurveMappingsReply *response = Arena::CreateMessage<SetCurveMappingsReply>(arena);
		ReplyHeader *header = Arena::CreateMessage<ReplyHeader>(arena);
		response->set_allocated_header(header);
		{
			std::unique_lock<std::shared_mutex> locked(lock_);
			auto &&iter = curve_mappers_.find(request->curve_group());
			SimpleCurveMapper *mapper = nullptr;
			if (iter == curve_mappers_.end()) {
				std::unique_ptr<SimpleCurveMapper> newmapper =
				    std::unique_ptr<SimpleCurveMapper>(new SimpleCurveMapper());
				mapper = newmapper.get();
				curve_mappers_.insert(std::pair<CurveGroup, std::unique_ptr<SimpleCurveMapper>>(
				    request->curve_group(), std::move(newmapper)));
			} else {
				mapper = iter->second.get();
			}

			mapper->clear_mappings();
			for (int i = 0; i < request->mappings_size(); i++) {
				auto const &mapping = request->mappings(i);
				PricingCurve from(mapping.from_id().type(), mapping.from_id().currency(),
						  mapping.from_id().index_family(), mapping.from_id().tenor());
				PricingCurve to(mapping.to_id().type(), mapping.to_id().currency(),
						mapping.to_id().index_family(), mapping.to_id().tenor());
				mapper->add_mapping(from, to);
			}
		}
		header->set_response_code(StandardResponseCode::SRC_OK);
		return response;
	}

	RegisterCurveDefinitionsReply *
	handle_register_curve_definitions_request(Arena *arena,
						  const RegisterCurveDefinitionsRequest *request) override final
	{
		RegisterCurveDefinitionsReply *response = Arena::CreateMessage<RegisterCurveDefinitionsReply>(arena);
		ReplyHeader *header = Arena::CreateMessage<ReplyHeader>(arena);
		response->set_allocated_header(header);
		{
			std::unique_lock<std::shared_mutex> locked(lock_);
			for (int i = 0; i < request->curve_definitions_size(); i++) {
				auto const &defn = request->curve_definitions(i);
				std::unique_ptr<IRCurveDefinition> copy_defn =
				    std::unique_ptr<IRCurveDefinition>(new IRCurveDefinition);
				copy_defn->CopyFrom(defn);
				int id = copy_defn->id();
				if (find_curve_definition(id)) {
					// Remove previous definition
					// fprintf(stderr, "Replacing curve definition %d\n", (int)id);
					curve_definitions_.erase(id);
				}
				curve_definitions_.insert(
				    std::pair<int, std::unique_ptr<IRCurveDefinition>>(id, std::move(copy_defn)));
			}
		}
		header->set_response_code(StandardResponseCode::SRC_OK);
		return response;
	}

	StatusCode check_par_sensitivities(const SetZeroCurvesRequest *request, const ZeroCurve &curve)
	{
		StatusCode status = StatusCode::kOk;
		bool found = false;
		for (int j = 0; j < request->par_sensitivities_size(); j++) {
			auto const &sens = request->par_sensitivities(j);
			if (sens.curve_definition_id() == curve.curve_definition_id()) {
				status = validate(curve, sens);
				found = true;
				break;
			}
		}
		if (!found)
			status = StatusCode::kVAL_MissingParSensitivitiesForCurve;
		return status;
	}

	SetZeroCurvesReply *handle_set_zero_curves_request(Arena *arena,
							   const SetZeroCurvesRequest *request) override final
	{
		SetZeroCurvesReply *response = Arena::CreateMessage<SetZeroCurvesReply>(arena);
		ReplyHeader *header = Arena::CreateMessage<ReplyHeader>(arena);
		response->set_allocated_header(header);
		header->set_response_code(StandardResponseCode::SRC_ERROR);
		StatusCode status = StatusCode::kOk;
		if (!is_valid_date(request->as_of_date())) {
			status = StatusCode::kVAL_BadBusinessDate;
		}
		// Right now we don't enforce the availability of par sensitivities
		// as sometimes when the curves were not bootstrapped by us the
		// sensitivities
		// are not available. However, if par sensitivies are supplied then
		// we check that all supplied curves have them
		for (int i = 0; status == StatusCode::kOk && i < request->forward_curves_size(); i++) {
			auto const &curve = request->forward_curves(i);
			status = validate(request->as_of_date(), curve);
			if (status == StatusCode::kOk && request->par_sensitivities_size() > 0) {
				status = check_par_sensitivities(request, curve);
			}
		}
		for (int i = 0; status == StatusCode::kOk && i < request->discount_curves_size(); i++) {
			auto const &curve = request->discount_curves(i);
			status = validate(request->as_of_date(), curve);
			if (status == StatusCode::kOk && request->par_sensitivities_size() > 0) {
				if (status == StatusCode::kOk && request->par_sensitivities_size() > 0) {
					status = check_par_sensitivities(request, curve);
				}
			}
		}
		if (status != StatusCode::kOk) {
			header->set_response_message(error_message(status));
			return response;
		}
		{
			std::unique_lock<std::shared_mutex> locked(lock_);
			for (int i = 0; status == StatusCode::kOk && i < request->forward_curves_size(); i++) {
				if (!find_curve_definition(request->forward_curves(i).curve_definition_id())) {
					status = StatusCode::kVAL_CurveDefinitionNotRegistered;
					break;
				}
			}
			for (int i = 0; status == StatusCode::kOk && i < request->discount_curves_size(); i++) {
				if (!find_curve_definition(request->discount_curves(i).curve_definition_id())) {
					status = StatusCode::kVAL_CurveDefinitionNotRegistered;
					break;
				}
			}

			if (status == StatusCode::kOk) {
				auto &&iter = curves_by_group_.find(request->curve_group());
				CurvesByGroup *group = nullptr;
				if (iter == curves_by_group_.end()) {
					std::unique_ptr<CurvesByGroup> newgroup =
					    std::unique_ptr<CurvesByGroup>(new CurvesByGroup(request->curve_group()));
					group = newgroup.get();
					curves_by_group_.insert(std::pair<CurveGroup, std::unique_ptr<CurvesByGroup>>(
					    request->curve_group(), std::move(newgroup)));
				} else {
					group = iter->second.get();
				}
				// We only need derivatives when generating base curves
				// Scenario curves are used for revalautions only
				int derive_order = request->scenario() == 0 ? 2 : 0;
				for (int i = 0; i < request->forward_curves_size(); i++) {
					auto const &curve = request->forward_curves(i);
					auto defn =
					    find_curve_definition(request->forward_curves(i).curve_definition_id());
					auto yieldcurve =
					    make_curve(request->as_of_date(), defn, curve, derive_order,
						       PRICING_CURVE_TYPE_FORWARD, request->qualifier(),
						       (short)request->cycle(), (short)request->scenario());
					group->add(std::move(yieldcurve), defn->id());
				}
				for (int i = 0; i < request->discount_curves_size(); i++) {
					auto const &curve = request->discount_curves(i);
					auto defn =
					    find_curve_definition(request->discount_curves(i).curve_definition_id());
					auto yieldcurve =
					    make_curve(request->as_of_date(), defn, curve, derive_order,
						       PRICING_CURVE_TYPE_DISCOUNT, request->qualifier(),
						       (short)request->cycle(), (short)request->scenario());
					group->add(std::move(yieldcurve), defn->id());
				}
				for (int i = 0; i < request->par_sensitivities_size(); i++) {
					auto const &sens = request->par_sensitivities(i);
					// Copy the zero curve par sensitivities to a dense col major matrix
					auto size = sens.num_instruments() * sens.num_maturities();
					void *p = GlobalAllocator.allocate(sizeof(Matrix) + sizeof(double) * size);
					auto matrix = std::unique_ptr<Matrix, Deleter<Matrix>>(
					    new (p) Matrix(sens.num_instruments(), sens.num_maturities()),
					    Deleter<Matrix>(&GlobalAllocator));
					for (auto &iter : sens.values()) {
						uint32_t n = iter.first;
						uint32_t col = n >> 16;
						uint32_t row = n & 0xFFFF;
						// column major matrix to allow matrix multiplication using BLAS
						uint32_t pos = col * sens.num_instruments() + row;
						assert((int)pos < size);
						matrix->data[pos] = iter.second;
					}
					CurveId curve_id = group->get_curve_id(sens.curve_definition_id());
					assert(curve_id); // because we just added the curve above
					debug("Added par sensitivities for curve %s\n",
					      curve_id_to_string(curve_id).c_str());
					// ravi_matrix_t A = { sens.num_instruments(), sens.num_maturities(),
					// matrix.get() };
					// ravi_matrix_dump_to(&A, "", stdout);
					group->add(curve_id, std::move(matrix));
				}
			}
		}
		if (status != StatusCode::kOk) {
			header->set_response_message(error_message(status));
		} else {
			header->set_response_code(StandardResponseCode::SRC_OK);
		}
		return response;
	}

	SetFixingsReply *handle_set_fixings_request(Arena *arena, const SetFixingsRequest *request) override final
	{
		SetFixingsReply *response = Arena::CreateMessage<SetFixingsReply>(arena);
		ReplyHeader *header = Arena::CreateMessage<ReplyHeader>(arena);
		response->set_allocated_header(header);
		header->set_response_code(StandardResponseCode::SRC_ERROR);
		std::vector<Value> values;
		auto &fixings = request->fixings_by_index_tenor();
		InterestRateIndex *index = index_service_->get_index(fixings.index(), fixings.tenor());
		if (index == nullptr) {
			header->set_response_message(error_message(StatusCode::kVAL_UnsupportedIndexOrTenor));
			return response;
		}
		StatusCode status = StatusCode::kOk;
		for (auto iter = std::begin(fixings.fixings());
		     status == StatusCode::kOk && iter != std::end(fixings.fixings()); iter++) {
			if (!is_valid_date(iter->first) /* || !index->is_valid_fixing_date(iter->first) */) {
				status = StatusCode::kVAL_BadFixingDate;
			} else if (iter->second < -0.05 || iter->second > 0.2) {
				status = StatusCode::kVAL_BadFixingValue;
			} else
				values.push_back(Value(iter->first, iter->second));
		}
		if (status != StatusCode::kOk) {
			header->set_response_message(error_message(status));
			return response;
		}
		{
			std::unique_lock<std::shared_mutex> locked(lock_);
			std::unique_ptr<TimeSeries> ts = std::unique_ptr<TimeSeries>(new TimeSeries());
			ts->add(values.size(), values.data());
			fixings_data_service_->set_fixings(index->id(), std::move(ts));
		}
		header->set_response_code(StandardResponseCode::SRC_OK);
		return response;
	}

	/**
	 * Performs valuations of cashflows within a pricing context.
	 */
	ValuationReply *handle_valuation_request(Arena *arena, const ValuationRequest *request) override final
	{
		ValuationReply *response = Arena::CreateMessage<ValuationReply>(arena);
		ReplyHeader *header = Arena::CreateMessage<ReplyHeader>(arena);
		response->set_allocated_header(header);
		header->set_response_code(StandardResponseCode::SRC_ERROR);
		if (!request->has_cashflows()) {
			header->set_response_message(error_message(StatusCode::kVAL_MissingCashflow));
			return response;
		}
		if (!request->has_pricing_context()) {
			header->set_response_message(error_message(StatusCode::kVAL_MissingPricingContext));
			return response;
		}
		if (!is_valid_date(request->pricing_context().as_of_date())) {
			header->set_response_message(error_message(StatusCode::kVAL_BadBusinessDate));
			return response;
		}
		{
			std::shared_lock<std::shared_mutex> locked(lock_);
			CurvesByGroup *curves_by_group = get_group(request->pricing_context().curve_group());
			if (curves_by_group == nullptr) {
				header->set_response_message(error_message(StatusCode::kVAL_CurveGroupNotFound));
				return response;
			}
			CurveMapper *mapper = get_mapper(request->pricing_context().curve_group());
			if (mapper == nullptr) {
				header->set_response_message(
				    error_message(StatusCode::kVAL_CurveMappingsForCurveGroupNotFound));
				return response;
			}

			// We would like to ensure that all the computations done
			// here use tightly controlled memory. We would like to meet several
			// goals.
			// a) Reuse memory allocated across requests - i.e. avoid repeated
			//   malloc() and free() calls. This means that some form of arena allocator
			//   is needed that is assigned to each request rather than created anew for
			//   each request.
			// b) Ensure that the allocated memory has locality - i.e. items that
			//   should be together in memory are indeed close together.
			// c) Taking the locality idea further even the curves should ideally be
			//   close to the cashflow data.
			// d) Next we want a solution that avoid repeated allocations in the leaf
			//   functions. Example of this was the multiply and divide functions in
			//   the autodiff implementation. These require temporary buffer space every time they
			//   are called.
			//   I recently (April 2017) changed this so that these routines now
			//   require caller to provide temp space so that the memory allocation concern is pushed
			//   down the chain. In general, ideally we want all the memory to be
			//   pre-allocated before any valuations are done, this require precomputing the memory
			//   required.
			// e) A detail of this is that we first create cashflow structures and
			//   then carry out valuations. Ideally these should be located close to each
			//   other in memory.
			// f) The allocator design should be such that we can allocate a bunch
			//   of memory for cashflow structures - and then alloc/free the memory required
			//   for valuations each time a valuation is done - reusing the same memory
			//   each time. That is the memory is organized so that first chunk is the cashflow
			//   structure, which is allocated once; this is followed by valuation data which
			//   grows and shrinks with each valuation but reuses the same memory.
			AllocatorSet *allocSet = get_threadspecific_allocators();
			allocSet->reset();
			// TODO this needs to be more speific
			int derivatives = request->pricing_context().from_scenario() == 0 ? 2 : 0;
			ValuationContextImpl ctx(request->pricing_context().as_of_date(), derivatives,
						 fixings_data_service_.get(),
						 request->pricing_context().payment_cutoff_date(),
						 request->pricing_context().is_todays_fixings_included());
			auto cashflows =
			    construct_cashflows(allocSet->cashflow_allocator, &request->cashflows(), ctx, mapper);
			if (!cashflows) {
				header->set_response_message(error_message(StatusCode::kVAL_CashflowConversionFailed));
				return response;
			}
			StatusCode status = StatusCode::kOk;
			Sensitivities sens(allocSet->sensitivities_allocator);
			auto result = Arena::CreateMessage<ValuationResult>(arena);
			response->set_allocated_result(result);
			for (int scenario = request->pricing_context().from_scenario();
			     status == StatusCode::kOk && scenario <= request->pricing_context().to_scenario();
			     scenario++) {
				CurveProviderImpl curve_provider(
				    request->pricing_context().as_of_date(), request->pricing_context().qualifier(),
				    request->pricing_context().cycle(), curves_by_group, scenario);
				int derivative_order = scenario == 0 ? 2 : 0;
				ctx.set_derivative_order(derivative_order);
				sens.reset();
				auto value = compute_present_value(allocSet->tempspace_allocator, ctx, cashflows,
								   &curve_provider, sens, status);
				if (status != StatusCode::kOk) {
					header->set_response_message(error_message(status));
				} else {
					if (scenario == 0) {
						add_irsensitivities(arena, curves_by_group, result, sens);
					}
					result->mutable_valuations()->insert(
					    MapPair<google::protobuf::int32, double>(scenario, value));
				}
			}
			if (status == StatusCode::kOk) {
				header->set_response_code(StandardResponseCode::SRC_OK);
			}
		}
		return response;
	}
};

std::unique_ptr<ValuationService> get_valuation_service(IndexService *index_service, CalendarService *calendar_service)
{
	return std::unique_ptr<ValuationServiceImpl>(new ValuationServiceImpl(index_service, calendar_service));
}
}
