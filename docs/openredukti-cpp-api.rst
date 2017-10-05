.. highlightlang:: cpp

===================
OpenRedukti C++ API
===================

This document described the C++ api that is available for programmers.

Basics
======
The OpenRedukti C++ api is designed to be fairly simple to use. While OpenRedukti uses C++ templates internally, these are not exposed at 
a user level. The api is presented as a simple set of modular classes or functions.

OpenRedukti uses a bunch of data types defined using Google Protocol Buffers. These reduce the need to manually maintain data types
and greatly improve the ability to make rapid changes to the codebase.

Namespace
=========
All components are put in the namespace ``redukti``.

Common Enums
============

::

   #include <enums.pb.h>

For reasons of efficiency all internal data structures use integer codes rather than 
strings for values. Most of the enums used are defined in the protocol buffers definition file 
`enums.proto <https://github.com/redukti/OpenRedukti/blob/master/proto/enums.proto>`_. 

redukti::Currency   
   Defines currency codes

redukti::IsdaIndex
   Defines commonly used ISDA Index names

redukti::IndexFamily
   Defines families of indices

redukti::DayCountFraction
   Defines supported ISDA Day Count Fractions

redukti::CompoundingMethod
   Defines compounding methods; compatible with FpML.

redukti::BusinessCenter
   Defines commonly used ISDA codes for business centers. These are used to derive holiday calendars.

redukti::BusinessDayConvention
   Defines ISDA specified business day conventions.

redukti::Tenor
   Defines a tenor period such as 1M or 1Y. The codes are designed so that larger terms have higher codes.

redukti::PeriodUnit
   Defines the unit in which a period is measured, e.g. Days or Months, etc.

redukti::RollConvention
   Defines ISDA specified roll conventions for calculating periods for interest rate streams.

redukti::JointCalendarRule
   Defines how multiple business centers are to be combined for the purposes of computing holidays.

redukti::InterpolatorType
   Defines the supported interpolation methods.

redukti::CurveGroup
   Defines curve group ids for the purposes of grouping curves used in pricing.

redukti::IRRateType
   Defines the type of value being used in a curve definition.

redukti::PricingCurveType
   Defines the usage of a curve in a pricing scenario.

redukti::MarketDataQualifier
   Defines the type of Market Data being used.

redukti::MaturityGenerationRule
   Defines the rule name for generating maturities of instruments in a curve.

Date types
==========

::

   #include <date.h>

Dates are used extensively in OpenRedukti. To make it efficient to manipulate dates an integer representation of a Date is
chosen. The date library is based upon:

* `Date Algorithms by Howard Hinnant <http://howardhinnant.github.io/date_algorithms.html>`_.
* The date implementation in QuantLib.

The following key types are defined:

::

   // This class provides a Period (length + TimeUnit) class
   // and implements a limited algebra.
   // Must be standard layout for C compatibility
   class Period
   {
   public:
      // Construct a period from length and unit
      Period(int n, PeriodUnit unit) noexcept;

      // Default constructor : 0D period
      Period() noexcept;

      int length() const noexcept;
      PeriodUnit units() const noexcept;

      bool operator==(const Period &p) const noexcept;

      // Normalisation converts weeks to days and
      // years to months
      Period normalised() const noexcept;

      // Converts a tenor to period representation
      // Must be updated if definition of Tenor changes.
      static Period tenor_to_period(Tenor tenor);
   };

::

   enum Weekday {
      Sunday = 0,
      Monday = 1,
      Tuesday = 2,
      Wednesday = 3,
      Thursday = 4,
      Friday = 5,
      Saturday = 6,
      Sun = 0,
      Mon = 1,
      Tue = 2,
      Wed = 3,
      Thu = 4,
      Fri = 5,
      Sat = 6
   };

::

   // Month names
   enum Month {
      January = 1,
      February = 2,
      March = 3,
      April = 4,
      May = 5,
      June = 6,
      July = 7,
      August = 8,
      September = 9,
      October = 10,
      November = 11,
      December = 12,
      Jan = 1,
      Feb = 2,
      Mar = 3,
      Apr = 4,
      Jun = 6,
      Jul = 7,
      Aug = 8,
      Sep = 9,
      Oct = 10,
      Nov = 11,
      Dec = 12
   };

::

   // Date type. Uses an int to
   // represent a serial number.
   // this implementation is immutable - hence
   // thread-safe.
   typedef int32_t Date;

   struct YearMonthDay {
      short y;
      unsigned char m;
      unsigned char d;
   };

::

   // Returns number of days since civil 1899-12-31.  Negative values indicate
   //    days prior to 1899-12-31.
   // Preconditions:  y-m-d represents a date in the civil (Gregorian) calendar
   //                 m is in [1, 12]
   //                 d is in [1, last_day_of_month(y, m)]
   //                 y is "approximately" in
   //                   [numeric_limits<Int>::min()/366,
   //                   numeric_limits<Int>::max()/366]
   //                 Exact range of validity is:
   //                 [civil_from_days(numeric_limits<Int>::min()),
   //                  civil_from_days(numeric_limits<Int>::max()-719468+25569)]
   // Notes: The original algorithm has been modified to make
   // the serial date match Excel dates. This is done by making the start
   // date 31/Dec/1899 rather than 1/Jan/1970.
   constexpr Date make_date(unsigned d, unsigned m, int y) noexcept;

   constexpr Date make_date(YearMonthDay ymd);

::

   // Returns year/month/day triple in civil calendar
   // Preconditions:  z is number of days since 1899-12-31 and is in the range:
   //                   [numeric_limits<Int>::min(),
   //                   numeric_limits<Int>::max()-719468+25569].
   // Notes: The original algorithm has been modified to make
   // the serial date match Excel dates. This is done by making the start
   // date 31/Dec/1899 rather than 1/Jan/1970.
   constexpr YearMonthDay date_components(Date z);

   // Day of the year, where Jan 1 is 1, Jan 2 is 2, Feb 1 is 32 and so on.
   constexpr int day_of_year(YearMonthDay ymd);

   // Returns day of week in civil calendar [0, 6] -> [Sun, Sat]
   // Preconditions:  z is number of days since 1899-12-31 and is in the range:
   //                   [numeric_limits<Int>::min(), numeric_limits<Int>::max()-4].
   // Notes: The original algorithm has been modified to make
   // the serial date match Excel dates. This is done by making the start
   // date 31/Dec/1899 rather than 1/Jan/1970.
   constexpr unsigned char weekday(Date z) noexcept;

   // Preconditions: m is in [1, 12]
   // Returns: The number of days in the month m of common year
   // The result is always in the range [28, 31].
   constexpr unsigned last_day_of_month_common_year(unsigned m) noexcept;

   // Returns: true if y is a leap year in the civil calendar, else false
   constexpr bool is_leap(int y) noexcept;

   // Preconditions: m is in [1, 12]
   // Returns: The number of days in the month m of year y
   // The result is always in the range [28, 31].
   constexpr unsigned last_day_of_month(int y, unsigned m) noexcept;

   // Add/subtract periods from dates
   extern Date add(Date date, const Period &) noexcept;
   extern Date sub(Date date, const Period &) noexcept;

   // Construct an end of month date for the
   // given year and month
   constexpr Date end_of_month(int y, unsigned m) noexcept;

   // Test whether given date is the calendar end of the month
   constexpr bool is_end_of_month(YearMonthDay ymd) noexcept;

   // Preconditions: x <= 6 && y <= 6
   // Returns: The number of days from the weekday y to the weekday x.
   // The result is always in the range [0, 6].
   constexpr unsigned weekday_difference(unsigned x, unsigned y) noexcept;

   // Preconditions: wd <= 6
   // Returns: The weekday following wd
   // The result is always in the range [0, 6].
   constexpr unsigned next_weekday(unsigned wd) noexcept;

   // Preconditions: wd <= 6
   // Returns: The weekday prior to wd
   // The result is always in the range [0, 6].
   inline constexpr unsigned prev_weekday(unsigned wd) noexcept;

   // next given weekday following or equal to the given date
   // E.g., the Friday following Tuesday, January 15th, 2002
   //   was January 18th, 2002.
   // see also http://www.cpearson.com/excel/DateTimeWS.htm
   constexpr Date next_weekday(Date d, Weekday desired_weekday) noexcept;

   // n-th given weekday in the given month and year
   // E.g., the 4th Thursday of March, 1998 was March 26th,
   // 1998.
   YearMonthDay nth_weekday(unsigned n, unsigned wd, unsigned month, int year);

   constexpr bool is_weekend(unsigned wd);

   // Min allowed date is Jan 1st 1901
   // This is imposed by OpenRedukti
   // This is helpful because then 0 can be used to represent an invalid date
   constexpr Date minimum_date() noexcept;

   // We limit the max date so that we can ensure date values
   // fit in 24 bits
   // Dec 31st, 2199 
   constexpr Date maximum_date() noexcept;

   // Parse a date
   // Returns true on success
   bool parse_date(const char *s, Date *d) noexcept;

   // We need to ensure that 0 is not a valid date as this
   // helps us with protobuf representation of dates as integers
   // where unspecified value is 0.
   // Another requirement is to limit the max date so that
   // date values can fit into 24 bits.
   bool is_valid_date(Date date) noexcept;


Holiday Calendars
=================

::

   #include <calendars.h>

OpenRedukti comes with predefined calendar implementations for following Business Centers:

* ``AUSY``
* ``EUTA``
* ``GBLO``
* ``USNY``
* ``JPTO``
* ``BRSP``

These implementations are derived from the QuantLib library.

The Calendar api is as described below.

::

   // The Calendar interface provides the means to determine whether
   // a given date is a holiday for a business center or not. Also
   // the interface provides methods for adjusting dates as per the 
   // holiday calendar.
   // Immutable for thread safety.
   class Calendar
   {
   public:
      virtual ~Calendar() noexcept;
      
      virtual int id() const noexcept = 0;

      // Returns all the ids - relevant for calendars made by combining
      // others
      virtual void get_ids(std::array<BusinessCenter, 4> &ids) const noexcept;

      virtual bool is_holiday(Date d) const noexcept = 0;

      bool is_businessday(Date d) const noexcept;

      bool is_end_of_month(Date d) const noexcept;

      // Adjust the given date to be the last business day of the month
      Date end_of_month(Date d) const noexcept;

      // Adjusts a non-business day to the appropriate near business day
      //  with respect to the given convention.
      Date adjust(Date date, BusinessDayConvention convention = BusinessDayConvention::FOLLOWING) const noexcept;

      //  Advances the given date of the given number of business days and
      //  returns the result. Note that if unit is Days then business day
      // convention and eom flags are not used as the date is move by the
      // specified business days. For other period units the date is moved as
      // per raw calendar and then adjusted if it falls on a holiday
      Date advance(Date date, int n, PeriodUnit unit,
              BusinessDayConvention convention = BusinessDayConvention::FOLLOWING, bool endOfMonth = false) const
          noexcept;

      //  Advances the given date as specified by the given period and
      //  returns the result.
      //  The input date is not modified.
      Date advance(Date date, const Period &period,
              BusinessDayConvention convention = BusinessDayConvention::FOLLOWING, bool endOfMonth = false) const
          noexcept;

      // Calculates the number of business days between two given
      // dates and returns the result.
      //
      int business_days_between(Date from, Date to, bool includeFirst = true, bool includeLast = false) const
          noexcept;
   };

   struct JointCalendarParameters {
      std::array<BusinessCenter, 4> centers;
      JointCalendarParameters(BusinessCenter center1, BusinessCenter center2,
               BusinessCenter center3 = BusinessCenter::BUSINESS_CENTER_UNSPECIFIED,
               BusinessCenter center4 = BusinessCenter::BUSINESS_CENTER_UNSPECIFIED);
   };

   // The Calendar Service manages calendar instances. It has to meet following requirements:
   // a) It must always return the same Calendar instance for a given business center. Clients
   //    can assume that the instance will not go away or change in any way as long as the
   //    service is live.
   // b) Ditto for joint calendar instances.
   // c) Calendar instances must be immutable.
   class CalendarService
   {
   public:
      virtual ~CalendarService() {}
      // Return the calendar specified. Memory is managed by the
      // CalendarFactory so the caller must not delete.
      virtual const Calendar *get_calendar(BusinessCenter id) noexcept = 0;

      // Set a calendar to given instance.
      // The service will take ownership of the instance
      // May fail if calendar instance already set and has been
      // accessed by a client - i.e. new calendars can only be set prior to
      // any use.
      virtual bool set_calendar(BusinessCenter id, std::unique_ptr<Calendar> calendar) noexcept = 0;

      // Create joint calendar
      // Note that the order in which the business centers are given
      // should not matter - i.e. the constituents must be sorted and then
      // combined so that for a given combination the returned instance is
      // always the same
      virtual Calendar *get_calendar(JointCalendarParameters calendars,
                      JointCalendarRule rule = JointCalendarRule::JOIN_HOLIDAYS) noexcept = 0;
   };

   // Get the calendar factory
   extern CalendarService *get_calendar_factory() noexcept;

   // Utility for constructing a joint calendar
   extern const Calendar *build_calendar(CalendarService *calendar_service,
                     const google::protobuf::RepeatedField<google::protobuf::int32> &values,
                     JointCalendarRule rule = JointCalendarRule::JOIN_HOLIDAYS);

   // Utility for constructing a joint calendar
   const Calendar *build_calendar(CalendarService *calendar_service, const std::vector<BusinessCenter> &values,
                   JointCalendarRule rule = JointCalendarRule::JOIN_HOLIDAYS);


Day Count Fractions
===================

::

   #include <dayfractions.h>


OpenRedukti comes with support for following DayCountFraction implementations:

* ``30/360``
* ``30E/360``
* ``30E/360.ISDA``
* ``ACT/360``
* ``ACT/365.FIXED``
* ``ACT/ACT.ISDA``
* ``ACT/ACT.ISMA``
* ``BUS/252``

The implementation is derived from QuantLib.

::

   // Compute the difference between dates as per Day Count Convention.
   // The difference is measured in factional units of a year, where one year 1.0.
   // Must be immutable and thread-safe.
   // Clients must be able to hold references to these for the lifetime of
   // the application.
   class DayFraction
   {
   public:
      virtual ~DayFraction() {}

      // Calculate the difference d2-d2 as per convention
      // for the DayFraction; value is a decimal expressed as a year fraction.
      // So 1.0 means 1 year.
      virtual double year_fraction(Date d1, Date d2) const = 0;

      // Only used for ThirtyE360ISDA (30E/360.ISDA)
      // The finalPeriod flag indicates whether this fraction is for the
      // final period - i.e. d2 is maturity date. So typically,
      // when calculating the last calc period in a swap, this flag must be
      // set to true.
      virtual double year_fraction(Date d1, Date d2, bool finalPeriod) const = 0;

      // Used only for ACT/ACT.ISMA
      // refStart - If regular period or front stub then adjusted end date
      //    minus calculation period frequency (roll convention NONE),
      //    else adjusted start date
      // refEnd - If regular period or front stub then adjusted end date,
      //    else adjusted start date minus calculation period
      //    frequency (roll convention NONE)
      virtual double year_fraction(Date d1, Date d2, Date refStart, Date refEnd) const = 0;

      // Returns the ISDA name
      virtual DayCountFraction id() const = 0;
   };

   // Get an instance of a DayFraction
   // Requirements:
   // a) There must only be one instance associated with a particular DayCountFraction
   // b) The DayCountFraction implementation must be immutable and hence thread-safe
   extern const DayFraction *get_day_fraction(DayCountFraction dfc);

   // The BUS252 day fraction requires a calendar.
   // Requirements:
   // a) There must only be one instance associated with a particular DayCountFraction
   // b) The DayCountFraction implementation must be immutable and hence thread-safe
   extern const DayFraction *get_bus_252(CalendarService *calendarService, BusinessCenter center);

Index Types
===========

::

   #include <index.h>

The ``IndexDefinition`` type captures essential information for working with indices.

::

   // Captures information about an interest rate 
   // index so that various operations associated with the
   // index can be performed.
   message IndexDefinition {
      // This is the FpML / ISDA name of the index
      // Excludes tenor
      IsdaIndex isda_index = 1;
      // Index family may be common across several indices
      IndexFamily index_family = 2;
      // Currency of the index
      Currency currency = 3;
      // Tenor of the index
      // May be left unspecified to use as the default configuration
      // for all tenors for the index family
      Tenor tenor = 4; 
      // Number of business days between a value date and fixing date
      int32 fixing_lag = 5;
      // Used to select the business day convention
      // tenors <= short_tenor_threshold use the short tenor convention
      Tenor short_tenor_threshold = 6;
      // Convention used if tenor is <= short_tenor_threshold
      BusinessDayConvention short_tenor_convention = 7;
      // Convention used if tenor is > than short_tenor_threshold
      BusinessDayConvention long_tenor_convention = 8;
      // Whether to apply EOM roll convention for tenors >= month
      bool eom = 9;
      // fixing calendars are used to move from accrual start date
      // to fixing date, and also to move from fixing date to
      // value date
      repeated BusinessCenter fixing_calendars = 10;
      // How to combine fixing calendars
      JointCalendarRule fixing_calendars_join_rule = 11;
      // The value date is checked against the value date
      // calendars and if it falls on a holiday then
      // an adjustment is applied using business day convention 
      repeated BusinessCenter value_date_calendars = 12;
      // How value date calendars are to be combined
      JointCalendarRule value_date_calendars_join_rule = 13;
      // The index calendars are used to calculate the maturity date
      repeated BusinessCenter index_calendars = 14;
      // How index date calendars are to be combined
      JointCalendarRule index_calendars_join_rule = 15;
      // The day count fraction associated with the index
      DayCountFraction day_count_fraction = 16;
      // Is this Isda Index default for the currency and index family?
      bool default_for_index_family = 17;
   }


Here is an example of how this is defined in C++ code for ``USD LIBOR 1D`` index:

::

   IndexDefinition definition;
   definition.set_isda_index(IsdaIndex::USD_LIBOR_BBA);
   definition.set_index_family(IndexFamily::LIBOR);
   definition.set_currency(Currency::USD);
   definition.set_tenor(TENOR_1D);
   definition.set_fixing_lag(0);
   definition.set_short_tenor_threshold(TENOR_2W);
   definition.set_short_tenor_convention(BusinessDayConvention::FOLLOWING);
   definition.set_long_tenor_convention(BusinessDayConvention::MODIFIED_FOLLOWING);
   definition.add_fixing_calendars(GBLO);
   definition.set_fixing_calendars_join_rule(JOIN_HOLIDAYS);
   definition.add_value_date_calendars(GBLO);
   definition.add_value_date_calendars(USNY);
   definition.set_value_date_calendars_join_rule(JOIN_HOLIDAYS);   
   definition.add_index_calendars(GBLO);
   definition.add_index_calendars(USNY);
   definition.set_index_calendars_join_rule(JOIN_HOLIDAYS);
   definition.set_day_count_fraction(DayCountFraction::ACT_360);   
   definition.set_default_for_index_family(true);
   definition.set_eom(false);

To support other tenors, one can simply take above and change folloowing:

::

   definition.set_tenor(TENOR_UNSPECIFIED);
   definition.set_fixing_lag(2);
   definition.set_eom(true);


The ``IndexDefinition`` acts as a template for creating instances of the ``InterestRateIndex`` type. 

The C++ api for working with indices is given below::

   // Unique identifier for an index 
   typedef uint32_t IndexId;

   // Makes a unique identifier from the give ISDA index identifier and
   // tenor
   IndexId make_index_id(IsdaIndex isda_index, Tenor tenor);

   class IndexDefinition;

   // Base type for all indices
   class Index
   {
   public:
      virtual ~Index() {}
      virtual IndexId id() const = 0;
   };

   // An interest rate index representation. A requirement of 
   // OpenRedukti is that an each unique IndexId should map to one
   // InterestRateIndex instance - as this allows the code to freely 
   // reference such instances without fear of the reference going away.
   // Additionally a requirement is that the instance is immutable.
   class InterestRateIndex : public Index
   {
   public:
      virtual ~InterestRateIndex();
      virtual Currency currency() const;
      virtual IndexFamily family() const;
      virtual Tenor tenor() const;
      virtual IsdaIndex isda_index() const;

      // Given a fixing date, calculate the value date 
      // by applying the calendars, day conventions associated
      // with the index
      virtual Date value_date(Date fixing_date) const;

      // Given a value date, calculate the fixing date 
      // by applying the calendars, day conventions associated
      // with the index
      virtual Date fixing_date(Date accrual_start_date) const;

      // Given a value date calculate the maturity date
      // Appropriate calendars, day conventions and EOM rules
      // must be applied
      virtual Date maturity_date(Date value_date) const;
      virtual bool is_valid_fixing_date(Date date) const;
      virtual const Calendar *fixing_calendar() const;
      virtual const DayFraction *day_fraction() const;
      virtual BusinessDayConvention day_convention() const;
   };

   // The IndexService is responsible for returning instances of InterestRateIndex.
   // Note that the index service must ensure the following:
   // a) There will only ever be one instance of an InterestRateIndex for a given
   //    IndexId.
   // b) Clients must be free to hold on to references to such instances without
   //    fear of them going out of scope. So essentially these instances can only be
   //    deleted at system shutdown.
   // c) An InterestRateIndex instance must be immutable.
   class IndexService
   {
   public:
      virtual ~IndexService() {}

      // Adds a definition for use as a template for generating instances of
      // InterestRateIndex
      virtual bool register_index(const IndexDefinition &definition) = 0;

      // Obtains an instance of IntrestRateIndex - must return an existing instance
      // if already defined 
      virtual InterestRateIndex *get_index(IsdaIndex isda_index, Tenor tenor) = 0;

      // Obtains an instance of IntrestRateIndex - must return an existing instance
      // if already defined 
      virtual InterestRateIndex *get_index(Currency currency, IndexFamily index_family, Tenor tenor) = 0;
   };

   extern IndexService *get_default_index_service();

Useful Conversions
==================

::

   #include <converters.h>


The api is as follows::

   class Converter
   {
   public:
      virtual ~Converter() {}
      virtual BusinessCenter business_center_from_string(const char *value) const;
      virtual BusinessDayConvention business_day_convention_from_string(const char *s) const;
      virtual PeriodUnit period_unit_from_string(const char *s) const;
      virtual bool period_from_string(const char *periodName, Period *p) const;
      virtual DayCountFraction day_count_fraction_from_string(const char *value) const;
      virtual Tenor tenor_from_period_unit_and_len(PeriodUnit unit, int value) const;
      virtual bool tenor_to_period_unit_and_multiplier(Tenor value, PeriodUnit *unit, int *mult) const;
      virtual std::string tenor_to_string(Tenor tenor) const;
      virtual RollConvention roll_convention_from_string(const char *s) const;
      virtual Currency currency_from_string(const char *s) const;
      virtual const char *currency_to_string(Currency value) const;
      virtual IsdaIndex isda_index_from_string(const char *s) const;
      virtual const char *isda_index_to_string(IsdaIndex value) const;
      virtual CompoundingMethod compounding_method_from_string(const char *value) const;
      virtual IndexFamily index_family_from_string(const char *value) const;
      virtual const char *index_family_to_string(IndexFamily value) const;
      virtual const char *period_unit_to_string(PeriodUnit period_unit) const;
      virtual int tenor_to_days(Tenor tenor) const;
      virtual InterpolatorType interpolator_type_from_string(const char *s) const;
      virtual PricingCurveType pricing_curve_type_from_string(const char *s) const;
      virtual IRRateType rate_type_from_string(const char *s) const;
      virtual CurveGroup curve_group_from_string(const char *value) const;
      virtual MaturityGenerationRule maturity_generation_rule_from_string(const char *value) const;
   };

   extern const Converter *get_default_converter();

Automatic Differentiation
=========================

::

   #include <autodiff.h>

OpenRedukti makes use of automatic differentiation techniques to compute derivatives. This approach enables 
computation of derivatives more accurately than would be possible using numeric differentation. On the other
hand, it is possible to implement fairly complex derivatives without having to construct the derivatives by
hand.

This approach does have the drawback that it is compute and memory intensive. Hence to improve performance 
special care is taken with regards to memory management.

The implementation of the adouble type is based on implementation of a vector-mode hyper-dual numbers
written by: Jeffrey A. Fike at Stanford University, Department of Aeronautics and Astronautics.

The core API is as follows::

   // WARNING
   //
   // This is a low level module that must be used with care.
   // In general this module requires the caller to allocate memory
   // correctly - as it assumes that all supplied arguments are
   // properly sized and allocated.

   /* autodiff variable */
   struct redukti_adouble_t {
      // derivative order
      uint32_t order_ : 2;
      // number of variables 
      uint32_t vars_;
      // data 
      double data_[1];

      redukti_adouble_t(const redukti_adouble_t &) = delete;
      redukti_adouble_t &operator=(const redukti_adouble_t &) = delete;
   };

   // Compute memory requirement for given number of variables and order
   // Supported orders are 0,1,2.
   size_t redukti_adouble_alloc_size(int vars, int order);

   // Initialize A; caller must have allocated memory of correct
   // size.
   void redukti_adouble_init(redukti_adouble_t *A, int n_vars, int order, int var, double v);

   // A = B
   // must be same size
   void redukti_adouble_assign(redukti_adouble_t *A, const redukti_adouble_t *B);

   // A = A + alpha*B
   void redukti_adouble_add(redukti_adouble_t *A, redukti_adouble_t *B, double alpha);

   // A = A*scalar
   void redukti_adouble_scalar_multiply(redukti_adouble_t *A, double alpha);

   // A = A*B
   // A = A*A also works
   // temp must be same size as A
   void redukti_adouble_multiply(redukti_adouble_t *A, redukti_adouble_t *B, redukti_adouble_t *temp);

   // A = A/B 
   // temp1, temp2 must be same size as A
   void redukti_adouble_divide(redukti_adouble_t *A, redukti_adouble_t *B, redukti_adouble_t *temp1,
                redukti_adouble_t *temp2);

   // A = exp(A)
   // temp must be same size as A
   void redukti_adouble_exp(redukti_adouble_t *A, redukti_adouble_t *temp);

   // A = log(A)
   // temp must be same size as A
   void redukti_adouble_log(redukti_adouble_t *A, redukti_adouble_t *temp);

   // A = A^p
   // temp must be same size as A
   void redukti_adouble_power(redukti_adouble_t *A, double p, redukti_adouble_t *temp);

   // A = abs(A)
   void redukti_adouble_abs(redukti_adouble_t *A);

   // A = sin(A)
   // temp must be same size as A
   void redukti_adouble_sin(redukti_adouble_t *A, redukti_adouble_t *temp);

   // A = cos(A)
   // temp must be same size as A
   void redukti_adouble_cos(redukti_adouble_t *A, redukti_adouble_t *temp);

   // A = tan(A)
   // temp must be same size as A
   void redukti_adouble_tan(redukti_adouble_t *A, redukti_adouble_t *temp);

   // Dumps contents of A
   void redukti_adouble_dump(redukti_adouble_t *v, FILE *out, const char *desc);

   // A = A + alpha 
   void redukti_adouble_scalar_add(redukti_adouble_t *A, double alpha);

   // Get A's value
   double redukti_adouble_get_value(redukti_adouble_t *A);

   // Get first derivative with respect to variable 
   double redukti_adouble_get_derivative1(redukti_adouble_t *A, int parameter);

   // Get second derivative with respect to variables
   double redukti_adouble_get_derivative2(redukti_adouble_t *A, int parameter1, int parameter2);

   // Set A's value
   void redukti_adouble_set_value(redukti_adouble_t *A, double v);

   // Set first derivative with respect to variable 
   void redukti_adouble_set_derivative1(redukti_adouble_t *A, int parameter, double v);

   // Set second derivative with respect to variables
   void redukti_adouble_set_derivative2(redukti_adouble_t *A, int parameter1, int parameter2, double v);


It is best to try to use the scriting API to get an understanding of how automatic differentiation works.

Calculation Schedules
=====================

::

   #include <schedule.h>

The following protocol buffer definitions capture parameters for creating calculation scheduled::

   enum StubLocation {
      STUB_TYPE_AUTO = 0;
      SHORT_FRONT_STUB = 1;
      LONG_FRONT_STUB = 2;
      SHORT_BACK_STUB = 3;
      LONG_BACK_STUB = 4;
   }

   message ScheduleParameters {
      int32 effective_date = 1;
      int32 termination_date = 2;
      RollConvention roll_convention = 3;
      int32 first_regular_period_start_date = 4;
      int32 last_regular_period_end_date = 5;
      Tenor term = 6;
      Tenor calculation_frequency = 7;
      Tenor payment_frequency = 8;
      StubLocation stub_location = 9;
      int32 first_payment_date = 10;
      int32 last_regular_payment_date = 11;
      BusinessDayConvention period_convention = 12;
      BusinessDayConvention payment_convention = 13;
      sint32 payment_lag = 14;
      repeated BusinessCenter period_calendars = 15;
      repeated BusinessCenter payment_calendars = 16;
   }

   message Schedule {
      bool has_front_stub = 1;
      bool has_back_stub = 2;
      repeated int32 adjusted_start_dates = 3;
      repeated int32 adjusted_end_dates = 4;
      // payment date may be zero if not applicable for a period
      repeated int32 adjusted_payment_dates = 5; 
   }


The C++ api to generate a schedule from given parameters is as defined below::

   class ScheduleParameters;
   class Schedule;

   // Build a schedule as per the schedule parameters
   // If succesful returns true
   extern bool build_schedule(ScheduleParameters &params, Schedule &schedule) noexcept;

   // Adjusts a date as per roll convention specified
   extern Date adjust_date(Date d, RollConvention rc) noexcept;


Memory Allocators
=================

::

   #include <allocators.h>

OpenRedukti uses a bunch of allocators that aim to reduce the overhead in allocating and releasing memory.
The general interface implemented by all the allocators are::

   // IMPORTANT
   //
   // The allocators defined below are NOT thread safe
   // You must ensure that an allocator (other than the
   // MallocAllocator to be accurate) is never shared across
   // threads
   //
   // Secondly these allocators are fine tuned to requirements
   // in this project and are not general purpose.

   // Generic allocator interface
   class Allocator
   {
   public:
      virtual ~Allocator() noexcept;

      // Allocate at least size bytes
      // A size of 0 will result in nullptr being returned
      virtual void *allocate(size_t size) noexcept;

      void *safe_allocate(size_t size) noexcept;

      // Depending upon the type of allocator a deallocate may
      // not do anything
      virtual void deallocate(void *address) noexcept;
   };

When objects are allocated and then captured in std::unique_ptr, it is necessary to provide a deleter object to 
correctly deallocate memory. For this, the following is provided::

   // Utility for associating a deleter with a
   // unique_ptr when memory was allocated using an allocator.
   //
   // Example:
   //  Allocator *A;
   //  std::unique_ptr<YieldCurve, Deleter<YieldCurve>>(
   // new (*A) YieldCurve(), Deleter<YieldCurve>(A));
   //
   template <typename T> class Deleter
   {
   public:
      Deleter(Allocator *A = nullptr) : A_(A);
      void operator()(T *p);
   };


An extension of the Allocator interface provides allocation strategies where all memory is released at once rather than
object at a time.::

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
      virtual ~RegionAllocator() noexcept;

      virtual void *allocate(size_t size) noexcept;

      // Deallocate does nothing
      void deallocate(void *address) noexcept override final {}

      // Resets the allocator so that all memory
      // is either freed and available for reuse
      virtual void release() noexcept;
   };


We have a FixedRegionAllocator that allocates from a predefined memory buffer.::

   // This is an allocator that returns memory from a fixed
   // sized memory buffer. The buffer may be externally provided or
   // owned. When the buffer is exhausted any allocation requests
   // will fail and allocate() will return nullptr.
   //
   // As it is a RegionAllocator, deallocate() is a no-op
   struct FixedRegionAllocator : public RegionAllocator {

      // memory externally supplied
      FixedRegionAllocator(char *start, char *end) noexcept;

      // memory externally supplied
      FixedRegionAllocator(void *start, size_t n) noexcept;

      // Acquire memory
      // Memory will be owned by this instance
      FixedRegionAllocator(size_t n) noexcept;

      // Current position
      size_t pos() const noexcept;

      // Sets current position
      // This is useful for scenarios where the user
      // wants to use the allocator in a stack like fashion
      // This is used by FixedRegionAllocatorGuard to
      // undo allocation upon destruction
      void pos(size_t i) noexcept;

   };

Since often memory can be allocated and deallocated in a stack like fashion, a FixedRegionAllocator can be used in 
combination with a guard to save/restore the allocation state, effectively releasing memory when the guard destructs.
For this we have::

   // This guard can be used to restore a FixedRegionAllocator to
   // its previous allocation state. It relies on the fact that
   // a FixedRegionAllocator is a bump the pointer allocator, and
   // can be restored by simply reseting the pointer to the previous
   // position
   class FixedRegionAllocatorGuard
   {
   public:
      FixedRegionAllocatorGuard(FixedRegionAllocator *A);
      ~FixedRegionAllocatorGuard();
   };


For scenarios where OpenRedukti is being used as a server, it is often the case that each request is served by a thread,
and while the thread executes it needs to allocate temporary memory for performing calculations. To faclitate this usage,
OpenRedukti provides some predefined thread specific allocators.::

   // Each thread is given a set of allocators to use
   // To obtain the thread specific allocator set call
   // get_threadspecific_allocators().
   struct AllocatorSet {
      RegionAllocator *cashflow_allocator;
      RegionAllocator *sensitivities_allocator;
      FixedRegionAllocator *tempspace_allocator;

      // Resets all the allocators
      // Use this after the thread has finished serving so that
      // the allocators are properly initialized for the next request
      void reset();
   };

   // Retrieves the thread specific allocator set.
   extern AllocatorSet *get_threadspecific_allocators();


Interpolators
=============

::

   #include <interpolators.h>

OpenRedukti supports the most common interpolators used in interest rate curves. The api for setting up interpolators is
described below.::

   struct InterpolationOptions;

   class Interpolator
   {
         public:
      virtual ~Interpolator() {}

      // Interpolate at x
      virtual double interpolate(double x) = 0;

      // Interpolate at x
      // And also compute sensitivities of value at x
      // to the various terms in the data set.
      // Both first order and second order sensitivies
      // can be computed depending upon how the
      // the interpolator was created.
      // Uses automatic differentiation
      virtual std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>
      interpolate_with_sensitivities(double x, FixedRegionAllocator *A) = 0;

      // Interpolate at x
      // And also compute sensitivities of value at x
      // to the various terms in the data set.
      // Both first order and second order sensitivies
      // can be computed depending upon how the
      // the interpolator was created.
      // Uses numeric differentiation
      virtual std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>
      interpolate_with_numeric_sensitivities(double x, FixedRegionAllocator *A) = 0;

      // If underlying values have changed, this
      // method can be called to reinitialise the
      // interpolator.
      virtual void update() = 0;

      // Only available on Monotone Convex interpolator as it is an
      // interest rate aware interpolator - for everything else
      // an exception will be thrown.
      virtual double forward(double x);

      // Return the interpolator type
      virtual InterpolatorType type() const = 0;

      // Returns 0 if derivatives are not enabled
      // Returns 1 if first order derivatives are enabled
      // Returns 2 if both first and second order derivatives are enabled
      virtual int order() const = 0;

      // Returns the options that are enabled
      virtual void get_options(InterpolationOptions &optons) const = 0;
   };

   struct InterpolationOptions {
      bool monotoneconvex_inputs_are_forwards;
      double cubic_left_condition_value;
      double cubic_right_condition_value;
      bool extrapolate;
      int differentiation_order;
   };

   // Return an Interpolator of the desired type.
   // The x and y arrays will be referenced by the Interpolator,
   // and therefore the caller must carefully manage
   // changes.
   extern std::unique_ptr<Interpolator, Deleter<Interpolator>>
   make_interpolator(InterpolatorType type, double *x, double *y, unsigned int size, Allocator *A,
           const InterpolationOptions &options = InterpolationOptions());


Interest Rate Curves
====================

::

   #include <curve.h>

OpenRedukti supports Zero Curves that are continuously compounded. Alternate representation using discount factors
is also supported.

There are a bunch of protocol buffers types related to curves.::

   // Curve configuration instance
   // For efficiency it is better to maintain
   // this separately from actual curve data
   // The definitions are static i.e. they do not change
   // from day to day
   message IRCurveDefinition {
      // All curve definitions must be given a unique id
      // This can be considered to be some sort of primary key
      // for the definition - i.e. no two curve definitions may
      // have the same id
      int32 id = 1;
      // The curve group is intended to allow the different
      // configurations of the same curve to be created for
      // different use cases, e.g. different interpolation methods
      // may be used for IM versus VM, or a reduced set of tenors
      // may be used for computing Liquidity Margin
      CurveGroup curve_group = 2;
      Currency currency = 3;
      IndexFamily index_family = 4;
      // Tenor is optional; if specified implies a tenor
      // specific curve
      Tenor tenor = 5;
      InterpolatorType interpolator_type = 6;
      // If interpolated_on is discount factors then it means
      // that the interpolator should operate on discount factors
      // rather than zero rates
      IRRateType interpolated_on = 7;
      // The maturity generation rule defines how the the bootstrapper
      // should generate the maturities of the curve
      MaturityGenerationRule maturity_generation_rule = 8;
      // If the curve is defined to have fixed maturity tenors
      // then a list of tenors is needed 
      // If the maturities are defined from input instruments then
      // tenors need not be defined
      repeated Tenor tenors = 9;
   }

   message ZeroCurve {
      int32 curve_definition_id = 1;
      repeated int32 maturities = 2;
      repeated double values = 3;
   }

   message ZeroCurveParSensitivities {
      int32 curve_definition_id = 1;
      int32 num_instruments = 2;
      int32 num_maturities = 3;
      // Map from <row,col> to value
      // The lower 16 bits represent the row index
      // The higher 16 bits represent the column index
      // We use this format as protobuf requires the map keys to be
      // integral type
      map<uint32, double> values = 4;
   }

The api for setting up and using curves is as follows::

   // Curve identifier
   typedef uint64_t CurveId;

   // Constructs a curve id by combining the constituents
   extern CurveId make_curve_id(PricingCurveType type, Currency ccy, IndexFamily index_family, Tenor tenor,
                 Date as_of_date, short int cycle = 0,
                 MarketDataQualifier qual = MarketDataQualifier::MDQ_NORMAL, short int scenario = 0);
   // Extracts the constituents from a curve id
   extern bool curve_id_components(CurveId id, PricingCurveType &type, Currency &ccy, IndexFamily &index_family,
               Tenor &tenor, Date &as_of_date, short int &cycle, MarketDataQualifier &qual,
               short int &scenario);

   // Gets a string representation of the curve Id,
   // note that this is an expensive operation so use only for
   // debugging
   extern std::string curve_id_to_string(CurveId id);

   class Curve
   {
   public:
      virtual ~Curve() noexcept;
      double time_from_reference(Date d) const noexcept;
      virtual const DayFraction &day_fraction() const noexcept;
      virtual Date as_of_date() const noexcept;
      virtual Date last_maturity() const noexcept;
      CurveId id() const noexcept;
      std::string name() const noexcept;
      virtual bool is_valid() const noexcept;
   };

   class YieldCurve : public Curve
   {
   public:
      virtual ~YieldCurve() noexcept;

      virtual double discount(double time) const noexcept;

      // Discount factors
      // These methods return the discount factor from a given date or time
      // to the reference date.  In the latter case, the time is calculated
      // as a fraction of year from the reference date.
      virtual double discount(Date d) const noexcept;

      // Zero-yield rates
      // These methods return the implied zero-yield rate for a
      // given date or time.  In the former case, the time is
      // calculated as a fraction of year from the reference date.
      virtual double zero_rate(Date d) const noexcept;

      // The resulting interest rate has the same day-counting rule
      // used by the term structure. The same rule should be used
      // for calculating the passed time t.
      virtual double zero_rate(double t) const noexcept;

      // Forward rates
      // These methods returns the forward interest rate between two dates
      // or times.  In the former case, times are calculated as fractions
      // of year from the reference date.
      // If both dates (times) are equal the instantaneous forward rate is
      // returned.
      virtual double forward_rate(Date d1, Date d2) const noexcept;

      // The resulting interest rate has the same day-counting rule
      // used by the term structure. The same rule should be used
      // for calculating the passed times t1 and t2.
      virtual double forward_rate(double t1, double t2) const noexcept;

      // Instantaneous forward rate
      virtual double forward(double t) const noexcept;

      // Gets the sensitivities to pillars using the underlying
      // interpolator.
      virtual std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>
      get_sensitivities(double x, FixedRegionAllocator *A) const noexcept;

      // The offset of the last pillar.
      // The first pillar is numbered 1.
      virtual int last_pillar() const noexcept;

      // Update the rates
      virtual void update_rates(const double *rates, size_t n) noexcept;

      // Value at pillar point
      virtual double value(int pillar) const noexcept;

      // maturity time from ref date
      virtual double maturity_time(int pillar) const;

      // maturity date for a pillar
      virtual Date maturity_date(int pillar) const;

      double last_maturity_time() const;

      virtual std::vector<std::unique_ptr<YieldCurve, Deleter<YieldCurve>>>
      get_bumped_curves(Allocator *A, double h = 0.00001) const noexcept;

      virtual std::unique_ptr<YieldCurve, Deleter<YieldCurve>> get_bumped_curve(Allocator *A, int pillar,
                                   double h = 0.00001) const
          noexcept;

      virtual void dump(FILE *fp = stderr) const noexcept;

      virtual InterpolatorType interpolator_type() const noexcept;

   };

   // When referencing a curve it is useful to have some
   // indirecton as this allows the curve to be modified without
   // affecting the client code. This is particularly needed when
   // bootstrapping curves. The CurveReference interface provides this
   // indirection.
   class CurveReference
   {
   public:
      virtual ~CurveReference() noexcept;
      virtual YieldCurve *get() const noexcept;
   };

   // Wraps a curve pointer
   class CurveWrapper : public CurveReference
   {
   public:
      CurveWrapper(YieldCurve *curve = nullptr) noexcept;
      virtual YieldCurve *get() const noexcept;
      void set(YieldCurve *c);
   };

   // Construct a curve
   // @param A - Memory allocator
   // @param id - ID of the curve
   // @param as_of_date - As of date
   // @param maturities - Curve pillar points
   // @param values - interpretation depends upon type below
   // @param n - Size of the arrays above
   // @param interpolator - Type of interpolator to be used
   // @param rateType - ZeroRate, DiscountFactor or FowardRate
   // @param derive_order - the order to which node sensitivities are to be
   // computed
   // @fraction - day count fraction
   //
   // Note that the curve object will copy the maturities and values arrays
   // so caller need not retain these arrays. Since the arrays are copied
   // changes to original values do not impact the curve. You can invoke
   // the method update_rates() to update the values after the curve is
   // created.
   extern std::unique_ptr<YieldCurve, Deleter<YieldCurve>>
   make_curve(Allocator *A, CurveId id, Date as_of_date, Date maturities[], double values[], size_t n,
         InterpolatorType interpolator, IRRateType type = IRRateType::ZERO_RATE, int deriv_order = 0,
         DayCountFraction fraction = DayCountFraction::ACT_365_FIXED) noexcept;

   class IRCurveDefinition;
   class ZeroCurve;

   extern std::unique_ptr<YieldCurve, Deleter<YieldCurve>>
   make_curve(Date as_of_date, const IRCurveDefinition *defn, const ZeroCurve &curve, int deriv_order,
         PricingCurveType type = PRICING_CURVE_TYPE_UNSPECIFIED, MarketDataQualifier mdq = MDQ_NORMAL,
         short int cycle = 0, short int scenario = 0);


Time Series / Fixings
=====================

::

   #include <timeseries.h>
   #include <fixings.h>


The timeseries type enables a set of date/value pairs to be managed. Values may be looked up by date. The
api is quite simple::

   class Value
   {
   public:
      Value();
      Value(Date d, double v = 0.0);
      Date date() const;
      double value() const;
   };

   class TimeSeries
   {
   public:
      TimeSeries();
      TimeSeries(size_t num_values, Value *data);
      void add(size_t num_values, Value *data); 
      ~TimeSeries();
      Value *begin();
      const Value *cbegin() const;
      Value *end();
      const Value *cend();
      bool find(Date d, double &value);
   };

The fixings service is a type that maintains timeseries data for a bunch of indices and allows values to be looked by
index.::

   class FixingDataService
   {
   public:
      FixingDataService();
      ~FixingDataService();

      void set_fixings(IndexId id, std::unique_ptr<TimeSeries> &&ts);
      TimeSeries *get_fixings(IndexId id);
   };

Cashflows
=========

::

   #include <cashflows.h>

OpenRedukti pricing approach is to convert instruments to cashflows and then price the cashflows. Once an instrument is
represented as a set of cashflows, OpenRedukti does not care what the original representation was. 

The client supplies cashflow data in the form of following protocol buffer types::

   // Simple cashflow (known amount)
   message CFSimple {
      Currency currency = 1;
      double amount = 2;
      int32 payment_date = 3;
      IsdaIndex trade_index = 4;
   }

   // A floating calculation period
   message CFFloatingPeriod {
      double notional = 1;
      double spread = 2;
      int32 accrual_start_date = 3;
      int32 accrual_end_date = 4;
      IsdaIndex index = 5;
      Tenor tenor = 6;
      IsdaIndex index2 = 7;
      Tenor tenor2 = 8;
   }

   // Floating Coupon 
   // May contain more than one calculation period
   message CFFloating { 
      Currency currency = 1;
      repeated CFFloatingPeriod floating_periods = 2;
      CompoundingMethod compounding_method = 3;
      DayCountFraction day_count_fraction = 4;
      int32 payment_date = 5;
   }

   // FRA Cashflow
   message CFFra {
      Currency currency = 1;
      double fixed_rate = 2;
      int32 payment_date = 3;
      CFFloatingPeriod floating_period = 4;
      DayCountFraction day_count_fraction = 5;
   }

   // OIS cashflow
   message CFOis {
      IsdaIndex index = 1;
      double notional = 2;
      int32 accrual_start_date = 3;
      int32 accrual_end_date = 4;
      int32 payment_date = 5;
      DayCountFraction day_count_fraction = 6;
   }

   // Single cashflow
   // This is a union type
   message CFSingle {
      oneof cashflow {
         CFSimple simple = 1;
         CFFloating floating = 2;
         CFFra fra = 3;
         CFOis ois = 4;
      }
   }

   // Cashflow stream
   message CFStream {
      repeated CFSingle cashflows = 1;
      // multiplicative factor
      // typically used to reverse direction using -1.0
      double factor = 2;
   }

   // Collection of cashflow streams
   message CFCollection {
      repeated CFStream streams = 1;
   }

The external cashflow definition must be converted to an internal representation for pricing. The api for performing this
conversion is described below.::

   // We separate out the concept of a Cashflow
   // definition (covered here) versus the valuation of
   // cashflows covered in cashflow_pricing.h.

   // The protobuf definition for a Cashflow Collection
   class CFCollection;

   // We need a way to refer to logical curve types
   // without having to reference real curves - the PricingCurve
   // helps us do that. Each PricingCurve instance represents
   // a logical identifier for a curve that will be resolved when
   // pricing via a CurveProvider implementation.
   class PricingCurve
   {

   public:
      // Defaults to 0 which is okay as it maps to unspecified
      // values component wise
      PricingCurve();
      PricingCurve(PricingCurveType type, Currency currency, IndexFamily index_family = INDEX_FAMILY_UNSPECIFIED,
              Tenor tenor = TENOR_UNSPECIFIED);
      explicit PricingCurve(uint32_t id);


      Currency currency() const;
      IndexFamily index_family() const;
      Tenor tenor() const;
      PricingCurveType curve_type() const;
      uint32_t id() const;
      bool is_valid() const;
      // Ordering is not meaningful - its purpose is to allow
      // insertion into containers
      bool operator<(const PricingCurve &c2) const;
      bool operator==(const PricingCurve &c2) const;
      bool operator!=(const PricingCurve &c2) const;
      // Get a string representation of the PricingCurve
      // Note that this is an expensive operation so use only for
      // debugging
      std::string name() const;
   };

   // Create a PricingCurve with specified type, and currency, index family
   // and tenor taken from the supplied curve Id.
   extern PricingCurve make_pricing_curve(PricingCurveType type, CurveId id);

   // When generating cashflows we do not know what actual curves will
   // be used - and whether the forward and discount curves map to the same
   // curve or different curves, or whether different tenor curves map to
   // different curves or the same curve. The CurveMapper allows the caller
   // to provide a mapping to the desired 'logical' curve. The mapping is
   // logical so that given a logical curve id, another function must obtain
   // an instance of the real curve.
   class CurveMapper
   {
   public:
      virtual ~CurveMapper();
      virtual PricingCurve map_index_tenor(PricingCurveType curve_type, Currency currency,
                       IndexFamily family = IndexFamily::INDEX_FAMILY_UNSPECIFIED,
                       Tenor tenor = Tenor::TENOR_UNSPECIFIED) const;
   };

   class ValuationContext
   {
   public:
      virtual ~ValuationContext();
      virtual Date evaluation_date() const;
      virtual Date payment_cutoff_date() const;
      virtual int derivative_order() const;
      // include today's fixing (e.g. eod)
      // If false then the curve will be used to determine the
      // rate. The bootstrapper requires this to be false so we set
      // the default value to false
      virtual bool include_todays_fixing() const;
      // Retrieve a fixing.
      // If the fixing date is < evaluation date then the absence of a fixing
      // will be an error reported via status. If the fixing date is ==
      // evaluation date then a missing fixing is not treated as error -
      // instead the method will return false; in all cases a true return
      // value indicates that the fixing was found and is set
      virtual bool get_fixing(IndexId fixing_key, Date fixing_date, double &fixing, StatusCode &status) const;
   };

   class Cashflows;

   // Converts the CFCollection to internal cashflow format
   extern Cashflows *construct_cashflows(RegionAllocator *A, const CFCollection *cfcollection, const ValuationContext &ctx,
                     const CurveMapper *curve_mapper);

Cashflow Pricing
================

::
   
   #include <cashflow_pricing.h>

Once cashflows are converted to internal format, and Zero Rate / Discount Factor curves are available, you can invoke the
cashflow pricing functions described below.::

   class Sensitivities;

   // First order sensitivities (i.e. delta)
   struct Sensitivities1D {
   public:
      Sensitivities1D(const CurveReference *curve, Allocator *A);
      ~Sensitivities1D();
      YieldCurve *curve() const;
      double at(size_t i) const;
      double &at(size_t i);
      int count() const;
   };

   // Second order sensitivities (i.e. gamma)
   struct Sensitivities2D {
   public:
      Sensitivities2D(const CurveReference *curve1, const CurveReference *curve2, Allocator *A);
      ~Sensitivities2D();
      YieldCurve *curve1() const;
      YieldCurve *curve2() const;
      double at(size_t i, size_t j);
      double &at(size_t i, size_t j);
      int count1() const;
      int count2() const;
   };

   class Sensitivities
   {
   public:
      Sensitivities(Allocator *A);
      ~Sensitivities();
      // Find or add
      Sensitivities1D *first_order_sensitivities(YieldCurve *curve);
      // Find or add
      Sensitivities2D *second_order_sensitivities(YieldCurve *curve1, YieldCurve *curve2);
      Sensitivities1D *find_first_order_sensitivities(YieldCurve *curve) const;
      Sensitivities1D *find_first_order_sensitivities(CurveId id) const;
      Sensitivities2D *find_second_order_sensitivities(YieldCurve *curve1, YieldCurve *curve2) const;
      Sensitivities2D *find_second_order_sensitivities(CurveId id1, CurveId id2) const;

      void reset();
      // Find or add
      const CurveReference *get(YieldCurve *curve);
      const CurveReference *find(YieldCurve *curve) const;
      void get_curve_ids(std::vector<CurveId> &ids) const;
   };

   class Cashflows;

   // Calculate sensitivities (delta and gamma) numerically
   // and store in supplied container
   extern void compute_sensitivity_numerically(FixedRegionAllocator *allocator, const Cashflows *flows,
                      const CurveReference *discount_curve, const CurveReference *forward_curve1,
                      const CurveReference *forward_curve2, Sensitivities *sensitivities,
                      StatusCode &status, double h = 0.00001);

   // When cashflows are defined, they reference logical curves via
   // PricingCurve identifiers. At the time of valuation these logical curves
   // must be mapped to physical instances of curves - the CurveProvider
   // interfaces= defines such a component.
   class CurveProvider
   {
   public:
      virtual ~CurveProvider() {}
      virtual const CurveReference *get_curve(PricingCurve curve) const = 0;
   };

   // Calculate PV and if ValuationContext.derivative_order > 0
   // then also delta and gamma
   extern double compute_present_value(FixedRegionAllocator *A, const ValuationContext &ctx, const Cashflows *flows,
                   const CurveProvider *mapping_provider, Sensitivities &sensitivities,
                   StatusCode &status);

   extern double compute_present_value(FixedRegionAllocator *A, const Cashflows *flows,
                   const CurveReference *discount_curve, const CurveReference *forward_curve1,
                   const CurveReference *forward_curve2, const ValuationContext &ctx,
                   Sensitivities &sensitivities, StatusCode &status);

Curve Building
==============

::

   #include <bootstrap.h>

Most of the data required to build curves is described in protocol buffers types.::

   message ParInstrument {
      // instrument type is used to decide the pricing algorithm to use
      string instrument_type = 1;
      // The instrument_key is a way to identify the instrument
      // for futures it must be MmmYY where MMM is the expiry month
      // for Fras its must be nnxnnF  
      string instrument_key = 2;
      // A reference to a curve within the owning set
      int32 discount_curve_definition_id = 3;
      // A reference to a curve within the owning set
      int32 forward_curve_definition_id = 4;
      // For instruments that reset on the floating side there
      // needs to be a floating reset frequency
      // Only required if different from the curve tenor
      Tenor floating_tenor = 5;
   };

   message ParRates {
      // We assume that all instrument definitions can be located by 
      // a numeric id - that is given the id there is a way to locate the
      // instrument, maybe by looking up in a database
      // Note that the boostrapper does not use these ids
      repeated int32 instrument_ids = 1;
      repeated double values = 2;
   }

   message ParCurve {
      int32 curve_definition_id = 1;
      repeated ParInstrument instruments = 2;
      ParRates par_rates = 3;
   }

   message ParCurveSet {
      // It is not clear that any meaningful values can be
      // assigned to ctycle, qualifier or scenario prior to bootstrapping
      // so these fields probably only make sense afterwards
      int32 as_of_date = 1;
      int32 cycle = 2;
      MarketDataQualifier qualifier = 3; 
      int32 scenario = 4;
      repeated ParCurve par_curves = 5;
   }

   enum SolverType {
      SOLVER_TYPE_LEVENBERG_MARQUARDT = 0;
      SOLVER_TYPE_LINEAR_LEAST_SQUARE = 1;
      SOLVER_TYPE_LINEAR_LUFACTOR = 2;
   }

   // The bootstrap request is self contained
   // i.e. all required data must be submitted so that
   // the request can be handled in a stateless manner
   message BootstrapCurvesRequest {
      int32 business_date = 1;
      repeated IRCurveDefinition curve_definitions = 2;  
      ParCurveSet par_curve_set = 3;
      // If true will attempt to generate par sensitvities  
      bool generate_par_sensitivities = 4;
      SolverType solver_type = 5;
      int32 max_solver_iterations = 6;
   }

   message BootstrapCurvesReply {
      ReplyHeader header = 1;
      repeated ZeroCurve curves = 2;
      // The sensitivity of zero rates to par rates
      repeated ZeroCurveParSensitivities par_sensitivities = 3;
   }


The api for invoking the curve builder is relatively simple.::

   class CurveBuilderService
   {
   public:
      virtual ~CurveBuilderService();
      virtual BootstrapCurvesReply *handle_bootstrap_request(google::protobuf::Arena *arena,
                               const BootstrapCurvesRequest *request);
   };

   std::unique_ptr<CurveBuilderService> get_curve_builder_service();

The curve building service uses Ravi scripting to define the cashflows for the instruments used in the curve.
To understand how this works, it is necessary to first understand the scripting interface, hence this subject will
be covered in that section.

Valuation Service
=================

::

   #include <valuation.h>

The Valuation Service brings together some of the other components of OpenRedukti. It enables deployment of OpenRedukti as 
a server. The service accepts all the market data via protocol buffers format messages, and then enables clients to invoke
pricing of cashflows.

The message definitions used by this service are as follows.::

   message PricingContext {
      int32 as_of_date = 1;
      MarketDataQualifier qualifier = 2;
      int32 cycle = 3;
      int32 payment_cutoff_date = 4;
      int32 derivative_order = 5;
      bool is_todays_fixings_included = 6;
      CurveGroup curve_group = 7;
      // Starting scenario; 0 is the current scenario,
      // historical scenarios start from 1 and go up.
      int32 from_scenario = 8;
      int32 to_scenario = 9;
   }

   enum SensitivityTypeCode {
      STC_ZERO = 0;
      STC_PAR = 1;
   }

   enum SensitivityOrderCode {
      SOC_DELTA = 0;
      SOC_GAMMA = 1;
   }

   enum SensitivityRiskCode {
      SRC_UNSPECIFIED = 0;
      SRC_FORWARD = 1;
      SRC_DISCOUNT = 2;
   }

   message IRCurveSensitivities {
      // Zeror or PAR sensitivities?
      SensitivityTypeCode sensitivity_type = 1;
      // Delta or Gamma ?
      SensitivityOrderCode order = 2;
      // For each dimension a curve identifier is needed
      // For delta curves there is only one dimension so only 
      // one curve will be present
      // Second order sensitivities have two dimensions
      // -1 if not applicable
      int32 curve_definition_id_1 = 3;
      int32 curve_definition_id_2 = 4;

      // For each dimension the type of risk being measured
      // is required
      SensitivityRiskCode risk_type_1 = 5;
      SensitivityRiskCode risk_type_2 = 6;

      // Map from <row,col> to value
      // The lower 32 bits represent the row index
      // The higher 32 bits represent the column index
      // We use this format as protobuf requires the map keys to be
      // integral type
      map<uint32, double> values = 7;
   }

   message ValuationRequest {
      PricingContext pricing_context = 1;
      CFCollection cashflows = 2;
   }

   message ValuationResult {
      // Valuations by scenario
      map<int32, double> valuations = 1;
      // Sensitivities for scenario 0 only
      repeated IRCurveSensitivities sensitivities = 2;
   }

   message ValuationReply {
      ReplyHeader header = 1; 
      ValuationResult result = 2;
   }

   message FixingsByIndexTenor {
      IsdaIndex index = 1;
      Tenor tenor = 2;
      // Map of fixing date to fixing value
      map<int32, double> fixings = 3;
   }

   // Publish fixings data to the backend
   message SetFixingsRequest {
      FixingsByIndexTenor fixings_by_index_tenor = 1;
   }

   message SetFixingsReply {
      ReplyHeader header = 1;
   }

   // Publish curve definitions to the backend
   message RegisterCurveDefinitionsRequest {
      repeated IRCurveDefinition curve_definitions = 1;
   }

   message RegisterCurveDefinitionsReply {
      ReplyHeader header = 1;
   }

   // Publish zero curves to the backend
   message SetZeroCurvesRequest {
      int32 as_of_date = 1;
      int32 cycle = 2;
      MarketDataQualifier qualifier = 3; 
      int32 scenario = 4;
      CurveGroup curve_group = 5;
      repeated ZeroCurve forward_curves = 6;
      repeated ZeroCurve discount_curves = 7;
      repeated ZeroCurveParSensitivities par_sensitivities = 8;
   }

   message SetZeroCurvesReply {
      ReplyHeader header = 1;
   }

   // A logical way of identifying a curve
   // Note that curves are assumed to belong to the same
   // context - i.e. business date, curve group,
   // scenario etc.
   message PricingCurveIdentifier {
      PricingCurveType type = 1; 
      Currency currency = 2;
      IndexFamily index_family = 3;
      Tenor tenor = 4;
   }

   // This mapping says that whenever the cashflow
   // would have looked for a curve with logical id
   // 'from_id' it should use logical curve with 'to_id'.
   // Note that the mapping is not recursive, i.e.
   // if 'to_id' was mapped also then that would not affect
   // the outcome of mapping 'from_id'
   message CurveMapping {
      PricingCurveIdentifier from_id = 1;
      PricingCurveIdentifier to_id = 2;
   }

   // Before any cashflow valuation can be done
   // one of the pre-requisites is to provide 
   // mappings for the logical curves. If a mapping
   // is not provided then the logical curve maps
   // to itself. The aim of the mapping is to allow
   // the cashflow pricer to be unaware of actual curve
   // assignments when performing valuations. Note that the
   // sensitivites are calculated against each logic curve
   // so the mappings affect PV and sensitivities.
   message SetCurveMappingsRequest {
      CurveGroup curve_group = 1;
      //PricingCurveType default_curve_type = 2;
      repeated CurveMapping mappings = 3;
   }

   message SetCurveMappingsReply {
      ReplyHeader header = 1;
   }

   message CurveInterpolationRequest {
      int32 business_date = 1;
      IRCurveDefinition definition = 2;
      ZeroCurve curve = 3;
      // Specify the rate type for which values are being
      // requested. If forward rate is requested then
      // forward_tenor attribute can be set to request a 
      // specific tenor
      IRRateType rate_type = 4;
      // If forward rates are requested then
      // Specify the forward tenor for which forward rates
      // should be returned; if not specified then the
      // tenor associated with the index will be returned
      Tenor forward_tenor = 5;
      repeated int32 dates = 6;
   }

   message CurveInterpolationReply {
      ReplyHeader header = 1; 
      repeated double values = 2;
   }

The api for interacting with the ValuationService is shown below.::

   class ValuationService
   {
   public:
      virtual ~ValuationService();
      virtual CurveInterpolationReply *
      handle_curve_interpolation_request(google::protobuf::Arena *arena,
                     const CurveInterpolationRequest *request);
      virtual SetCurveMappingsReply *handle_set_curve_mappings_request(google::protobuf::Arena *arena,
                               const SetCurveMappingsRequest *request);
      virtual SetZeroCurvesReply *handle_set_zero_curves_request(google::protobuf::Arena *arena,
                              const SetZeroCurvesRequest *request);
      virtual RegisterCurveDefinitionsReply *
      handle_register_curve_definitions_request(google::protobuf::Arena *arena,
                       const RegisterCurveDefinitionsRequest *request);
      virtual SetFixingsReply *handle_set_fixings_request(google::protobuf::Arena *arena,
                            const SetFixingsRequest *request);
      virtual ValuationReply *handle_valuation_request(google::protobuf::Arena *arena,
                         const ValuationRequest *request);
   };


