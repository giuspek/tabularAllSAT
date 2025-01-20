#include <cinttypes>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstring>
#include <fstream>

#include <algorithm>
#include <vector>
#include <set>
#include <gmpxx.h>
#include <iostream>

#include "arena.hpp"
#include "compact.hpp"
#include "heap.hpp"
#include "shared.hpp"

// Core data structure for clauses.
bool debugging = false;

struct Clause
{
  // We use bit fields to save space and accordingly have to limit the
  // maximum glue to fit the number of bits ('27') for the 'glue' field.

  static const unsigned MAXIMUM_GLUE = (1u << 26) - 1;

  unsigned size;          // Actual size of 'literals'.
  unsigned glue : 26;     // Glucose level (LBD).
  bool garbage : 1;       // Marked as garbage clause.
  bool reason : 1;        // Active reason clause.
  bool redundant : 1;     // Learned redundant clause.
  bool backtrue : 1;      // Learned redundant clause.
  unsigned char used : 2; // Used since last reduction.
  unsigned searched;      // Last replacement search position.
  unsigned *position;     // Position in the clause for implied_w.

#ifdef LOGGING
  uint64_t id; // For debugging only.
#endif

  unsigned literals[]; // Actually 'literals[size]'.

  // Provides simple iteration of the literals of the clause.

  unsigned *begin() { return literals; }
  unsigned *end() { return literals + size; }

  // Support for the low-level clause arena memory allocators.

  Clause *next() { return (Clause *)(literals + size); }

  // The number of words has to match the number of unsigned words 'n'
  // needed for allocation as given to the arena 'allocate ()' function.
  // It is determined by first computing the number of needed bytes.

  static size_t bytes(unsigned size)
  {
    assert(sizeof(unsigned) == 4);
    return 4 * (size_t)size + (size_t)(((Clause *)0)->literals);
  }

  static size_t words(unsigned size)
  {
    size_t bytes = Clause::bytes(size);
    assert(!(bytes & 3));
    return bytes / 4;
  }
};

typedef class Arena<Clause> Clauses;

// The reference to a clause (or a literal).

typedef unsigned Reference;
typedef Compact<Reference> References;

// Unit and decision reasons are references with specific constants.  We
// need to distinguish unit from decision reasons due to chronological
// backtracking which might assign units on higher decision levels with of
// course still assignment level zero.  We use 'INVALID' aka 'UINT_MAX' and
// one less which are checked never ever to point to the start of a clause.

static const Reference DECISION_REASON = INVALID;
static const Reference UNIT_REASON = INVALID - 1;
static const Reference BACKTRUE_REASON = INVALID - 2;

// The largest position of a clause in the arena tagged with zero needs to
// be smaller than the 'UNIT_REASON = INVALID-1 = 0xfffffffe' thus the
// actual position smaller than '(INVALID-1)/2 = 0x7fffffff'.

static const unsigned MAXIMUM_POSITION_OF_CLAUSE = 0x7ffffffd;

// Variable flags.

struct Flags
{
  bool minimizable; // Cache that this variable is minimizable.
  bool poison;      // Cache that this variable is not minimizable.
  bool seen;        // Variable is resolved or in the learned clause.

  Flags() : minimizable(false), poison(false), seen(false) {}
};

struct Phase
{
  // Most consistent assignment reset and reused in best phase rephasing.

  signed char best = 0;

  // Saved phase during assignment is fall back in stable mode and the main
  // phase decision heuristic in focused mode.

  signed char saved = 0;

  // Maximum consistent assignment since last rephase takes priority in
  // deciding the phase (value) during stable mode over saved phase.

  signed char target = 0;
};

// Unbiased exponential moving averages (check the ADAM paper for details).

struct Average
{
  double value, alpha, beta, biased, exp;

  Average(double a) : value(0), alpha(a), beta(1 - a), biased(0), exp(1) {}

  operator double() { return value; }

  void operator+=(double y)
  {
    biased += alpha * (y - biased);
    double new_exp = exp * beta;
    exp = (new_exp == exp ? 0 : new_exp);
    value = exp ? biased / (1 - exp) : biased;
  }
};

// Common statistics for profiling and scheduling functions.

struct Statistics
{

  uint64_t added_clauses;
  uint64_t original_clauses;
  uint64_t analyzed;
  uint64_t chronological;
  uint64_t conflicts;
  uint64_t decisions;
  uint64_t deduced_literals;
  unsigned fixed;
  uint64_t learned_literals;
  uint64_t propagations;
  uint64_t recycled_clauses;
  uint64_t reductions;
  uint64_t redundant_clauses;
  uint64_t rephased;
  uint64_t reported;
  uint64_t restarts;
  uint64_t simplifications;
  mpz_class n_assignments_partial;

  Statistics() { memset(this, 0, sizeof *this); }

  double average(double a, double b) { return b ? a / b : 0; }
  double percent(double a, double b) { return average(100 * a, b); }

  void print(double t)
  {

    uint64_t minimized_literals = deduced_literals - learned_literals;

    // Print statistics in alphabetic order.
    gmp_printf("c %-20s %Zd\n", "n-partial-assignments", n_assignments_partial);
    printf("c %-20s %13" PRIu64 " %13.2f %% conflicts\n", "analyzed:", analyzed, percent(analyzed, conflicts));
    printf("c %-20s %13" PRIu64 " %13.2f %% analyzed\n", "chronological:", chronological, percent(chronological, analyzed));
    printf("c %-20s %13" PRIu64 " %13.2f per second\n", "conflicts:", conflicts, average(conflicts, t));
    printf("c %-20s %13" PRIu64 " %13.2f per conflict\n", "decisions:", decisions, average(decisions, conflicts));
    printf("c %-20s %13.2f\n", "learned-clause-size:", average(learned_literals, conflicts));
    printf("c %-20s %13.2f %%\n", "minimized-literals:", percent(minimized_literals, deduced_literals));
    printf("c %-20s %13" PRIu64 " %13.2f per second\n", "propagations:", propagations, average(propagations, t));
    printf("c %-20s %13" PRIu64 " %13.2f %% added\n", "recycled-clauses:", recycled_clauses, percent(recycled_clauses, added_clauses));
    printf("c %-20s %13" PRIu64 " %13.2f interval\n", "reductions:", reductions, average(conflicts, reductions));
    printf("c %-20s %13" PRIu64 " %13.2f interval\n", "rephased:", rephased, average(conflicts, rephased));
    printf("c %-20s %13" PRIu64 " %13.2f interval\n", "restarts:", restarts, average(conflicts, restarts));
    printf("c %-20s %13" PRIu64 " %13.2f %% reductions\n", "simplifications:", simplifications, percent(simplifications, reductions));
    printf("c %-20s %13" PRIu64 " %13.2f %% conflicts\n", "units:", reported, percent(reported, conflicts));
  }
};

// Limits for scheduling functions.

struct Limits
{
  uint64_t reduce = 0;
  uint64_t rephase = 0;
  uint64_t restart = 0;
};

// Internal list of alphabetically ordered options.

struct SolverOptions : Options
{
  Option *begin() { return &backjump_limit; }
  Option backjump_limit = Option("backjump-limit", 0, 0, INT_MAX, "chronological backjumping limit");
  Option initial_phase = Option("initial-phase", 0, 0, 1, "default initial phase assigned to variables");
  Option reduce_fraction = Option("reduce-fraction", 75, 0, 100, "percentage of reduced redundant clause");
  Option reduce_interval = Option("reduce-interval", 1000, 1, INT_MAX, "conflicts between reduction");
  Option reduce = Option("reduce", 1, 0, 1, "reduction (forgetting) of learned clauses");
  Option rephase_interval = Option("rephase-interval", 1000, 10, 100000, "conflicts between resetting phases");
  Option rephase = Option("rephase", 1, 0, 1, "reset phases in regular intervals");
  Option restart_interval = Option("restart-interval", INT_MAX, 1, INT_MAX, "conflicts between restarts");
  Option restart_margin = Option("restart-margin", 100, 0, 100, "margin in percent fast-over-slow glue average for restart");
  Option restart = Option("restart", 0, 0, 1, "schedule restarts in regular intervals");
  Option tier1_max_glue = Option("tier1-max-glue", 2, 0, INT_MAX, "tier1 glucose level for always kept clauses");
  Option tier2_max_glue = Option("tier2-max-glue", 6, 0, INT_MAX, "tier2 glucose level for clauses with second chance");
  Option vsids_decay = Option("vsids-decay", 50, 0, 500, "VSIDS scores decay (alpha) in per mille");
  Option output_file = Option("output-file", 0, 0, 1, "set if printing models as stdout (0) or in output.txt (1)");
  Option *end() { return 1 + &output_file; }
};

// Decision control stack frame.

struct Frame
{
  unsigned decision; // Decision literal of frame.
  unsigned trail;    // Trail height when entering frame.
  bool pulled;       // Needed to determine glue during analysis.

  Frame(unsigned d = 0, unsigned t = 0) : decision(d), trail(t), pulled(false) {}
};

// The actual solver code starts here.

class Solver : public Shared
{
  bool inconsistent = false; // Found or learned empty clause.
  bool initialized = false;  // Limits need to be initialized.
  bool iterating = false;    // Learned unit clause.

  bool enumerate_total = false;

  std::vector<Frame> control; // One frame per decision level.

  unsigned propagated = 0; // Propagated literals in 'trail'.
  std::vector<unsigned> trail;

  unsigned best_assigned = 0;
  unsigned target_assigned = 0;
  unsigned aggressive_level_limit = 0;
  unsigned lit_to_flip;

  Clauses clauses;                // The clause data memory arena.
  std::vector<References> matrix; // Maps literals to references.

  std::vector<Reference> reasons; // Maps variables to reasons.
  std::vector<Phase> phases;      // Maps variables to phases.
  std::vector<Flags> flags;       // Maps variables to flags.
  std::vector<int> position;
  
  Heap heap; // Decision priority heap.

  std::vector<unsigned> analyzed; // Variables analyzed.
  std::vector<unsigned> pulled;   // Decision levels analyzes.

  std::vector<unsigned> minimized; // Minimizable variables.
  std::vector<unsigned> poisoned;  // Poisoned variables.

  Limits limits;
  Statistics statistics;
  SolverOptions options;

  Average slow_glue_average = Average(1e-5);
  Average fast_glue_average = Average(3e-2);

#ifdef LOGGING
  void debug_clause(Clause *, const char *, ...)
      __attribute__((format(printf, 3, 4)));
  void debug_binary(unsigned, unsigned, const char *, ...)
      __attribute__((format(printf, 4, 5)));
  void debug_clauses();
#else
#define debug_clause(...) \
  do                      \
  {                       \
  } while (0)
#define debug_binary(...) \
  do                      \
  {                       \
  } while (0)
#endif

  Reference new_binary(bool backreason = false);
  void bump_clause(Clause *);
  Clause *dereference(Reference ref);
  Reference new_clause(bool redundant = false, unsigned glue = 0, bool backreason = false);
  void watch_clause(unsigned lit, unsigned other, Reference);
  void watch_binary(unsigned lit, unsigned other);
  bool satisfied_clause(Clause *);
  void delete_clause(Clause *);

  void assign(unsigned lit, Reference reason);
  void unassign(unsigned lit);

  unsigned decide_phase(unsigned);
  unsigned decide_variable();
  void decide();

  bool propagate(unsigned &falsified, Reference &conflict);
  void backtrack(unsigned new_level);

  void minimize();
  bool minimize(unsigned lit, unsigned depth = 0);
  void analyze_literal(unsigned lit,
                       unsigned &open, unsigned &glue, unsigned &jump);
  bool analyze_conflict(unsigned propagated_literal, Reference conflict);
  void bump_analyzed_variables();

  void report(char type);

  unsigned first_redundant = INVALID;
  unsigned fixed_at_reduce = 0;
  void mark_reason_clauses(bool);
  void gather_reduce_candidates(std::vector<unsigned> &);
  unsigned mark_garbage(std::vector<unsigned> &);
  void mark_satisfied_clauses_as_garbage_too(unsigned &);
  void flush_garbage_watches(unsigned first_garbage, bool simplify,
                             std::vector<unsigned> &);
  void recycle_garbage_clauses(unsigned first_garbage,
                               std::vector<unsigned> &);
  void map_reasons(unsigned first_garbage, std::vector<unsigned> &);
  unsigned gather_and_mark_reduced_clauses();
  void recycle_clauses_and_flush_watches(unsigned first_garbage);
  bool reducing();
  void reduce();

  struct
  {
    uint64_t u, v;
  } reluctant = {1, 1};
  void set_restart_limit();
  bool restarting();
  void restart();

  bool rephasing();
  void rephase_best();
  void rephase_inverted();
  void rephase_original();
  void update_phases();
  void rephase();

  State cdcl();
  bool block_chrono(unsigned M);
  mpz_class get_models_covered_by_assignment(unsigned length);
  int to_dimacs(unsigned lit);
  unsigned get_implicant_lifting();
  unsigned get_aggressive_implicant_lifting();
  unsigned get_aggressive_projected_implicant_lifting();
  unsigned get_aggressive_projected_implicant_lifting_bis();
  unsigned get_projection_lifting();
  void initialize_internal_solving();

  friend struct less_relevant;

  // The six virtual functions of 'Interface' called from 'Shared'
  // functions are here.  The other functions above are all private.

  const char *one_word_name() { return "TabularAllSAT+"; }
  const char *one_line_description()
  {
    return "TabularAllSAT+: CB-CDCL AllSAT Solver without blocking clauses";
  }
  void adjust_to_increased_variables();
  void add_simplified_irredundant_clause();
  State solve_internally();
  void print_statistics(double t) { statistics.print(t); }
  Options &iterate_options() { return options; }

public:
  Solver();
};

Solver::Solver() : inconsistent(false),
                   initialized(false),
                   iterating(false),
                   propagated(0),
                   first_redundant(INVALID),
                   fixed_at_reduce(0)
{
  control.push_back(Frame(INVALID, 0));
}

void Solver::adjust_to_increased_variables()
{
  assert(flags.size() < variables);
  phases.resize(variables);
  flags.resize(variables);
  matrix.resize(2 * variables);
  reasons.resize(variables);
  position.resize(variables);
  assert(literals == matrix.size());
  assert(variables == reasons.size());
}

static inline bool tagged(Reference ref) { return ref & 1; }
static inline unsigned untag(Reference ref) { return ref >> 1; }

static inline Reference tag_binary_or_large(unsigned d, unsigned t)
{
  assert(t == 0 || t == 1);
  Reference res = (d << 1) | t;
  assert(tagged(res) == t);
  assert(untag(res) == d);
  assert(res != UNIT_REASON);
  assert(res != DECISION_REASON);
  return res;
}

static inline Reference tag_binary(unsigned other)
{
  assert(other <= MAXIMUM_INTERNAL_LITERAL);
  return tag_binary_or_large(other, 1);
}

static inline Reference tag_large(unsigned pos)
{
  assert(pos <= MAXIMUM_POSITION_OF_CLAUSE);
  return tag_binary_or_large(pos, 0);
}

inline Clause *Solver::dereference(Reference ref)
{
  assert(!tagged(ref));
  unsigned i = untag(ref);
  return clauses[i];
}

void Solver::delete_clause(Clause *clause)
{
  Statistics &s = statistics;
  if (clause->redundant)
    assert(s.redundant_clauses), s.redundant_clauses--;
  debug_clause(clause, "delete");
}

/*------------------------------------------------------------------------*/

#ifdef LOGGING

void Solver::debug_clause(Clause *clause, const char *fmt, ...)
{
  assert(clause);
  if (!logging())
    return;
  debug_prefix();
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
  if (clause->redundant)
    printf(" redundant glue %u", clause->glue);
  else
    fputs(" irredundant", stdout);
  printf(" size %u clause[%" PRIu64 "]", clause->size, clause->id);
  for (auto lit : *clause)
    printf(" %s", debug_literal(lit));
  debug_suffix();
}

void Solver::debug_binary(unsigned lit, unsigned other,
                          const char *fmt, ...)
{
  if (!logging())
    return;
  debug_prefix();
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
  printf(" binary clause %s %s",
         debug_literal(lit), debug_literal(other));
  debug_suffix();
}

void Solver::debug_clauses()
{
  for (auto clause : clauses)
    debug_clause(clause, "clauses[%u]", clauses.position(clause));
}

#endif

#define debug_reference(REF, ...) \
  debug_clause(dereference(REF), __VA_ARGS__)

void Solver::assign(unsigned lit, Reference reason)
{

  // With chronological backtracking enabled variables are assigned at an
  // assignment level computed from the reason which might be smaller than
  // the current decision level (and thus the assignment out-of-order). In
  // the following part we compute the assignment level explicitly.

  unsigned assignment_level;
  {
    if (reason == DECISION_REASON)
      assignment_level = level;
    else if (reason == BACKTRUE_REASON)
      assignment_level = level;
    else if (reason == UNIT_REASON)
      assignment_level = 0;
    else if (!tagged(reason))
    {
      assignment_level = 0;
      Clause *clause = dereference(reason);
      for (auto other : *clause)
      {
        if (other == lit)
          continue;
        unsigned other_level = levels[index(other)];
        if (other_level > assignment_level)
          assignment_level = other_level;
      }
    }
    else
    {
      unsigned other = untag(reason);
      assignment_level = levels[index(other)];
    }
  }

  // Here comes the usual assignment logic.

  unsigned idx = index(lit);
  levels[idx] = assignment_level;
  if (!assignment_level)
  {
    statistics.fixed++;
    reason = UNIT_REASON;
  }
  reasons[idx] = reason;
  phases[idx].saved = sign(lit) ? -1 : 1;

  unsigned not_lit = negate(lit);
  assert(!values[lit]);
  assert(!values[not_lit]);
  values[lit] = 1;
  values[not_lit] = -1;
  trail.push_back(lit);

  position[idx] = trail.size() - 1;

#ifdef LOGGING
  if (logging())
  {
    const char *s = debug_literal(lit);
    if (reason == DECISION_REASON)
      debug("assign %s as decision", s);
    else if (reason == UNIT_REASON)
      debug("assign %s as unit", s);
    else if (reason == BACKTRUE_REASON)
      debug("assign %s with virtual backtrue reason", s);
    else if (tagged(reason))
      debug_binary(lit, untag(reason), "assign %s with reason", s);
    else
      debug_reference(reason, "assign %s reason", s);
    if (assignment_level != level)
      debug("out-of-order assigned %s", s);
  }
#endif
}

void Solver::watch_clause(unsigned lit, unsigned blocking, Reference ref)
{
  debug_reference(ref, "watching %s with blocking literal %s in",
                  debug_literal(lit), debug_literal(blocking));
  matrix[lit].push_back(tag_large(blocking));
  matrix[lit].push_back(ref);
}

void Solver::watch_binary(unsigned lit, unsigned blocking)
{
  debug_binary(lit, blocking, "watching %s with blocking %s in",
               debug_literal(lit), debug_literal(blocking));
  matrix[lit].push_back(tag_binary(blocking));
}

Reference Solver::new_binary(bool backreason)
{
  assert(clause.size() == 2);
  unsigned lit = clause[0], other = clause[1];
  debug_binary(lit, other, "new binary clause");

  if (!backreason)
  {
    watch_binary(lit, other);
    watch_binary(other, lit);
  }
  return tag_binary(other);
}

void Solver::bump_clause(Clause *c)
{
  if (!c->redundant)
    return;
  c->used = 1 + (c->glue <= options.tier2_max_glue);
  debug_clause(c, "bumped");
}

Reference Solver::new_clause(bool redundant, unsigned glue, bool backreason)
{

  if (!redundant)
    added_original_non_binary_clauses+=1;

  size_t size = clause.size();
  assert(3 <= size);
  assert(size <= MAXIMUM_INTERNAL_VARIABLE);

  size_t words = Clause::words(size);
  assert(words >= 4);

  size_t pos = clauses.size();
  if (clauses.size() > MAXIMUM_POSITION_OF_CLAUSE)
    error("maximum clause arena size exceeded");

  if (redundant && first_redundant == INVALID)
    first_redundant = pos;

  Clause *c = clauses.allocate(words);
  Reference ref = tag_large(pos);

  c->size = size;
  c->glue = glue < Clause::MAXIMUM_GLUE ? glue : Clause::MAXIMUM_GLUE;
  assert(c->glue == glue);
  c->searched = 2;

#ifdef LOGGING
  c->id =
#endif
      statistics.added_clauses++;
  if (redundant)
    statistics.redundant_clauses++;

  c->garbage = false;
  c->reason = false;
  c->redundant = redundant;
  c->used = 0;

  memcpy(c->literals, &clause[0], size * sizeof(unsigned));

  debug_clause(c, "new");
  if (redundant)
    bump_clause(c);

  // Watch the first two literals of the clause.
  if (!backreason)
  {
    watch_clause(clause[0], clause[1], ref);
    watch_clause(clause[1], clause[0], ref);
  }

  assert(dereference(ref) == c);

  return ref;
}

void Solver::add_simplified_irredundant_clause()
{
  size_t size = clause.size();
  if (size > 1){
    unsigned lit = clause[0];
    unsigned other = clause[1];
    int idx1 = index(lit) + 1;
    if (lit % 2 == 1)
      idx1 = -idx1;
    int idx2 = index(other) + 1;
    if (other % 2 == 1)
      idx2 = -idx2;
  }
  if (!size)
  {
    debug("parsed empty clause");
    inconsistent = true;
  }
  else if (size == 1)
  {
    unsigned unit = clause[0];
    debug("parsed unit clause %s", debug_literal(unit));

    int idx = index(unit) + 1;
    if (unit % 2 == 1)
      idx = -idx;
    assert(!level);
    assign(unit, UNIT_REASON);
  }
  else if (size == 2)
    new_binary();
  else
    new_clause();
}

bool Solver::propagate(unsigned &falsified, Reference &conflict)
{

  bool res = true;

  uint64_t propagations = 0;

  while (res && propagated != trail.size())
  {

    propagations++;
    unsigned lit = trail[propagated++];
    unsigned not_lit = negate(lit);
    debug("propagating %s", debug_literal(lit));

    auto &not_lit_watches = matrix[not_lit];
    auto begin_not_lit_watches = not_lit_watches.begin();
    auto end_not_lit_watches = not_lit_watches.end();
    auto p = begin_not_lit_watches, q = p;

    while (p != end_not_lit_watches)
    {

      Reference watch = *q++ = *p++;
      unsigned blocking = untag(watch);
      signed char blocking_value = values[blocking];

      if (tagged(watch))
      {

        if (blocking_value > 0)
          continue;

        Reference reason = tag_binary(not_lit);

        if (blocking_value < 0)
        {
          debug_binary(not_lit, blocking, "conflicting");
          falsified = blocking;
          conflict = reason;
          res = false;
          break;
        }

        assign(blocking, reason);
        continue;
      }

      Reference ref = *q++ = *p++;
      Clause *c = dereference(ref);
      if (blocking_value > 0)
        continue;

      unsigned *literals = c->literals;
      unsigned other = literals[0] ^ literals[1] ^ not_lit;
      signed char other_value = values[other];
      if (other_value > 0)
      {
        q[-2] = tag_large(other);
        continue;
      }

      unsigned *end_literals = literals + c->size;
      unsigned *searched = literals + c->searched;
      unsigned *l, replacement;

      for (l = searched; l != end_literals; l++)
        if (values[replacement = *l] >= 0)
        {
        FOUND_REPLACEMENT:
          debug_clause(c, "unwatching %s in", debug_literal(not_lit));
          literals[0] = other, literals[1] = replacement, *l = not_lit;
          watch_clause(replacement, other, ref);
          if (!initial_watched[index(replacement)]){
            initial_watched[index(replacement)] = true;
            heap.update_watched(index(replacement));
          }
          c->searched = l - literals;
          q -= 2;
          goto CONTINUE_WITH_NEXT_REFERENCE;
        }

      for (l = literals + 2; l != searched; l++)
        if (values[replacement = *l] >= 0)
          goto FOUND_REPLACEMENT;

      if (other_value < 0)
      {
        res = false;
        conflict = ref;
        debug_reference(conflict, "conflicting");
        break;
      }

      assign(other, ref);
    CONTINUE_WITH_NEXT_REFERENCE:;
    }

    while (p != end_not_lit_watches)
      *q++ = *p++;

    if (q != end_not_lit_watches)
      not_lit_watches.resize(q - begin_not_lit_watches);
  }

  if (!res)
    statistics.conflicts++;
  statistics.propagations += propagations;

  return res;
}

unsigned Solver::decide_variable()
{
  unsigned decision;
  while (values[decision = literal(*heap)])
    heap.pop();
  debug("highest scored unassigned variable %s with score %g",
        debug_literal(decision), heap.score(index(decision)));
  return decision;
}

unsigned Solver::decide_phase(unsigned decision)
{
  unsigned idx = index(decision);
  if (positive_polarity[idx] == 1 && negative_polarity[idx] == 0)
    return decision;
  if (positive_polarity[idx] == 0 && negative_polarity[idx] == 1)
    return negate(decision);
  signed char value = phases[idx].saved;
  if (!value)
    value = options.initial_phase ? 1 : -1;
  if (value < 0)
    decision = negate(decision);
  return decision;
}

void Solver::decide()
{
  level++;
  statistics.decisions++;

  unsigned decision = decide_variable();

  decision = decide_phase(decision);
  debug("decision %s", debug_literal(decision));

  assert(level == control.size());
  control.push_back(Frame(decision, trail.size()));

  assign(decision, DECISION_REASON);
}

void Solver::unassign(unsigned lit)
{
  debug("unassign %s", debug_literal(lit));
  unsigned not_lit = negate(lit);
  assert(values[lit] > 0);
  assert(values[not_lit] < 0);
  values[lit] = values[not_lit] = 0;
  unsigned idx = index(lit);
  if (!heap.contains(idx))
    heap.push(idx);
}

void Solver::backtrack(unsigned jump)
{
  assert(jump < level);
  assert(control.size() == level + 1);

  if (jump + 1 == level)
    debug("chronological backtracking to level %u", jump);
  else
    debug("non-chronological backjumping to level %u", jump);

  size_t begin = control[jump + 1].trail;
  size_t end = trail.size();
  size_t i = begin, j = i;

  while (i != end)
  {
    unsigned lit = trail[i++];
    unsigned lit_level = levels[index(lit)];
    if (lit_level > jump)
      unassign(lit);
    else
    {
      trail[j++] = lit;
      debug("keeping out-of-order assigned %s", debug_literal(lit));
    }
  }

  trail.resize(j);
  propagated = begin;
  control.resize(jump + 1);
  level = jump;
}

bool Solver::minimize(unsigned lit, unsigned depth)
{
  if (depth > 1000)
    return false;
  unsigned idx = index(lit);
  Flags &f = flags[idx];
  if (depth && f.seen)
    return true;
  if (f.poison)
    return false;
  if (f.minimizable)
    return true;
  Reference reason = reasons[idx];
  if (reason == DECISION_REASON)
    return false;
  if (reason == UNIT_REASON)
    return true;
  if (reason == BACKTRUE_REASON)
    return false;
  bool res = true;
  if (tagged(reason))
    res = minimize(untag(reason), depth + 1);
  else
  {
    Clause *clause = dereference(reason);
    for (auto other : *clause)
      if (other != lit && !minimize(other, depth + 1))
      {
        res = false;
        break;
      }
  }
  if (res)
  {
    f.minimizable = true;
    minimized.push_back(idx);
  }
  else
  {
    f.poison = true;
    poisoned.push_back(idx);
  }
  return res;
}

void Solver::minimize()
{
  statistics.deduced_literals += clause.size();
  auto begin = clause.begin();
  auto end = clause.end();
  auto p = begin, q = p;
  for (p = begin; p != end; p++)
    if (!minimize(*p))
      *q++ = *p;
  clause.resize(q - begin);
  statistics.learned_literals += clause.size();
  for (auto idx : minimized)
    flags[idx].minimizable = false;
  for (auto idx : poisoned)
    flags[idx].poison = false;
  minimized.clear();
  poisoned.clear();
}

struct larger_score
{
  const Heap &heap;
  larger_score(const Heap &h) : heap(h) {}
  bool operator()(unsigned a, unsigned b)
  {
    return (
      (heap.is_important(a) > heap.is_important(b)) ||
      (heap.is_important(a) == heap.is_important(b) && heap.score(a) > heap.score(b)) ||
      (heap.is_important(a) == heap.is_important(b) && heap.score(a) == heap.score(b) && heap.is_watched(a) && !heap.is_watched(b)) ||
      (heap.is_important(a) == heap.is_important(b) && heap.score(a) == heap.score(b) && heap.is_watched(a) == heap.is_watched(b) && b > a)
    );
  }
};

void Solver::bump_analyzed_variables()
{
  std::sort(analyzed.begin(), analyzed.end(), larger_score(heap));
  for (auto idx : analyzed)
  {
    heap.bump(idx);
    debug("bumped %s with new score %g",
          debug_variable(idx), heap.score(idx));
  }
  heap.bump();
}

// We factored out this code which analyzes a literal.  Without tagged
// reasons (and virtual binary clauses) we could keep it in the
// 'analyze_conflict' function below, but we want to reuse it for the tagged
// binary clause reason and also for analyzing the propagated literal.

inline void
Solver::analyze_literal(unsigned lit,
                        unsigned &open, unsigned &glue, unsigned &jump)
{
  unsigned idx = index(lit);
  Flags &f = flags[idx];
  if (f.seen)
    return;
  unsigned l = levels[idx];
  if (!l)
    return;
  debug("analyzing literal %s", debug_literal(lit));
  assert(values[lit] < 0);
  analyzed.push_back(idx);
  f.seen = true;
  if (l == level)
  {
    open++;
    return;
  }
  clause.push_back(lit);
  if (l > jump)
  {
    jump = l;
    if (clause.size() > 2)
      std::swap(clause[1], clause.back());
  }
  Frame &frame = control[l];
  if (frame.pulled)
    return;
  frame.pulled = true;
  pulled.push_back(l);
  glue++;
}

bool Solver::analyze_conflict(unsigned falsified, Reference reason)
{
  assert(clause.empty());
  assert(analyzed.empty());
  assert(reason != DECISION_REASON);
  assert(reason != UNIT_REASON);
  assert(reason != BACKTRUE_REASON);
  unsigned conflict_level;

 
  {
    // Compute maximum assignment level of all conflicting literals and
    // backtrack to that decision level first.  This is necessary to support
    // chronological backtracking and thus out-of-order assignments.
    // With chronological backtracking and out-of-order assignments a
    // conflicting clause could become forcing.

    unsigned forced_unit = INVALID;

    if (tagged(reason))
    {
      unsigned other = untag(reason);
      unsigned other_level = levels[index(other)];
      unsigned falsified_level = levels[index(falsified)];
      conflict_level = std::max(other_level, falsified_level);
      if (other_level < conflict_level)
        forced_unit = falsified;
    }
    else
    {
      conflict_level = 0;
      unsigned literals_on_conflict_level = 0;
      Clause *clause = dereference(reason);
      for (auto lit : *clause)
      {
        unsigned lit_level = levels[index(lit)];
        if (lit_level > conflict_level)
        {
          conflict_level = lit_level;
          literals_on_conflict_level = 1;
          forced_unit = lit;
        }
        else if (lit_level == conflict_level)
          literals_on_conflict_level++;
      }
      if (literals_on_conflict_level > 1)
        forced_unit = INVALID;
      else
        assert(forced_unit != INVALID);
    }

    if (conflict_level == level)
      debug("conflict level %u matches decision level", conflict_level);
    else
    {
      assert(conflict_level < level);
      debug("conflict level %u smaller decision level", conflict_level);
      backtrack(conflict_level);
    }
  }

  // Here starts the standard conflict analysis which is independent of the
  // chronological backtracking code above.

  if (!level)
  { // Conflict on root level.
    inconsistent = true;
    trace('a');
    return false;
  }

  statistics.analyzed++;
  
  clause.push_back(0); // Make room for 'not_uip'.

  unsigned open = 0; // Remaining on current level.
  unsigned glue = 0; // Number of pulled in levels.
  unsigned jump = 0; // Backjump level.

  size_t t = trail.size(); // Walk trail backwards.
  unsigned uip = negate(falsified);

  if (level > 1)
    update_phases();

  if (tagged(reason))
    analyze_literal(falsified, open, glue, jump);

  unsigned decision_uip = control[conflict_level].decision;
  for (;;)
  {

    assert(reason != DECISION_REASON);
    assert(reason != UNIT_REASON);

    if (reason == BACKTRUE_REASON)
    {
      unsigned i;
      for (i = 1; i <= conflict_level; i++)
      {
        debug("analyzing control %d vs %d ", control[i].decision, negate(control[i].decision));
        analyze_literal(negate(control[i].decision), open, glue, jump);
      }
    }
    else if (tagged(reason))
    {
      unsigned other = untag(reason);
      debug_binary(uip, other, "analyzing");
      analyze_literal(other, open, glue, jump);
    }
    else
    {
      Clause *clause = dereference(reason);
      debug_clause(clause, "analyzing");
      bump_clause(clause);
      for (auto lit : *clause)
        analyze_literal(lit, open, glue, jump);
    }

    unsigned uip_idx;
    do
    {
      assert(t);
      uip = trail[--t];
      uip_idx = index(uip);
    } while (!flags[uip_idx].seen || levels[uip_idx] != level);

    assert(open);

    if (uip == decision_uip)
      break;

    reason = reasons[uip_idx];
  }

  assert(glue == pulled.size());

  
  debug("decision UIP %s", debug_literal(uip));
  debug("glucose level %u (LBD/glue)", glue);
  debug("original jump level %u", jump);

  slow_glue_average += glue;
  fast_glue_average += glue;

  assert(uip != INVALID);
  unsigned not_uip = negate(uip);
  clause[0] = not_uip;
  debug_vector(clause, "last UIP");

  minimize();

  size_t size = clause.size();
  iterating = (size == 1);
  assert(size > 0);

  assert(level > jump);
  if (level - jump > options.backjump_limit)
  {
#ifdef LOGGING
    if (options.backjump_limit)
      debug("forcing chronological backtracking "
            "instead of jumping over %u level",
            level - jump);
    else
      debug("always forcing chronological backtracking");
#endif
    jump = level - 1;
  }

  if (jump + 1 == level)
    statistics.chronological++;
  std::vector<unsigned> reinserting_lits;
  std::vector<unsigned> reinserting_reasons;
  if (size == 1){
    for (auto lit: trail){
      if (index(lit) != index(not_uip) && important[index(lit)] != 0 && (reasons[index(lit)] == DECISION_REASON || (reasons[index(lit)] == BACKTRUE_REASON && levels[index(lit)] < levels[index(not_uip)]))){
        reinserting_lits.push_back(lit);
        reinserting_reasons.push_back(reasons[index(lit)]);
      }
    }
    backtrack(0);
  }
  else  
    backtrack(jump);

  trace('a');
  Reference learned;
  if (size == 1)
    learned = UNIT_REASON;
  else if (size == 2)
    learned = new_binary();
  else
    learned = new_clause(true, glue);

  assign(not_uip, learned);


  if (size == 1){
    for (unsigned i = 0; i < reinserting_lits.size(); i++){
      if (reinserting_reasons[i] == DECISION_REASON){
        level++;
        control.push_back(Frame(reinserting_lits[i], trail.size()));
      }
      assign(reinserting_lits[i], reinserting_reasons[i]);
    }
  }

  aggressive_level_limit = level;
  
  bump_analyzed_variables();

  for (auto idx : analyzed)
    flags[idx].seen = false;
  analyzed.clear();

  for (auto l : pulled)
    control[l].pulled = false;
  pulled.clear();

  clause.clear(); 

  return true;
}

void Solver::report(char type)
{
  return;
  if (verbosity < 0)
    return;
  if (!(statistics.reported++ & 15))
    printf("c\n"
           "c   seconds,MB,reductions,restarts,conflicts,redundant,glue,remaining variables\n"
           "c\n");
  unsigned remaining = variables - statistics.fixed;
  Statistics &s = statistics;
  double MB = current_resident_set_size() / (double)(1 << 20);
  printf("c %c %7.2f %2.0f %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64
         " %.2f %u %.0f%%\n",
         type, process_time(), MB,
         s.reductions, s.restarts, s.conflicts, s.redundant_clauses,
         (double)slow_glue_average, remaining,
         variables ? 100.0 * remaining / variables : 0);
  fflush(stdout);
}

void Solver::mark_reason_clauses(bool mark)
{
  size_t begin = control[1].trail;
  size_t end = trail.size();
  for (size_t i = begin; i != end; i++)
  {
    unsigned lit = trail[i];
    unsigned idx = index(lit);
    Reference reason = reasons[idx];
    if (reason == DECISION_REASON)
      continue;
    if (reason == UNIT_REASON)
      continue;
    if (reason == BACKTRUE_REASON)
      continue;
    if (tagged(reason))
      continue;
    Clause *clause = dereference(reason);
    assert(clause->reason != mark);
    clause->reason = mark;
  }
}

void Solver::map_reasons(unsigned first_garbage,
                         std::vector<unsigned> &map)
{
  size_t begin = control[1].trail;
  size_t end = trail.size();
  for (size_t i = begin; i != end; i++)
  {
    unsigned lit = trail[i];
    unsigned idx = index(lit);
    Reference reason = reasons[idx];
    if (reason == DECISION_REASON)
      continue;
    if (reason == UNIT_REASON)
      continue;
    if (tagged(reason))
      continue;
    unsigned pos = untag(reason);
    if (pos < first_garbage)
      continue;
    unsigned delta = pos - first_garbage;
    assert(delta < map.size());
    unsigned mapped = map[delta];
    if (mapped != INVALID)
      reasons[idx] = tag_large(mapped);
  }
}

void Solver::gather_reduce_candidates(std::vector<unsigned> &candidates)
{
  mark_reason_clauses(true);
  Clauses::iterator begin = clauses[first_redundant], end = clauses.end();
  for (Clauses::iterator i = begin; i != end; ++i)
  {
    Clause *clause = *i;
    if (clause->reason)
      continue;
    if (!clause->redundant)
      continue;
    if (clause->glue <= options.tier1_max_glue)
      continue;
    if (clause->used)
    {
      clause->used--;
      continue;
    }
    unsigned pos = clauses.position(clause);
    candidates.push_back(pos);
    debug_clause(clause, "gathered reduce candidate");
  }
  mark_reason_clauses(false);
  message(1, "gathered %zu redundant reduce candidates", candidates.size());
}

struct less_relevant
{
  Solver *solver;
  less_relevant(Solver *s) : solver(s) {}
  bool operator()(unsigned p, unsigned q)
  {
    Clause *c = solver->clauses[p];
    Clause *d = solver->clauses[q];
    if (c->glue < d->glue)
      return false;
    if (c->glue > d->glue)
      return true;
    if (c->size < d->size)
      return false;
    if (c->size > d->size)
      return true;
    return p < q;
  }
};

unsigned Solver::mark_garbage(std::vector<unsigned> &candidates)
{
  unsigned target = options.reduce_fraction / 100.0 * candidates.size();
  unsigned marked = 0, first_garbage = INVALID;
  for (auto pos : candidates)
  {
    if (marked == target)
      break;
    auto clause = clauses[pos];
    if (pos < first_garbage)
      first_garbage = pos;
    debug_clause(clause, "marked garbage less relevant");
    clause->garbage = true;
    marked++;
  }
  message(1, "marked %u less relevant clauses as garbage", marked);
  return first_garbage;
}

bool Solver::satisfied_clause(Clause *clause)
{
  for (auto lit : *clause)
    if (fixed(lit) > 0)
      return true;
  return false;
}

void Solver::mark_satisfied_clauses_as_garbage_too(unsigned &first_garbage)
{
  unsigned marked = 0;
  for (auto clause : clauses)
  {
    if (clause->garbage)
      continue;
    if (!satisfied_clause(clause))
      continue;
    unsigned pos = clauses.position(clause);
    if (pos < first_garbage)
      first_garbage = pos;
    debug_clause(clause, "marked satisfied");
    clause->garbage = true;
    marked++;
  }
  message(1, "marked %u satisfied clauses as garbage", marked);
}

void Solver::flush_garbage_watches(unsigned first_garbage, bool simplify,
                                   std::vector<unsigned> &map)
{
  size_t flushed = 0;

  for (auto lit : literals)
  {
    auto &watches = matrix[lit];
    auto begin = watches.begin();
    auto end = watches.end();
    auto p = begin, q = p;
    signed char lit_fixed = simplify ? fixed(lit) : 0;
    while (p != end)
    {
      Reference first = *p++;
      if (!tagged(first))
      {
        Reference second = *p++;
        unsigned pos = untag(second);
        if (first_garbage <= pos)
        {
          unsigned delta = pos - first_garbage;
          assert(delta <= map.size());
          unsigned mapped = map[delta];
          if (mapped == INVALID)
            continue;
          *q++ = first, *q++ = tag_large(mapped);
        }
        else
          *q++ = first, *q++ = second;
      }
      else if (simplify)
      {
        unsigned other = untag(first);
        auto other_fixed = fixed(other);
        if (lit_fixed <= 0 && other_fixed <= 0)
          *q++ = first;
        else if (lit > other)
        {
          debug_binary(lit, other, "flushing satisfied");
          trace('d', lit, other);
          flushed++;
        }
      }
      else
        *q++ = first;
    }
    watches.resize(q - begin);
    watches.shrink_to_fit();
  }

  if (simplify)
    message(1, "simplified and flushed %zu binary clauses", flushed);
  else
  {
    assert(!flushed);
    message(2, "no binary clause simplified");
  }
}

void Solver::recycle_garbage_clauses(unsigned first_garbage,
                                     std::vector<unsigned> &map)
{
  Clause *begin = clauses[first_garbage];
  Clause *end = *clauses.end();

  assert(clauses.size() <= UINT_MAX);
  assert((unsigned *)end - (unsigned *)begin <= UINT_MAX);

  unsigned traversed = 0, mapped = first_garbage;
  unsigned copied = 0, recycled = 0;

  Clause *src = begin, *dst = begin;

  if (first_garbage < first_redundant)
    first_redundant = INVALID;

  while (src != end)
  {

    assert(clauses[first_garbage + traversed] == src);

    unsigned size = src->size;
    size_t bytes = Clause::bytes(size);
    unsigned words = bytes / 4;
    assert(Clause::words(size) == words);

    while (map.size() < traversed)
      map.push_back(INVALID);

    assert(map.size() == traversed);

    if (src->garbage)
    {
      trace('d', size, src->literals);
      statistics.recycled_clauses++;
      delete_clause(src);
      map.push_back(INVALID);
      assert(map[traversed] == INVALID);
      recycled++;
    }
    else
    {
      memmove(dst, src, bytes);
      map.push_back(mapped);
      assert(map[traversed] == mapped);
      debug_clause(dst, "copied and keeping");
      dst = (Clause *)((char *)dst + bytes);
      if (first_redundant == INVALID && dst->redundant)
        first_redundant = mapped;
      mapped += words;
      copied++;
    }

    src = (Clause *)((char *)src + bytes);
    traversed += words;
  }

  clauses.resize(mapped);

  message(2, "copied %u large non-garbage clauses", copied);
  message(1, "recycled %u large clauses", recycled);
}

bool Solver::reducing()
{
  if (!options.reduce)
    return false;
  return limits.reduce < statistics.conflicts;
}

unsigned Solver::gather_and_mark_reduced_clauses()
{
  std::vector<unsigned> candidates;

  if (first_redundant != INVALID)
    gather_reduce_candidates(candidates);

  unsigned first_garbage = INVALID;

  if (!candidates.empty())
  {
    std::sort(candidates.begin(),
              candidates.end(), less_relevant(this));
    first_garbage = mark_garbage(candidates);
  }

  if (fixed_at_reduce < statistics.fixed)
  {
    statistics.simplifications++;
    mark_satisfied_clauses_as_garbage_too(first_garbage);
  }

  debug("first garbage clause at position %u", first_garbage);
  return first_garbage;
}

void Solver::recycle_clauses_and_flush_watches(unsigned first_garbage)
{
  std::vector<unsigned> map;
  recycle_garbage_clauses(first_garbage, map);
  bool simplify = fixed_at_reduce < statistics.fixed;
  flush_garbage_watches(first_garbage, simplify, map);
  map_reasons(first_garbage, map);
}

void Solver::reduce()
{
  start("reduce");
  Statistics &s = statistics;
  s.reductions++;

  message(1, "clause data-base reduction %" PRIu64 "", s.reductions);
  unsigned first_garbage = gather_and_mark_reduced_clauses();
  if (first_garbage != INVALID)
    recycle_clauses_and_flush_watches(first_garbage);

  fixed_at_reduce = s.fixed;
  uint64_t base = options.reduce_interval;
  uint64_t scaled = base * s.reductions;
  limits.reduce = s.conflicts + scaled;
  message(2, "new reduce limit at %" PRIu64 " after %" PRIu64 " conflicts", limits.reduce, scaled);

  stop("reduce");
  report('-');
}

bool Solver::restarting()
{
  if (!options.restart)
    return false;
  if (!level)
    return false;
  Limits &l = limits;
  Statistics &s = statistics;
  if (l.restart >= s.conflicts)
    return false;
  unsigned percent = options.restart_margin;
  double scale = 1.0 + percent / 100.0;
  double slow = slow_glue_average;
  double fast = fast_glue_average;
  double padded = scale * slow;
  bool res = padded < fast;
  message(3, "%srestarting at conflict %" PRIu64 " as %.2f = %.2f * slow %.2f %s fast glue %.2f",
          res ? "" : "not-", s.conflicts, padded, scale, slow,
          res ? "<" : ">=", fast);
  return res;
}

void Solver::set_restart_limit()
{
  Statistics &s = statistics;
  uint64_t u = reluctant.u, v = reluctant.v;
  if ((u & -u) == v)
    u++, v = 1;
  else
    v *= 2;
  uint64_t interval = v * options.restart_interval;
  reluctant.u = u, reluctant.v = v;
  limits.restart = s.conflicts + interval;
  message(2, "new restart limit at %" PRIu64 " after %" PRIu64 " conflicts", limits.restart, interval);
}

void Solver::restart()
{
  Statistics &s = statistics;
  s.restarts++;
  message(2, "restart %" PRIu64 " triggered after %" PRIu64 " conflicts",
          s.restarts, s.conflicts);
  backtrack(0);
  set_restart_limit();
}

void Solver::update_phases(void)
{
  assert(level);
  size_t consistently_assigned = control[level].trail;
  bool update_best_assignment = (consistently_assigned > best_assigned);
  bool update_target_assignment = (consistently_assigned > target_assigned);
  if (!update_best_assignment && !update_target_assignment)
    return;
  for (size_t i = control[1].trail; i != consistently_assigned; i++)
  {
    unsigned lit = trail[i];
    Phase &phase = phases[index(lit)];
    signed char value = sign(lit) ? -1 : 1;
    if (update_best_assignment)
      phase.best = value;
    if (update_target_assignment)
      phase.target = value;
  }
  if (update_best_assignment)
  {
    best_assigned = consistently_assigned;
    message(2, "updated best assigned to %u", best_assigned);
  }
  if (update_target_assignment)
  {
    target_assigned = consistently_assigned;
    message(2, "updated target assigned to %u", target_assigned);
  }
}

bool Solver::rephasing()
{
  if (!options.rephase)
    return false;
  return statistics.conflicts > limits.rephase;
}

void Solver::rephase_original()
{
  signed char original_phase = options.initial_phase ? 1 : -1;
  for (auto idx : variables)
    phases[idx].target = phases[idx].saved = original_phase;
  report('O');
}

void Solver::rephase_inverted()
{
  signed char inverted_phase = options.initial_phase ? -1 : 1;
  for (auto idx : variables)
    phases[idx].target = phases[idx].saved = inverted_phase;
  report('I');
}

void Solver::rephase_best()
{
  signed char best;
  for (auto idx : variables)
    if ((best = phases[idx].best))
      phases[idx].target = phases[idx].saved = best;
  best_assigned = 0;
  report('B');
}

static double nlogpow(uint64_t count, unsigned pow)
{
  assert(count), assert(pow);
  double res = count, factor = log10(count + 9);
  for (unsigned i = 0; i != pow; i++)
    res *= factor;
  return res;
}

void Solver::rephase()
{
  Statistics &s = statistics;
  s.rephased++;
  message(1, "rephasing %" PRIu64 " by resetting phases", s.rephased);

  switch (s.rephased % 4)
  {
  case 0:
    rephase_original();
    break;
  case 1:
    rephase_best();
    break;
  case 2:
    rephase_inverted();
    break;
  case 3:
    rephase_best();
    break;
  }

  target_assigned = 0;

  uint64_t base = options.rephase_interval;
  uint64_t scaled = nlogpow(s.rephased, 3) * base;
  limits.rephase = s.conflicts + scaled;
  message(2, "new rephase limit at %" PRIu64 " after %" PRIu64 " conflicts", limits.reduce, scaled);
}


bool Solver::block_chrono(unsigned M)
{
  if (level == 0 || M == 0)
  {
    return true;
  }

  unsigned jump = M - 1;
  unsigned lit = control[M].decision;

  backtrack(jump);
  aggressive_level_limit = level;
  
  if (lit != lit_to_flip){
    assign(negate(lit), BACKTRUE_REASON);
    // assign(negate(lit_to_flip), BACKTRUE_REASON);
  }
  else
    assign(negate(lit), BACKTRUE_REASON);
  return false;
}


unsigned Solver::get_aggressive_projected_implicant_lifting_bis(){
  std::vector<unsigned> copy_trail;
  std::vector<unsigned> copy_trail2;

  if (level == 0)
  {
    return level;
  }

  for (auto lit : trail){
    copy_trail.push_back(lit);
    copy_trail2.push_back(lit);
  }

  std::set<unsigned> trail_minimal_set;

  std::vector<References> W(literals);
  std::vector<unsigned> N(added_original_non_binary_clauses);
  
  unsigned size = trail.size();

  unsigned base = 0;

  for (unsigned i = 0; i < added_original_non_binary_clauses; i++){
    Clause *c = clauses[base];
    if (c->redundant)
      continue;
    for (int j = 0; j < c->size; j++){
      unsigned lit = c->literals[j];
      if (values[lit] > 0){
        W[lit].push_back(i);
        N[i]++;
      }
    }
    base += Clause::bytes(c->size) / 4;
  }

  // TRY TO REMOVE ONLY IMPORTANT
  for(int idx = size - 1; idx >= 0; idx--){

    unsigned lit = trail[idx];

    if (levels[index(lit)] <= aggressive_level_limit){
      trail_minimal_set.insert(lit);
      continue;
    }

    if (projected && important[index(lit)] == 0){
      trail_minimal_set.insert(lit);
      continue;
    }

    if (reasons[index(lit)] == BACKTRUE_REASON){
      trail_minimal_set.insert(lit);
      continue;
    }

    auto &not_lit_watches = matrix[lit];
    auto begin_not_lit_watches = not_lit_watches.begin();
    auto end_not_lit_watches = not_lit_watches.end();
    auto p = begin_not_lit_watches, q = p;
    bool one_binary_require_lit = false;

    while (p != end_not_lit_watches)
    {

      // Get the watched clause and the blocking literal, plus its value
      Reference watch = *q++ = *p++;

      if (tagged(watch))
      {
        // Binary clause: other is for sure in implicant
        unsigned other = untag(watch);
        signed char other_value = values[other];
        if (other_value <= 0 || std::find(copy_trail.begin(), copy_trail.end(), other) == copy_trail.end())
        {
          trail_minimal_set.insert(lit);
          one_binary_require_lit = true;
          continue;
        }
      }
      else
      {
        // We skip it since we use the De Harbe algorithm for normal clauses
        *q++;
        *p++;
      }
    }

    if(one_binary_require_lit){
      continue;
    }
    
    // TODO: GET ERROR HERE!
    bool one_clause_require_lit = false;
    // De Harbe algorithm
    for (auto ref : W[lit]){
      if (N[ref] == 1){
        trail_minimal_set.insert(lit);
        one_clause_require_lit = true;
        break;
      } 
    }

    if (!one_clause_require_lit){
      copy_trail.erase(copy_trail.begin() + idx);
      for (auto ref : W[lit]){
        N[ref]--;
      }
    }

  }

  backtrack(aggressive_level_limit);

  unsigned falsified;
  Reference conflict;
  for (unsigned i = 0; i < copy_trail2.size(); i++){
    unsigned lit = copy_trail2[i];
    if (levels[index(lit)] <= aggressive_level_limit){
      continue;
    }
    if (trail_minimal_set.find(lit) != trail_minimal_set.end() && !values[lit] && !values[negate(lit)]){
      if (projected && important[index(lit)] == 0)
        continue;
      if (reasons[index(lit)] == BACKTRUE_REASON){
        assign(lit, BACKTRUE_REASON);
        continue;
      }
      level++;
      control.push_back(Frame(copy_trail2[i], trail.size()));
      assign(lit, DECISION_REASON);
      bool t = propagate(falsified, conflict);
    }
  }

  return level;
}



unsigned Solver::get_aggressive_implicant_lifting(){
  std::vector<unsigned> copy_trail;
  std::vector<unsigned> copy_trail2;
  std::vector<Reference> copy_reasons;

  for (auto lit : trail){
    copy_trail.push_back(lit);
    copy_trail2.push_back(lit);
    copy_reasons.push_back(reasons[index(lit)]);
  }

  std::set<unsigned> trail_minimal_set;

  std::vector<References> W(literals);
  std::vector<unsigned> N(added_original_non_binary_clauses);
  
  unsigned size = trail.size();

  unsigned base = 0;

  for (unsigned i = 0; i < added_original_non_binary_clauses; i++){
    Clause *c = clauses[base];
    if (c->redundant)
      continue;
    for (int j = 0; j < c->size; j++){
      unsigned lit = c->literals[j];
      if (values[lit] > 0){
        W[lit].push_back(i);
        N[i]++;
      }
    }
    base += Clause::bytes(c->size) / 4;
  }

  // TRY TO REMOVE ONLY IMPORTANT
  for(int idx = size - 1; idx >= 0; idx--){

    unsigned lit = trail[idx];

    copy_trail.erase(copy_trail.begin() + idx);

    auto &not_lit_watches = matrix[lit];
    auto begin_not_lit_watches = not_lit_watches.begin();
    auto end_not_lit_watches = not_lit_watches.end();
    auto p = begin_not_lit_watches, q = p;
    bool one_binary_require_lit = false;

    while (p != end_not_lit_watches)
    {

      // Get the watched clause and the blocking literal, plus its value
      Reference watch = *q++ = *p++;

      if (tagged(watch))
      {
        // Binary clause: other is for sure in implicant
        unsigned other = untag(watch);
        signed char other_value = values[other];
        if (other_value <= 0 || std::find(copy_trail.begin(), copy_trail.end(), other) == copy_trail.end())
        {
          trail_minimal_set.insert(lit);
          one_binary_require_lit = true;
          continue;
        }
      }
      else
      {
        // We skip it since we use the De Harbe algorithm for normal clauses
        *q++;
        *p++;
      }
    }

    if(one_binary_require_lit){
      continue;
    }
    
    // TODO: GET ERROR HERE!
    bool one_clause_require_lit = false;
    // De Harbe algorithm
    for (auto ref : W[lit]){
      if (N[ref] == 1){
        trail_minimal_set.insert(lit);
        one_clause_require_lit = true;
        break;
      } 
    }

    if (!one_clause_require_lit){
      for (auto ref : W[lit]){
        N[ref]--;
      }
    }

  }

  // TODO: now remove element from trail that do not belong to trail_minimal_set
  backtrack(0);

  unsigned falsified;
  Reference conflict;
  for (unsigned i = 0; i < copy_trail2.size(); i++){
    unsigned lit = copy_trail2[i];
    if (trail_minimal_set.find(lit) != trail_minimal_set.end() && !values[lit] && !values[negate(lit)]){
      if (projected && important[index(lit)] == 0)
        continue;
      level++;
      #ifdef MATHSAT
        if (has_event) { msat_dpll_callback_notify_new_level(event); }
      #endif
      control.push_back(Frame(copy_trail2[i], trail.size()));
      assign(lit, DECISION_REASON);
      bool t = propagate(falsified, conflict);
    }
  }

  return level;
}


unsigned Solver::get_implicant_lifting()
{
  // std::set<unsigned> removed;
  unsigned M = 0;

  unsigned size = trail.size();
  for(int idx = size - 1; idx >= 0; idx--)
  {
    unsigned lit = trail[idx];

    if (reasons[index(lit)] != DECISION_REASON)
    {
      if (levels[index(lit)] > M)
        if (!projected || important[index(lit)] == 1)
          M = levels[index(lit)];
      continue;
    }
    
    if (matrix[lit].size() == 0)
      continue;
    
    auto &not_lit_watches = matrix[lit];
    auto begin_not_lit_watches = not_lit_watches.begin();
    auto end_not_lit_watches = not_lit_watches.end();
    auto p = begin_not_lit_watches, q = p;
    while (p != end_not_lit_watches)
    {
      
      // Get the watched clause and the blocking literal, plus its value
      Reference watch = *q++ = *p++;
      
      if (tagged(watch))
      {
        // Binary clause: other is for sure in implicant
        unsigned other = untag(watch);
        signed char other_value = values[other];
        if (other_value <= 0 || position[index(other)] > idx)
        {
          if (levels[index(lit)] > M)
            if (!projected || important[index(lit)] == 1)
              M = levels[index(lit)];
            // break;
        }
        else
        {
          if (levels[index(other)] > M)
            if (!projected || important[index(other)] == 1)
              M = levels[index(other)];
            continue;
        }
      }
      else
      {
        unsigned blocking = untag(*(p - 1));
        Reference ref = *q++ = *p++;
        Clause *c = dereference(ref);
        signed char blocking_value = values[blocking];
        
        if (c->redundant)
          continue;

        unsigned other = c->literals[0] == lit ? c->literals[1] : c->literals[0];

        if (values[other] > 0)
        {
          q[-2] = tag_large(other);
          blocking = other;
          blocking_value = values[blocking];
        }
          
        if (blocking_value <= 0 || position[index(blocking)] > idx)
        {
          if (levels[index(lit)] > M)
            if (!projected || important[index(lit)] == 1)
              M = levels[index(lit)];
            // break;
        }
        else
        {
          if (levels[index(blocking)] > M)
            if (!projected || important[index(blocking)] == 1)
              M = levels[index(blocking)];
            continue;
        }
      }
    }
  }

  return M;
}


unsigned Solver::get_aggressive_projected_implicant_lifting(){
  std::vector<unsigned> copy_trail;
  std::vector<unsigned> copy_trail2;

  if (level == 0)
  {
    return level;
  }

  for (auto lit : trail){
    copy_trail.push_back(lit);
    copy_trail2.push_back(lit);
  }

  std::set<unsigned> trail_minimal_set;

  std::vector<References> W(literals);
  std::vector<unsigned> N(added_original_non_binary_clauses);
  
  unsigned size = trail.size();

  unsigned base = 0;

  for (unsigned i = 0; i < added_original_non_binary_clauses; i++){
    Clause *c = clauses[base];
    if (c->redundant)
      continue;
    for (int j = 0; j < c->size; j++){
      unsigned lit = c->literals[j];
      if (values[lit] > 0){
        W[lit].push_back(i);
        N[i]++;
      }
    }
    base += Clause::bytes(c->size) / 4;
  }

  // TRY TO REMOVE ONLY IMPORTANT
  for(int idx = size - 1; idx >= 0; idx--){

    unsigned lit = trail[idx];

    if (important[index(lit)] == 0){
      trail_minimal_set.insert(lit);
      continue;
    }

    copy_trail.erase(copy_trail.begin() + idx);

    auto &not_lit_watches = matrix[lit];
    auto begin_not_lit_watches = not_lit_watches.begin();
    auto end_not_lit_watches = not_lit_watches.end();
    auto p = begin_not_lit_watches, q = p;
    bool one_binary_require_lit = false;

    while (p != end_not_lit_watches)
    {

      // Get the watched clause and the blocking literal, plus its value
      Reference watch = *q++ = *p++;

      if (tagged(watch))
      {
        // Binary clause: other is for sure in implicant
        unsigned other = untag(watch);
        signed char other_value = values[other];
        if (other_value <= 0 || std::find(copy_trail.begin(), copy_trail.end(), other) == copy_trail.end())
        {
          trail_minimal_set.insert(lit);
          one_binary_require_lit = true;
          continue;
        }
      }
      else
      {
        // We skip it since we use the De Harbe algorithm for normal clauses
        *q++;
        *p++;
      }
    }

    if(one_binary_require_lit){
      continue;
    }
    
    // TODO: GET ERROR HERE!
    bool one_clause_require_lit = false;
    // De Harbe algorithm
    for (auto ref : W[lit]){
      if (N[ref] == 1){
        trail_minimal_set.insert(lit);
        one_clause_require_lit = true;
        break;
      } 
    }

    if (!one_clause_require_lit){
      for (auto ref : W[lit]){
        N[ref]--;
      }
    }

  }

  // TODO: now remove element from trail that do not belong to trail_minimal_set
  backtrack(0);

  unsigned falsified;
  Reference conflict;
  for (unsigned i = 0; i < copy_trail2.size(); i++){
    unsigned lit = copy_trail2[i];
    if (trail_minimal_set.find(lit) != trail_minimal_set.end() && !values[lit] && !values[negate(lit)]){
      if (projected && important[index(lit)] == 0)
        continue;
      level++;
      #ifdef MATHSAT
        if (has_event) { msat_dpll_callback_notify_new_level(event); }
      #endif
      control.push_back(Frame(copy_trail2[i], trail.size()));
      assign(lit, DECISION_REASON);
      bool t = propagate(falsified, conflict);
    }
  }

  return level;
}


unsigned Solver::get_projection_lifting()
{
  unsigned M = 0;

  unsigned size = trail.size();
  for(int idx = size - 1; idx >= 0; idx--)
  {
    unsigned lit = trail[idx];

    // if index(lit) is not in important
    if (projected && important[index(lit)] == 0)
    {
      continue;
    }
    else{
      if (levels[index(lit)] > M){
        M = levels[index(lit)];
        lit_to_flip = lit;
      }
      continue;
    }
    
  }
  return M;
}


mpz_class Solver::get_models_covered_by_assignment(unsigned length)
{
  mpz_class models_per_assignment = 0;
  if (!projected)
    mpz_ui_pow_ui(models_per_assignment.get_mpz_t(), 2, variables - length);
  else
    mpz_ui_pow_ui(models_per_assignment.get_mpz_t(), 2, projected_count - length);
  return models_per_assignment;
}

int Solver::to_dimacs(unsigned lit){
  int idx = index(lit) + 1;
  if (lit % 2 == 1)
    idx = -idx;
  return idx;
}

State Solver::cdcl()
{
  options.backjump_limit = 0;
  mpz_class model_count = 0;
  mpz_class assignments_count = 0;

  std::ostream* outfile = nullptr; // Pointer to the output stream

  if (options.output_file.value == 0) {
    // If output_file is empty, print to stdout
    outfile = &std::cout;
  } else {
    // Otherwise, open the specified file and write to it
    std::ofstream* file_stream = new std::ofstream("output.txt");
    if (file_stream->is_open()) {
      outfile = file_stream;
    } else {
      std::cerr << "Error: Could not open file " << options.output_file.value << " for writing." << std::endl;
      exit(-2);
    }
  }


  // std::ofstream outfile("../../../../data/datasets/output_projected.txt");
  //std::ofstream outfile("output_projected.txt");


  for (;;)
  {
    unsigned falsified;
    Reference conflict;
    if (!propagate(falsified, conflict))
    {
      if (!analyze_conflict(falsified, conflict))
      {
        if (true)
          {
            printf("s MODEL COUNT\n");
            gmp_printf("%Zd\n", model_count);
          }
        return UNSATISFIABLE;
      }
    }
    else
    {
      if (iterating)
        iterating = false, report('i');
      if (trail.size() == variables)
      {
        debug("Found model");
        options.restart = 0;
        assignments_count++;
        statistics.n_assignments_partial++;

        unsigned M;
        bool endTrue;

        if (!enumerate_total){
          M = get_aggressive_projected_implicant_lifting_bis();
          M = get_projection_lifting();
        }
        else{
          if (projected){
            M = get_projection_lifting();
          }
          else{
            M = level;
          }
        }

        debug("CURRENT M: %d", M);
        
        if (level > M)
          backtrack(M);
        
        if (!model_count_flag && !assignment_count_flag)
        {
          for (unsigned i = 0; i < trail.size() ; i++)
          {
            unsigned lit = trail[i];
            int idx = to_dimacs(lit);

            // If not important in projected, skip it
            if (projected && important[index(lit)] == 0)
              continue;

            (*outfile) << idx;
            (*outfile) << " ";
          }
          (*outfile) << "\n";
          outfile->flush();
        }

        if(!projected)
          model_count += get_models_covered_by_assignment(trail.size());
        else
        {
          unsigned projected_counter = 0;

          for (unsigned i = 0; i < trail.size(); i++)
          {
            unsigned lit = trail[i];
            if (important[index(lit)] == 1)
              projected_counter += 1;
          }

          model_count += get_models_covered_by_assignment(projected_counter);  
        }    
        endTrue = block_chrono(M);
        if (inconsistent || endTrue)
        {
          if (true)
          {
            printf("s MODEL COUNT\n");
            gmp_printf("%Zd\n", model_count);
          }
          return UNSATISFIABLE;
        }
        continue;
      }
      else if (terminated)
        return UNKNOWN;
      else if (restarting())
        restart();
      else if (rephasing())
        rephase();
      else if (reducing())
        reduce();
      else
        decide();
    }
  }
}

void Solver::initialize_internal_solving()
{
  assert(!initialized);

  limits.reduce = options.reduce_interval;
  message(1, "initial reduce limit of %" PRIu64 " conflicts",
          limits.reduce);

  set_restart_limit();

  limits.rephase = options.rephase_interval;
  message(1, "initial rephase limit of %" PRIu64 " conflicts",
          limits.rephase);

  double alpha = options.vsids_decay / 1000.0, decay = 1 - alpha;
  if (!added_original_clauses)
  {
    printf("s MODEL COUNT\n");
    gmp_printf("%Zd\n", get_models_covered_by_assignment(0));
    exit(20);
  }
  heap.resize(variables, DLCS, initial_watched, important);
  heap.set_decay(decay);
  message(1, "setting VSIDS decay to %g "
             "(alpha specified as %g = %u/1000)",
          decay, alpha, (unsigned)options.vsids_decay);        

  initialized = true;
}

State Solver::solve_internally()
{
  State res;
  if (inconsistent)
  {
    printf("s MODEL COUNT\n0\n");
    return UNSATISFIABLE;
  }
  else
  {
    if (!initialized)
      initialize_internal_solving();
    report('*');
    start("search");
    res = cdcl();
    stop("search");
    if (verbosity >= 0)
    {
      report(res == 10 ? '1' : res == 20 ? '0'
                                         : '?');
      fputs("c\n", stdout), fflush(stdout);
    }
  }
  return res;
}

int main(int argc, char **argv)
{
  Solver solver;
  int res = solver.main(argc, argv);
  return res;
}
