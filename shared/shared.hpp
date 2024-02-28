#ifndef _Shared_hpp_included
#define _Shared_hpp_included

/*------------------------------------------------------------------------*/

#include <cassert>
#include <cstdio>
#include <climits>

/*------------------------------------------------------------------------*/

#include <vector>
#include <set>
#include "heap.hpp"

/*------------------------------------------------------------------------*/

// The 'Shared' class below contains code shared among all solver variants
// and before we have helper data-structures including a generic interface.

/*------------------------------------------------------------------------*/

// These constants are used as internal state encoding as well as result
// values for the 'solve' function of a solver.

enum State {
  READY = 0,            // Initial solver state.
  ADDING = 1,           // Adding clause.
  UNKNOWN = 0,          // Unknown solver result.
  SATISFIABLE = 10,     // All clauses are satisfied.
  UNSATISFIABLE = 20,   // Determined unsatisfiability.
};




// Common invalid literal, invalid variable index, decision level etc.

const unsigned INVALID = UINT_MAX;

// We reserve one bit for bit-stuffing and thus in principle only support
// 2^30 variables (instead of 2^31).  We further have to reserve one
// additional value above the literal range to distinguish tagged literals
// from the decision reason constant 'UINT_MAX = 0xfffffff' and the unit
// reason constant '0xffffffe' one less.  Therefore we require
//
//     (((MAXIMUM_INTERNAL_LITERAL << 1)) | 1)  <  0xffffffe
//
// which gives (note that the tag does not really matter here)
//
//     MAXIMUM_INTERNAL_LITERAL < 0x7ffffff
//
// thus
//
//     MAXIMUM_INTERNAL_LITERAL == 0x7fffffd
//
// as the maximum literal is negated and thus odd.

const unsigned MAXIMUM_INTERNAL_LITERAL = (1u << 31) - 3; // 0x7fffffd

// The following two constants are a consequence of the argument above.

const unsigned MAXIMUM_INTERNAL_VARIABLE = (1u << 30) - 2;
const int MAXIMUM_EXTERNAL_VARIABLE_INDEX = (1 << 30) - 1;

// The individual solver variants only have to implement the following
// abstract virtual member functions. We factored out this interface for
// documentation purposes only.  The parent class of an actual solver
// variant should instead be the 'Shared' class defined below.

struct Interface
{
  virtual ~Interface () { }

  // Short name and description of the solver variant.
  //
  virtual const char * one_word_name () = 0;
  virtual const char * one_line_description () = 0;

  // Notify the internal solver variant to adjust its own internal
  // data-structures to the new increased number of 'variables' and
  // 'literals' triggered by 'import_literal'.
  //
  virtual void adjust_to_increased_variables () = 0;

  // Add the temporary simplified clause as irredundant clause. This
  // function is virtual to keep the code for importing literals resizing
  // data-structures as well as simplifying the clause code shared among all
  // solver variants even though 'resizing' needs the 'resize' call back.
  //
  virtual void add_simplified_irredundant_clause () = 0;


  // Solve and return '10=SATISFIED' and '20=UNSATISFIED' as above.
  //
  virtual State solve_internally () = 0;

  // At the end of solving 'main' wants to print solver specific statistics.
  //
  virtual void print_statistics (double current_process_time) { }

  // Allows to traverse over all options.
  //
  virtual struct Options & iterate_options () = 0;
};

/*------------------------------------------------------------------------*/

// This 'Range' class allows to easily iterate over all variables and all
// literals with the following range-based for-loop idioms:
//
//    for (auto idx : variables)
//      do_something_with (idx);
//
//    for (auto lit : literals) {
//      do_something_with (lit);
//
// We only use it for the 'variables' and 'literals' fields.

class Range
{
  class iterator
  {
    unsigned element;

  public:

    iterator (unsigned i) : element (i) { }
    void operator++ () { assert (element != UINT_MAX); element++; }
    const unsigned & operator * () const { return element; }
    friend bool operator != (const iterator & a, const iterator & b) {
      return a.element != b.element;
    }
  };

  unsigned bound;               // The range is from '0' to 'bound-1'.

public:

  Range (unsigned b = 0) : bound (b) { }

  operator unsigned () const { return bound; }
  void operator++ () { assert (bound != UINT_MAX); bound++; }

  void operator+= (unsigned increment) {
    assert (bound <= UINT_MAX - increment);
    bound += increment;
  }

  iterator end () const { return iterator (bound); }
  iterator begin () const { return iterator (0); }
};

/*------------------------------------------------------------------------*/

// Code for internal performance profiling (used as 'Solver::profiles').

class Profiles
{
  struct Profile
  {
    const char * name;
    double time;

    Profile (const char * n) : name (n), time (0) { }
  };

  struct Profiling
  {
    size_t idx;
    double start;

    Profiling (size_t i, double s) : idx (i), start (s) { }
  };

  std::vector<Profile> profiles;
  std::vector<Profiling> started;

  void flush (Profiling &, double);
  size_t find (const char *);

public:
  
  void start (const char *);
  double stop (const char *);

  void print ();
};

/*------------------------------------------------------------------------*/

struct Option
{
  unsigned value, initial, minimum, maximum;
  const char * name, * description;
  Option (const char * n, unsigned i, unsigned l, unsigned h,
          const char * d) :
    value (i), initial (i), minimum (l), maximum (h),
    name (n), description (d) {
    assert (minimum <= initial);
    assert (initial <= maximum);
    assert (maximum <= (unsigned) INT_MAX);
  }
  operator unsigned () { return value; }
  void operator = (unsigned new_value) {
    if (new_value < minimum) new_value = minimum;
    if (new_value > maximum) new_value = maximum;
    value = new_value;
  }

  struct iterator
  {
    virtual Option * begin () = 0;
    virtual Option * end () = 0;
  };
};

// Used to provide an iterator over all solver options.

struct Options {
  virtual Option * begin () = 0;
  virtual Option * end () = 0;
};

/*------------------------------------------------------------------------*/

#ifdef LOGGING

// Fixed non-resizable buffers for logging literals and variables.

struct Buffer
{
  char chars[256];
  size_t size = 0;
};

struct Buffers
{
  Buffer buffers[4];
  size_t current = 0;

  Buffer & next () {
    Buffer & res = buffers[current++];
    if (current == sizeof buffers / sizeof *buffers) current = 0;
    return res;
  }
};

#endif

/*------------------------------------------------------------------------*/

// This is the actual parent class for all solver variants.

class Shared : Interface
{
public:

  static unsigned sign (unsigned lit) { return lit & 1; }
  static unsigned index (unsigned lit) { return lit / 2; }
  static unsigned negate (unsigned lit) { return lit ^ 1; }
  static unsigned literal (unsigned idx) { return 2*idx; }

protected:                      // Code used by each solver variant.

  State state = READY;		// Actually of type 'State'.
  bool terminated = false;

  static const char * decode_state (int);
  void transition (State new_state);

  Range variables;              // Number of variables.
  Range literals;               // Number of literals (= 2*variables).

  unsigned level = 0;               // Current decision level.

  std::vector<unsigned> levels;     // Maps variables to assignment level.
  std::vector<signed char> marks;   // Mark variables for parsing.
  std::vector<signed char> values;  // Maps literals to values (-1,0,1).

  signed char fixed (int lit) {
    signed char res = values[lit];
    if (!res) return 0;
    if (levels [index (lit)]) return 0;
    return res;
  }

#ifndef NDEBUG
  std::vector<int> original;    // Original zero separated clause literals.
#endif

  void resize (unsigned new_variables);

  // Import added literal in 'add (int)' in the range '[1,-1,2,-2,...]' to
  // internal literals in the range '[0,1,2,3,...]' and if necessary resize
  // data-structures of 'Shared' and through 'resize' of the solver variant.

  unsigned import_literal (int external_literal);

  // Maps internal literals '[0,1,2,3,...]' to external '[1,-1,2,-2,...]'.
  //
  int export_literal (unsigned internal_literal) const;

  std::vector<unsigned> clause;         // Simplified clause.
  std::vector<unsigned> imported;       // Unsimplified clause.

  bool simplify_clause ();              // Returns 'false' if tautological.

  size_t added_original_clauses = 0;    // Number of added original clauses.

  int verbosity = 0;		// Negative if 'quiet'.
  bool model_count_flag = false;
  char debug_info = 0;
  bool assignment_count_flag = false;
  bool projected = false;
  std::vector<unsigned> DLCS;
  std::vector<bool> initial_watched;
  std::vector<bool> positive_polarity;
  std::vector<bool> negative_polarity;
  char counter = 0;

  // Proof tracing part.

  bool ascii = false;		// Use ASCII version of DRAT proof format.
  FILE * proof = 0;		// Trace to this file if non-zero.

  void trace (size_t size, const unsigned *);

  void trace (char type, size_t size, const unsigned * literals) {
    assert (type == 'a' || type == 'd');
    if (!proof) return;
    if (!ascii) fputc (type, proof);
    else if (type == 'd') fputs ("d ", proof);
    trace (size, literals);
  }

  void trace (char type, const std::vector<unsigned> & v) {
    trace (type, v.size (), &v[0]);
  }

  void trace (char type, unsigned lit, unsigned other) {
    unsigned literals[2] = { lit, other };
    trace (type, 2, literals);
  }

  void trace (char type) { trace (type, clause); }

  // Resource usage part.

  static double process_time ();
  static size_t maximum_resident_set_size ();
  static size_t current_resident_set_size ();

  // Internal profiling.

  Profiles profiles;
  void start (const char * name) { profiles.start (name); }
  void stop (const char * name) { profiles.stop (name); }
  friend class Profiles;

  bool parse_long_option (const char * arg);

#ifdef LOGGING

  Buffers buffers;

  bool logging () const { return verbosity == INT_MAX; }

  void debug_prefix ();         // Start a debugging line.
  void debug_suffix ();         // Finishes debugging line with new-line etc.

  // Generic 'printf' style debugging/logging function.  The code is only
  // included if 'LOGGING' is defined (with './configure -l' - the default
  // for './configure -g') and still needs to be enabled at run-time.

  // The '__attribute__ ...' checks correct format strings at compile-time.

  void debug (const char * fmt, ...) __attribute__ ((format (printf, 2, 3)));

  // The following two functions give variable and literal strings with
  // their assigned value and decision level (if assigned) and uses
  // 'values', 'levels' and 'buffer'.

  const char * debug_literal (unsigned lit);
  const char * debug_variable (unsigned idx);

  void debug_vector (std::vector<unsigned> &, const char * fmt, ...)
                     __attribute__ ((format (printf, 3, 4)));
#else

  // The debugging code is not compiled in if 'LOGGING' is undefined.  This
  // allows to make sure that it is completely ignored if not enabled.

#define debug(...) do { } while (0)
#define debug_vector(...) do { } while (0)

#endif

  // Generic messages which print some header and finish the message with a
  // new line.  The '__attribute__ ...' parts make it easy to find incorrect
  // format strings at compile time, not matching the given arguments.  The
  // 'level' for message specifies the requested 'verbosity' level.

  void message (int level, const char * fmt, ...)
    __attribute__ ((format (printf, 3, 4)));

  // Same as message exists the process.
  //
  void error (const char * fmt, ...) __attribute__ ((format (printf, 2, 3)));

  char parse_error[128];
  friend struct Parser;

public:                                 // Common API for all solvers.

  Shared ();
  ~Shared ();

  void set_verbosity (int verbosity);	// Change verbosity (-1=quiet).
  void force_ascii_proof_format ();	// Enforce ASCII DRAT format.
  void set_proof (FILE * file);		// Start proof tracing.

  static const char * version ();       // Global version number.
  static const char * gitid ();         // Return Git-ID.
  static const char * build ();		// Return compilation flags.

  void add (int elit, bool backtrue = false);                       // Import and add literal.
                                        // Clauses are terminated by '0'.

  int solve ();				// External 'solve' function,
  					// returns 10=SATISFIABLE or
					// 20=UNSATISFIABLE through
					// calling 'solve_internally'.

  // There is a way to asynchronously stop solving, which for instance is
  // used to flag time out to the solver if the time limit is reached.
  //
  virtual void terminate_solving () { terminated = true; }


  int value (int lit);                  // Returns 'lit' if assigned true,
                                        // and '-lit' if assigned false,
                                        // and '0' if unassigned.
#ifndef NDEBUG
  void check_model ();                  // Check model to satisfy clauses.
#endif
  void print_model (int max_var,
                    FILE *);            // Print in competition format.

  const char * parse (FILE *,           // Parse DIMACS file and add its
                      int & max_var);   // clauses.  Returns parse error.

  void print_all_option_ranges ();
  void print_all_option_usage ();
  void print_all_option_values ();

  void print_resources ();
  void print_full_statistics ();

  int main (int argc, char ** argv);    // Main of stand-alone solver.

private:

  void catch_signals ();
  void reset_signals ();

  void set_alarm (unsigned seconds);
  void reset_alarm ();

  static volatile Shared * singleton;	// To catch signals.
  static void catch_signal (int sig);	// Actual signal catcher.
  static void catch_alarm (int sig);	// Actual alarm signal catcher.
};

/*------------------------------------------------------------------------*/

#endif
