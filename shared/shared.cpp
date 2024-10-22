const char * basic_usage =
"usage: %s [ <option> ... ] [ <dimacs> [ <proof> ] ]\n"
"\n"
"with the list of basic options\n"
"\n"
"  -a | --ascii       use ASCII variant of DRAT proof format\n"
"  -h | --help        print this command line option summary\n"
#ifdef LOGGING
"  -l | --logging     print very verbose logging information\n"
#endif
"  -n | --no-witness  do not print witness if satisfiable\n"
"  -c                 count number of solutions\n"
"  -e                 print number of assignments retrieved\n"
"  -q | --quiet       disable all messages\n"
"  -v | --verbose     increase verbosity\n"
"\n"
"  -t <seconds>       set time limit\n"
"\n"
"  --all-options      print full list of all options\n"
"  --option-ranges    print option value ranges\n"
"  --version          print version number\n"
"\n"
"and '<dimacs>' the input file in DIMACS format, optionally compressed with\n"
"'gzip', 'bzip2' or 'xz'.  If no input file is specified the solver reads\n"
"from '<stdin>'.  If the '<proof>' is given a DRAT (actually DRUP) proof is\n"
"written to that file unless the path name '-' is used in which case it is\n"
"written to '<stdout>'.  In this case the output format is the ASCII version\n"
"of DRAT, which can also be enforced by specifying '--ascii'.\n";

#include "shared.hpp"

#include <cstring>

const char * all_options_suffix =
"Without '=' and no argument the value of the internal option is set to '1'\n"
"and to '0' if the name is preceded with 'no'.  Further 'true' and 'false'\n"
"are treated as synonyms for '1' and '0'.\n";

void Shared::print_all_option_ranges () {
  for (auto o : iterate_options ())
    printf ("%s %u %u %u\n", o.name, o.minimum, o.initial, o.maximum);
}

void Shared::print_all_option_usage () {
  printf (basic_usage, one_word_name ());
  printf ("\nThe full list of internal options is as follows:\n\n");
  size_t len = 0;
  for (int i = 0; i != 2; i++) {
    for (auto option : iterate_options ()) {

      char type[32];
      if (option.minimum == option.maximum)
	sprintf (type, "%u", option.minimum);
      else if (option.minimum == 0 && option.maximum == 1)
	strcpy (type, "bool");
      else if (option.maximum == INT_MAX)
	sprintf (type, "%u..INT_MAX", option.minimum);
      else
	sprintf (type, "%u..%u", option.minimum, option.maximum);

      char name[80];
      sprintf (name, "%s=%s", option.name, type);

      if (i) {
	fputs ("  --", stdout);
	fputs (name, stdout);
	size_t tmp = strlen (name);
	assert (tmp <= len);
	for (size_t i = tmp; i != len; i++) fputc (' ', stdout);
	fputs ("  ", stdout);
	fputs (option.description, stdout);
	fputs (" (default ", stdout);
	char value[16];
	if (option.initial == INT_MAX) strcpy (value, "INT_MAX");
	else if (option.maximum == 1 && option.initial) strcpy (value, "true");
	else if (option.maximum == 1 && !option.initial) strcpy (value, "false");
	else sprintf (value, "%u", option.value);
	fputs (value, stdout);
	fputs (")\n", stdout);
      } else {
	size_t tmp = strlen (name);
	if (tmp > len) len = tmp;
      }
    }
  }
  fputc ('\n', stdout);
  fputs (all_options_suffix, stdout);
}

void Shared::print_all_option_values () {

  size_t printed_option_lines = 0;

  for (Option & option : iterate_options ()) {
    if (option.value != option.initial) {
      if (!printed_option_lines++) fputs ("c\n", stdout);
      printf ("c using '--%s=%u' (different from default %u)\n",
	       option.name, option.value, option.initial);
    } else if (verbosity) {
      printed_option_lines++;
      printf ("c using '--%s=%u' (same as default %u)\n",
	       option.name, option.value, option.initial);
    }
  }

#if 0
  if (printed_option_lines)
    fputs ("c\n", stdout), fflush (stdout);
#endif
}


#include <cctype>
#include <cstdarg>
#include <cstdlib>

Shared::Shared () { }
Shared::~Shared () { }

void Shared::set_verbosity (int new_verbosity) { verbosity = new_verbosity; }
void Shared::force_ascii_proof_format () { ascii = true; }

#ifdef LOGGING // Start of section with debugging / logging code.

void Shared::debug_prefix () {
  assert (logging ());
  printf ("c DEBUG %u ", level);
}

void Shared::debug_suffix () {
  assert (logging ());
  fputc ('\n', stdout);
  fflush (stdout);
}

void Shared::debug (const char * fmt, ...) {
  if (!logging ())
    return;
  debug_prefix ();
  va_list ap;
  va_start (ap, fmt);
  vprintf (fmt, ap);
  va_end (ap);
  debug_suffix ();
}

const char * Shared::debug_literal (unsigned lit) {
  assert (lit < literals);
  Buffer & buffer = buffers.next ();
  char * res = buffer.chars;
  int exported = export_literal (lit);
  sprintf (res, "%u(%d)", lit, exported);
  int value = values[lit];
  if (value) {
    char * tmp = res + strlen (res);
    unsigned idx = index (lit);
    unsigned level = levels[idx];
    if (level == INVALID) sprintf (tmp, "=%d", value);
    else sprintf (tmp, "@%u=%d", level, value);
  }
  assert (strlen (res) < sizeof buffer.chars);
  return res;
}

const char * Shared::debug_variable (unsigned idx) {
  unsigned lit = literal (idx);
  const char * tmp = debug_literal (lit);
  Buffer & buffer = buffers.next ();
  char * res = buffer.chars;
  sprintf (res, "variable %u literal %s", idx, tmp);
  assert (strlen (res) < sizeof buffer.chars);
  return res;
  return res;
}

void Shared::debug_vector (std::vector<unsigned> & v,
                           const char * fmt, ...) {
  if (!logging ())
    return;
  debug_prefix ();
  va_list ap;
  va_start (ap, fmt);
  vprintf (fmt, ap);
  va_end (ap);
  printf (" size %zu clause", v.size ());
  for (auto lit : v)
    printf (" %s", debug_literal (lit));
  debug_suffix ();
}

#endif // End of section with debugging / logging code.

const char * Shared::decode_state (int state) {
  switch (state) {
    case READY: return "READY";
    case ADDING: return "ADDING";
    case SATISFIABLE: return "SATISFIABLE";
    case UNSATISFIABLE: return "UNSATISFIABLE";
  }
  return "INVALID";
}

void Shared::transition (State new_state) {
  debug ("transition to state '%s' from state '%s'",
         decode_state (new_state), decode_state (state));
  state = new_state;
}

void Shared::resize (unsigned new_variables) {
  if (variables >= new_variables) return;
  variables = new_variables;
  literals = 2*variables;
  levels.resize (variables);
  marks.resize (variables);
  values.resize (literals);
  adjust_to_increased_variables ();
  debug ("resized to %u variables", (unsigned) variables);
}

unsigned Shared::import_literal (int elit) {
  assert (elit);
  assert (elit != INT_MIN);
  int eidx = abs (elit);
  unsigned iidx = eidx - 1;
  if (iidx >= variables) resize (eidx);
  unsigned ilit = literal (iidx);
  if (elit < 0)
    ilit = negate (ilit);
  debug ("imported external literal %d as internal literal %u", elit, ilit);
  return ilit;
}

int Shared::export_literal (unsigned internal_literal) const {
  assert (internal_literal < literals);
  unsigned internal_index = index (internal_literal);
  assert (internal_index < (unsigned) INT_MAX);
  int external_literal = internal_index + 1;
  if (sign (internal_literal))
    external_literal = -external_literal;
  return external_literal;
}

void Shared::message (int level, const char * fmt, ...) {
  if (verbosity < level)
    return;
  fputs ("c ", stdout);
  va_list ap;
  va_start (ap, fmt);
  vprintf (fmt, ap);
  va_end (ap);
  fputc ('\n', stdout);
  fflush (stdout);
}

void Shared::error (const char * fmt, ...) {
  fputs (one_word_name (), stderr);
  fputs (": error: ", stderr);
  va_list ap;
  va_start (ap, fmt);
  vfprintf (stderr, fmt, ap);
  va_end (ap);
  fputc ('\n', stderr);
  fflush (stderr);
  exit (1);
}

void Shared::set_proof (FILE * file) { 
  if (!file) error ("zero proof file argument in 'set_proof");
  if (proof) error ("can not call 'set_proof' twice");
  if (state != READY)
    error ("'set_proof' called in '%s' state", decode_state (state));
  if (added_original_clauses)
    error ("'set_proof' called after adding clauses");
  proof = file;
}

void Shared::trace (size_t size, const unsigned * literals) {
  assert (proof);
  if (ascii) {
    for (const unsigned * p = literals, * end = p + size; p != end; p++)
      fprintf (proof, "%d ", export_literal (*p));
    fputs ("0\n", proof);
  } else {
    for (const unsigned * p = literals, * end = p + size; p != end; p++) {
      unsigned tmp = *p + 2;
      while (tmp & ~0x7f)
	fputc ((tmp & 0x7f) | 0x80, proof), tmp >>= 7;
      fputc (tmp, proof);
    }
    fputc (0, proof);
  }
}

void Shared::add (int elit, bool backtrue) {
  if (elit == INT_MIN || abs (elit) > MAXIMUM_EXTERNAL_VARIABLE_INDEX)
    error ("invalid literal in 'add (%d)'", elit);
  if (state != READY && state != ADDING)
    error ("can not call 'add (%d)' in '%s' state",
           elit, decode_state (state));

  if (state != ADDING) transition (ADDING);
#ifndef NDEBUG
  original.push_back (elit);
#endif
  if (elit) {
    unsigned ilit = import_literal (elit);
    imported.push_back (ilit);
    if(elit > 0)
      positive_polarity[abs(elit) - 1] = true;
    else
      negative_polarity[abs(elit) - 1] = true;
    if (counter++ < 2)
      initial_watched[abs(elit) - 1] = true;
  } else {
    for(auto l: imported)
      debug ("imported literal %s", debug_literal (l));

    counter = 0;
    added_original_clauses++;
    assert (clause.empty ());
    bool tautological = false;
    for (auto ilit : imported) {
      signed char value = values[ilit];
      if (value > 0)
        {
          debug_vector (imported, "original[%zu] %s satisfied unsimplified",
                        added_original_clauses, debug_literal (ilit));
          tautological = true;
          break;
        }
      if (value < 0 && !backtrue)
        continue;
      unsigned iidx = index (ilit);
      DLCS[iidx] += 1;
      signed char mark = marks[iidx];
      if (sign (ilit))
        mark = -mark;
      if (mark < 0)
        {
          debug_vector (imported,
	                "original[%zu] tautological unsimplified",
			added_original_clauses);
          tautological = true;
          break;
        }
      if (mark > 0)
        continue;
      marks[iidx] = sign (ilit) ? -1 : 1;
      clause.push_back (ilit);
    }
    for (auto ilit : clause)
      marks[index (ilit)] = 0;
    if (tautological) trace ('d', imported);
    else {
      debug_vector (imported, "original[%zu] imported",
                    added_original_clauses);
      add_simplified_irredundant_clause ();
      if (proof && clause.size () != imported.size ())
        trace ('a'), trace ('d', imported);
    }
    imported.clear ();
    clause.clear ();
    transition (READY);
  }
}

int Shared::solve () {
  if (state == ADDING)
    error ("can not call 'solve ()' without terminating clause with '0'");
  if (state != READY && state != SATISFIABLE)
    error ("invalid '%s' state calling 'solve ()'", decode_state (state));
  State res = solve_internally ();
  if (res != state) transition (res);
  debug ("'solve ()' call returns %d ('%s')", state, decode_state (state));
  terminated = false;
  return res;
}

int Shared::value (int elit) {
  if (!elit || elit == INT_MIN ||
      abs (elit) > MAXIMUM_EXTERNAL_VARIABLE_INDEX)
    error ("invalid literal in 'value (%d)'", elit);
  int eidx = abs (elit);
  unsigned iidx = eidx - 1;
  if (iidx >= variables) return 0;
  unsigned ilit = literal (iidx);
  if (elit < 0) ilit = negate (ilit);
  assert (ilit < literals);
  signed char value = values[ilit];
  return !value ? 0 : value < 0 ? -elit : elit;
}

struct Parser
{
  Shared & shared;
  FILE * file;
  size_t lineno = 0, litlineno = 0;
  
  int & variables, clauses = 0, important_size = 0, last = '\n';
  Parser (Shared & s, FILE * f, int & m) :
    shared (s), file (f), variables (m) { }
  int get_char () {
    if (last == '\n') lineno++;
    int res = getc (file);
    if (res == '\r') {
      res = getc (file);
      if (res != '\n') res = '\r';
    }
    return last = res;
  }
  bool get_int (int & res, int & ch) {
    litlineno = lineno;
    int sign = (ch == '-' ? -1 : 1);
    if (ch == '-' && (ch = get_char ()) == '0') return false;
    if (!isdigit (ch)) return false;
    res = ch - '0';
    while (isdigit (ch = get_char ())) {
      if (INT_MAX/10 < res) return false;
      res *= 10;
      int digit = ch - '0';
      if (INT_MAX - digit < res) return false;
      res += digit;
    }
    res *= sign;
    return true;
  }

  const char * error (const char * fmt) {
    sprintf (shared.parse_error, "line %zu: %s", lineno, fmt);
    assert (strlen (shared.parse_error) < sizeof shared.parse_error);
    return shared.parse_error;
  }
  const char * invalid_literal () {
    assert (litlineno), lineno = litlineno;
    return error ("invalid literal");
  }
  const char * end_of_file_in_comment () {
    return error ("end-of-file in comment");
  }
  const char * no_new_line_after_carriage_return () {
    return error ("carriage return without new line");
  }
  const char * parse ();
};

static bool string_to_unsigned (const char * str, unsigned & res) {
  const char * p = str;
  char ch = *p++;
  if (!isdigit (ch)) return false;
  unsigned tmp = ch - '0';
  ch = *p++;
  if (ch) {
    if (!isdigit (ch)) return false;
    if (!tmp) return false;
    tmp = 10*tmp + (ch - '0');
    while ((ch = *p++)) {
      if (!isdigit (ch)) return false;
      if (UINT_MAX/10 < tmp) return false;
      tmp *= 10;
      unsigned digit = ch - '0';
      if (UINT_MAX - digit < tmp) return false;
      tmp += digit;
    }
  }
  res = tmp;
  return true;
}

static bool match_long_option (const char * arg,
                               const char * name, unsigned & res) {
  if (arg[0] != '-' || arg[1] != '-') return false;
  if (arg[2] == 'n' && arg[3] == 'o' && arg[4] == '-' &&
      !strcmp (arg + 5, name)) { res = 0; return true; }
  const char * a = arg + 2;
  const char * n = name;
  while (*a && *a == *n)
    a++, n++;
  if (*n) return false;
  if (!*a) { res = 1; return true; }
  if (*a++ != '=') return false;
  if (!strcmp (a, "true")) { res = 1; return true; }
  if (!strcmp (a, "false")) { res = 0; return true; }
  return string_to_unsigned (a, res);
}

bool Shared::parse_long_option (const char * arg) {
  for (Option & option : iterate_options ()) {
    unsigned res;
    if (!match_long_option (arg, option.name, res)) continue;
    if (res < option.minimum) return false;
    if (res > option.maximum) return false;
    option = res;
    return true;
  }
  return false;
}

const char * Parser::parse () {
  int ch, comments = 0;
#ifndef NDEBUG
  std::vector<char> buffer;
#endif
  while ((ch = get_char ()) == 'c') {
    comments++;
#ifdef NDEBUG
    while ((ch = get_char ()) != '\n')
      if (ch == '\r') return no_new_line_after_carriage_return ();
      else if (ch == EOF) return end_of_file_in_comment ();
#else
    assert (buffer.empty ());
    while ((ch = get_char ()) != '\n') {
      if (ch == '\r') return no_new_line_after_carriage_return ();
      else if (ch == EOF) return end_of_file_in_comment ();
      if (ch != ' ' && ch != '\t') buffer.push_back (ch);
    }
    buffer.push_back (0);
    const char * arg = &buffer[0];
    shared.parse_long_option (arg);
    buffer.clear ();
#endif
  }
  if (ch == EOF) {
    if (comments) {
      assert (lineno > 1), lineno--;
      return error ("expected header after comments");
    } else return error ("empty input file");
  }
  if (ch != 'p') return error ("expected header or comment");
  for (const char * p = " cnf "; *p; p++)
    if (get_char () != *p) return error ("invalid header");
  if (!get_int (variables, ch = get_char ()) || variables < 0)
    return error ("expected number of variables");
  if (ch != ' ') return error ("expected space after variables");
  if (variables > MAXIMUM_EXTERNAL_VARIABLE_INDEX)
    return error ("maximum number of variables exceeded");
  if (!get_int (clauses, ch = get_char ()) || clauses < 0)
    return error ("expected number of clauses");
  if (ch == ' '){
    if (!get_int (important_size, ch = get_char ()) || important_size < 0 || important_size > variables)
      //return error ("expected number of important atoms (projected enumeration)");
      {}
    else
      shared.projected = true; 
  } 
  if (ch != '\n')
    do {
      if (ch == EOF)
	return error ("end-of-file in header before end-of-line");
      else if (ch == '\r') return no_new_line_after_carriage_return ();
      else if (ch != ' ' && ch != '\t')
	return error ("invalid character after header before end-of-line");
      ch = get_char ();
    } while (ch != '\n');
  if (!shared.projected)
    shared.message (0, "parsed header 'p cnf %d %d'", variables, clauses);
  else
    shared.message (0, "parsed header 'p cnf %d %d %d'", variables, clauses, important_size);
  shared.resize (variables);
  shared.DLCS.resize(variables, 0);
  shared.positive_polarity.resize(variables);
  shared.negative_polarity.resize(variables);
  shared.initial_watched.resize(variables);
  shared.important.resize(variables);

  if (shared.projected == true)
    shared.important.resize(variables, false);
  else
    shared.important.resize(variables, true);
  
  // projected extension
  while (shared.projected && (ch = get_char()) == 'c') {
    for (const char * p = " p show "; *p; p++){
      if (get_char () != *p){
        // This is a regular comment line, skip it
        while (ch != '\n' && ch != EOF) ch = get_char();
        if (ch != '\n') return error("invalid character in comment line");
      }
    }
    while (ch != '\n' && ch != EOF) {
      int var;
      if (!get_int(var, ch = get_char()) || var <= 0){
        return error("expected positive variable number in 'c p show' line");
      }
      if (var > variables)
        return error("too many important variables");

      shared.important[var - 1] = true;
      shared.projected_count++;
    }
    if (ch != '\n') return error("invalid character in 'c p show' line");
  } 

  int elit = 0, parsed = 0;
  while (ch != EOF) {
    if (shared.terminated) return 0;
    if (ch == 'c') {
SKIP_COMMENT:
      while ((ch = get_char ()) != '\n')
	if (ch == '\r') return no_new_line_after_carriage_return ();
	else if (ch == EOF) return end_of_file_in_comment ();
      ch = get_char ();
      continue;
    }
    if (ch == '\r') return no_new_line_after_carriage_return ();
    if (ch == ' ' || ch == '\t' || ch == '\n'){ ch = get_char (); continue;}
    if (ch != '-' && !isdigit (ch))
      return error ("expected literal or comment");
    if (!get_int (elit, ch)) return invalid_literal ();
    if (ch != 'c' && ch != ' ' && ch != '\t' && ch != '\n')
      return error ("expected white space after literal");
    if (!clauses) return error ("no clause expected");
    if (parsed == 1 && clauses == 1)
      return error ("more than one clause");
    if (parsed == clauses) return error ("too many clauses");
    assert (elit != INT_MIN);
    if (abs (elit) > variables) return invalid_literal ();
    parsed += !elit;
    shared.add (elit);
    if (ch == 'c') goto SKIP_COMMENT;

    ch = get_char ();
  }
  if (elit) {
    assert (litlineno), lineno = litlineno;
    return error ("zero missing after last literal");
  }
  if (parsed != clauses) return error ("clause missing");
  if (parsed == 1) shared.message (0, "parsed single clause");
  else shared.message (0, "parsed all %d clauses", parsed);

  return 0;
}

const char *
Shared::parse (FILE * file, int & maximum_variable) {
  Parser parser (*this, file, maximum_variable);
  return parser.parse ();
}

extern "C" {
#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
};

double Shared::process_time () {
  struct rusage u;
  double res;
  if (getrusage (RUSAGE_SELF, &u)) return 0;
  res = u.ru_utime.tv_sec + 1e-6 * u.ru_utime.tv_usec;
  res += u.ru_stime.tv_sec + 1e-6 * u.ru_stime.tv_usec;
  return res;
}

size_t Shared::maximum_resident_set_size (void) {
  struct rusage u;
  if (getrusage (RUSAGE_SELF, &u)) return 0;
  return ((size_t) u.ru_maxrss) << 10;
}

size_t Shared::current_resident_set_size (void) {
  char path[48];
  sprintf (path, "/proc/%d/statm", (int) getpid ());
  FILE *file = fopen (path, "r");
  if (!file) return 0;
  size_t dummy, rss;
  int scanned = fscanf (file, "%zu %zu", &dummy, &rss);
  fclose (file);
  return scanned == 2 ? rss * sysconf (_SC_PAGESIZE) : 0;
}

void Shared::print_resources () {
  printf ("c %-34s %13.2f seconds\n",
          "process-time:", process_time ());
  double mb = maximum_resident_set_size () / (double) (1<<20);
  printf ("c %-34s %13.2f MB\n",
          "maximum-resident-set-size:", mb);
}

void Profiles::flush (Profiles::Profiling & profiling, double now) {
  assert (profiling.start <= now);
  double delta = now - profiling.start;
  profiling.start = now;
  profiles[profiling.idx].time += delta;
}

size_t Profiles::find (const char * name) {
  size_t size = profiles.size ();
  for (size_t i = 0; i != size; i++)
    if (!strcmp (profiles[i].name, name))
      return i;
  profiles.push_back (Profile (name));
  return size;
}

void Profiles::start (const char * name) {
  size_t idx = find (name);
  double now = Shared::process_time ();
  started.push_back (Profiling (idx, now));
}

double Profiles::stop (const char * name) {
  assert (!started.empty ());
  Profiling & profiling = started.back ();
#ifndef NDEBUG
  assert (profiling.idx < profiles.size ());
  Profile & profile = profiles[profiling.idx];
  assert (!strcmp (profile.name, name));
#else
  (void) name;
#endif
  double now = Shared::process_time ();
  flush (profiling, now);
  started.pop_back ();
  return now;
}

void Profiles::print () {
  double now = Shared::process_time ();
  for (auto & profiling : started) flush (profiling, now);
  size_t size = profiles.size ();
  for (unsigned i = 0; i != size; i++)
    for (unsigned j = i+1; j != size; j++)
      if (profiles[i].time < profiles[j].time ||
          (profiles[i].time == profiles[j].time &&
	   strcmp (profiles[i].name, profiles[j].name) < 0))
	std::swap (profiles[i], profiles[j]);
  for (auto & profile : profiles)
    printf ("c %11.2f  %6.2f%%  %s\n",
            profile.time, now ? 100*profile.time / now : 0, profile.name);
  fputs ("c -----------------------------------\n", stdout);
  printf ("c %11.2f  100.00%%  total\n", now);
}

#ifndef NDEBUG

// Checking the model on the original formula is extremely useful for
// testing and debugging.  This 'checker' aborts if an unsatisfied clause is
// found and prints the clause on '<stderr>' for debugging purposes.

void Shared::check_model () {
  debug ("checking model");
  bool satisfied = false;
  size_t start = 0, count = 1;
  for (size_t i = 0; i != original.size (); i++) {
    int lit = original[i];
    if (lit) {
      if (!satisfied && value (lit) == lit) satisfied = true;
    } else if (!satisfied) {
      fprintf (stderr, "%s: unsatisfied clause[%zu]:\n",
               one_word_name (), count);
      for (size_t j = start; j != i; j++)
        fprintf (stderr, "%d ", original[j]);
      fputs ("0\n", stderr);
      fflush (stderr);
      abort ();
      exit (1);
    } else {
      count++;
      satisfied = false;
      start = i + 1;
    }
  }
}

#endif

class Line
{
  char buffer[80];
  size_t size = 0;

public:

  void put_char (char ch) {
    assert (size < sizeof buffer);
    buffer[size++] = ch;
  }

  void put_str (const char * s) {
    while (*s) put_char (*s++);
  }

  void put_int (int i) {
    char tmp[16];
    sprintf (tmp, " %d", i);
    size_t len = strlen (tmp);
    if (size + len > 76) {
      flush ();
      put_char ('v');
    }
    put_str (tmp);
  }

  void flush () {
    if (!size) return;
    for (size_t i = 0; i != size; i++)
      fputc (buffer[i], stdout);
    fputc ('\n', stdout);
    size = 0;
  }
};

// Printing the model in the format of the SAT competition, e.g.,
//
//   v -1 2 3 0
//
// Always prints a full assignments even if not all values are set.

void Shared::print_model (int maximum_variable, FILE * file) {
  Line line;
  line.put_char ('v');
  int idx = 0;
  assert (0 <= maximum_variable);
  while (idx++ != maximum_variable) {
    int tmp = value (idx);
    if (tmp) line.put_int (tmp == idx ? idx : -idx);
  }
  line.put_str (" 0");
  line.flush ();
}

static bool has_suffix (const char * a, const char * b) {
  size_t l = strlen (a), k = strlen (b);
  return l >= k && !strcmp (a + l - k, b);
}

static bool looks_like_compressed (const char * name) {
  if (has_suffix (name, ".gz")) return true;
  if (has_suffix (name, ".bz2")) return true;
  if (has_suffix (name, ".xz")) return true;
  return false;
}

static bool looks_like_cnf (const char * name) {
  if (has_suffix (name, ".cnf")) return true;
  if (has_suffix (name, ".dimacs")) return true;
  return false;
}

static FILE * read_from_pipe (const char * fmt, const char * path) {
  char * cmd = (char*) malloc (strlen (fmt) + strlen (path));
  sprintf (cmd, fmt, path);
  FILE * res = popen (cmd, "r");
  free (cmd);
  return res;
}

int Shared::main (int argc, char ** argv) {

  unsigned time_limit = 0;
  bool witness_printing_enabled = true;

  const char * proof_name = 0, * dimacs_name = 0;
  int close_dimacs_file = 2, close_proof_file = 0;
  FILE * dimacs_file = 0, * proof_file = 0;

  for (int i = 1; i != argc; i++) {
    const char * arg = argv[i];
    if (!strcmp (arg, "-a") || !strcmp (arg, "--ascii"))
      ascii = true;
    else if (!strcmp (arg, "--all-options")) {
      print_all_option_usage (); return 0;
    } else if (!strcmp (arg, "--option-ranges")) {
      print_all_option_ranges (); return 0;
    } else if (!strcmp (arg, "-h") || !strcmp (arg, "--help")) {
      printf (basic_usage, one_word_name ()); return 0;
    } else if (!strcmp (arg, "-l") || !strcmp (arg, "--logging"))
#ifdef LOGGING
      verbosity = INT_MAX;
#else
      error ("compiled without logging (use './configure --logging')");
#endif
    else if (!strcmp (arg, "-n") || !strcmp (arg, "--no-witness"))
      witness_printing_enabled = false;
    else if (!strcmp (arg, "-q") || !strcmp (arg, "--quiet"))
      verbosity = -1;
    else if (!strcmp (arg, "-v") || !strcmp (arg, "--verbose"))
      verbosity += (verbosity >= 0 && verbosity != INT_MAX);
    else if (!strcmp (arg, "-c") || !strcmp (arg, "--count"))
      model_count_flag = true;
    else if (!strcmp (arg, "--total-length"))
      debug_info = 1;
    else if (!strcmp (arg, "-e"))
      assignment_count_flag = true;
    else if (!strcmp (arg, "--version")) {
      printf ("%s\n", version ()); return 0;
    } else if (!strcmp (arg, "-t")) {
      if (++i == argc) error ("argument to '-t' missing");
      arg = argv[i];
      if (!string_to_unsigned (arg, time_limit) || !time_limit)
	error ("expected positive number in '-t %s'", arg);
    } else if (arg[0] == '-' && arg[1] == '-') {
      if (!parse_long_option (arg))
	error ("invalid long option '%s'", arg);
    } else if (arg[0] == '-' && arg[1])
      error ("invalid option '%s' (try '-h')", arg);
    else if (proof_name)
      error ("too many arguments '%s', '%s' and '%s' (try '-h')",
	     dimacs_name, proof_name, arg);
    else if (dimacs_name) proof_name = arg;
    else dimacs_name = arg;
  }

  if (proof_name) {
    if (!strcmp (proof_name, "-")) {
      proof_name = "<stdout>";
      proof_file = stdout;
      ascii = true;
    } else if (looks_like_compressed (proof_name))
      error ("can not write compressed proof file '%s'", proof_name);
    else if (looks_like_cnf (proof_name))
      error ("will not '%s' with proof which looks like a CNF file",
              proof_name);
    else if (!(proof_file = fopen (proof_name, "w")))
      error ("can not write to proof file '%s'", proof_name);
    else close_proof_file = 1;
  }

  if (!dimacs_name || !strcmp (dimacs_name, "-")) {
    dimacs_name = "<stdin>";
    close_dimacs_file = 0;
    dimacs_file = stdin;
  } else if (has_suffix (dimacs_name, ".gz"))
   dimacs_file = read_from_pipe ("gunzip -c -d %s", dimacs_name);
  else if (has_suffix (dimacs_name, ".bz2"))
   dimacs_file = read_from_pipe ("bunzip2 -c -d %s", dimacs_name);
  else if (has_suffix (dimacs_name, ".xz"))
   dimacs_file = read_from_pipe ("xz -c -d %s", dimacs_name);
  else if (!(dimacs_file = fopen (dimacs_name, "r")))
    error ("could not open and read '%s'", dimacs_name);
  else close_dimacs_file = 1;

  if (verbosity >= 0)
    printf ("c %s\nc Version %s %s\nc Compiled with '%s'\nc\n",
            one_line_description (), version (), gitid (), build ());

  {
    size_t printed_option_lines = 0;

    if (proof_file) {
      set_proof (proof_file),
      message (0, "enabled proof tracing to '%s'", proof_name);
      printed_option_lines++;
    }

    catch_signals ();
    if (time_limit) {
      message (0, "setting time limit of %u seconds", time_limit);
      set_alarm (time_limit);
      printed_option_lines++;
    }

    if (verbosity >= 0 && printed_option_lines) fputs ("c\n", stdout);
  }

  message (0, "reading from '%s'", dimacs_name);
  int maximum_variable;
  const char * failed = parse (dimacs_file, maximum_variable);
  if (close_dimacs_file == 1) fclose (dimacs_file);
  if (close_dimacs_file == 2) pclose (dimacs_file);
  if (failed) {
    fprintf (stderr, "%s: parse error in '%s': %s\n",
             one_word_name (), dimacs_name, failed);
    fflush (stderr);
    exit (1);
  }

  int res = UNKNOWN;
  if (terminated)
    message (0, "time limit hit during parsing");
  else {
    if (verbosity >= 0) print_all_option_values ();
    res = solve ();
  }

  if (close_proof_file) fclose (proof_file);

  if (res == 10) {
#ifndef NDEBUG
    check_model ();
#endif
    printf ("s SATISFIABLE\n");
    fflush (stdout);
    if (witness_printing_enabled)
      print_model (maximum_variable, stdout);
  }
  else if (res == 20) {}
    // printf ("s UNSATISFIABLE\n");
    
  fflush (stdout);

  if (time_limit) reset_alarm ();
  reset_signals ();

  if (verbosity >= 0) {
    fputs ("c\n", stdout);
    print_full_statistics ();
    printf ("c\nc exit code %d\n", res);
    fflush (stdout);
  }

  return res;
}

void Shared::print_full_statistics () {
  double p = process_time ();
  profiles.print ();
  fputs ("c\n", stdout);
  print_statistics (p);
  fputs ("c\n", stdout);
  print_resources ();
}

// Print statistics even if killed after receiving a signal.  This is not
// completely safe as we lazily use standard non-reentrant functions such as
// 'printf'.  Occasionally these might fail (with the result that the
// signal handling code hangs and/or does not produce any output).  We also
// need access to a (single) global solver instances, which is register with
// the signal handlers through 'catch_signals' below.

#include <csignal>

volatile Shared * Shared::singleton;

#define SIGNALS \
SIGNAL(SIGABRT) \
SIGNAL(SIGBUS) \
SIGNAL(SIGILL) \
SIGNAL(SIGINT) \
SIGNAL(SIGSEGV) \
SIGNAL(SIGTERM)

#define SIGNAL(SIG) \
static void (*saved_ ## SIG ## _handler)(int);
SIGNALS

void Shared::reset_signals () {
  assert (singleton == this);
  singleton = 0;
#undef SIGNAL
#define SIGNAL(SIG) \
  signal (SIG, saved_ ## SIG ## _handler);
  SIGNALS
}

void Shared::catch_signal (int sig) {
  Shared * solver = (Shared*) singleton;
  solver->reset_signals ();
  if (solver && solver->verbosity >= 0) {
    printf ("c\nc caught signal %d\nc\n", sig);
    solver->print_full_statistics ();
    printf ("c\nc raising signal %d\n", sig);
    fflush (stdout);
  }
  raise (sig);
}

static void (*saved_SIGALRM_handler)(int);

void Shared::reset_alarm () {
  assert (singleton == this);
  signal (SIGALRM, saved_SIGALRM_handler);
}

void Shared::catch_alarm (int sig) {
  assert (sig == SIGALRM);
  Shared * solver = (Shared*) singleton;
  solver->reset_alarm ();
  solver->terminate_solving ();
  (void) sig;
  return;
}

void Shared::catch_signals () {
  assert (!singleton);
  singleton = this;
#undef SIGNAL
#define SIGNAL(SIG) \
  saved_ ## SIG ## _handler = signal (SIG, catch_signal);
  SIGNALS
}

void Shared::set_alarm (unsigned seconds) {
  assert (singleton);
  assert (this == singleton);
  saved_SIGALRM_handler = signal (SIGALRM, catch_alarm);
  alarm (seconds);
}
