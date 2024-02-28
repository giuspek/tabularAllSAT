#ifndef _arena_hpp_INCLUDED
#define _arena_hpp_INCLUDED

#include <cassert>
#include <cstdlib>

// Arena for allocating clauses consecutively and referencing respectively
// accessing them through 32-bit (4 byte) unsigned indices.  We keep the
// arena somewhat generic as we want to factor it out into a shared header
// file and at the same make it independent of the actual clause structure.
// The client clause class 'C' need to support 'next ()' which should point
// to the end of the the clause.  This also needs to be consistent with 'n'
// the number of words given as argument to 'allocate (size_t n)' which in
// turn should be identical to the two 'Clause::words (...)' functions.

template<class C> class Arena
{
  unsigned * _begin = 0, * _end = 0, * _allocated = 0;

public:

  class iterator
  {
    C * clause;
  public:
    iterator (C * c) : clause (c) { }
    void operator++ () { clause = clause->next (); }
    C * operator * () const { return clause; }
    friend bool operator != (const iterator & a, const iterator & b) {
      return a.clause != b.clause;
    }
  };

  ~Arena () { free (_begin); }

  size_t size () { return _end - _begin; }
  size_t capacity () { return _allocated - _begin; }

  C * allocate (size_t n) {
    size_t old_size = size ();
    size_t new_size = old_size + n;
    size_t old_capacity = capacity ();
    if (new_size > old_capacity) {
      size_t new_capacity = old_capacity ? 2*old_capacity : 1;
      while (new_size > new_capacity)
	new_capacity *= 2;
      size_t new_bytes = new_capacity * sizeof *_begin;
      _begin = (unsigned*) realloc (_begin, new_bytes);
      _end = _begin + old_size;
      _allocated = _begin + new_capacity;
    }
    C * res = (C *) _end;
    _end += n;
    return res;
  }

  unsigned position (C * c) {
    unsigned * w = (unsigned *) c;
    assert (_begin <= w), assert (w <= _end);
    return w - _begin;
  }

  void resize (size_t new_size) {
    assert (new_size <= size ());
    _end = _begin + new_size;
  }

  // This is not the clause number 'pos' in the arena but the clause which
  // starts at the word position 'pos' (and extends 'word_clause' words).

  C * operator [] (size_t pos) {
    unsigned * w = _begin + pos;
    C * res = (C*) w;
    assert (res->next () <= (C*) _end);
    return res;
  }

  iterator begin () { return iterator ((C*) _begin); }
  iterator end () { return iterator ((C*) _end); }
};

#endif
