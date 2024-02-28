#ifndef _compact_hpp_INCLUDED
#define _compact_hpp_INCLUDED

#include <cassert>
#include <climits>
#include <cstdlib>

// Like Vector but more compact and with less maximal capacity. The main
// usage is to reduce the size of watcher stacks from 24 to 16 bytes which
// at the same time makes it access faster and in general reducing cache
// access pressure during propagation.

template<class T> class Compact
{
  T * _begin = 0;
  unsigned _size = 0, _capacity = 0;

  bool full () const { return _size == _capacity; }

  void enlarge () {
    _capacity = _capacity ? 2*_capacity : 1;
    _begin = (T*) realloc (_begin, _capacity * sizeof *_begin);
  }

public:

  Compact () { }
  ~Compact () { free (_begin); }

  Compact (const Compact & other) {
    assert (other.size () <= UINT_MAX);
    _size = other.size ();
    _capacity = other.capacity ();
    _begin = (T *) malloc (_capacity * sizeof *_begin);
    T * p = _begin;
    for (const auto & ref : other)
      *p++ = ref;
  }

  bool empty () const { return !_size; }

  size_t size () const { return _size; }

  size_t capacity () const { return _capacity; }

  void push_back (T ref) {
    if (full ()) enlarge ();
    _begin[_size++] = ref;
  }

  T & operator [] (unsigned i) { return _begin[i]; }

  void resize (size_t new_size) {
    assert (new_size <= _size);
    _size = new_size;
  }

  void shrink_to_fit () {
    if (!_capacity) return;
    while (_capacity/2 > _size) _capacity /= 2;
    _begin = (T*) realloc (_begin, _capacity * sizeof *_begin);
  }

  const T * begin () const { return _begin; }
  const T * end () const { return _begin + _size; }

  T * begin () { return _begin; }
  T * end () { return _begin + _size; }
};

#endif
