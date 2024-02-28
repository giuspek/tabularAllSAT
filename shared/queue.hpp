#ifndef _queue_hpp_INCLUDED
#define _queue_hpp_INCLUDED

#include <climits>

// Variable-Move-To-Front (VMTF) decision variable queue with time stamping.

class Queue
{
  static const unsigned INVALID = UINT_MAX;

  struct Link
  {
    unsigned stamp;
    unsigned prev, next;
  };

  unsigned stamp = 0, search = INVALID;
  unsigned first = INVALID, last = INVALID;
  std::vector<Link> links;

  void restamp () {
    stamp = 0;
    for (auto & link : links)
      link.stamp = ++stamp;
  }

public:

  unsigned stamped (unsigned idx) const { return links[idx].stamp; }

  unsigned operator * () const { return search; }
  void pop () { search = links[search].prev; }

  void enqueue (unsigned idx, bool update = false) {
    assert (idx != INVALID);
    if (last == INVALID) first = idx;
    else links[last].next = idx;
    Link & link = links[idx];
    link.prev = last;
    link.next = INVALID;
    link.stamp = ++stamp;
    if (!stamp) restamp ();
    last = idx;
    if (update || search == INVALID) search = idx;
  }

  void dequeue (unsigned idx) {
    Link & link = links[idx];
    if (link.next == INVALID) last = link.prev;
    else links[link.next].prev = link.prev;
    if (link.prev == INVALID) first = link.next;
    else links[link.prev].next = link.next;
  }

  void resize (unsigned variables) {
    unsigned idx = links.size ();
    links.resize (variables);
    while (idx != variables)
      enqueue (idx++);
    if (variables) search = variables - 1;
  }

  void update (unsigned idx) {
    if (search == INVALID || links[search].stamp < links[idx].stamp)
      search = idx;
  }
};

#endif
