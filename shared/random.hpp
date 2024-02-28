#ifndef _random_hpp_INCLUDED
#define _random_hpp_INCLUDED

class Random
{
  uint64_t prev = 0;
public:
  uint64_t next_uint64 () {
    prev *= 6364136223846793005ul;
    prev += 1442695040888963407ul;
    return prev;
  }
  unsigned next_unsigned () { return next_uint64 () >> 32; }
  double next_double () { return next_unsigned () / 4294967296.0; }
};

#endif

