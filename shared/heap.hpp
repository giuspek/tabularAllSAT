#ifndef _heap_h_INCLUDED
#define _heap_h_INCLUDED

#include <cassert>
#include <climits>
#include <set>

class Heap
{
  static const unsigned INVALID = UINT_MAX;

  struct Node
  {
    double score;
    double DLCS;
    unsigned watched;
    unsigned important;
    unsigned cc_value;
    unsigned child = INVALID, prev = INVALID, next = INVALID;
    Node(double s = 0, double d = 0, unsigned b = 0, unsigned i = 0) : score(s), DLCS(d), watched(b), important(i)
    {
    }
  };

  const double limit = 1e150;

  double factor = 1.05;
  double increment = 1.0;
  std::vector<Node> nodes;
  std::vector<unsigned> DLCS_scores;
  unsigned root = INVALID;

  unsigned merge(unsigned a, unsigned b)
  {
    if (a == INVALID)
      return b;
    if (b == INVALID)
      return a;
    unsigned parent, child;
    if (
        (nodes[b].important > nodes[a].important) ||
        ((nodes[b].important == nodes[a].important) && (nodes[b].score + 0.5 * nodes[b].DLCS) > (nodes[a].score + 0.5 * nodes[a].DLCS)) || 
        ((nodes[b].important == nodes[a].important) && (nodes[b].score + 0.5 * nodes[b].DLCS) == (nodes[a].score + 0.5 * nodes[a].DLCS) && nodes[b].watched > nodes[a].watched) ||
        ((nodes[b].important == nodes[a].important) && (nodes[b].score + 0.5 * nodes[b].DLCS) == (nodes[a].score + 0.5 * nodes[a].DLCS) && nodes[b].watched == nodes[a].watched && b > a)
      ){
      parent = b, child = a;
      }
    else
      parent = a, child = b;
    auto &node_parent = nodes[parent];
    unsigned parent_child = node_parent.child;
    auto &node_child = nodes[child];
    node_child.next = parent_child;
    if (parent_child != INVALID)
      nodes[parent_child].prev = child;
    node_child.prev = parent;
    node_parent.child = child;
    node_parent.prev = node_parent.next = INVALID;
    return parent;
  }

  unsigned collapse(unsigned idx)
  {
    if (idx == INVALID)
      return INVALID;
    unsigned next = idx, tail = INVALID;
    do
    {
      unsigned a = next;
      assert(a != INVALID);
      auto &node_a = nodes[a];
      unsigned b = node_a.next;
      if (b == INVALID)
      {
        node_a.prev = tail;
        tail = a;
        break;
      }
      next = nodes[b].next;
      unsigned tmp = merge(a, b);
      assert(tmp != INVALID);
      nodes[tmp].prev = tail;
      tail = tmp;
    } while (next != INVALID);

    unsigned res = INVALID;
    while (tail != INVALID)
    {
      unsigned prev = nodes[tail].prev;
      res = merge(res, tail);
      tail = prev;
    }
    return res;
  }

  void remove(unsigned idx)
  {
    assert(idx != INVALID);
    auto &node_idx = nodes[idx];
    unsigned prev = node_idx.prev;
    unsigned next = node_idx.next;
    assert(prev != INVALID);
    auto &node_prev = nodes[prev];
    node_idx.prev = INVALID;
    if (node_prev.child == idx)
      node_prev.child = next;
    else
      node_prev.next = next;
    if (next != INVALID)
      nodes[next].prev = prev;
  }

  void rescore(double scale)
  {
    assert(scale > 0);
    double inverse = 1 / scale;
    for (auto &node : nodes)
      node.score *= inverse;
    increment *= inverse;
  }

  void update(unsigned idx, double new_score)
  {
    auto &node_idx = nodes[idx];
    double old_score = node_idx.score;
    assert(old_score <= new_score);
    if (old_score == new_score)
      return;
    node_idx.score = new_score;
    unsigned tmp = root;
    if (tmp == idx)
      return;
    if (node_idx.prev == INVALID)
      return;
    remove(idx);
    root = merge(tmp, idx);
  }

public:
  void set_decay(double d)
  {
    assert(0.5 <= d), assert(d <= 1.0);
    factor = 1 / d;
  }

  void update_watched(unsigned idx)
  {
    auto &node_idx = nodes[idx];
    node_idx.watched = true;
    unsigned tmp = root;
    if (tmp == idx)
      return;
    if (node_idx.prev == INVALID)
      return;
    remove(idx);
    root = merge(tmp, idx);
  }


  unsigned operator*() const { return root; }

  bool empty() const { return root == INVALID; }

  double score(unsigned idx) const { return nodes[idx].score + 1.5 * nodes[idx].DLCS; }

  char is_watched(unsigned idx) const { return nodes[idx].watched; }

  char is_important(unsigned idx) const { return nodes[idx].important; }

  bool contains(unsigned idx) const
  {
    return root == idx || nodes[idx].prev != INVALID;
  }

  void push(unsigned idx)
  {
    assert(idx != INVALID);
    assert(!contains(idx));
    nodes[idx].child = INVALID;
    root = merge(root, idx);
    assert(contains(idx));
  }

  void resize(unsigned variables, std::vector<unsigned> DLCS_score, std::vector<bool> watched, std::vector<unsigned> important)
  {
    unsigned old_size = nodes.size(), idx = old_size;

    while (idx != variables)
    {
      Node n = Node(1.0, DLCS_score[idx], watched[idx], important[idx]);
      nodes.push_back(n);
      idx++;  
    }
    while (idx != old_size)
    {
      push(--idx);
    }
  }

  void pop()
  {
    unsigned res = root;
    assert(res != INVALID);
    unsigned child = nodes[res].child;
    root = collapse(child);
    assert(!contains(res));
  }

  void bump(unsigned idx)
  {
    auto &node = nodes[idx];
    double new_score = node.score + increment;
    update(idx, new_score);
    if (new_score > limit)
      rescore(new_score);
  }

  void bump()
  {
    increment *= factor;
    if (increment > limit)
      rescore(increment);
  }

  unsigned get_size(){
    return nodes.size();
  }

  void get_DLCS_scores(){
    // printf all scores
    printf("DLCS SCORES\n");
    for (unsigned i = 0; i < nodes.size(); i++)
    {
      printf("NODE %d: %d (watched? %d)\n", i, nodes[i].DLCS, nodes[i].watched);
    }
  }
};

#endif
