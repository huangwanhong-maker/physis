// Frames of reference as a graph: the frame groupoid of a world, and the holonomy of its loops.
//
// The frames of a world are the objects of a groupoid whose morphisms are the changes of
// coordinates between them: time-dependent rigid motions (kinematics.hpp). It is presented
// by a graph: each frame but the root is given by one edge to a parent frame, so that the
// frames form a tree, and the transformation between any two is the composition along the
// unique path between them through their lowest common ancestor. This is how NAIF's SPICE
// toolkit (Frames Required Reading) and Orekit (Maisonobe, Pommier and Parraud 2010, "Orekit:
// an open source library for operational flight dynamics") organise frames, and it makes
// the groupoid consistent by construction: on a tree every loop is trivial.
//
// A model of the world often gives a second route between two frames: the celestial and
// terrestrial frames of the Earth are joined both through the celestial intermediate origin
// and through the equinox (IERS Conventions 2010, ch. 5); a satellite's orbit frame through
// two models of the Earth's orientation. Such a route is a *relation*, an extra edge that
// closes a loop with the tree. Going round the loop gives the transformation
//   H = T_tree(b → a) ∘ E(a → b),
// the holonomy of the loop: the identity exactly when the two routes agree. A graph of
// frames with values in SE(3) is a lattice gauge field on the graph (physis `gauge`), a
// change of the frames' origins and axes is a gauge transformation, under which H is
// conjugated; the deviation of H from the identity, the angle of its rotation and the slide
// along its screw axis (kinematics.hpp, Chasles), is gauge invariant. Measured, it is the
// disagreement of the two models; the formalism of Cosmopoiesis (§6) reports it as a check of
// every conversion.
//
// The graph is generic in the time argument: Time is whatever the edges read (a time scale
// of the metrology domain, seconds of a simulation). Edges are functions of it and may throw
// (std::domain_error) where their model or data do not reach.

#pragma once

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <physis/geometry/kinematics.hpp>

namespace physis::geometry {

template <class Time>
class FrameGraph {
 public:
  // A change of coordinates from one frame to another at a time.
  using Edge = std::function<StateTransform(const Time&)>;

  struct Frame {
    std::string name, about;
    int parent = -1;  // -1 for the root
    Edge to_parent;   // coordinates in this frame → coordinates in the parent
    int depth = 0;
  };
  struct Relation {
    std::string name, about;
    int from = -1, to = -1;
    Edge transform;  // coordinates in `from` → coordinates in `to`, by a second route
  };

  // The root: the frame every other is related to (for the Solar System, the ICRS).
  int root(std::string name, std::string about = "") {
    if (!frames_.empty()) throw std::logic_error("FrameGraph: the root is the first frame");
    frames_.push_back({std::move(name), std::move(about), -1, nullptr, 0});
    return 0;
  }
  // A frame given by its change of coordinates to a parent.
  int add(std::string name, int parent, Edge to_parent, std::string about = "") {
    check(parent);
    if (find(name) >= 0) throw std::invalid_argument("FrameGraph: a frame named " + name + " exists");
    frames_.push_back({std::move(name), std::move(about), parent, std::move(to_parent), frames_[parent].depth + 1});
    return int(frames_.size()) - 1;
  }
  int add(std::string name, std::string_view parent, Edge to_parent, std::string about = "") {
    return add(std::move(name), require(parent), std::move(to_parent), std::move(about));
  }
  // A second route between two frames, closing a loop with the tree.
  int relate(std::string name, int from, int to, Edge transform, std::string about = "") {
    check(from);
    check(to);
    relations_.push_back({std::move(name), std::move(about), from, to, std::move(transform)});
    return int(relations_.size()) - 1;
  }

  int find(std::string_view name) const {
    for (std::size_t i = 0; i < frames_.size(); i++)
      if (frames_[i].name == name) return int(i);
    return -1;
  }
  int require(std::string_view name) const {
    const int i = find(name);
    if (i < 0) throw std::invalid_argument("FrameGraph: no frame named " + std::string(name));
    return i;
  }
  const std::vector<Frame>& frames() const { return frames_; }
  const std::vector<Relation>& relations() const { return relations_; }
  const Frame& frame(int i) const {
    check(i);
    return frames_[i];
  }

  // The frames from `from` up to the common ancestor and down to `to`, inclusive.
  std::vector<int> path(int from, int to) const {
    check(from);
    check(to);
    std::vector<int> up{from}, down{to};
    int a = from, b = to;
    while (frames_[a].depth > frames_[b].depth) up.push_back(a = frames_[a].parent);
    while (frames_[b].depth > frames_[a].depth) down.push_back(b = frames_[b].parent);
    while (a != b) {
      up.push_back(a = frames_[a].parent);
      down.push_back(b = frames_[b].parent);
    }
    down.pop_back();  // the common ancestor is already the last of `up`
    up.insert(up.end(), down.rbegin(), down.rend());
    return up;
  }

  // The change of coordinates from frame `from` to frame `to` at a time, along the tree.
  StateTransform transform(int from, int to, const Time& t) const {
    check(from);
    check(to);
    StateTransform up, down;  // from → ancestor, to → ancestor
    int a = from, b = to;
    while (frames_[a].depth > frames_[b].depth) {
      up = frames_[a].to_parent(t) * up;
      a = frames_[a].parent;
    }
    while (frames_[b].depth > frames_[a].depth) {
      down = frames_[b].to_parent(t) * down;
      b = frames_[b].parent;
    }
    while (a != b) {
      up = frames_[a].to_parent(t) * up;
      a = frames_[a].parent;
      down = frames_[b].to_parent(t) * down;
      b = frames_[b].parent;
    }
    return down.inverse() * up;
  }
  StateTransform transform(std::string_view from, std::string_view to, const Time& t) const {
    return transform(require(from), require(to), t);
  }

  // The holonomy of the loop a relation closes: the relation's route from `from` to `to`,
  // then the tree's route back. The identity when the two routes agree.
  StateTransform holonomy(int relation, const Time& t) const {
    if (relation < 0 || relation >= int(relations_.size())) throw std::out_of_range("FrameGraph: no such relation");
    const Relation& r = relations_[relation];
    return transform(r.to, r.from, t) * r.transform(t);
  }

 private:
  void check(int i) const {
    if (i < 0 || i >= int(frames_.size())) throw std::out_of_range("FrameGraph: no such frame");
  }
  std::vector<Frame> frames_;
  std::vector<Relation> relations_;
};

}  // namespace physis::geometry
