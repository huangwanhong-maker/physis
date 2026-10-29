// The committed versions of an evolving field (decision 0012; docs/theory/fields.md §7; formalism §5.6).
//
// A field that a law advances is published as a sequence of versions, each a field of its own, committed
// at the instant it describes and never changed afterwards. A query at t reads the latest version
// committed at or before t; until the first is committed there is none (missing). A reader holds the
// version it read (a shared pointer to a constant field), so that a version committed while it reads, by
// the solver's next step, does not change what it sees: the observers released at one instant read one
// state (formalism Proposition 8.6).

#pragma once

#include <cstddef>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>
#include <vector>

#include <physis/fields/field.hpp>

namespace physis::fields {

template <PointEvaluable F>
class History {
 public:
  using value_type = typename F::value_type;
  using point_type = typename F::point_type;

  explicit History(Descriptor d) : d_(std::move(d)) {}
  History(const History& o) : d_(o.d_) {
    std::lock_guard<std::mutex> lock(o.mutex_);
    commits_ = o.commits_;
  }

  const Descriptor& descriptor() const { return d_; }

  // Commits a version at t; the instants of the commits do not decrease.
  void commit(double t, F version) {
    auto v = std::make_shared<const F>(std::move(version));
    std::lock_guard<std::mutex> lock(mutex_);
    if (!commits_.empty() && t < commits_.back().first) throw std::invalid_argument("physis::fields::History: a version committed before the last");
    commits_.emplace_back(t, std::move(v));
  }

  std::size_t versions() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return commits_.size();
  }

  // The latest version committed at or before t, with its instant; null before the first.
  std::pair<double, std::shared_ptr<const F>> version_at(double t) const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::pair<double, std::shared_ptr<const F>> found{0.0, nullptr};
    for (const auto& c : commits_) {
      if (c.first > t) break;
      found = c;
    }
    return found;
  }

  Evaluation<value_type> at(const Query<point_type>& q) const {
    const auto [when, v] = version_at(q.t);
    if (!v) return refused<value_type>(Quality::missing, "no version of the field was committed by this instant");
    auto e = v->at(q);
    if (e.ok() && e.note.empty() && when != q.t) e.note = "the latest version committed by the instant";
    return e;
  }

 private:
  Descriptor d_;
  mutable std::mutex mutex_;
  std::vector<std::pair<double, std::shared_ptr<const F>>> commits_;
};

}  // namespace physis::fields
