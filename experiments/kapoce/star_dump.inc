
// ---- appended by experiments/kapoce/build_lbounds.sh (not part of KaPoCE) ----
// The root star packing of star_bound(), unchanged, with its stars written to
// `out`, one per line, centre first; returns the same bound as star_bound().
int star_bound_dump(const Instance &inst, int limit, std::ostream &out) {
    auto bound = star_bound_packing(inst, limit);
    for (auto &star : bound.stars_in_random_order()) {
        for (size_t i = 0; i < star.nodes.size(); ++i)
            out << (i ? " " : "") << star.nodes[i];
        out << "\n";
    }
    return bound.potential.get_bound();
}
