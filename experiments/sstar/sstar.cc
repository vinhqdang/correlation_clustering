// sstar: star-packing lower bound for cluster editing on the sparse support.
//
// A star (c; L) is a centre c and at least two leaves L with every centre-leaf
// pair an edge and every leaf-leaf pair a non-edge.  Stars that share no pair
// give the lower bound sum (|L| - 1).  Every pair of a star is an edge or at
// distance two, so the packing lives on the distance-two support P and never
// needs the n x n matrices of a dense implementation.
//
// The heuristic follows the star bound described for the KaPoCE solver
// (Blasius et al., SEA 2022): a greedy start that colours the free
// neighbourhood of every vertex, then local search over the stars in random
// order with three moves: merging two stars with the same centre (+1),
// removing one leaf and inserting two stars or leaves through the freed pairs
// (+1), and a plateau move that replaces the removed leaf by another
// candidate.  This is an independent implementation on hash maps over the
// used pairs; memory is O(n + m + used pairs).
//
// usage: sstar STARS_OUT TIME_LIMIT SEED [MIN_TIME] [MAX_UNCHANGED] < graph.gr
//   graph.gr   PACE format (p cep n m, 1-indexed edges)
//   STARS_OUT  one star per line: centre, then the leaves (0-indexed)
//   stops when rounds > MAX_UNCHANGED * improving rounds (default 5), but not
//   before MIN_TIME seconds (default 0), and at the latest after TIME_LIMIT.
// prints "init VALUE SECONDS" and "star VALUE SECONDS".
// environment: SSTAR_KMAX caps the leaves per star (the pairs of a star grow
// quadratically; only needed on graphs with very high degrees).
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace std;
using u64 = uint64_t;

// open addressing map u64 -> int32 with backward-shift deletion
struct PairMap {
    vector<u64> keys;
    vector<int32_t> vals;
    size_t mask = 0, count = 0;
    static constexpr u64 EMPTY = ~0ull;
    explicit PairMap(size_t cap = 1024) { init(cap); }
    void init(size_t cap) {
        size_t c = 16;
        while (c < 2 * cap) c <<= 1;
        keys.assign(c, EMPTY);
        vals.assign(c, 0);
        mask = c - 1;
        count = 0;
    }
    static inline size_t h(u64 k) {
        k ^= k >> 33; k *= 0xff51afd7ed558ccdull; k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53ull; k ^= k >> 33;
        return (size_t)k;
    }
    void grow() {
        vector<u64> ok = move(keys);
        vector<int32_t> ov = move(vals);
        init(ok.size());
        for (size_t i = 0; i < ok.size(); ++i)
            if (ok[i] != EMPTY) put(ok[i], ov[i]);
    }
    inline int32_t get(u64 k, int32_t dflt = -1) const {
        size_t i = h(k) & mask;
        while (keys[i] != EMPTY) {
            if (keys[i] == k) return vals[i];
            i = (i + 1) & mask;
        }
        return dflt;
    }
    void put(u64 k, int32_t v) {
        if (2 * (count + 1) > keys.size()) grow();
        size_t i = h(k) & mask;
        while (keys[i] != EMPTY) {
            if (keys[i] == k) { vals[i] = v; return; }
            i = (i + 1) & mask;
        }
        keys[i] = k; vals[i] = v; ++count;
    }
    void erase(u64 k) {
        size_t i = h(k) & mask;
        while (keys[i] != EMPTY && keys[i] != k) i = (i + 1) & mask;
        if (keys[i] == EMPTY) return;
        size_t j = i;
        for (;;) {
            j = (j + 1) & mask;
            if (keys[j] == EMPTY) break;
            size_t home = h(keys[j]) & mask;
            // move keys[j] into the hole at i if its home is not in (i, j]
            bool between = (i <= j) ? (i < home && home <= j) : (i < home || home <= j);
            if (!between) { keys[i] = keys[j]; vals[i] = vals[j]; i = j; }
        }
        keys[i] = EMPTY; --count;
    }
};

struct Graph {
    int n = 0;
    vector<int64_t> off;
    vector<int> adj;
    int deg(int u) const { return (int)(off[u + 1] - off[u]); }
    bool is_edge(int u, int v) const {
        if (deg(u) > deg(v)) swap(u, v);
        return binary_search(adj.begin() + off[u], adj.begin() + off[u + 1], v);
    }
    int common(int u, int v) const {
        int c = 0;
        int64_t i = off[u], j = off[v];
        while (i < off[u + 1] && j < off[v + 1]) {
            if (adj[i] < adj[j]) ++i;
            else if (adj[i] > adj[j]) ++j;
            else { ++c; ++i; ++j; }
        }
        return c;
    }
};

static Graph read_pace(FILE *f) {
    Graph g;
    char line[512];
    long long m = 0;
    vector<pair<int, int>> e;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'c' || line[0] == '\n') continue;
        if (line[0] == 'p') {
            char a[16], b[16];
            sscanf(line, "%15s %15s %d %lld", a, b, &g.n, &m);
            e.reserve(m);
            continue;
        }
        int u, v;
        if (sscanf(line, "%d %d", &u, &v) != 2) continue;
        --u; --v;
        if (u == v) continue;
        e.emplace_back(u, v);
    }
    vector<int64_t> d(g.n + 1, 0);
    for (auto [u, v] : e) { d[u + 1]++; d[v + 1]++; }
    for (int i = 0; i < g.n; ++i) d[i + 1] += d[i];
    g.off = d;
    g.adj.assign(d[g.n], 0);
    vector<int64_t> pos(d.begin(), d.end() - 1);
    for (auto [u, v] : e) { g.adj[pos[u]++] = v; g.adj[pos[v]++] = u; }
    // sort and remove parallel edges
    vector<int64_t> noff(g.n + 1, 0);
    int64_t w = 0;
    for (int u = 0; u < g.n; ++u) {
        auto b = g.adj.begin() + g.off[u], en = g.adj.begin() + g.off[u + 1];
        sort(b, en);
        auto last = unique(b, en);
        noff[u] = w;
        for (auto it = b; it != last; ++it) g.adj[w++] = *it;
    }
    noff[g.n] = w;
    g.adj.resize(w);
    g.off = noff;
    return g;
}

struct Star {
    int c = -1;
    vector<int> L;  // leaves, unordered
    int version = 0;
    bool alive = false;
};

struct Cand {
    // type 0: new star (c; a, b); type 1: leaf v added to star sid
    int type, c, a, b, sid, ver, v, pair;  // pair: index of the freed pair
};

struct Packing {
    const Graph &g;
    int n;
    PairMap owner;  // used pair -> star id
    vector<Star> stars;
    vector<int> free_ids;
    vector<vector<int>> by_center;
    long long value = 0;
    mt19937_64 gen;
    int cand_cap = 16, total_cap = 256;
    size_t kmax = (size_t)-1;  // at most this many leaves per star (memory)
    PairMap p3cache;

    Packing(const Graph &g_, u64 seed) : g(g_), n(g_.n), owner(1 << 16), by_center(g_.n),
                                          gen(seed), p3cache(1 << 12) {}

    inline u64 key(int a, int b) const {
        if (a > b) swap(a, b);
        return (u64)a * (u64)n + (u64)b;
    }
    inline bool pfree(int a, int b) const { return owner.get(key(a, b)) < 0; }

    int new_id() {
        if (!free_ids.empty()) { int i = free_ids.back(); free_ids.pop_back(); return i; }
        stars.emplace_back();
        return (int)stars.size() - 1;
    }
    int add_star(int c, const vector<int> &L) {
        int id = new_id();
        Star &s = stars[id];
        s.c = c; s.L = L; s.alive = true; s.version++;
        for (size_t i = 0; i < L.size(); ++i) {
            owner.put(key(c, L[i]), id);
            for (size_t j = i + 1; j < L.size(); ++j) owner.put(key(L[i], L[j]), id);
        }
        by_center[c].push_back(id);
        value += (long long)L.size() - 1;
        return id;
    }
    void remove_star(int id) {
        Star &s = stars[id];
        for (size_t i = 0; i < s.L.size(); ++i) {
            owner.erase(key(s.c, s.L[i]));
            for (size_t j = i + 1; j < s.L.size(); ++j) owner.erase(key(s.L[i], s.L[j]));
        }
        auto &bc = by_center[s.c];
        bc.erase(find(bc.begin(), bc.end(), id));
        value -= (long long)s.L.size() - 1;
        s.alive = false; s.version++; s.L.clear();
        free_ids.push_back(id);
    }
    void add_leaf(int id, int v) {
        Star &s = stars[id];
        owner.put(key(s.c, v), id);
        for (int x : s.L) owner.put(key(v, x), id);
        s.L.push_back(v); s.version++;
        value += 1;
    }
    void remove_leaf(int id, int v) {
        Star &s = stars[id];
        auto it = find(s.L.begin(), s.L.end(), v);
        s.L.erase(it);
        owner.erase(key(s.c, v));
        for (int x : s.L) owner.erase(key(v, x));
        s.version++;
        value -= 1;
    }
    // v can join the leaves of s: edge c-v free, v non-adjacent to and free with every leaf
    bool leaf_ok(const Star &s, int v) const {
        if (s.L.size() >= kmax || v == s.c || !pfree(s.c, v)) return false;
        for (int x : s.L)
            if (x == v || g.is_edge(v, x) || !pfree(v, x)) return false;
        return true;
    }
    bool can_add(const Cand &k) const {
        if (k.type == 0) return pfree(k.c, k.a) && pfree(k.c, k.b) && pfree(k.a, k.b);
        // checked on the current leaves: the star may have changed since
        const Star &s = stars[k.sid];
        return s.alive && s.c == k.c && leaf_ok(s, k.v);
    }
    // apply; returns the id of a new star (type 0) or the star id (type 1)
    int apply(const Cand &k) {
        if (k.type == 0) return add_star(k.c, {k.a, k.b});
        add_leaf(k.sid, k.v);
        return k.sid;
    }
    void undo(const Cand &k, int id) {
        if (k.type == 0) remove_star(id);
        else remove_leaf(id, k.v);
    }

    int p3count(int a, int b) {
        u64 kk = key(a, b);
        int c = p3cache.get(kk, -1);
        if (c >= 0) return c;
        int cm = g.common(a, b);
        c = g.is_edge(a, b) ? g.deg(a) + g.deg(b) - 2 - 2 * cm : cm;
        p3cache.put(kk, c);
        return c;
    }
    int score(const Cand &k) {
        if (k.type == 0) return p3count(k.c, k.a) + p3count(k.c, k.b) + p3count(k.a, k.b);
        const Star &s = stars[k.sid];
        int t = p3count(k.v, s.c);
        for (int x : s.L) t += p3count(k.v, x);
        return t;
    }

    // candidates that use the free pair (a, b); at most cand_cap, from a random offset
    void candidates(int a, int b, int pidx, vector<Cand> &out) {
        size_t start = out.size();
        auto full = [&]() { return out.size() - start >= (size_t)cand_cap; };
        if (g.is_edge(a, b)) {
            for (int t = 0; t < 2 && !full(); ++t) {
                int x = t ? b : a, y = t ? a : b;  // centre x, leaf y
                int d = g.deg(x);
                if (d > 1) {
                    int r = (int)(gen() % d);
                    for (int q = 0; q < d && !full(); ++q) {
                        int z = g.adj[g.off[x] + (q + r) % d];
                        if (z == y || !pfree(x, z) || g.is_edge(y, z) || !pfree(y, z)) continue;
                        out.push_back({0, x, y, z, -1, 0, -1, pidx});
                    }
                }
                for (int sid : by_center[x]) {
                    if (full()) break;
                    const Star &s = stars[sid];
                    bool ok = s.L.size() < kmax;
                    for (int l : s.L)
                        if (l == y || g.is_edge(y, l) || !pfree(y, l)) { ok = false; break; }
                    if (ok) out.push_back({1, x, -1, -1, sid, s.version, y, pidx});
                }
            }
        } else {
            // new star (z; a, b) with a common neighbour z
            {
                int64_t i = g.off[a], j = g.off[b];
                vector<int> zs;
                while (i < g.off[a + 1] && j < g.off[b + 1]) {
                    if (g.adj[i] < g.adj[j]) ++i;
                    else if (g.adj[i] > g.adj[j]) ++j;
                    else { zs.push_back(g.adj[i]); ++i; ++j; }
                }
                int d = (int)zs.size();
                if (d) {
                    int r = (int)(gen() % d);
                    for (int q = 0; q < d && !full(); ++q) {
                        int z = zs[(q + r) % d];
                        if (pfree(z, a) && pfree(z, b)) out.push_back({0, z, a, b, -1, 0, -1, pidx});
                    }
                }
            }
            // leaf e added to a star with centre c that has leaf o
            for (int t = 0; t < 2 && !full(); ++t) {
                int e = t ? b : a, o = t ? a : b;
                int d = g.deg(e);
                if (!d) continue;
                int r = (int)(gen() % d);
                for (int q = 0; q < d && !full(); ++q) {
                    int c = g.adj[g.off[e] + (q + r) % d];
                    if (!pfree(c, e)) continue;
                    int sid = owner.get(key(c, o));
                    if (sid < 0) continue;
                    const Star &s = stars[sid];
                    if (s.c != c) continue;
                    if (leaf_ok(s, e)) out.push_back({1, c, -1, -1, sid, s.version, e, pidx});
                }
            }
        }
    }

    // ------------------------------------------------------------------ start
    void greedy_start() {
        vector<int> order(n);
        for (int i = 0; i < n; ++i) order[i] = i;
        vector<u64> tie(n);
        for (auto &t : tie) t = gen();
        sort(order.begin(), order.end(), [&](int a, int b) {
            return g.deg(a) != g.deg(b) ? g.deg(a) > g.deg(b) : tie[a] < tie[b];
        });
        const int CAPD = 3000;
        vector<int> F, color, cdeg, ord;
        vector<vector<int>> cg;
        for (int u : order) {
            F.clear();
            for (int64_t p = g.off[u]; p < g.off[u + 1]; ++p)
                if (pfree(u, g.adj[p])) F.push_back(g.adj[p]);
            if (F.size() < 2) continue;
            if ((int)F.size() > CAPD) {
                shuffle(F.begin(), F.end(), gen);
                F.resize(CAPD);
            }
            int f = (int)F.size();
            // conflict graph: adjacent, or the non-edge is already used
            cg.assign(f, {});
            for (int i = 0; i < f; ++i)
                for (int j = i + 1; j < f; ++j)
                    if (g.is_edge(F[i], F[j]) || !pfree(F[i], F[j])) {
                        cg[i].push_back(j); cg[j].push_back(i);
                    }
            // degeneracy order (smallest degree first), coloured in reverse
            cdeg.assign(f, 0);
            int maxd = 0;
            for (int i = 0; i < f; ++i) { cdeg[i] = (int)cg[i].size(); maxd = max(maxd, cdeg[i]); }
            vector<vector<int>> bucket(maxd + 1);
            for (int i = 0; i < f; ++i) bucket[cdeg[i]].push_back(i);
            vector<char> done(f, 0);
            ord.clear();
            int dcur = 0;
            while ((int)ord.size() < f) {
                if (dcur > maxd) dcur = 0;
                bool found = false;
                for (int d = 0; d <= maxd && !found; ++d) {
                    while (!bucket[d].empty()) {
                        int x = bucket[d].back(); bucket[d].pop_back();
                        if (done[x] || cdeg[x] != d) continue;
                        done[x] = 1; ord.push_back(x); found = true;
                        for (int y : cg[x])
                            if (!done[y]) { cdeg[y]--; bucket[cdeg[y]].push_back(y); }
                        break;
                    }
                }
                if (!found) break;
            }
            reverse(ord.begin(), ord.end());
            color.assign(f, -1);
            int ncol = 0;
            vector<char> used;
            for (int x : ord) {
                used.assign(ncol + 1, 0);
                for (int y : cg[x]) if (color[y] >= 0) used[color[y]] = 1;
                int c = 0;
                while (used[c]) ++c;
                color[x] = c;
                ncol = max(ncol, c + 1);
            }
            vector<vector<int>> cls(ncol);
            for (int i = 0; i < f; ++i) cls[color[i]].push_back(F[i]);
            for (auto &L : cls)
                for (size_t a = 0; a + 1 < L.size();) {
                    size_t b = L.size() - a > kmax ? a + kmax : L.size();
                    if (b - a >= 2) add_star(u, vector<int>(L.begin() + a, L.begin() + b));
                    a = b;
                }
        }
    }

    // --------------------------------------------------------- local search
    bool try_merge(int id) {
        Star &s = stars[id];
        for (int t : vector<int>(by_center[s.c])) {
            if (t == id) continue;
            const Star &o = stars[t];
            bool ok = s.L.size() + o.L.size() <= kmax;
            for (int x : s.L) {
                for (int y : o.L)
                    if (x == y || g.is_edge(x, y) || !pfree(x, y)) { ok = false; break; }
                if (!ok) break;
            }
            if (!ok) continue;
            vector<int> leaves = o.L;
            remove_star(t);
            for (int y : leaves) add_leaf(id, y);
            return true;
        }
        return false;
    }

    void try_improve(int id) {
        if (try_merge(id)) return;
        vector<int> leaves = stars[id].L;
        shuffle(leaves.begin(), leaves.end(), gen);
        vector<Cand> C;
        vector<pair<int, int>> freed;
        for (int v : leaves) {
            Star &s = stars[id];
            if (!s.alive || find(s.L.begin(), s.L.end(), v) == s.L.end()) break;
            bool p3 = s.L.size() == 2;
            int c = s.c, a = s.L[0], b = s.L[1];
            freed.clear();
            if (p3) {
                remove_star(id);
                freed = {{c, a}, {c, b}, {a, b}};
            } else {
                remove_leaf(id, v);
                freed.push_back({c, v});
                for (int x : stars[id].L) freed.push_back({v, x});
                if ((int)freed.size() > total_cap / 4) {
                    shuffle(freed.begin() + 1, freed.end(), gen);
                    freed.resize(total_cap / 4);
                }
            }
            C.clear();
            for (int p = 0; p < (int)freed.size() && (int)C.size() < total_cap; ++p)
                candidates(freed[p].first, freed[p].second, p, C);
            // two candidates through different freed pairs: +1
            bool two = false;
            for (size_t i = 0; i < C.size() && !two; ++i) {
                if (!can_add(C[i])) continue;
                int ida = apply(C[i]);
                for (size_t j = 0; j < C.size(); ++j) {
                    if (C[j].pair == C[i].pair || !can_add(C[j])) continue;
                    apply(C[j]);
                    two = true;
                }
                if (!two) undo(C[i], ida);
            }
            if (!two) {
                vector<size_t> ok;
                for (size_t i = 0; i < C.size(); ++i) if (can_add(C[i])) ok.push_back(i);
                if (!ok.empty()) {
                    size_t pick;
                    if (uniform_real_distribution<double>(0, 1)(gen) < 0.8) {
                        pick = ok[0];
                        int best = score(C[pick]);
                        for (size_t q = 1; q < ok.size(); ++q) {
                            int sc = score(C[ok[q]]);
                            if (sc < best) { best = sc; pick = ok[q]; }
                        }
                    } else {
                        pick = ok[gen() % ok.size()];
                    }
                    apply(C[pick]);
                } else if (p3) {
                    id = add_star(c, {a, b});
                } else {
                    add_leaf(id, v);
                }
            }
            if (two || p3) return;
        }
    }

    bool check() const {
        long long val = 0;
        PairMap seen(1 << 12);
        for (size_t id = 0; id < stars.size(); ++id) {
            const Star &s = stars[id];
            if (!s.alive) continue;
            if (s.L.size() < 2) return false;
            val += (long long)s.L.size() - 1;
            for (size_t i = 0; i < s.L.size(); ++i) {
                if (!g.is_edge(s.c, s.L[i])) return false;
                if (seen.get(key(s.c, s.L[i])) >= 0) return false;
                seen.put(key(s.c, s.L[i]), 1);
                for (size_t j = i + 1; j < s.L.size(); ++j) {
                    if (s.L[i] == s.L[j] || g.is_edge(s.L[i], s.L[j])) return false;
                    if (seen.get(key(s.L[i], s.L[j])) >= 0) return false;
                    seen.put(key(s.L[i], s.L[j]), 1);
                }
            }
        }
        return val == value;
    }
};

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: sstar STARS_OUT TIME_LIMIT SEED [MIN_TIME] [MAX_UNCHANGED] < g.gr\n");
        return 2;
    }
    auto t0 = chrono::steady_clock::now();
    auto el = [&]() { return chrono::duration<double>(chrono::steady_clock::now() - t0).count(); };
    double T = atof(argv[2]);
    u64 seed = strtoull(argv[3], nullptr, 10);
    double min_time = argc > 4 ? atof(argv[4]) : 0.0;
    int max_unchanged = argc > 5 ? atoi(argv[5]) : 5;
    Graph g = read_pace(stdin);
    Packing P(g, seed);
    if (getenv("SSTAR_CAND")) P.cand_cap = atoi(getenv("SSTAR_CAND"));
    if (getenv("SSTAR_TOTAL")) P.total_cap = atoi(getenv("SSTAR_TOTAL"));
    if (getenv("SSTAR_KMAX")) P.kmax = (size_t)atol(getenv("SSTAR_KMAX"));
    P.greedy_start();
    printf("init %lld %.3f\n", P.value, el());
    fflush(stdout);
    long long rounds = 0, improvements = 0;
    while (el() < T && ((long long)max_unchanged * improvements >= rounds || el() < min_time)) {
        long long old = P.value;
        vector<int> ids;
        for (size_t i = 0; i < P.stars.size(); ++i) if (P.stars[i].alive) ids.push_back((int)i);
        shuffle(ids.begin(), ids.end(), P.gen);
        for (size_t q = 0; q < ids.size(); ++q) {
            if (P.stars[ids[q]].alive) P.try_improve(ids[q]);
            if ((q & 255) == 0 && el() >= T) break;
        }
        ++rounds;
        improvements += P.value > old;
        fprintf(stderr, "round %lld value %lld time %.1f\n", rounds, P.value, el());
    }
    if (!P.check()) {
        fprintf(stderr, "internal check failed\n");
        return 1;
    }
    FILE *out = fopen(argv[1], "w");
    for (auto &s : P.stars) {
        if (!s.alive) continue;
        vector<int> L = s.L;
        sort(L.begin(), L.end());
        fprintf(out, "%d", s.c);
        for (int x : L) fprintf(out, " %d", x);
        fputc('\n', out);
    }
    fclose(out);
    printf("star %lld %.3f\n", P.value, el());
    printf("rounds %lld %lld\n", rounds, improvements);
    return 0;
}
