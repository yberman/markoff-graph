/*
 * libmarkoff.c
 *
 * Small single-threaded C backend for the graph of solutions modulo a prime p to
 *
 *     x^2 + y^2 + z^2 = x*y*z + A*x + B*y + C*z + D  (mod p).
 *
 * The Python wrapper exposes triples only.  The C ABI returns encoded vertices
 *
 *     V = x + p*y + p*p*z,     0 <= x,y,z < p,
 *
 * as uint32_t arrays.  We reject p > 1621 so every possible V fits in uint32_t.
 *
 * Internal lookup design:
 *
 *   - A solution node stores x,y,z as uint16_t and its union-find parent index.
 *   - Three dense two-slot tables map coordinate pairs to node indices:
 *
 *         xy[(x,y)] -> up to two nodes, differing in z
 *         xz[(x,z)] -> up to two nodes, differing in y
 *         yz[(y,z)] -> up to two nodes, differing in x
 *
 *     This matches the Vieta moves: each move keeps two coordinates fixed and
 *     replaces the third root of the corresponding quadratic.
 *
 * Exported ABI:
 *
 *     int build_graph(int A, int B, int C, int D, int prime,
 *                     uint32_t capacity, int *V_len,
 *                     uint32_t *Vs, uint32_t *root)
 *
 * The caller owns Vs and root, each of length at least capacity.
 * On success, Vs[0..V_len) are encoded vertices and root[0..V_len) are encoded
 * component roots parallel to Vs.  The return value is the component count.
 *
 * Return values:
 *   >=0 connected component count
 *   -1  invalid input pointer, modulus <= 1, or capacity too large for int
 *   -2  too many solutions for caller-provided capacity
 *   -3  allocation failure
 *   -4  internal error: a Vieta neighbor was not found among solutions
 *   -5  modulus > 1621, so uint32 vertex encoding is not guaranteed safe
 *   -6  modulus is composite; this builder requires a prime field
 *   -7  internal error: more than two solutions share a fixed coordinate pair
 */

#include <stdint.h>
#include <stdlib.h>
#include <limits.h>

#if defined(_WIN32) || defined(__CYGWIN__)
#  define MARKOFF_EXPORT __declspec(dllexport)
#else
#  define MARKOFF_EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MARKOFF_MAX_PRIME_U32 1621
#define U32_MISSING UINT32_MAX

typedef struct {
    uint16_t x;
    uint16_t y;
    uint16_t z;
    uint32_t parent;
} SolutionNode;

static int is_prime_small(int n) {
    if (n < 2) return 0;
    if (n == 2 || n == 3) return 1;
    if ((n % 2) == 0) return 0;
    for (int d = 3; d * d <= n; d += 2) {
        if ((n % d) == 0) return 0;
    }
    return 1;
}

static int modp_ll(long long a, int p) {
    long long r = a % (long long)p;
    if (r < 0) r += p;
    return (int)r;
}

static uint32_t encode_v_u32(int x, int y, int z, int p) {
    uint64_t pp = (uint64_t)(uint32_t)p * (uint64_t)(uint32_t)p;
    uint64_t v = (uint64_t)(uint32_t)x
               + (uint64_t)(uint32_t)p * (uint64_t)(uint32_t)y
               + pp * (uint64_t)(uint32_t)z;
    return (uint32_t)v;
}

static size_t pair_offset(int a, int b, int p) {
    return 2u * ((size_t)(uint32_t)a * (size_t)(uint32_t)p + (size_t)(uint32_t)b);
}

static void fill_missing(uint32_t *table, size_t n) {
    for (size_t i = 0; i < n; ++i) table[i] = U32_MISSING;
}

static int table_insert(uint32_t *table, int p, int a, int b, uint32_t node_index) {
    size_t off = pair_offset(a, b, p);
    if (table[off] == U32_MISSING) {
        table[off] = node_index;
        return 1;
    }
    if (table[off + 1u] == U32_MISSING) {
        table[off + 1u] = node_index;
        return 1;
    }
    return 0;
}

static uint32_t find_by_x(const uint32_t *yz, const SolutionNode *nodes,
                          int p, int y, int z, int wanted_x) {
    size_t off = pair_offset(y, z, p);
    for (int k = 0; k < 2; ++k) {
        uint32_t ix = yz[off + (size_t)k];
        if (ix != U32_MISSING && nodes[ix].x == (uint16_t)wanted_x) return ix;
    }
    return U32_MISSING;
}

static uint32_t find_by_y(const uint32_t *xz, const SolutionNode *nodes,
                          int p, int x, int z, int wanted_y) {
    size_t off = pair_offset(x, z, p);
    for (int k = 0; k < 2; ++k) {
        uint32_t ix = xz[off + (size_t)k];
        if (ix != U32_MISSING && nodes[ix].y == (uint16_t)wanted_y) return ix;
    }
    return U32_MISSING;
}

static uint32_t find_by_z(const uint32_t *xy, const SolutionNode *nodes,
                          int p, int x, int y, int wanted_z) {
    size_t off = pair_offset(x, y, p);
    for (int k = 0; k < 2; ++k) {
        uint32_t ix = xy[off + (size_t)k];
        if (ix != U32_MISSING && nodes[ix].z == (uint16_t)wanted_z) return ix;
    }
    return U32_MISSING;
}

static uint32_t uf_find(SolutionNode *nodes, uint32_t a) {
    uint32_t r = a;
    while (nodes[r].parent != r) r = nodes[r].parent;

    while (nodes[a].parent != a) {
        uint32_t b = nodes[a].parent;
        nodes[a].parent = r;
        a = b;
    }
    return r;
}

/* Deterministic union: the smaller root index becomes the root. */
static void uf_union_min(SolutionNode *nodes, uint32_t a, uint32_t b) {
    uint32_t ra = uf_find(nodes, a);
    uint32_t rb = uf_find(nodes, b);
    if (ra == rb) return;
    if (ra < rb) nodes[rb].parent = ra;
    else nodes[ra].parent = rb;
}

static int add_solution(
    SolutionNode *nodes,
    uint32_t *Vs,
    uint32_t *n,
    uint32_t capacity,
    uint32_t *xy,
    uint32_t *xz,
    uint32_t *yz,
    int x,
    int y,
    int z,
    int p
) {
    if (*n >= capacity) return -2;

    uint32_t i = *n;
    nodes[i].x = (uint16_t)x;
    nodes[i].y = (uint16_t)y;
    nodes[i].z = (uint16_t)z;
    nodes[i].parent = i;
    Vs[i] = encode_v_u32(x, y, z, p);

    if (!table_insert(xy, p, x, y, i)) return -7;
    if (!table_insert(xz, p, x, z, i)) return -7;
    if (!table_insert(yz, p, y, z, i)) return -7;

    *n = i + 1u;
    return 1;
}

static void free_all(SolutionNode *nodes, uint32_t *xy, uint32_t *xz, uint32_t *yz,
                     int *sqrt_count, int *sqrt_roots) {
    free(nodes);
    free(xy);
    free(xz);
    free(yz);
    free(sqrt_count);
    free(sqrt_roots);
}

MARKOFF_EXPORT int build_graph(
    int A, int B, int C, int D,
    int prime,
    uint32_t capacity,
    int *V_len,
    uint32_t *Vs,
    uint32_t *root
) {
    if (prime <= 1 || capacity > (uint32_t)INT_MAX ||
        V_len == NULL || Vs == NULL || root == NULL) {
        return -1;
    }
    if (prime > MARKOFF_MAX_PRIME_U32) return -5;
    if (!is_prime_small(prime)) return -6;

    const int p = prime;
    const int Am = modp_ll(A, p);
    const int Bm = modp_ll(B, p);
    const int Cm = modp_ll(C, p);
    const int Dm = modp_ll(D, p);

    const size_t pair_slots = 2u * (size_t)(uint32_t)p * (size_t)(uint32_t)p;
    SolutionNode *nodes = (SolutionNode *)malloc((size_t)capacity * sizeof(SolutionNode));
    uint32_t *xy = (uint32_t *)malloc(pair_slots * sizeof(uint32_t));
    uint32_t *xz = (uint32_t *)malloc(pair_slots * sizeof(uint32_t));
    uint32_t *yz = (uint32_t *)malloc(pair_slots * sizeof(uint32_t));
    if (nodes == NULL || xy == NULL || xz == NULL || yz == NULL) {
        free_all(nodes, xy, xz, yz, NULL, NULL);
        return -3;
    }
    fill_missing(xy, pair_slots);
    fill_missing(xz, pair_slots);
    fill_missing(yz, pair_slots);

    uint32_t n = 0;

    if (p == 2) {
        for (int x = 0; x < p; ++x) {
            for (int y = 0; y < p; ++y) {
                for (int z = 0; z < p; ++z) {
                    int lhs = modp_ll((long long)x*x + (long long)y*y + (long long)z*z, p);
                    int rhs = modp_ll(1LL*x*y*z + (long long)Am*x + (long long)Bm*y +
                                      (long long)Cm*z + Dm, p);
                    if (lhs == rhs) {
                        int ok = add_solution(nodes, Vs, &n, capacity, xy, xz, yz, x, y, z, p);
                        if (ok < 0) {
                            *V_len = (int)n;
                            free_all(nodes, xy, xz, yz, NULL, NULL);
                            return ok;
                        }
                    }
                }
            }
        }
    } else {
        int *sqrt_count = (int *)calloc((size_t)p, sizeof(int));
        int *sqrt_roots = (int *)malloc((size_t)2 * (size_t)p * sizeof(int));
        if (sqrt_count == NULL || sqrt_roots == NULL) {
            free_all(nodes, xy, xz, yz, sqrt_count, sqrt_roots);
            return -3;
        }
        for (int i = 0; i < 2 * p; ++i) sqrt_roots[i] = -1;

        for (int t = 0; t < p; ++t) {
            int s = modp_ll((long long)t * t, p);
            int c = sqrt_count[s];
            if (c < 2) {
                sqrt_roots[2*s + c] = t;
                sqrt_count[s] = c + 1;
            }
        }

        const int inv2 = (p + 1) / 2;

        for (int x = 0; x < p; ++x) {
            for (int y = 0; y < p; ++y) {
                /* z^2 - (x*y + C)*z + (x^2 + y^2 - A*x - B*y - D) = 0 */
                int S = modp_ll(1LL*x*y + Cm, p);
                int T = modp_ll((long long)x*x + (long long)y*y
                              - (long long)Am*x - (long long)Bm*y - Dm, p);
                int disc = modp_ll((long long)S*S - 4LL*T, p);
                int cnt = sqrt_count[disc];

                for (int k = 0; k < cnt; ++k) {
                    int r = sqrt_roots[2*disc + k];
                    int z = modp_ll((long long)(S + r) * inv2, p);
                    int ok = add_solution(nodes, Vs, &n, capacity, xy, xz, yz, x, y, z, p);
                    if (ok < 0) {
                        *V_len = (int)n;
                        free_all(nodes, xy, xz, yz, sqrt_count, sqrt_roots);
                        return ok;
                    }
                }
            }
        }

        free(sqrt_count);
        free(sqrt_roots);
    }

    *V_len = (int)n;
    if (n == 0) {
        free_all(nodes, xy, xz, yz, NULL, NULL);
        return 0;
    }

    for (uint32_t i = 0; i < n; ++i) {
        int x = (int)nodes[i].x;
        int y = (int)nodes[i].y;
        int z = (int)nodes[i].z;

        int sx = modp_ll(1LL*y*z + Am - x, p);
        int sy = modp_ll(1LL*x*z + Bm - y, p);
        int sz = modp_ll(1LL*x*y + Cm - z, p);

        uint32_t ix = find_by_x(yz, nodes, p, y, z, sx);
        uint32_t iy = find_by_y(xz, nodes, p, x, z, sy);
        uint32_t iz = find_by_z(xy, nodes, p, x, y, sz);
        if (ix == U32_MISSING || iy == U32_MISSING || iz == U32_MISSING) {
            free_all(nodes, xy, xz, yz, NULL, NULL);
            return -4;
        }

        uf_union_min(nodes, i, ix);
        uf_union_min(nodes, i, iy);
        uf_union_min(nodes, i, iz);
    }

    for (uint32_t i = 0; i < n; ++i) {
        nodes[i].parent = uf_find(nodes, i);
    }

    int components = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (nodes[i].parent == i) components += 1;
        SolutionNode *r = &nodes[nodes[i].parent];
        root[i] = encode_v_u32((int)r->x, (int)r->y, (int)r->z, p);
    }

    free_all(nodes, xy, xz, yz, NULL, NULL);
    return components;
}

#ifdef __cplusplus
}
#endif
