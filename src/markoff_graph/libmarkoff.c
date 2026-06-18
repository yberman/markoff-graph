/*
 * libmarkoff.c
 *
 * Single-threaded C backend for connected components of the Markoff-type
 * surface over a prime field:
 *
 *     x^2 + y^2 + z^2 = x*y*z + A*x + B*y + C*z + D  (mod p).
 *
 * Public ABI shape:
 *
 *     markoff_build(...) returns an opaque MarkoffGraph* handle.
 *     markoff_nodes(graph) returns a pointer to graph->nodes.
 *     markoff_free(graph) releases all C-owned memory.
 *
 * The permanent graph storage is one natural node array.  Each node stores its
 * coordinates and the index of the root node of its connected component.
 *
 * Temporary pair lookup tables are used only while building components.  They
 * are freed before the graph handle is returned.
 */

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#if defined(_WIN32) || defined(__CYGWIN__)
#  define MARKOFF_EXPORT __declspec(dllexport)
#else
#  define MARKOFF_EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum {
    MARKOFF_OK = 0,
    MARKOFF_ERR_BAD_ARGUMENT = -1,
    MARKOFF_ERR_TOO_MANY_NODES = -2,
    MARKOFF_ERR_ALLOC = -3,
    MARKOFF_ERR_NEIGHBOR_NOT_FOUND = -4,
    MARKOFF_ERR_STORAGE_BOUND = -5,
    MARKOFF_ERR_COMPOSITE_MODULUS = -6,
    MARKOFF_ERR_PAIR_OVERFLOW = -7,
    MARKOFF_ERR_EIGEN_FAILED = -8
};

#define U32_MISSING UINT32_MAX
#define MARKOFF_MAX_PRIME_BY_NODE_COUNT 46340u

typedef struct {
    uint16_t x;
    uint16_t y;
    uint16_t z;
    uint16_t pad;
    uint32_t root;
    uint32_t neighbor0;
    uint32_t neighbor1;
    uint32_t neighbor2;
    double eigenvector;
} MarkoffNode;

typedef struct {
    const MarkoffNode *root;
    uint32_t root_index;
    uint32_t size;
    double eigenvalue;
} MarkoffComponent;

typedef struct {
    int32_t A;
    int32_t B;
    int32_t C;
    int32_t D;
    uint32_t prime;
    uint32_t node_count;
    uint32_t component_count;
    MarkoffNode *nodes;
    MarkoffComponent *components;
} MarkoffGraph;

typedef struct {
    uint32_t size;
    uint32_t nnz;
    uint32_t root_index;
    double *data;
    int *indices;
    int *indptr;
    uint32_t *nodes;
} MarkoffCSR;

typedef struct {
    uint32_t *xy;
    uint32_t *xz;
    uint32_t *yz;
    size_t slots;
} PairTables;

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

static int prime_fits_storage_bound(int p) {
    if (p <= 1) return 0;
    if ((uint32_t)p > (uint32_t)UINT16_MAX) return 0;
    return 2ull * (uint64_t)(uint32_t)p * (uint64_t)(uint32_t)p <= (uint64_t)UINT32_MAX;
}

static uint32_t default_capacity(int p) {
    if (p == 2) return 8u;
    return (uint32_t)(2ull * (uint64_t)(uint32_t)p * (uint64_t)(uint32_t)p);
}

static size_t pair_offset(int a, int b, int p) {
    return 2u * ((size_t)(uint32_t)a * (size_t)(uint32_t)p + (size_t)(uint32_t)b);
}

static void fill_missing(uint32_t *table, size_t n) {
    for (size_t i = 0; i < n; ++i) table[i] = U32_MISSING;
}

static void pair_tables_free(PairTables *tables) {
    if (tables == NULL) return;
    free(tables->xy);
    free(tables->xz);
    free(tables->yz);
    tables->xy = NULL;
    tables->xz = NULL;
    tables->yz = NULL;
    tables->slots = 0;
}

static int pair_tables_init(PairTables *tables, int p) {
    if (tables == NULL) return MARKOFF_ERR_BAD_ARGUMENT;

    tables->xy = NULL;
    tables->xz = NULL;
    tables->yz = NULL;
    tables->slots = 0;

    const uint64_t slots64 = 2ull * (uint64_t)(uint32_t)p * (uint64_t)(uint32_t)p;
    if (slots64 > (uint64_t)SIZE_MAX / sizeof(uint32_t)) return MARKOFF_ERR_BAD_ARGUMENT;

    const size_t slots = (size_t)slots64;
    tables->xy = (uint32_t *)malloc(slots * sizeof(uint32_t));
    tables->xz = (uint32_t *)malloc(slots * sizeof(uint32_t));
    tables->yz = (uint32_t *)malloc(slots * sizeof(uint32_t));
    if (tables->xy == NULL || tables->xz == NULL || tables->yz == NULL) {
        pair_tables_free(tables);
        return MARKOFF_ERR_ALLOC;
    }

    tables->slots = slots;
    fill_missing(tables->xy, slots);
    fill_missing(tables->xz, slots);
    fill_missing(tables->yz, slots);
    return MARKOFF_OK;
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

static uint32_t find_by_x(const uint32_t *yz, const MarkoffNode *nodes,
                          int p, int y, int z, int wanted_x) {
    size_t off = pair_offset(y, z, p);
    for (int k = 0; k < 2; ++k) {
        uint32_t ix = yz[off + (size_t)k];
        if (ix != U32_MISSING && nodes[ix].x == (uint16_t)wanted_x) return ix;
    }
    return U32_MISSING;
}

static uint32_t find_by_y(const uint32_t *xz, const MarkoffNode *nodes,
                          int p, int x, int z, int wanted_y) {
    size_t off = pair_offset(x, z, p);
    for (int k = 0; k < 2; ++k) {
        uint32_t ix = xz[off + (size_t)k];
        if (ix != U32_MISSING && nodes[ix].y == (uint16_t)wanted_y) return ix;
    }
    return U32_MISSING;
}

static uint32_t find_by_z(const uint32_t *xy, const MarkoffNode *nodes,
                          int p, int x, int y, int wanted_z) {
    size_t off = pair_offset(x, y, p);
    for (int k = 0; k < 2; ++k) {
        uint32_t ix = xy[off + (size_t)k];
        if (ix != U32_MISSING && nodes[ix].z == (uint16_t)wanted_z) return ix;
    }
    return U32_MISSING;
}

static uint32_t uf_find(MarkoffNode *nodes, uint32_t a) {
    uint32_t r = a;
    while (nodes[r].root != r) r = nodes[r].root;

    while (nodes[a].root != a) {
        uint32_t b = nodes[a].root;
        nodes[a].root = r;
        a = b;
    }
    return r;
}

/* Deterministic union: the smaller root index becomes the root. */
static void uf_union_min(MarkoffNode *nodes, uint32_t a, uint32_t b) {
    uint32_t ra = uf_find(nodes, a);
    uint32_t rb = uf_find(nodes, b);
    if (ra == rb) return;
    if (ra < rb) nodes[rb].root = ra;
    else nodes[ra].root = rb;
}

static int add_solution(
    MarkoffNode *nodes,
    uint32_t *n,
    uint32_t capacity,
    PairTables *tables,
    int x,
    int y,
    int z,
    int p
) {
    if (*n >= capacity) return MARKOFF_ERR_TOO_MANY_NODES;

    uint32_t i = *n;
    nodes[i].x = (uint16_t)x;
    nodes[i].y = (uint16_t)y;
    nodes[i].z = (uint16_t)z;
    nodes[i].pad = 0;
    nodes[i].root = i;
    nodes[i].neighbor0 = U32_MISSING;
    nodes[i].neighbor1 = U32_MISSING;
    nodes[i].neighbor2 = U32_MISSING;
    nodes[i].eigenvector = 0.0;

    if (!table_insert(tables->xy, p, x, y, i)) return MARKOFF_ERR_PAIR_OVERFLOW;
    if (!table_insert(tables->xz, p, x, z, i)) return MARKOFF_ERR_PAIR_OVERFLOW;
    if (!table_insert(tables->yz, p, y, z, i)) return MARKOFF_ERR_PAIR_OVERFLOW;

    *n = i + 1u;
    return MARKOFF_OK;
}

static int enumerate_solutions(
    MarkoffNode *nodes,
    uint32_t *n,
    uint32_t capacity,
    PairTables *tables,
    int A,
    int B,
    int C,
    int D,
    int p
) {
    if (p == 2) {
        for (int x = 0; x < p; ++x) {
            for (int y = 0; y < p; ++y) {
                for (int z = 0; z < p; ++z) {
                    int lhs = modp_ll((long long)x*x + (long long)y*y + (long long)z*z, p);
                    int rhs = modp_ll(1LL*x*y*z + (long long)A*x + (long long)B*y +
                                      (long long)C*z + D, p);
                    if (lhs == rhs) {
                        int status = add_solution(nodes, n, capacity, tables, x, y, z, p);
                        if (status != MARKOFF_OK) return status;
                    }
                }
            }
        }
        return MARKOFF_OK;
    }

    int *sqrt_count = (int *)calloc((size_t)p, sizeof(int));
    int *sqrt_roots = (int *)malloc((size_t)2 * (size_t)p * sizeof(int));
    if (sqrt_count == NULL || sqrt_roots == NULL) {
        free(sqrt_count);
        free(sqrt_roots);
        return MARKOFF_ERR_ALLOC;
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
            int S = modp_ll(1LL*x*y + C, p);
            int T = modp_ll((long long)x*x + (long long)y*y
                          - (long long)A*x - (long long)B*y - D, p);
            int disc = modp_ll((long long)S*S - 4LL*T, p);
            int cnt = sqrt_count[disc];

            for (int k = 0; k < cnt; ++k) {
                int r = sqrt_roots[2*disc + k];
                int z = modp_ll((long long)(S + r) * inv2, p);
                int status = add_solution(nodes, n, capacity, tables, x, y, z, p);
                if (status != MARKOFF_OK) {
                    free(sqrt_count);
                    free(sqrt_roots);
                    return status;
                }
            }
        }
    }

    free(sqrt_count);
    free(sqrt_roots);
    return MARKOFF_OK;
}

static int compute_components(MarkoffNode *nodes, uint32_t n, const PairTables *tables,
                              int A, int B, int C, int p, uint32_t *component_count) {
    for (uint32_t i = 0; i < n; ++i) {
        int x = (int)nodes[i].x;
        int y = (int)nodes[i].y;
        int z = (int)nodes[i].z;

        int sx = modp_ll(1LL*y*z + A - x, p);
        int sy = modp_ll(1LL*x*z + B - y, p);
        int sz = modp_ll(1LL*x*y + C - z, p);

        uint32_t ix = find_by_x(tables->yz, nodes, p, y, z, sx);
        uint32_t iy = find_by_y(tables->xz, nodes, p, x, z, sy);
        uint32_t iz = find_by_z(tables->xy, nodes, p, x, y, sz);
        if (ix == U32_MISSING || iy == U32_MISSING || iz == U32_MISSING) {
            return MARKOFF_ERR_NEIGHBOR_NOT_FOUND;
        }

        nodes[i].neighbor0 = ix;
        nodes[i].neighbor1 = iy;
        nodes[i].neighbor2 = iz;

        uf_union_min(nodes, i, ix);
        uf_union_min(nodes, i, iy);
        uf_union_min(nodes, i, iz);
    }

    uint32_t components = 0;
    for (uint32_t i = 0; i < n; ++i) {
        nodes[i].root = uf_find(nodes, i);
        if (nodes[i].root == i) components += 1u;
    }

    *component_count = components;
    return MARKOFF_OK;
}



static int build_component_data(MarkoffGraph *graph) {
    if (graph == NULL) return MARKOFF_ERR_BAD_ARGUMENT;

    uint32_t n = graph->node_count;
    uint32_t component_count = graph->component_count;

    if (component_count == 0) {
        graph->components = NULL;
        return MARKOFF_OK;
    }

    MarkoffComponent *components = (MarkoffComponent *)calloc((size_t)component_count, sizeof(MarkoffComponent));
    if (components == NULL) return MARKOFF_ERR_ALLOC;

    uint32_t c = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (graph->nodes[i].root == i) {
            if (c >= component_count) {
                free(components);
                return MARKOFF_ERR_EIGEN_FAILED;
            }
            components[c].root = &graph->nodes[i];
            components[c].root_index = i;
            components[c].size = 0;
            components[c].eigenvalue = 0.0;
            c += 1u;
        }
    }
    if (c != component_count) {
        free(components);
        return MARKOFF_ERR_EIGEN_FAILED;
    }

    for (uint32_t ci = 0; ci < component_count; ++ci) {
        uint32_t root = components[ci].root_index;
        uint32_t size = 0;
        for (uint32_t i = 0; i < n; ++i) {
            if (graph->nodes[i].root == root) size += 1u;
        }
        components[ci].size = size;
    }

    graph->components = components;
    return MARKOFF_OK;
}

static uint32_t component_size_for_root(const MarkoffGraph *graph, uint32_t root_index) {
    if (graph == NULL || graph->components == NULL) return 0;
    for (uint32_t i = 0; i < graph->component_count; ++i) {
        if (graph->components[i].root_index == root_index) return graph->components[i].size;
    }
    return 0;
}

static void markoff_csr_free_internal(MarkoffCSR *csr) {
    if (csr == NULL) return;
    free(csr->data);
    free(csr->indices);
    free(csr->indptr);
    free(csr->nodes);
    free(csr);
}

MARKOFF_EXPORT int markoff_component_csr(
    const MarkoffGraph *graph,
    uint32_t root_index,
    MarkoffCSR **out_csr
) {
    if (out_csr == NULL) return MARKOFF_ERR_BAD_ARGUMENT;
    *out_csr = NULL;
    if (graph == NULL) return MARKOFF_ERR_BAD_ARGUMENT;
    if (root_index >= graph->node_count) return MARKOFF_ERR_BAD_ARGUMENT;
    if (graph->nodes[root_index].root != root_index) return MARKOFF_ERR_BAD_ARGUMENT;

    const uint32_t n = graph->node_count;
    const uint32_t m = component_size_for_root(graph, root_index);
    if (m == 0) return MARKOFF_ERR_BAD_ARGUMENT;
    if (m > (uint32_t)INT_MAX) return MARKOFF_ERR_BAD_ARGUMENT;
    if (3ull * (uint64_t)m > (uint64_t)INT_MAX) return MARKOFF_ERR_BAD_ARGUMENT;

    MarkoffCSR *csr = (MarkoffCSR *)calloc(1u, sizeof(MarkoffCSR));
    uint32_t *local_pos = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    if (csr == NULL || local_pos == NULL) {
        free(csr);
        free(local_pos);
        return MARKOFF_ERR_ALLOC;
    }

    csr->size = m;
    csr->nnz = 3u * m;
    csr->root_index = root_index;
    csr->data = (double *)malloc((size_t)csr->nnz * sizeof(double));
    csr->indices = (int *)malloc((size_t)csr->nnz * sizeof(int));
    csr->indptr = (int *)malloc((size_t)(m + 1u) * sizeof(int));
    csr->nodes = (uint32_t *)malloc((size_t)m * sizeof(uint32_t));

    if (csr->data == NULL || csr->indices == NULL || csr->indptr == NULL || csr->nodes == NULL) {
        free(local_pos);
        markoff_csr_free_internal(csr);
        return MARKOFF_ERR_ALLOC;
    }

    for (uint32_t i = 0; i < n; ++i) local_pos[i] = U32_MISSING;

    uint32_t j = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (graph->nodes[i].root == root_index) {
            if (j >= m) {
                free(local_pos);
                markoff_csr_free_internal(csr);
                return MARKOFF_ERR_EIGEN_FAILED;
            }
            local_pos[i] = j;
            csr->nodes[j] = i;
            j += 1u;
        }
    }
    if (j != m) {
        free(local_pos);
        markoff_csr_free_internal(csr);
        return MARKOFF_ERR_EIGEN_FAILED;
    }

    for (uint32_t row = 0; row < m; ++row) {
        uint32_t gi = csr->nodes[row];
        uint32_t nbs[3] = {
            graph->nodes[gi].neighbor0,
            graph->nodes[gi].neighbor1,
            graph->nodes[gi].neighbor2,
        };
        csr->indptr[row] = (int)(3u * row);
        for (uint32_t k = 0; k < 3u; ++k) {
            uint32_t nb = nbs[k];
            if (nb == U32_MISSING || nb >= n || local_pos[nb] == U32_MISSING) {
                free(local_pos);
                markoff_csr_free_internal(csr);
                return MARKOFF_ERR_NEIGHBOR_NOT_FOUND;
            }
            uint32_t pos = 3u * row + k;
            csr->data[pos] = 1.0;
            csr->indices[pos] = (int)local_pos[nb];
        }
    }
    csr->indptr[m] = (int)csr->nnz;

    free(local_pos);
    *out_csr = csr;
    return MARKOFF_OK;
}

MARKOFF_EXPORT void markoff_csr_free(MarkoffCSR *csr) {
    markoff_csr_free_internal(csr);
}

MARKOFF_EXPORT int markoff_build(
    int A,
    int B,
    int C,
    int D,
    int prime,
    MarkoffGraph **out_graph
) {
    if (out_graph == NULL) return MARKOFF_ERR_BAD_ARGUMENT;
    *out_graph = NULL;

    if (prime <= 1) return MARKOFF_ERR_BAD_ARGUMENT;
    if (!prime_fits_storage_bound(prime)) return MARKOFF_ERR_STORAGE_BOUND;
    if (!is_prime_small(prime)) return MARKOFF_ERR_COMPOSITE_MODULUS;

    const int p = prime;
    const int Am = modp_ll(A, p);
    const int Bm = modp_ll(B, p);
    const int Cm = modp_ll(C, p);
    const int Dm = modp_ll(D, p);
    const uint32_t capacity = default_capacity(p);

    if ((size_t)capacity > SIZE_MAX / sizeof(MarkoffNode)) return MARKOFF_ERR_BAD_ARGUMENT;

    MarkoffNode *nodes = (MarkoffNode *)malloc((size_t)capacity * sizeof(MarkoffNode));
    if (nodes == NULL) return MARKOFF_ERR_ALLOC;

    PairTables tables;
    int status = pair_tables_init(&tables, p);
    if (status != MARKOFF_OK) {
        free(nodes);
        return status;
    }

    uint32_t n = 0;
    status = enumerate_solutions(nodes, &n, capacity, &tables, Am, Bm, Cm, Dm, p);
    if (status == MARKOFF_OK) {
        uint32_t component_count = 0;
        status = compute_components(nodes, n, &tables, Am, Bm, Cm, p, &component_count);
        if (status == MARKOFF_OK) {
            if (n == 0) {
                free(nodes);
                nodes = NULL;
            } else {
                MarkoffNode *shrunk = (MarkoffNode *)realloc(nodes, (size_t)n * sizeof(MarkoffNode));
                if (shrunk != NULL) nodes = shrunk;
            }

            MarkoffGraph *graph = (MarkoffGraph *)malloc(sizeof(MarkoffGraph));
            if (graph == NULL) {
                free(nodes);
                return MARKOFF_ERR_ALLOC;
            }

            graph->A = (int32_t)A;
            graph->B = (int32_t)B;
            graph->C = (int32_t)C;
            graph->D = (int32_t)D;
            graph->prime = (uint32_t)p;
            graph->node_count = n;
            graph->component_count = component_count;
            graph->nodes = nodes;
            graph->components = NULL;

            status = build_component_data(graph);
            if (status != MARKOFF_OK) {
                free(graph->components);
                free(graph->nodes);
                free(graph);
                pair_tables_free(&tables);
                return status;
            }

            pair_tables_free(&tables);
            *out_graph = graph;
            return MARKOFF_OK;
        }
    }

    pair_tables_free(&tables);
    free(nodes);
    return status;
}

MARKOFF_EXPORT void markoff_free(MarkoffGraph *graph) {
    if (graph == NULL) return;
    free(graph->components);
    graph->components = NULL;
    free(graph->nodes);
    graph->nodes = NULL;
    free(graph);
}

MARKOFF_EXPORT uint32_t markoff_node_count(const MarkoffGraph *graph) {
    if (graph == NULL) return 0;
    return graph->node_count;
}

MARKOFF_EXPORT uint32_t markoff_component_count(const MarkoffGraph *graph) {
    if (graph == NULL) return 0;
    return graph->component_count;
}

MARKOFF_EXPORT const MarkoffNode *markoff_nodes(const MarkoffGraph *graph) {
    if (graph == NULL) return NULL;
    return graph->nodes;
}

MARKOFF_EXPORT const MarkoffComponent *markoff_components(const MarkoffGraph *graph) {
    if (graph == NULL) return NULL;
    return graph->components;
}

#ifdef __cplusplus
}
#endif
