// chiProposal46_4ary_fast.c — 4-ary, 3x4 block, R = 11/12 (FAST)
// Giữ: Dither D (seam-safe), EA-2V (U1,U2), intra cost, seam trái/phải + seam với hàng trên.
// Bỏ: snapshot/logging, DP theo cột, forward/backward nhiều lượt, refine 2×2.

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stddef.h>

/* =========================== Geometry =========================== */
#define BH 3
#define BW 4
#define N_PER_BLOCK (BH*BW)      /* 12 */
#define SYMS_PER_BLOCK 11        /* 11 data/block (R = 11/12) */

/* ===================== Vị trí data & index ====================== */
#define INDEX_POS (1*BW + 1)     /* (r=1,c=1) */
static int g_index_pos = INDEX_POS;
static int g_data_pos[SYMS_PER_BLOCK];

/* ===================== Trọng số cost ===================== */
static const int W_EXTREME   = 4;
static const int W_DIFF      = 2;
static const int W_RUN_ROW   = 2;
static const int W_RUN_COL   = 2;
static const int W_SEAM      = 10;
static const int W_SEAM_XTRM = 6;
static const int W_DC        = 4;

/* ================== Lân cận nội khối (17 cạnh) ================== */
struct Edge { int u, v; };
static struct Edge EDGES[32];
static int E_CNT = 0;

/* =================== Mặt nạ (seam-safe) & Dither =================== */
/* Dither D: 0 ở biên trái/trên; checker nội thất; index=0 */
static int g_mask_D[N_PER_BLOCK];
static int g_mask_U1[N_PER_BLOCK];
static int g_mask_U2[N_PER_BLOCK];

/* =================== Init vị trí, cạnh, mặt nạ =================== */
static void init_positions_once(void) {
    static int inited = 0;
    if (inited) return;
    int t = 0;
    for (int p = 0; p < N_PER_BLOCK; ++p) {
        if (p == g_index_pos) continue;
        g_data_pos[t++] = p;
    }
    inited = 1;
}

static void init_masks_once(void) {
    static int inited = 0;
    if (inited) return;
    for (int r = 0; r < BH; ++r) {
        for (int c = 0; c < BW; ++c) {
            int p = r * BW + c;
            int D = 0, U1 = 0, U2 = 0;
            /* D: nội thất checker, 0 ở mép trên/trái & index */
            if (!(r == 0 || c == 0) && p != g_index_pos) D = (((r + c) & 1) == 0) ? 1 : 3;
            /* U1: theo cột (nội thất), xoá ở hàng trên & index */
            if (c == 1) U1 = 1; else if (c == BW - 2) U1 = 3;
            if (r == 0 || p == g_index_pos) U1 = 0;
            /* U2: theo hàng (nội thất), xoá ở cột trái & index */
            if (r == 1) U2 = 1; else if (r == BH - 2) U2 = 3;
            if (c == 0 || p == g_index_pos) U2 = 0;

            g_mask_D[p]  = D  & 3;
            g_mask_U1[p] = U1 & 3;
            g_mask_U2[p] = U2 & 3;
        }
    }
    inited = 1;
}

static void init_edges_once(void) {
    if (E_CNT) return;
    int e = 0;
    for (int r = 0; r < BH; ++r)
        for (int c = 0; c < BW - 1; ++c) {
            int u = r * BW + c;
            EDGES[e].u = u; EDGES[e].v = u + 1; ++e;
        }
    for (int r = 0; r < BH - 1; ++r)
        for (int c = 0; c < BW; ++c) {
            int u = r * BW + c;
            EDGES[e].u = u; EDGES[e].v = u + BW; ++e;
        }
    E_CNT = e; /* 17 */
    init_positions_once();
    init_masks_once();
}

/* ===================== Helpers & cost ===================== */
static inline int addmod4(int a, int b) { return (a + b) & 3; }
static inline int submod4(int a, int b) { return (a - b) & 3; }
static inline int is_extreme(int x)     { return (x == 0 || x == 3); }
static inline int is_bad_pair(int a, int b) { return ((a == 0 && b == 3) || (a == 3 && b == 0)); }

static int cost_intrablock(const int g[N_PER_BLOCK]) {
    int cost = 0;
    if (W_EXTREME) {
        int cnt = 0; for (int i = 0; i < N_PER_BLOCK; ++i) cnt += is_extreme(g[i]);
        cost += W_EXTREME * cnt;
    }
    if (W_DIFF) {
        int s = 0; for (int i = 0; i < E_CNT; ++i) {
            int d = g[EDGES[i].u] - g[EDGES[i].v]; if (d < 0) d = -d; s += d;
        }
        cost += W_DIFF * s;
    }
    if (W_RUN_ROW) {
        for (int r = 0; r < BH; ++r) {
            int run = 0;
            for (int c = 0; c < BW; ++c) {
                int x = g[r * BW + c];
                if (is_extreme(x)) { ++run; if (run > 1) cost += W_RUN_ROW * (run - 1); } else run = 0;
            }
        }
    }
    if (W_RUN_COL) {
        for (int c = 0; c < BW; ++c) {
            int run = 0;
            for (int r = 0; r < BH; ++r) {
                int x = g[r * BW + c];
                if (is_extreme(x)) { ++run; if (run > 1) cost += W_RUN_COL * (run - 1); } else run = 0;
            }
        }
    }
    if (W_DC) {
        int sum = 0; for (int i = 0; i < N_PER_BLOCK; ++i) sum += g[i];
        int dev = sum - 18; if (dev < 0) dev = -dev;
        cost += W_DC * dev;
    }
    return cost;
}

static int cost_seam_cols(const int cl[BH], const int cr[BH]) {
    int cost = 0;
    for (int r = 0; r < BH; ++r) {
        int a = cl[r], b = cr[r];
        int d = a - b; if (d < 0) d = -d;
        cost += W_SEAM * d;
        if (W_SEAM_XTRM && is_bad_pair(a, b)) cost += W_SEAM_XTRM;
        if (W_RUN_ROW && is_extreme(a) && is_extreme(b) && (a == b)) cost += W_RUN_ROW; /* run ngang băng seam */
    }
    return cost;
}

static int cost_seam_rows(const int rt[BW], const int rb[BW]) {
    int cost = 0;
    for (int c = 0; c < BW; ++c) {
        int a = rt[c], b = rb[c];
        int d = a - b; if (d < 0) d = -d;
        cost += W_SEAM * d;
        if (W_SEAM_XTRM && is_bad_pair(a, b)) cost += W_SEAM_XTRM;
        if (W_RUN_COL && is_extreme(a) && is_extreme(b) && (a == b)) cost += W_RUN_COL; /* run dọc băng seam */
    }
    return cost;
}

/* ===================== EA-2V block masks (2-bit) ===================== */
static inline int apply_U_fwd(int tf, int pos, int x) {
    int m = 0;
    switch (tf & 3) {
        case 0: m = 0; break;
        case 1: m = g_mask_U1[pos]; break;
        case 2: m = g_mask_U2[pos]; break;
        default: m = (g_mask_U1[pos] + g_mask_U2[pos]) & 3; break;
    }
    return addmod4(x & 3, m);
}
static inline int apply_U_inv(int tf, int pos, int y) {
    int m = 0;
    switch (tf & 3) {
        case 0: m = 0; break;
        case 1: m = g_mask_U1[pos]; break;
        case 2: m = g_mask_U2[pos]; break;
        default: m = (g_mask_U1[pos] + g_mask_U2[pos]) & 3; break;
    }
    return submod4(y & 3, m);
}

/* ========== Candidate builder (tạo grid + biên) ========== */
typedef struct {
    int grid[N_PER_BLOCK];
    int left[BH], right[BH];
    int top[BW], bottom[BW];
    int cost_intra;
} Cand;

static void build_base_grid_from_syms(const int* in_syms, int base[N_PER_BLOCK]) {
    for (int p = 0; p < N_PER_BLOCK; ++p) base[p] = 0;
    for (int k = 0; k < SYMS_PER_BLOCK; ++k) {
        int pos = g_data_pos[k];
        base[pos] = in_syms[k] & 3;
    }
    base[g_index_pos] = 0;
    for (int p = 0; p < N_PER_BLOCK; ++p) {
        if (p == g_index_pos) continue;
        base[p] = addmod4(base[p], g_mask_D[p]);
    }
}

static void make_cand_from_base(const int base[N_PER_BLOCK], int tf, Cand* out) {
    for (int p = 0; p < N_PER_BLOCK; ++p) {
        int v = base[p];
        if (p != g_index_pos) v = apply_U_fwd(tf, p, v);
        out->grid[p] = v & 3;
    }
    out->grid[g_index_pos] = tf & 3;
    for (int r = 0; r < BH; ++r) {
        out->left[r]  = out->grid[r * BW + 0] & 3;
        out->right[r] = out->grid[r * BW + (BW - 1)] & 3;
    }
    for (int c = 0; c < BW; ++c) {
        out->top[c]    = out->grid[0 * BW + c] & 3;
        out->bottom[c] = out->grid[(BH - 1) * BW + c] & 3;
    }
    out->cost_intra = cost_intrablock(out->grid);
}

/* ===================== I/O tiện ích biên ===================== */
static void write_block_to_page_from_cand(int** PAGE, int r0, int c0, const Cand* cand) {
    for (int r = 0; r < BH; ++r)
        for (int c = 0; c < BW; ++c)
            PAGE[r0 + r][c0 + c] = cand->grid[r * BW + c] & 3;
}
static void read_page_top_row(int** PAGE, int r0, int c0, int out[BW]) {
    for (int c = 0; c < BW; ++c) out[c] = PAGE[r0 - 1][c0 + c] & 3;
}

/* ========== Lấp mép còn hở (nếu Page_Size không chia hết) ========== */
static void fill_uncovered_edges_after_encode(int** PAGE, int Page_Size) {
    const int full_cols = (Page_Size / BW) * BW;
    const int full_rows = (Page_Size / BH) * BH;
    if (full_cols < Page_Size && full_cols > 0) {
        int csrc = full_cols - 1;
        for (int r = 0; r < full_rows; ++r) {
            int v = PAGE[r][csrc] & 3;
            for (int c = full_cols; c < Page_Size; ++c) PAGE[r][c] = v;
        }
    }
    if (full_rows < Page_Size && full_rows > 0) {
        int rsrc = full_rows - 1;
        for (int r = full_rows; r < Page_Size; ++r)
            for (int c = 0; c < Page_Size; ++c)
                PAGE[r][c] = PAGE[rsrc][c] & 3;
    }
}

/* ======================  PUBLIC API: FAST ENCODER  ====================== */
void Encode_khanhProposal1112_4ary(int** PAGE, int* input_1Ddata, int Page_Size)
{
    init_edges_once(); /* cũng khởi tạo positions & masks */
    const int BR = Page_Size / BH;
    const int BC = Page_Size / BW;
    if (BR <= 0 || BC <= 0) return;

    /* Workspace cho 1 hàng */
    Cand* C = (Cand*)malloc(sizeof(Cand) * (size_t)BC * 4);
    int* U  = (int*)malloc(sizeof(int) * (size_t)BC * 4);   /* unary = intra + top seam */
    int* back = (int*)malloc(sizeof(int) * (size_t)BC * 4);
    int  dp_prev[4], dp_curr[4];

    /* Quét theo hàng: mỗi hàng tối ưu Viterbi trái->phải 1 lần */
    for (int i = 0; i < BR; ++i) {
        const int r0 = i * BH;

        /* 1) Build 4 ứng viên cho từng block trong hàng + tính unary (intra + top seam) */
        for (int j = 0; j < BC; ++j) {
            int in[SYMS_PER_BLOCK];
            for (int k = 0; k < SYMS_PER_BLOCK; ++k)
                in[k] = input_1Ddata[(i * BC + j) * SYMS_PER_BLOCK + k] & 3;

            int base[N_PER_BLOCK];
            build_base_grid_from_syms(in, base);

            for (int t = 0; t < 4; ++t) {
                Cand* cj = &C[j * 4 + t];
                make_cand_from_base(base, t, cj);
                int u = cj->cost_intra;
                if (i > 0) { /* seam với hàng trên -> unary */
                    int neigh_top[BW];
                    read_page_top_row(PAGE, r0, j * BW, neigh_top);
                    u += cost_seam_rows(cj->top, neigh_top);
                }
                U[j * 4 + t] = u;
            }
        }

        /* 2) Viterbi trên hàng: pairwise = seam trái/phải giữa (j-1) và (j) */
        for (int s = 0; s < 4; ++s) { dp_prev[s] = U[0 * 4 + s]; back[0 * 4 + s] = -1; }

        for (int j = 1; j < BC; ++j) {
            for (int s = 0; s < 4; ++s) {
                int best = INT_MAX, bestt = 0;
                for (int t = 0; t < 4; ++t) {
                    int pair = cost_seam_cols(C[j * 4 + s].left, C[(j - 1) * 4 + t].right);
                    int val = dp_prev[t] + pair;
                    if (val < best) { best = val; bestt = t; }
                }
                dp_curr[s] = best + U[j * 4 + s];
                back[j * 4 + s] = bestt;
            }
            for (int s = 0; s < 4; ++s) dp_prev[s] = dp_curr[s];
        }

        /* 3) Truy vết & ghi PAGE cho hàng i */
        int best_s = 0, bestE = dp_prev[0];
        for (int s = 1; s < 4; ++s) if (dp_prev[s] < bestE) { bestE = dp_prev[s]; best_s = s; }

        int* states = (int*)malloc(sizeof(int) * (size_t)BC);
        int s = best_s;
        for (int j = BC - 1; j >= 0; --j) {
            states[j] = s;
            s = back[j * 4 + s];
            if (s < 0 && j > 0) s = 0;
        }
        for (int j = 0; j < BC; ++j)
            write_block_to_page_from_cand(PAGE, r0, j * BW, &C[j * 4 + states[j]]);
        free(states);
    }

    fill_uncovered_edges_after_encode(PAGE, Page_Size);

    free(back); free(U); free(C);
}

/* ======================  PUBLIC API: DECODER  ====================== */
void Decode_khanhProposal1112_4ary(int* output_1Ddata, double** PAGE, int Page_Size)
{
    init_edges_once();
    const int BLOCK_ROWS = Page_Size / BH;
    const int BLOCK_COLS = Page_Size / BW;

    int grid12[N_PER_BLOCK];
    int out_block[SYMS_PER_BLOCK];

    for (int i = 0; i < BLOCK_ROWS; ++i) {
        for (int j = 0; j < BLOCK_COLS; ++j) {
            int r0 = i * BH, c0 = j * BW;

            for (int r = 0; r < BH; ++r)
                for (int c = 0; c < BW; ++c)
                    grid12[r * BW + c] = ((int)PAGE[r0 + r][c0 + c]) & 3;

            int tf = grid12[g_index_pos] & 3;

            for (int p = 0; p < N_PER_BLOCK; ++p) {
                if (p == g_index_pos) continue;
                grid12[p] = apply_U_inv(tf, p, grid12[p]) & 3;
            }
            for (int p = 0; p < N_PER_BLOCK; ++p) {
                if (p == g_index_pos) continue;
                grid12[p] = submod4(grid12[p], g_mask_D[p]) & 3;
            }
            for (int k = 0; k < SYMS_PER_BLOCK; ++k) {
                int pos = g_data_pos[k];
                out_block[k] = grid12[pos] & 3;
            }
            for (int k = 0; k < SYMS_PER_BLOCK; ++k)
                output_1Ddata[(i * BLOCK_COLS + j) * SYMS_PER_BLOCK + k] = out_block[k] & 3;
        }
    }
}
