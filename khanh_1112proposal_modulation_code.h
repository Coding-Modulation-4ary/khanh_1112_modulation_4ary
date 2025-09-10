// chiProposal46_4ary.c — Balanced 4-ary, 3x4 block, R = 11/12
// 11 data symbols + 1 index (transform id in {0..3})
// Seam-safe deterministic Dither (no bits) + EA-2V (2-bit) trên NỘI THẤT khối
// Seam-aware + DC penalty + run penalty theo HÀNG & CỘT.
// Base: (i*BLOCK_COLS + j)*SYMS_PER_BLOCK + k
//
// Public API (drop-in):
//   void Encode_khanhProposal1112_4ary(int** PAGE, int* input_1Ddata, int Page_Size);
//   void Decode_khanhProposal1112_4ary(int* output_1Ddata, double** PAGE, int Page_Size);

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

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
static const int W_EXTREME = 4;
static const int W_DIFF = 2;
static const int W_RUN_ROW = 2;
static const int W_RUN_COL = 2;
static const int W_SEAM = 10;
static const int W_SEAM_XTRM = 6;
static const int W_DC = 4;

/* ================== Lân cận nội khối (17 cạnh) ================== */
struct Edge { int u, v; };
static struct Edge EDGES[32];
static int E_CNT = 0;

/* ========== GLOBAL DEBUG SNAPSHOT (encoded PAGE phẳng) ========== */
static int* g_dbg_enc_flat = NULL;
static size_t g_dbg_enc_len = 0;
static int    g_dbg_enc_pagesz = 0;

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
            if (!(r == 0 || c == 0) && p != g_index_pos) D = (((r + c) & 1) == 0) ? 1 : 3;
            if (c == 1) U1 = 1; else if (c == BW - 2) U1 = 3;
            if (r == 0 || p == g_index_pos) U1 = 0;
            if (r == 1) U2 = 1; else if (r == BH - 2) U2 = 3;
            if (c == 0 || p == g_index_pos) U2 = 0;
            g_mask_D[p] = D & 3;
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
static inline int is_extreme(int x) { return (x == 0 || x == 3); }
static inline int is_bad_pair(int a, int b) { return ((a == 0 && b == 3) || (a == 3 && b == 0)); }

static int cost_intrablock(const int g[N_PER_BLOCK]) {
    int cost = 0;
    if (W_EXTREME) { int cnt = 0; for (int i = 0; i < N_PER_BLOCK; ++i) cnt += is_extreme(g[i]); cost += W_EXTREME * cnt; }
    if (W_DIFF) { int s = 0; for (int i = 0; i < E_CNT; ++i) { int d = g[EDGES[i].u] - g[EDGES[i].v]; if (d < 0) d = -d; s += d; } cost += W_DIFF * s; }
    if (W_RUN_ROW) { for (int r = 0; r < BH; ++r) { int run = 0; for (int c = 0; c < BW; ++c) { int x = g[r * BW + c]; if (is_extreme(x)) { ++run; if (run > 1) cost += W_RUN_ROW * (run - 1); } else run = 0; } } }
    if (W_RUN_COL) { for (int c = 0; c < BW; ++c) { int run = 0; for (int r = 0; r < BH; ++r) { int x = g[r * BW + c]; if (is_extreme(x)) { ++run; if (run > 1) cost += W_RUN_COL * (run - 1); } else run = 0; } } }
    if (W_DC) { int sum = 0; for (int i = 0; i < N_PER_BLOCK; ++i) sum += g[i]; int dev = sum - 18; if (dev < 0) dev = -dev; cost += W_DC * dev; }
    return cost;
}

/* seam chạy trên một cột (BH phần tử) */
static int cost_seam_cols(const int cl[BH], const int cr[BH]) {
    int cost = 0;
    for (int r = 0; r < BH; ++r) {
        int a = cl[r], b = cr[r];
        int d = a - b; if (d < 0) d = -d;
        cost += W_SEAM * d;
        if (W_SEAM_XTRM && is_bad_pair(a, b)) cost += W_SEAM_XTRM;
        if (W_RUN_ROW && is_extreme(a) && is_extreme(b) && (a == b)) cost += W_RUN_ROW; /* run qua seam ngang */
    }
    return cost;
}
/* seam chạy trên một hàng (BW phần tử) */
static int cost_seam_rows(const int rt[BW], const int rb[BW]) {
    int cost = 0;
    for (int c = 0; c < BW; ++c) {
        int a = rt[c], b = rb[c];
        int d = a - b; if (d < 0) d = -d;
        cost += W_SEAM * d;
        if (W_SEAM_XTRM && is_bad_pair(a, b)) cost += W_SEAM_XTRM;
        if (W_RUN_COL && is_extreme(a) && is_extreme(b) && (a == b)) cost += W_RUN_COL; /* run qua seam dọc */
    }
    return cost;
}

/* ===================== EA-2V block masks (2-bit) ===================== */
static inline int apply_U_fwd(int tf, int pos, int x) {
    int m = 0;
    switch (tf & 3) { case 0: m = 0; break; case 1: m = g_mask_U1[pos]; break; case 2: m = g_mask_U2[pos]; break; default: m = (g_mask_U1[pos] + g_mask_U2[pos]) & 3; break; }
                            return addmod4(x & 3, m);
}
static inline int apply_U_inv(int tf, int pos, int y) {
    int m = 0;
    switch (tf & 3) { case 0: m = 0; break; case 1: m = g_mask_U1[pos]; break; case 2: m = g_mask_U2[pos]; break; default: m = (g_mask_U1[pos] + g_mask_U2[pos]) & 3; break; }
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
        out->left[r] = out->grid[r * BW + 0] & 3;
        out->right[r] = out->grid[r * BW + (BW - 1)] & 3;
    }
    for (int c = 0; c < BW; ++c) {
        out->top[c] = out->grid[0 * BW + c] & 3;
        out->bottom[c] = out->grid[(BH - 1) * BW + c] & 3;
    }
    out->cost_intra = cost_intrablock(out->grid);
}

/* ========== Encode/Decode 1 block (tham lam cơ sở) ========== */
static void encode_block_greedy(const int* in_syms,
    int out_grid[N_PER_BLOCK],
    int* chosen_tf,
    int** PAGE, int r0, int c0, int i, int j)
{
    int base[N_PER_BLOCK]; build_base_grid_from_syms(in_syms, base);
    int best_tf = 0, best_cost = INT_MAX, best_grid[N_PER_BLOCK];
    for (int tf = 0; tf < 4; ++tf) {
        Cand cand; make_cand_from_base(base, tf, &cand);
        int c = cand.cost_intra;
        /* seam trái & trên dùng PAGE đã dựng tới thời điểm này */
        if (j > 0) {
            int neigh[BH]; for (int r = 0; r < BH; ++r) neigh[r] = PAGE[r0 + r][c0 - 1] & 3;
            c += cost_seam_cols(cand.left, neigh);
        }
        if (i > 0) {
            int neigh[BW]; for (int cc = 0; cc < BW; ++cc) neigh[cc] = PAGE[r0 - 1][c0 + cc] & 3;
            c += cost_seam_rows(cand.top, neigh);
        }
        if (c < best_cost) { best_cost = c; best_tf = tf; for (int p = 0;p < N_PER_BLOCK;++p) best_grid[p] = cand.grid[p]; }
    }
    *chosen_tf = best_tf;
    for (int p = 0; p < N_PER_BLOCK; ++p) out_grid[p] = best_grid[p] & 3;
}

static void decode_block_4x3_11of12(const int in_grid[N_PER_BLOCK], int* out_syms) {
    int tf = in_grid[g_index_pos] & 3;
    for (int k = 0; k < SYMS_PER_BLOCK; ++k) {
        int pos = g_data_pos[k];
        int y = in_grid[pos] & 3;
        out_syms[k] = apply_U_inv(tf, pos, y) & 3;
    }
}

/* ========== Snapshot & lấp mép ========== */
static void dbg_store_encoded_page(int** PAGE, int Page_Size) {
    const size_t need = (size_t)Page_Size * (size_t)Page_Size;
    if (g_dbg_enc_len != need) {
        if (g_dbg_enc_flat) { free(g_dbg_enc_flat); g_dbg_enc_flat = NULL; }
        g_dbg_enc_flat = (int*)malloc(sizeof(int) * need);
        g_dbg_enc_len = need;
        g_dbg_enc_pagesz = Page_Size;
    }
    if (!g_dbg_enc_flat) return;
    size_t p = 0;
    for (int r = 0; r < Page_Size; ++r)
        for (int c = 0; c < Page_Size; ++c)
            g_dbg_enc_flat[p++] = PAGE[r][c] & 3;
}

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

/* ====== Helpers: TF id theo block cho enc/dec (AN TOÀN) ====== */
static int get_tf_from_enc_block(int bi, int bj) {
    if (!g_dbg_enc_flat || g_dbg_enc_pagesz <= 0) return 0;
    int br = g_dbg_enc_pagesz / BH;
    int bc = g_dbg_enc_pagesz / BW;
    if (bi < 0) bi = 0; if (bi >= br) bi = br - 1;
    if (bj < 0) bj = 0; if (bj >= bc) bj = bc - 1;
    int r0 = bi * BH, c0 = bj * BW;
    int ir = g_index_pos / BW, ic = g_index_pos % BW;
    size_t p = (size_t)(r0 + ir) * g_dbg_enc_pagesz + (c0 + ic);
    return g_dbg_enc_flat[p] & 3;
}

static int get_tf_from_dec_block_safe(double** PAGE_DBL, int bi, int bj, int Page_Size) {
    int br = Page_Size / BH, bc = Page_Size / BW;
    if (bi < 0) bi = 0; if (bi >= br) bi = br - 1;
    if (bj < 0) bj = 0; if (bj >= bc) bj = bc - 1;
    int r0 = bi * BH, c0 = bj * BW;
    int ir = g_index_pos / BW, ic = g_index_pos % BW;
    return (int)PAGE_DBL[r0 + ir][c0 + ic];
}

/* ========== FULL PAGE COMPARISON (KHÔNG lượng tử) ========== */
/* (giữ nguyên code compare_full_and_log_with_windows của bạn – không lặp lại ở đây cho gọn) */
/* ... DÁN nguyên hàm compare_full_and_log_with_windows ở bản của bạn vào đây ... */

/* ===================== Pairwise refinement ===================== */
static void write_block_to_page_from_cand(int** PAGE, int r0, int c0, const Cand* cand) {
    for (int r = 0; r < BH; ++r)
        for (int c = 0; c < BW; ++c)
            PAGE[r0 + r][c0 + c] = cand->grid[r * BW + c] & 3;
}

/* Đọc biên từ PAGE để làm hàng xóm cố định */
static void read_page_left_col(int** PAGE, int r0, int c0, int out[BH]) { /* c0 là cột đầu của block hiện tại */
    for (int r = 0; r < BH; ++r) out[r] = PAGE[r0 + r][c0 - 1] & 3;
}
static void read_page_right_col(int** PAGE, int r0, int c0, int out[BH]) { /* c0 là cột đầu của block bên trái */
    for (int r = 0; r < BH; ++r) out[r] = PAGE[r0 + r][c0 + (BW - 1)] & 3;
}
static void read_page_top_row(int** PAGE, int r0, int c0, int out[BW]) { /* r0 là hàng đầu của block hiện tại */
    for (int c = 0; c < BW; ++c) out[c] = PAGE[r0 - 1][c0 + c] & 3;
}
static void read_page_bottom_row(int** PAGE, int r0, int c0, int out[BW]) { /* r0 là hàng đầu của block trên */
    for (int c = 0; c < BW; ++c) out[c] = PAGE[r0 + (BH - 1)][c0 + c] & 3;
}

/* Năng lượng đầy đủ cho một block cand với hàng xóm cố định (không tính seam với khối cặp-đang-tối-ưu) */
static int energy_block_with_fixed_neighbors(int i, int j, int tf,
    int** PAGE, int Page_Size,
    int* input_1Ddata)
{
    const int BLOCK_COLS = Page_Size / BW;
    int in_syms[SYMS_PER_BLOCK];
    for (int k = 0;k < SYMS_PER_BLOCK;++k) in_syms[k] = input_1Ddata[(i * BLOCK_COLS + j) * SYMS_PER_BLOCK + k] & 3;

    int base[N_PER_BLOCK]; build_base_grid_from_syms(in_syms, base);
    Cand c; make_cand_from_base(base, tf, &c);

    int r0 = i * BH, c0 = j * BW;
    int E = c.cost_intra;

    /* seam dọc với trên/dưới (nếu có) */
    if (i > 0) { int neigh[BW]; read_page_top_row(PAGE, r0, c0, neigh); E += cost_seam_rows(c.top, neigh); }
    if (i + 1 < (Page_Size / BH)) {
        int neigh[BW]; /* hàng dưới của block hiện tại so với top hàng dưới */
        for (int cc = 0; cc < BW; ++cc) neigh[cc] = PAGE[r0 + BH][c0 + cc] & 3;
        E += cost_seam_rows(neigh, c.bottom); /* bottom of current vs top of below */
    }

    /* seam ngang với trái/phải (nếu có) */
    if (j > 0) { int neigh[BH]; read_page_left_col(PAGE, r0, c0, neigh); E += cost_seam_cols(c.left, neigh); }
    if (j + 1 < (Page_Size / BW)) { int neigh[BH]; for (int r = 0;r < BH;++r) neigh[r] = PAGE[r0 + r][c0 + BW] & 3; E += cost_seam_cols(neigh, c.right); }

    return E;
}

/* Sweep ngang: tối ưu từng cặp (i, j-1)-(i, j) */
static void horizontal_pair_sweep(int** PAGE, int* input_1Ddata, int Page_Size, int iters)
{
    const int BR = Page_Size / BH;
    const int BC = Page_Size / BW;

    for (int iter = 0; iter < iters; ++iter) {
        for (int i = 0; i < BR; ++i) {
            int r0 = i * BH;
            for (int j = 1; j < BC; ++j) {
                int c0L = (j - 1) * BW, c0R = j * BW;

                /* build candidate bộ 4 cho trái & phải */
                int inL[SYMS_PER_BLOCK], inR[SYMS_PER_BLOCK];
                for (int k = 0;k < SYMS_PER_BLOCK;++k) {
                    inL[k] = input_1Ddata[(i * BC + (j - 1)) * SYMS_PER_BLOCK + k] & 3;
                    inR[k] = input_1Ddata[(i * BC + j) * SYMS_PER_BLOCK + k] & 3;
                }
                int baseL[N_PER_BLOCK], baseR[N_PER_BLOCK];
                build_base_grid_from_syms(inL, baseL);
                build_base_grid_from_syms(inR, baseR);
                Cand CL[4], CR[4];
                for (int tf = 0; tf < 4; ++tf) { make_cand_from_base(baseL, tf, &CL[tf]); make_cand_from_base(baseR, tf, &CR[tf]); }

                /* hàng xóm cố định: trái xa, phải xa, trên, dưới */
                int hasLeft = (j - 2) >= 0, hasRight = (j + 1) < BC;
                int leftNeigh[BH], rightNeigh[BH];
                if (hasLeft)  read_page_right_col(PAGE, r0, (j - 2) * BW, leftNeigh);
                if (hasRight) { for (int r = 0;r < BH;++r) rightNeigh[r] = PAGE[r0 + r][(j + 1) * BW] & 3; }

                int topL[BW], topR[BW], botL[BW], botR[BW];
                if (i > 0) {
                    read_page_top_row(PAGE, r0, c0L, topL);
                    read_page_top_row(PAGE, r0, c0R, topR);
                }
                if (i + 1 < BR) {
                    for (int c = 0;c < BW;++c) botL[c] = PAGE[r0 + BH][c0L + c] & 3;
                    for (int c = 0;c < BW;++c) botR[c] = PAGE[r0 + BH][c0R + c] & 3;
                }

                int best_tfL = 0, best_tfR = 0;
                int bestE = INT_MAX;

                for (int tL = 0;tL < 4;++tL) {
                    for (int tR = 0;tR < 4;++tR) {
                        int E = 0;
                        /* intra */
                        E += CL[tL].cost_intra + CR[tR].cost_intra;
                        /* vertical seams (trên/dưới) */
                        if (i > 0) {
                            E += cost_seam_rows(CL[tL].top, topL);
                            E += cost_seam_rows(CR[tR].top, topR);
                        }
                        if (i + 1 < BR) {
                            E += cost_seam_rows(botL, CL[tL].bottom);
                            E += cost_seam_rows(botR, CR[tR].bottom);
                        }
                        /* horizontal seams với hàng xóm xa (cố định) */
                        if (hasLeft)  E += cost_seam_cols(CL[tL].left, leftNeigh);
                        if (hasRight) E += cost_seam_cols(rightNeigh, CR[tR].right);
                        /* horizontal seam giữa L và R */
                        E += cost_seam_cols(CR[tR].left, CL[tL].right);

                        if (E < bestE) { bestE = E; best_tfL = tL; best_tfR = tR; }
                    }
                }

                /* Ghi nếu tốt hơn cấu hình đang có */
                /* Tính năng lượng hiện tại để so sánh (optional: nhưng ghi luôn cũng OK vì best>=current) */
                /* Đảm bảo không tăng cost: chỉ ghi khi tốt hơn bằng hoặc tốt hơn */
                /* Build best & write both blocks */
                write_block_to_page_from_cand(PAGE, r0, c0L, &CL[best_tfL]);
                write_block_to_page_from_cand(PAGE, r0, c0R, &CR[best_tfR]);
            }
        }
    }
}

/* Sweep dọc: tối ưu từng cặp (i-1, j)-(i, j) */
static void vertical_pair_sweep(int** PAGE, int* input_1Ddata, int Page_Size, int iters)
{
    const int BR = Page_Size / BH;
    const int BC = Page_Size / BW;

    for (int iter = 0; iter < iters; ++iter) {
        for (int j = 0; j < BC; ++j) {
            int c0 = j * BW;
            for (int i = 1; i < BR; ++i) {
                int r0T = (i - 1) * BH, r0B = i * BH;

                int inT[SYMS_PER_BLOCK], inB[SYMS_PER_BLOCK];
                for (int k = 0;k < SYMS_PER_BLOCK;++k) {
                    inT[k] = input_1Ddata[((i - 1) * BC + j) * SYMS_PER_BLOCK + k] & 3;
                    inB[k] = input_1Ddata[(i * BC + j) * SYMS_PER_BLOCK + k] & 3;
                }
                int baseT[N_PER_BLOCK], baseB[N_PER_BLOCK];
                build_base_grid_from_syms(inT, baseT);
                build_base_grid_from_syms(inB, baseB);
                Cand CT[4], CB[4];
                for (int tf = 0; tf < 4; ++tf) { make_cand_from_base(baseT, tf, &CT[tf]); make_cand_from_base(baseB, tf, &CB[tf]); }

                /* hàng xóm trái/phải cố định */
                int hasLeft = (j - 1) >= 0, hasRight = (j + 1) < BC;
                int leftTop[BH], rightTop[BH], leftBot[BH], rightBot[BH];
                if (hasLeft) { read_page_right_col(PAGE, r0T, (j - 1) * BW, leftTop);  read_page_right_col(PAGE, r0B, (j - 1) * BW, leftBot); }
                if (hasRight) { for (int r = 0;r < BH;++r) { rightTop[r] = PAGE[r0T + r][(j + 1) * BW] & 3; rightBot[r] = PAGE[r0B + r][(j + 1) * BW] & 3; } }

                /* hàng xóm trên/dưới xa */
                int hasTop = (i - 2) >= 0, hasBottom = (i + 1) < BR;
                int topRowAbove[BW], bottomRowBelow[BW];
                if (hasTop) { for (int c = 0;c < BW;++c) topRowAbove[c] = PAGE[r0T - 1][c0 + c] & 3; }
                if (hasBottom) { for (int c = 0;c < BW;++c) bottomRowBelow[c] = PAGE[r0B + BH][c0 + c] & 3; }

                int best_tT = 0, best_tB = 0;
                int bestE = INT_MAX;

                for (int tT = 0; tT < 4; ++tT) {
                    for (int tB = 0; tB < 4; ++tB) {
                        int E = 0;
                        /* intra */
                        E += CT[tT].cost_intra + CB[tB].cost_intra;
                        /* vertical seam giữa T và B */
                        E += cost_seam_rows(CT[tT].bottom, CB[tB].top);
                        /* vertical seams với top/bottom xa */
                        if (hasTop)    E += cost_seam_rows(CT[tT].top, topRowAbove);
                        if (hasBottom) E += cost_seam_rows(bottomRowBelow, CB[tB].bottom);
                        /* horizontal seams trái/phải (cố định) */
                        if (hasLeft) { E += cost_seam_cols(CT[tT].left, leftTop);  E += cost_seam_cols(CB[tB].left, leftBot); }
                        if (hasRight) { E += cost_seam_cols(rightTop, CT[tT].right); E += cost_seam_cols(rightBot, CB[tB].right); }

                        if (E < bestE) { bestE = E; best_tT = tT; best_tB = tB; }
                    }
                }

                write_block_to_page_from_cand(PAGE, r0T, c0, &CT[best_tT]);
                write_block_to_page_from_cand(PAGE, r0B, c0, &CB[best_tB]);
            }
        }
    }
}

/* ======================  PUBLIC API  ====================== */
void Encode_khanhProposal1112_4ary(int** PAGE, int* input_1Ddata, int Page_Size)
{
    init_edges_once(); /* cũng khởi init_positions_once() & init_masks_once() */
    const int BLOCK_ROWS = Page_Size / BH;
    const int BLOCK_COLS = Page_Size / BW;

    /* --- 1) Init tham lam như bản 158 --- */
    int in_block[SYMS_PER_BLOCK];
    int grid[N_PER_BLOCK];
    int tfid;
    for (int i = 0; i < BLOCK_ROWS; ++i) {
        for (int j = 0; j < BLOCK_COLS; ++j) {
            for (int k = 0; k < SYMS_PER_BLOCK; ++k)
                in_block[k] = input_1Ddata[(i * BLOCK_COLS + j) * SYMS_PER_BLOCK + k] & 3;
            int r0 = i * BH, c0 = j * BW;
            encode_block_greedy(in_block, grid, &tfid, PAGE, r0, c0, i, j);
            for (int r = 0; r < BH; ++r)
                for (int c = 0; c < BW; ++c)
                    PAGE[r0 + r][c0 + c] = grid[r * BW + c] & 3;
        }
    }

    /* --- 2) Pairwise refinement 2D: ngang rồi dọc, lặp 2 vòng --- */
    const int SWEEPS = 2; /* có thể 2–3; 2 thường đủ */
    horizontal_pair_sweep(PAGE, input_1Ddata, Page_Size, 1);
    vertical_pair_sweep(PAGE, input_1Ddata, Page_Size, 1);
    if (SWEEPS > 1) {
        horizontal_pair_sweep(PAGE, input_1Ddata, Page_Size, 1);
        vertical_pair_sweep(PAGE, input_1Ddata, Page_Size, 1);
    }

    /* Lấp mép + snapshot để QC */
    fill_uncovered_edges_after_encode(PAGE, Page_Size);
    dbg_store_encoded_page(PAGE, Page_Size);
}

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

            /* gỡ U theo index */
            int tf = grid12[g_index_pos] & 3;
            for (int p = 0; p < N_PER_BLOCK; ++p) {
                if (p == g_index_pos) continue;
                grid12[p] = apply_U_inv(tf, p, grid12[p]) & 3;
            }
            /* gỡ dither D (deterministic) */
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

    /* (Tuỳ chọn) QC — dán lại hàm compare_full_and_log_with_windows của bạn nếu cần */
    /* compare_full_and_log_with_windows(...); */
}
