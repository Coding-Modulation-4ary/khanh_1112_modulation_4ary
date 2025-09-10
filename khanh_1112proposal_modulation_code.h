// chiProposal46_4ary.c — Balanced 4-ary, 3x4 block, R = 11/12
// 11 data symbols + 1 index (transform id in {0..3})
// Seam-aware + DC penalty + run penalty theo HÀNG & CỘT.
// Base: (i*BLOCK_COLS + j)*SYMS_PER_BLOCK + k
//
// Public API (drop-in):
//   void Encode_khanhProposal1112_4ary(int** PAGE, int* input_1Ddata, int Page_Size);
//   void Decode_khanhProposal1112_4ary(int* output_1Ddata, double** PAGE, int Page_Size);

#include <stdio.h>
#include <stdlib.h>

/* =========================== Geometry =========================== */
#define BH 3
#define BW 4
#define N_PER_BLOCK (BH*BW)      /* 12 */
#define SYMS_PER_BLOCK 11        /* 11 data/block (R = 11/12) */

/* ===================== Vị trí data & index ====================== */
/* Đặt index VÀO GIỮA cho bền hơn khi detect theo y */
#define INDEX_POS (1*BW + 1)      /* (r=1,c=1) */

static int g_index_pos = INDEX_POS;
static int g_data_pos[SYMS_PER_BLOCK];

/* ================= Biến đổi đảo được trên {0,1,2,3} ================ */
static inline int tf_identity(int x) { return x & 3; }
static inline int tf_reflect(int x) { return (3 - x) & 3; }
static inline int tf_rotp1(int x) { return (x + 1) & 3; }
static inline int tf_rotm1(int x) { return (x + 3) & 3; }

typedef int (*sym_tf)(int);
static sym_tf TFWD[4] = { tf_identity, tf_reflect, tf_rotp1, tf_rotm1 };
static sym_tf TINV[4] = { tf_identity, tf_reflect, tf_rotm1, tf_rotp1 };

/* ===================== Trọng số cost ===================== */
static const int W_EXTREME = 4;  /* phạt {0,3} */
static const int W_DIFF = 2;  /* tổng |Δ| nội khối */
static const int W_RUN_ROW = 2;  /* chuỗi {0,3} theo HÀNG */
static const int W_RUN_COL = 2;  /* chuỗi {0,3} theo CỘT */
static const int W_SEAM = 10; /* biên trái/trên với khối đã ghi */
static const int W_SEAM_XTRM = 6;  /* + nếu là cặp 0-3/3-0 ở seam */
static const int W_DC = 4;  /* lệch tổng so với 18 */

/* ================== Lân cận nội khối (17 cạnh) ================== */
struct Edge { int u, v; };
static struct Edge EDGES[32];
static int E_CNT = 0;

/* ========== GLOBAL DEBUG SNAPSHOT (encoded PAGE phẳng) ========== */
static int* g_dbg_enc_flat = NULL;
static size_t g_dbg_enc_len = 0;
static int    g_dbg_enc_pagesz = 0;

/* =================== Init vị trí & cạnh =================== */
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

static void init_edges_once(void) {
    if (E_CNT) return;
    int e = 0;
    /* cạnh ngang */
    for (int r = 0; r < BH; ++r)
        for (int c = 0; c < BW - 1; ++c) {
            int u = r * BW + c;
            EDGES[e].u = u;
            EDGES[e].v = u + 1;
            ++e;
        }
    /* cạnh dọc */
    for (int r = 0; r < BH - 1; ++r)
        for (int c = 0; c < BW; ++c) {
            int u = r * BW + c;
            EDGES[e].u = u;
            EDGES[e].v = u + BW;
            ++e;
        }
    E_CNT = e; /* = 17 cho khối 3x4 */
    init_positions_once();
}

/* ===================== Cost helpers ===================== */
static inline int is_extreme(int x) { return (x == 0 || x == 3); }
static inline int is_bad_pair(int a, int b) { return ((a == 0 && b == 3) || (a == 3 && b == 0)); }

/* cost nội khối */
static int cost_intrablock(const int g[N_PER_BLOCK]) {
    int cost = 0;

    /* phạt {0,3} */
    if (W_EXTREME) {
        int cnt = 0; for (int i = 0; i < N_PER_BLOCK; ++i) cnt += is_extreme(g[i]);
        cost += W_EXTREME * cnt;
    }
    /* tổng |Δ| theo 17 cạnh */
    if (W_DIFF) {
        int s = 0; for (int i = 0; i < E_CNT; ++i) {
            int d = g[EDGES[i].u] - g[EDGES[i].v]; if (d < 0) d = -d; s += d;
        }
        cost += W_DIFF * s;
    }
    /* run penalty theo hàng */
    if (W_RUN_ROW) {
        for (int r = 0; r < BH; ++r) {
            int run = 0;
            for (int c = 0; c < BW; ++c) {
                int x = g[r * BW + c];
                if (is_extreme(x)) { ++run; if (run > 1) cost += W_RUN_ROW * (run - 1); }
                else run = 0;
            }
        }
    }
    /* run penalty theo cột */
    if (W_RUN_COL) {
        for (int c = 0; c < BW; ++c) {
            int run = 0;
            for (int r = 0; r < BH; ++r) {
                int x = g[r * BW + c];
                if (is_extreme(x)) { ++run; if (run > 1) cost += W_RUN_COL * (run - 1); }
                else run = 0;
            }
        }
    }
    /* DC về 18 (điểm giữa 12 ô với alphabet 0..3) */
    if (W_DC) {
        int sum = 0; for (int i = 0; i < N_PER_BLOCK; ++i) sum += g[i];
        int dev = sum - 18; if (dev < 0) dev = -dev;
        return cost + W_DC * dev;
    }
    return cost;
}

/* seam cost so với PAGE (đã có khối trái/trên) */
static int cost_seam(const int cand[N_PER_BLOCK], int** PAGE, int r0, int c0, int i, int j) {
    int cost = 0;
    /* biên trái */
    if (j > 0) {
        for (int r = 0; r < BH; ++r) {
            int a = cand[r * BW + 0];
            int b = PAGE[r0 + r][c0 - 1] & 3;
            int d = a - b; if (d < 0) d = -d;
            cost += W_SEAM * d;
            if (W_SEAM_XTRM && is_bad_pair(a, b)) cost += W_SEAM_XTRM;
            if (W_RUN_ROW && is_extreme(a) && is_extreme(b) && (a == b)) cost += W_RUN_ROW; /* run qua seam */
        }
    }
    /* biên trên */
    if (i > 0) {
        for (int c = 0; c < BW; ++c) {
            int a = cand[0 * BW + c];
            int b = PAGE[r0 - 1][c0 + c] & 3;
            int d = a - b; if (d < 0) d = -d;
            cost += W_SEAM * d;
            if (W_SEAM_XTRM && is_bad_pair(a, b)) cost += W_SEAM_XTRM;
            if (W_RUN_COL && is_extreme(a) && is_extreme(b) && (a == b)) cost += W_RUN_COL; /* run qua seam */
        }
    }
    return cost;
}

/* ========== Encode/Decode 1 block ========== */
static void encode_block_4x3_11of12(const int* in_syms,
    int out_grid[N_PER_BLOCK],
    int* chosen_tf,
    int** PAGE, int r0, int c0, int i, int j)
{
    int best_tf = 0, best_cost = 2147483647;
    int cand[N_PER_BLOCK], best[N_PER_BLOCK];
    for (int t = 0; t < N_PER_BLOCK; ++t) best[t] = 0;

    for (int tf = 0; tf < 4; ++tf) {
        for (int t = 0; t < N_PER_BLOCK; ++t) cand[t] = 0;

        /* ghi 11 symbol đã biến đổi */
        for (int k = 0; k < SYMS_PER_BLOCK; ++k) {
            int pos = g_data_pos[k];
            cand[pos] = TFWD[tf](in_syms[k] & 3);
        }
        /* ghi index */
        cand[g_index_pos] = (tf & 3);

        /* cost */
        int c = cost_intrablock(cand) + cost_seam(cand, PAGE, r0, c0, i, j);

        if (c < best_cost) {
            best_cost = c; best_tf = tf;
            for (int t = 0; t < N_PER_BLOCK; ++t) best[t] = cand[t];
        }
    }
    *chosen_tf = best_tf;
    for (int t = 0; t < N_PER_BLOCK; ++t) out_grid[t] = best[t];
}

static void decode_block_4x3_11of12(const int in_grid[N_PER_BLOCK], int* out_syms) {
    int tf = in_grid[g_index_pos] & 3;
    for (int k = 0; k < SYMS_PER_BLOCK; ++k) {
        int raw = in_grid[g_data_pos[k]] & 3;
        out_syms[k] = TINV[tf](raw) & 3;
    }
}

/* ========== Snapshot trang encode ========== */
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

/* ====== Lấp mép chưa phủ bởi lưới block (để compare/heatmap ổn định) ====== */
static void fill_uncovered_edges_after_encode(int** PAGE, int Page_Size) {
    const int full_cols = (Page_Size / BW) * BW;  /* số cột đã encode */
    const int full_rows = (Page_Size / BH) * BH;  /* số hàng đã encode */

    /* Lấp cột đuôi bằng cột mã hoá cuối */
    if (full_cols < Page_Size && full_cols > 0) {
        int csrc = full_cols - 1;
        for (int r = 0; r < full_rows; ++r) {
            int v = PAGE[r][csrc] & 3;
            for (int c = full_cols; c < Page_Size; ++c) PAGE[r][c] = v;
        }
    }
    /* (Phòng xa) nếu chiều dọc không chia hết cho BH */
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
    if (bj < 0) bj = 0; if (bj >= bc) bj = bc - 1;  /* CLAMP */

    int r0 = bi * BH, c0 = bj * BW;
    int ir = g_index_pos / BW, ic = g_index_pos % BW;
    size_t p = (size_t)(r0 + ir) * g_dbg_enc_pagesz + (c0 + ic);
    return g_dbg_enc_flat[p] & 3;
}

static int get_tf_from_dec_block_safe(double** PAGE_DBL, int bi, int bj, int Page_Size) {
    int br = Page_Size / BH;
    int bc = Page_Size / BW;
    if (bi < 0) bi = 0; if (bi >= br) bi = br - 1;
    if (bj < 0) bj = 0; if (bj >= bc) bj = bc - 1;

    int r0 = bi * BH, c0 = bj * BW;
    int ir = g_index_pos / BW, ic = g_index_pos % BW;
    /* Detector đã lượng hoá về {0..3} → chỉ cần ép kiểu */
    return (int)PAGE_DBL[r0 + ir][c0 + ic];
}

/* ========== FULL PAGE COMPARISON + NEIGHBOR WINDOWS (KHÔNG lượng tử) ========== */
static long long compare_full_and_log_with_windows(
    const char* summary_file,
    const char* coords_file,
    const char* heatmap_file,
    const char* perblock_csv,
    const char* windows_csv,       /* CSV chi tiết lân cận */
    const char* byrow_csv,         /* CSV mismatch theo hàng */
    double** PAGE_DBL, int Page_Size,
    int bh, int bw,
    int win                          /* bán kính cửa sổ: 1 → 3x3, 2 → 5x5 */
) {
    if (!g_dbg_enc_flat || g_dbg_enc_pagesz != Page_Size) {
        FILE* fsum = fopen(summary_file, "w");
        if (fsum) { fprintf(fsum, "NO_ENCODE_SNAPSHOT\n"); fclose(fsum); }
        return -1;
    }

    const int ENC_H = (Page_Size / bh) * bh;
    const int ENC_W = (Page_Size / bw) * bw;
    const int BR = (bh > 0) ? (ENC_H / bh) : 0;
    const int BC = (bw > 0) ? (ENC_W / bw) : 0;

    long long mism = 0;

    FILE* fcoords = fopen(coords_file, "w");
    FILE* fheat = fopen(heatmap_file, "w");
    FILE* fsum = fopen(summary_file, "w");
    FILE* fblk = fopen(perblock_csv, "w");
    FILE* fwin = fopen(windows_csv, "w");
    FILE* frow = fopen(byrow_csv, "w");

    long long* blk_cnt = NULL;
    if (BR > 0 && BC > 0) {
        blk_cnt = (long long*)malloc(sizeof(long long) * (size_t)BR * (size_t)BC);
        if (blk_cnt) for (int i = 0; i < BR * BC; ++i) blk_cnt[i] = 0;
    }
    long long* row_cnt = (long long*)malloc(sizeof(long long) * (size_t)Page_Size);
    if (row_cnt) for (int r = 0; r < Page_Size; ++r) row_cnt[r] = 0;

    if (fheat) fprintf(fheat, "# page_diff_full '.':match  'X':mismatch  ' ':outside encoded_area\n");

    if (fwin) {
        fprintf(fwin, "r,c,bi,bj,ri,ci,is_index,enc_tf,dec_tf,enc_sym,dec_sym,orig_sym,dec_data_sym,raw_center");
        const int winSz = (2 * win + 1) * (2 * win + 1);
        for (int t = 0; t < winSz; ++t) fprintf(fwin, ",enc_w%d", t);
        for (int t = 0; t < winSz; ++t) fprintf(fwin, ",dec_w%d", t);
        fputc('\n', fwin);
    }
    if (frow) fprintf(frow, "row,mismatch_count\n");

    size_t p = 0;
    for (int r = 0; r < Page_Size; ++r) {
        for (int c = 0; c < Page_Size; ++c, ++p) {
            const int in_encoded = (r < ENC_H && c < ENC_W);
            if (!in_encoded) { if (fheat) fputc(' ', fheat); continue; }

            int enc = g_dbg_enc_flat[p] & 3;
            /* Detector đã kẹp: chỉ cast */
            int dec = (int)PAGE_DBL[r][c];
            int diff = (enc != dec);

            if (diff) {
                ++mism;
                if (row_cnt) row_cnt[r]++;

                if (fcoords) fprintf(fcoords, "%d %d %d %d\n", r, c, enc, dec);

                if (blk_cnt) {
                    int bi = r / bh, bj = c / bw;
                    if (bi >= 0 && bi < BR && bj >= 0 && bj < BC) blk_cnt[bi * BC + bj] += 1;
                }

                if (fwin) {
                    int bi = r / bh, bj = c / bw;
                    int ri = r % bh, ci = c % bw;
                    int enc_tf = get_tf_from_enc_block(bi, bj);
                    int dec_tf = get_tf_from_dec_block_safe(PAGE_DBL, bi, bj, Page_Size);
                    int is_index = (ri == (g_index_pos / BW) && ci == (g_index_pos % BW));

                    int orig_sym = -1, dec_data_sym = -1;
                    if (!is_index) {
                        /* orig_sym = TINV[enc_tf](enc) */
                        switch (enc_tf & 3) {
                        case 0: orig_sym = enc & 3; break;
                        case 1: orig_sym = (3 - enc) & 3; break;
                        case 2: orig_sym = (enc + 3) & 3; break;
                        case 3: orig_sym = (enc + 1) & 3; break;
                        }
                        /* dec_data_sym = TINV[dec_tf](dec) */
                        switch (dec_tf & 3) {
                        case 0: dec_data_sym = dec & 3; break;
                        case 1: dec_data_sym = (3 - dec) & 3; break;
                        case 2: dec_data_sym = (dec + 3) & 3; break;
                        case 3: dec_data_sym = (dec + 1) & 3; break;
                        }
                    }

                    double raw_center = PAGE_DBL[r][c];

                    fprintf(fwin, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.10f",
                        r, c, bi, bj, ri, ci, is_index, enc_tf, dec_tf,
                        enc, dec, orig_sym, dec_data_sym, raw_center);

                    /* enc window (int) */
                    for (int dr = -win; dr <= win; ++dr) {
                        for (int dc = -win; dc <= win; ++dc) {
                            int rr = r + dr, cc = c + dc; int v = -1;
                            if (rr >= 0 && rr < Page_Size && cc >= 0 && cc < Page_Size) {
                                if (rr < ENC_H && cc < ENC_W) {
                                    size_t pp = (size_t)rr * g_dbg_enc_pagesz + cc;
                                    v = g_dbg_enc_flat[pp] & 3;
                                }
                                else v = -1;
                            }
                            fprintf(fwin, ",%d", v);
                        }
                    }
                    /* dec window (int) */
                    for (int dr = -win; dr <= win; ++dr) {
                        for (int dc = -win; dc <= win; ++dc) {
                            int rr = r + dr, cc = c + dc; int v = -1;
                            if (rr >= 0 && rr < Page_Size && cc >= 0 && cc < Page_Size) {
                                v = (int)PAGE_DBL[rr][cc];
                            }
                            fprintf(fwin, ",%d", v);
                        }
                    }
                    fputc('\n', fwin);
                }
            }

            if (fheat) fputc(diff ? 'X' : '.', fheat);
        }
        if (fheat) fputc('\n', fheat);
    }

    if (fsum) {
        double rate = (ENC_H > 0 && ENC_W > 0) ? ((double)mism / (double)(ENC_H * ENC_W)) : 0.0;
        fprintf(fsum, "page_size %d\n", Page_Size);
        fprintf(fsum, "encoded_area %d x %d  (BH=%d,BW=%d)\n", ENC_H, ENC_W, bh, bw);
        fprintf(fsum, "blocks %d x %d\n", BR, BC);
        fprintf(fsum, "mismatches %lld\n", mism);
        fprintf(fsum, "mismatch_rate_in_encoded_area %.9f\n", rate);
        fclose(fsum);
    }

    if (fcoords) fclose(fcoords);
    if (fheat)   fclose(fheat);

    if (fblk) {
        if (blk_cnt) {
            fprintf(fblk, "i,j,count\n");
            for (int bi = 0; bi < BR; ++bi)
                for (int bj = 0; bj < BC; ++bj) {
                    long long cnt = blk_cnt[bi * BC + bj];
                    if (cnt > 0) fprintf(fblk, "%d,%d,%lld\n", bi, bj, cnt);
                }
        }
        else {
            fprintf(fblk, "NO_BLOCK_GRID\n");
        }
        fclose(fblk);
    }

    if (frow) {
        if (row_cnt) {
            for (int r = 0; r < Page_Size; ++r) fprintf(frow, "%d,%lld\n", r, row_cnt[r]);
        }
        fclose(frow);
    }

    if (fwin) fclose(fwin);
    free(blk_cnt);
    free(row_cnt);
    return mism;
}

/* ======================  PUBLIC API  ====================== */
void Encode_khanhProposal1112_4ary(int** PAGE, int* input_1Ddata, int Page_Size)
{
    init_edges_once(); /* cũng khởi init_positions_once() */

    const int BLOCK_ROWS = Page_Size / BH;
    const int BLOCK_COLS = Page_Size / BW;

    int in_block[SYMS_PER_BLOCK];
    int grid[N_PER_BLOCK];
    int tfid;

    for (int i = 0; i < BLOCK_ROWS; ++i) {
        for (int j = 0; j < BLOCK_COLS; ++j) {

            for (int k = 0; k < SYMS_PER_BLOCK; ++k)
                in_block[k] = input_1Ddata[(i * BLOCK_COLS + j) * SYMS_PER_BLOCK + k] & 3;

            int r0 = i * BH, c0 = j * BW;

            encode_block_4x3_11of12(in_block, grid, &tfid, PAGE, r0, c0, i, j);

            for (int r = 0; r < BH; ++r)
                for (int c = 0; c < BW; ++c)
                    PAGE[r0 + r][c0 + c] = grid[r * BW + c] & 3;
        }
    }

    fill_uncovered_edges_after_encode(PAGE, Page_Size);
    dbg_store_encoded_page(PAGE, Page_Size);
}

void Decode_khanhProposal1112_4ary(int* output_1Ddata, double** PAGE, int Page_Size)
{
    init_edges_once();

    const int BLOCK_ROWS = Page_Size / BH;
    const int BLOCK_COLS = Page_Size / BW;

    int grid[N_PER_BLOCK];
    int out_block[SYMS_PER_BLOCK];

    for (int i = 0; i < BLOCK_ROWS; ++i) {
        for (int j = 0; j < BLOCK_COLS; ++j) {
            int r0 = i * BH, c0 = j * BW;

            for (int r = 0; r < BH; ++r) {
                for (int c = 0; c < BW; ++c) {
                    /* Detector đã lượng hoá về {0..3} → chỉ cast */
                    int v = (int)PAGE[r0 + r][c0 + c];
                    grid[r * BW + c] = v;
                }
            }

            decode_block_4x3_11of12(grid, out_block);

            for (int k = 0; k < SYMS_PER_BLOCK; ++k)
                output_1Ddata[(i * BLOCK_COLS + j) * SYMS_PER_BLOCK + k] = out_block[k] & 3;
        }
    }

    /* Giữ lại khối so sánh full-page (CSV, heatmap) nếu bạn vẫn dùng để QC */
    compare_full_and_log_with_windows(
        "compare_summary.txt",
        "mismatch_coords_full.txt",
        "page_diff_full.txt",
        "mismatch_per_block.csv",
        "mismatch_windows.csv",
        "mismatch_by_row.csv",
        PAGE, Page_Size, BH, BW,
        1
    );
}
