#ifndef SITEGA_TRAIN_H
#define SITEGA_TRAIN_H

#ifdef __cplusplus
extern "C" {
#endif

/* Parameters for SiteGA model training.
 * Mirrors the CLI arguments of andy05cell.cpp.
 * fg_path / bg_path / out_path: directories ending with '/' (Linux) or '\\' (Windows).
 * fg_file / bg_file: foreground and background FASTA file names (relative to fg_path/bg_path).
 * log_file: output log file name (relative to out_path).
 * seed: RNG seed (0 means time(NULL) for backward compatibility).
 */
struct TrainParams {
    const char* fg_path;    /* arg 1: path to fasta files */
    const char* fg_file;    /* arg 2: foreground fasta filename */
    const char* bg_file;    /* arg 3: background fasta filename */
    int max_lpd;            /* arg 4: maximal LPD length (default 6) */
    int motif_len;          /* arg 5: motif length in nucleotides */
    int size;               /* arg 6: number of LPDs */
    int olig_bg;            /* arg 7: k-mer length for background bias (default 6) */
    int infc;               /* arg 8: 1 = use info content in fitness, 0 = disregard */
    const char* bg_path;    /* optional: dir for bg_file (NULL = use fg_path) */
    const char* out_path;   /* arg 9: path to output files */
    int max_peak_len;       /* arg 10: maximal peak length (default 3000) */
    const char* log_file;   /* arg 11: output log filename */
    unsigned long seed;     /* optional: RNG seed (0 = time(NULL)) */
    int verbose;            /* optional: 0 = summary only, 1 = per-individual detail */
    int pop_size;           /* optional: population size (0 = default 100, max 500) */
    int num_motifs;         /* optional: number of best motifs to write (0 = default 20, capped by pop_size) */
};

/* Result of SiteGA model training.
 * mat_path / loc_path: full paths to the produced .mat and _loc files.
 * best_fit: fitness value of the best individual.
 * status: 0 on success, non-zero on error.
 */
struct TrainResult {
    char mat_path[512];
    char loc_path[512];
    double best_fit;
    int status;
};

/* Train a SiteGA model.
 * Returns 0 on success, non-zero on error.
 * On success, r->mat_path and r->loc_path contain the output file paths,
 * and r->best_fit contains the best fitness achieved.
 */
int sitega_train(const TrainParams* p, TrainResult* r);

#ifdef __cplusplus
}
#endif

#endif /* SITEGA_TRAIN_H */