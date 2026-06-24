// Stage 7: Python bindings for SiteGA training
// Build with: pip install . (uses pyproject.toml + setuptools + pybind11)
// Also exposes C ABI via sitega_train() for ctypes usage

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "sitega_train.h"
#include <cstring>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

// C++ wrapper for TrainParams that uses std::string for Python compatibility
struct PyTrainParams {
    std::string fg_path;       // full path to foreground FASTA (e.g. "data/PEAKS.fa")
    std::string bg_path;       // full path to background FASTA (e.g. "data/PEAKS_bg.fa")
    int max_lpd = 6;
    int motif_len = 12;
    int size = 28;
    int olig_bg = 6;
    int infc = 1;
    std::string out_path;
    int max_peak_len = 3000;
    std::string log_file;
    unsigned long seed = 0;
    int verbose = 0;
    int num_threads = 0;  // 0 = use OMP_NUM_THREADS env var
    int pop_size = 0;     // 0 = default 100, max 500 (MEGE)
    int num_motifs = 0;   // 0 = default 20, capped by pop_size

    // Internal storage for split paths (populated by to_c(), kept alive for C pointers)
    mutable std::string _fg_dir, _fg_file, _bg_dir, _bg_file;

    TrainParams to_c() const {
        // Split full paths into directory + filename for the C API
        // e.g. "data/PEAKS.fa" → dir="data/", file="PEAKS.fa"
        auto split_path = [](const std::string& full, std::string& dir, std::string& file) {
            size_t pos = full.find_last_of('/');
            if (pos == std::string::npos) {
                dir = "";
                file = full;
            } else {
                dir = full.substr(0, pos + 1);  // include trailing '/'
                file = full.substr(pos + 1);
            }
        };
        std::string fg_dir, fg_file, bg_dir, bg_file;
        split_path(fg_path, fg_dir, fg_file);
        split_path(bg_path, bg_dir, bg_file);
        _fg_dir = fg_dir; _fg_file = fg_file;
        _bg_dir = bg_dir; _bg_file = bg_file;

        TrainParams p;
        std::memset(&p, 0, sizeof(p));
        p.fg_path      = _fg_dir.c_str();
        p.fg_file      = _fg_file.c_str();
        p.bg_file      = _bg_file.c_str();
        p.bg_path      = _bg_dir.c_str();   // separate bg directory
        p.max_lpd      = max_lpd;
        p.motif_len    = motif_len;
        p.size         = size;
        p.olig_bg      = olig_bg;
        p.infc         = infc;
        p.out_path     = out_path.c_str();
        p.max_peak_len = max_peak_len;
        p.log_file     = log_file.c_str();
        p.seed         = seed;
        p.verbose      = verbose;
        p.pop_size     = pop_size;
        p.num_motifs   = num_motifs;
        return p;
    }
};

PYBIND11_MODULE(sitega, m) {
    m.doc() = "SiteGA: de novo motif search via genetic algorithm over LPDs";

    // Bind PyTrainParams with readable/writable attributes
    py::class_<PyTrainParams>(m, "TrainParams")
        .def(py::init<>())
        .def_readwrite("fg_path",      &PyTrainParams::fg_path)
        .def_readwrite("bg_path",      &PyTrainParams::bg_path)
        .def_readwrite("max_lpd",     &PyTrainParams::max_lpd)
        .def_readwrite("motif_len",   &PyTrainParams::motif_len)
        .def_readwrite("size",        &PyTrainParams::size)
        .def_readwrite("olig_bg",     &PyTrainParams::olig_bg)
        .def_readwrite("infc",        &PyTrainParams::infc)
        .def_readwrite("out_path",    &PyTrainParams::out_path)
        .def_readwrite("max_peak_len",&PyTrainParams::max_peak_len)
        .def_readwrite("log_file",    &PyTrainParams::log_file)
        .def_readwrite("seed",        &PyTrainParams::seed)
        .def_readwrite("verbose",     &PyTrainParams::verbose)
        .def_readwrite("num_threads",  &PyTrainParams::num_threads)
        .def_readwrite("pop_size",     &PyTrainParams::pop_size)
        .def_readwrite("num_motifs",   &PyTrainParams::num_motifs);

    // Convenience wrapper: train using TrainParams object
    m.def("train", [](const PyTrainParams& params) {
        TrainParams p = params.to_c();
#ifdef _OPENMP
        if (params.num_threads > 0)
            omp_set_num_threads(params.num_threads);
#endif
        TrainResult r;
        std::memset(&r, 0, sizeof(r));
        int rc = sitega_train(&p, &r);
        return py::make_tuple(
            rc,
            py::str(r.mat_path),
            py::str(r.loc_path),
            r.best_fit
        );
    }, py::arg("params"),
       "Train a SiteGA model from a TrainParams object. "
       "Returns (rc, mat_path, loc_path, best_fit). rc=0 on success.");

    // Convenience wrapper: train using keyword arguments
    m.def("train", [](
        const std::string& fg_path,
        const std::string& bg_path,
        int max_lpd, int motif_len, int size, int olig_bg, int infc,
        const std::string& out_path, int max_peak_len,
        const std::string& log_file,
        unsigned long seed, int verbose,
        int num_threads, int pop_size, int num_motifs
    ) {
        PyTrainParams p;
        p.fg_path = fg_path;
        p.bg_path = bg_path;
        p.max_lpd = max_lpd;
        p.motif_len = motif_len;
        p.size = size;
        p.olig_bg = olig_bg;
        p.infc = infc;
        p.out_path = out_path;
        p.max_peak_len = max_peak_len;
        p.log_file = log_file;
        p.seed = seed;
        p.verbose = verbose;
        p.num_threads = num_threads;
        p.pop_size = pop_size;
        p.num_motifs = num_motifs;
        TrainParams cp = p.to_c();
#ifdef _OPENMP
        if (num_threads > 0)
            omp_set_num_threads(num_threads);
#endif
        TrainResult r;
        std::memset(&r, 0, sizeof(r));
        int rc = sitega_train(&cp, &r);
        return py::make_tuple(
            rc,
            py::str(r.mat_path),
            py::str(r.loc_path),
            r.best_fit
        );
    }, py::arg("fg_path"), py::arg("bg_path"),
       py::arg("max_lpd"), py::arg("motif_len"), py::arg("size"),
       py::arg("olig_bg"), py::arg("infc"), py::arg("out_path"),
       py::arg("max_peak_len"), py::arg("log_file"),
       py::arg("seed") = 0UL, py::arg("verbose") = 0,
       py::arg("num_threads") = 0,
       py::arg("pop_size") = 0,
       py::arg("num_motifs") = 0,
       "Train a SiteGA model with keyword arguments. "
       "Returns (rc, mat_path, loc_path, best_fit). rc=0 on success.");
}