#include "bindings.h"
#include <torch/extension.h>

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
    m.def("isect_tiles", &isect_tiles_tensor);
    m.def("isect_offset_encode", &isect_offset_encode_tensor);

    m.def("rasterize_to_pixels_fwd", &rasterize_to_pixels_fwd_tensor);
    m.def("rasterize_to_pixels_bwd", &rasterize_to_pixels_bwd_tensor);
}