#include <c10/cuda/CUDAGuard.h>
#include <torch/extension.h>
#include <tuple>

#define N_THREADS 256

#define CHECK_CUDA(x) TORCH_CHECK(x.is_cuda(), #x " must be a CUDA tensor")
#define CHECK_CONTIGUOUS(x) TORCH_CHECK(x.is_contiguous(), #x " must be contiguous")
#define CHECK_INPUT(x)                                                                 \
    CHECK_CUDA(x);                                                                     \
    CHECK_CONTIGUOUS(x)
#define DEVICE_GUARD(_ten)                                                             \
    const at::cuda::OptionalCUDAGuard device_guard(device_of(_ten));

// https://github.com/pytorch/pytorch/blob/233305a852e1cd7f319b15b5137074c9eac455f6/aten/src/ATen/cuda/cub.cuh#L38-L46
#define CUB_WRAPPER(func, ...)                                                         \
    do {                                                                               \
        size_t temp_storage_bytes = 0;                                                 \
        func(nullptr, temp_storage_bytes, __VA_ARGS__);                                \
        auto &caching_allocator = *::c10::cuda::CUDACachingAllocator::get();           \
        auto temp_storage = caching_allocator.allocate(temp_storage_bytes);            \
        func(temp_storage.get(), temp_storage_bytes, __VA_ARGS__);                     \
    } while (false)

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
isect_tiles_tensor(const torch::Tensor &means2d, // [N, 2]
                   const torch::Tensor &radii,   // [N]
                   const uint32_t tile_size,
                   const uint32_t tile_width, const uint32_t tile_height,
                   const bool sort, const bool double_buffer);

torch::Tensor isect_offset_encode_tensor(const torch::Tensor &isect_ids, // [n_isects]
                                         const uint32_t tile_width,
                                         const uint32_t tile_height);

torch::Tensor rasterize_to_pixels_fwd_tensor(
    // Gaussian parameters
    const torch::Tensor &means2d,                   // [N, 2]
    const torch::Tensor &conics,                    // [N, 3]
    const torch::Tensor &amplitude,                 // [N]
    const torch::Tensor &beta,                      // [N]
    // image size
    const uint32_t image_width, const uint32_t image_height, const uint32_t tile_size,
    // intersections
    const torch::Tensor &tile_offsets, // [C, tile_height, tile_width]
    const torch::Tensor &flatten_ids   // [n_isects]
);

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor>
rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    const torch::Tensor &means2d,                   // [N, 2]
    const torch::Tensor &conics,                    // [N, 3]
    const torch::Tensor &amplitude,                 // [N]
    const torch::Tensor &beta,                      // [N]
    // image size
    const uint32_t image_width, const uint32_t image_height, const uint32_t tile_size,
    // intersections
    const torch::Tensor &tile_offsets, // [tile_height, tile_width]
    const torch::Tensor &flatten_ids,  // [n_isects]

    // gradients of outputs
    const torch::Tensor &grad_render // [image_height, image_width]
    );
