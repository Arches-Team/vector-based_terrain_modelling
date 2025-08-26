#include "bindings.h"
#include <cooperative_groups.h>
#include <cub/cub.cuh>
#include <cuda_runtime.h>

namespace cg = cooperative_groups;

/****************************************************************************
 * Gaussian Tile Intersection
 ****************************************************************************/

__global__ void isect_tiles(
    // parallelize over N
    const uint32_t N,
    // data
    const float2 *__restrict__ means2d,              // [N, 2]
    const int32_t *__restrict__ radii,               // [N]
    const int64_t *__restrict__ cum_tiles_per_gauss, // [N]
    const uint32_t tile_size, const uint32_t tile_width, const uint32_t tile_height,
    int32_t *__restrict__ tiles_per_gauss, // [N]
    int64_t *__restrict__ isect_ids,       // [n_isects]
    int32_t *__restrict__ flatten_ids      // [n_isects]
) {
    // parallelize over N.
    uint32_t idx = cg::this_grid().thread_rank();
    bool first_pass = cum_tiles_per_gauss == nullptr;
    if (idx >=  N)
        return;
    if (radii[idx] <= 0) {
        if (first_pass)
            tiles_per_gauss[idx] = 0;
        return;
    }

    float tile_radius = radii[idx] / static_cast<float>(tile_size);
    float tile_x = means2d[idx].x / tile_size;
    float tile_y = means2d[idx].y / tile_size;

    // tile_min is inclusive, tile_max is exclusive
    uint2 tile_min, tile_max;
    tile_min.x = min(max(0, (uint32_t)floor(tile_x - tile_radius)), tile_width);
    tile_min.y = min(max(0, (uint32_t)floor(tile_y - tile_radius)), tile_height);
    tile_max.x = min(max(0, (uint32_t)ceil(tile_x + tile_radius)), tile_width);
    tile_max.y = min(max(0, (uint32_t)ceil(tile_y + tile_radius)), tile_height);

    if (first_pass) {
        // first pass only writes out tiles_per_gauss
        tiles_per_gauss[idx] =
            static_cast<int32_t>((tile_max.y - tile_min.y) * (tile_max.x - tile_min.x));
        return;
    }

    int64_t cur_idx = (idx == 0) ? 0 : cum_tiles_per_gauss[idx - 1];
    for (int32_t i = tile_min.y; i < tile_max.y; ++i) {
        for (int32_t j = tile_min.x; j < tile_max.x; ++j) {
            int64_t tile_id = i * tile_width + j;
            isect_ids[cur_idx] = tile_id;
            // the flatten index in [N]
            flatten_ids[cur_idx] = static_cast<int32_t>(idx);
            ++cur_idx;
        }
    }
}

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
isect_tiles_tensor(const torch::Tensor &means2d, // [N, 2]
                   const torch::Tensor &radii,   // [N]
                   const uint32_t tile_size,
                   const uint32_t tile_width, const uint32_t tile_height,
                   const bool sort, const bool double_buffer) {
    DEVICE_GUARD(means2d);
    CHECK_INPUT(means2d);
    CHECK_INPUT(radii);

    uint32_t N, nnz, total_elems;

    N = means2d.size(0); // number of gaussians
    total_elems = N;

    uint32_t n_tiles = tile_width * tile_height;
    at::cuda::CUDAStream stream = at::cuda::getCurrentCUDAStream();

    // first pass: compute number of tiles per gaussian
    torch::Tensor tiles_per_gauss =
        torch::empty_like(radii, radii.options().dtype(torch::kInt32));

    int64_t n_isects;
    torch::Tensor cum_tiles_per_gauss;
    if (total_elems) {
        isect_tiles<<<(total_elems + N_THREADS - 1) / N_THREADS, N_THREADS, 0,
                      stream>>>(
            N,
            (float2 *)means2d.data_ptr<float>(), radii.data_ptr<int32_t>(),
            nullptr, tile_size, tile_width, tile_height,
            tiles_per_gauss.data_ptr<int32_t>(), nullptr, nullptr);
        cum_tiles_per_gauss = torch::cumsum(tiles_per_gauss.view({-1}), 0);
        n_isects = cum_tiles_per_gauss[-1].item<int64_t>();
    } else {
        n_isects = 0;
    }

    // second pass: compute isect_ids and flatten_ids as a packed tensor
    torch::Tensor isect_ids =
        torch::empty({n_isects}, radii.options().dtype(torch::kInt64));
    torch::Tensor flatten_ids =
        torch::empty({n_isects}, radii.options().dtype(torch::kInt32));
    if (n_isects) {
        isect_tiles<<<(total_elems + N_THREADS - 1) / N_THREADS, N_THREADS, 0,
                      stream>>>(
            N,
            (float2 *)means2d.data_ptr<float>(), radii.data_ptr<int32_t>(),
            cum_tiles_per_gauss.data_ptr<int64_t>(),
            tile_size, tile_width, tile_height, nullptr,
            isect_ids.data_ptr<int64_t>(), flatten_ids.data_ptr<int32_t>());
    }

    // optionally sort the Gaussians by isect_ids
    if (n_isects && sort) {
        torch::Tensor isect_ids_sorted = torch::empty_like(isect_ids);
        torch::Tensor flatten_ids_sorted = torch::empty_like(flatten_ids);

        // https://nvidia.github.io/cccl/cub/api/structcub_1_1DeviceRadixSort.html
        // DoubleBuffer reduce the auxiliary memory usage from O(N+P) to O(P)
        if (double_buffer) {
            // Create a set of DoubleBuffers to wrap pairs of device pointers
            cub::DoubleBuffer<int64_t> d_keys(isect_ids.data_ptr<int64_t>(),
                                              isect_ids_sorted.data_ptr<int64_t>());
            cub::DoubleBuffer<int32_t> d_values(flatten_ids.data_ptr<int32_t>(),
                                                flatten_ids_sorted.data_ptr<int32_t>());
            CUB_WRAPPER(cub::DeviceRadixSort::SortPairs, d_keys, d_values, n_isects, 0, 64,
                        stream);
            switch (d_keys.selector) {
            case 0: // sorted items are stored in isect_ids
                isect_ids_sorted = isect_ids;
                break;
            case 1: // sorted items are stored in isect_ids_sorted
                break;
            }
            switch (d_values.selector) {
            case 0: // sorted items are stored in flatten_ids
                flatten_ids_sorted = flatten_ids;
                break;
            case 1: // sorted items are stored in flatten_ids_sorted
                break;
            }
            // printf("DoubleBuffer d_keys selector: %d\n", d_keys.selector);
            // printf("DoubleBuffer d_values selector: %d\n", d_values.selector);
        } else {
            CUB_WRAPPER(cub::DeviceRadixSort::SortPairs, isect_ids.data_ptr<int64_t>(),
                        isect_ids_sorted.data_ptr<int64_t>(),
                        flatten_ids.data_ptr<int32_t>(),
                        flatten_ids_sorted.data_ptr<int32_t>(), n_isects, 0, 64, stream);
        }
        return std::make_tuple(tiles_per_gauss, isect_ids_sorted, flatten_ids_sorted);
    } else {
        return std::make_tuple(tiles_per_gauss, isect_ids, flatten_ids);
    }
}

__global__ void isect_offset_encode(const uint32_t n_isects,
                                    const int64_t *__restrict__ isect_ids,
                                    const uint32_t n_tiles,
                                    int32_t *__restrict__ offsets // [n_tiles]
) {
    // e.g., ids: [1, 1, 1, 3, 3], n_tiles = 6
    // counts: [0, 3, 0, 2, 0, 0]
    // cumsum: [0, 3, 3, 5, 5, 5]
    // offsets: [0, 0, 3, 3, 5, 5]
    uint32_t idx = cg::this_grid().thread_rank();
    if (idx >= n_isects)
        return;

    int64_t id_curr = isect_ids[idx];

    if (idx == 0) {
        // write out the offsets until the first valid tile (inclusive)
        for (uint32_t i = 0; i < id_curr + 1; ++i)
            offsets[i] = static_cast<int32_t>(idx);
    }
    if (idx == n_isects - 1) {
        // write out the rest of the offsets
        for (uint32_t i = id_curr + 1; i < n_tiles; ++i)
            offsets[i] = static_cast<int32_t>(n_isects);
    }

    if (idx > 0) {
        // visit the current and previous isect_id and check if the (cid, tile_id)
        // pair changes.
        int64_t isect_id_prev = isect_ids[idx - 1];
        if (isect_id_prev == id_curr)
            return;

        // write out the offsets between the previous and current tiles
        int64_t id_prev = isect_id_prev;
        for (uint32_t i = id_prev + 1; i < id_curr + 1; ++i)
            offsets[i] = static_cast<int32_t>(idx);
    }
}

torch::Tensor isect_offset_encode_tensor(const torch::Tensor &isect_ids, // [n_isects]
                                         const uint32_t tile_width,
                                         const uint32_t tile_height) {
    DEVICE_GUARD(isect_ids);
    CHECK_INPUT(isect_ids);

    uint32_t n_isects = isect_ids.size(0);
    torch::Tensor offsets = torch::empty({tile_height, tile_width},
                                         isect_ids.options().dtype(torch::kInt32));
    if (n_isects) {
        uint32_t n_tiles = tile_width * tile_height;
        at::cuda::CUDAStream stream = at::cuda::getCurrentCUDAStream();
        isect_offset_encode<<<(n_isects + N_THREADS - 1) / N_THREADS, N_THREADS, 0,
                              stream>>>(n_isects, isect_ids.data_ptr<int64_t>(),
                                        n_tiles, offsets.data_ptr<int32_t>());
    } else {
        offsets.fill_(0);
    }
    return offsets;
}

/****************************************************************************
 * Rasterization
 ****************************************************************************/

__global__ void rasterize_to_pixels_fwd_kernel(
    const uint32_t N, const uint32_t n_isects,
    const float2 *__restrict__ means2d,    // [N, 2]
    const float3 *__restrict__ conics,     // [N, 3]
    const float *__restrict__ amplitude,   // [N]
    const float *__restrict__ beta,        // [N]
    const uint32_t image_width, const uint32_t image_height, const uint32_t tile_size,
    const uint32_t tile_width, const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    float *__restrict__ render_colors // [image_height, image_width]
) {
    // each thread draws one pixel, but also timeshares caching gaussians in a
    // shared tile

    auto block = cg::this_thread_block();
    int32_t tile_id = block.group_index().y * tile_width + block.group_index().z;
    uint32_t i = block.group_index().y * tile_size + block.thread_index().y;
    uint32_t j = block.group_index().z * tile_size + block.thread_index().x;

    float px = (float)j + 0.5f;
    float py = (float)i + 0.5f;
    int32_t pix_id = i * image_width + j;

    // return if out of bounds
    // keep not rasterizing threads around for reading data
    bool inside = (i < image_height && j < image_width);
    bool done = !inside;

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    const uint32_t block_size = block.size();
    uint32_t num_batches = (range_end - range_start + block_size - 1) / block_size;

    extern __shared__ int s[];
    int32_t *id_batch = (int32_t *)s;                              // [block_size]
    float2 *xy_batch = (float2 *)&id_batch[block_size];    // [block_size]
    float3 *conic_batch = (float3 *)&xy_batch[block_size]; // [block_size]
    float *amplitude_batch = (float *)&conic_batch[block_size]; // [block_size]
    float *beta_batch = (float *)&amplitude_batch[block_size]; // [block_size]

    // index of most recent gaussian to write to this thread's pixel
    uint32_t cur_idx = 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing its
    // designated pixel
    uint32_t tr = block.thread_rank();

    float pix_out = 0.f;
    for (uint32_t b = 0; b < num_batches; ++b) {
        // resync all threads before beginning next batch
        // end early if entire tile is done
        if (__syncthreads_count(done) >= block_size) {
            break;
        }

        // each thread fetch 1 gaussian from front to back
        // index of gaussian to load
        uint32_t batch_start = range_start + block_size * b;
        uint32_t idx = batch_start + tr;
        if (idx < range_end) {
            int32_t g = flatten_ids[idx]; // flatten index in [N]
            id_batch[tr] = g;
            xy_batch[tr] = means2d[g];
            conic_batch[tr] = conics[g];
            amplitude_batch[tr] = amplitude[g];
            beta_batch[tr] = beta[g];
        }

        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        uint32_t batch_size = min(block_size, range_end - batch_start);
        for (uint32_t t = 0; (t < batch_size) && !done; ++t) {
            const float3 conic = conic_batch[t];
            const float2 xy = xy_batch[t];
            const float2 delta = {xy.x - px, xy.y - py};
            const float sigma =
                0.5f * (conic.x * delta.x * delta.x + conic.z * delta.y * delta.y) +
                conic.y * delta.x * delta.y;
            float alpha = __expf(-__powf(sigma, beta_batch[t])) * amplitude_batch[t];
            if (sigma < 0.f) {
                continue;
            }

            int32_t g = id_batch[t];
            const float vis = alpha;
            pix_out += vis;

            cur_idx = batch_start + t;
        }
    }

    if (inside) {
        render_colors[pix_id] = pix_out;
    }
}

torch::Tensor rasterize_to_pixels_fwd_tensor(
    // Gaussian parameters
    const torch::Tensor &means2d,   // [N, 2]
    const torch::Tensor &conics,    // [N, 3]
    const torch::Tensor &amplitude, // [N]
    const torch::Tensor &beta,      // [N]
    // image size
    const uint32_t image_width, const uint32_t image_height, const uint32_t tile_size,
    // intersections
    const torch::Tensor &tile_offsets, // [tile_height, tile_width]
    const torch::Tensor &flatten_ids   // [n_isects]
) {
    DEVICE_GUARD(means2d);
    CHECK_INPUT(means2d);
    CHECK_INPUT(conics);
    CHECK_INPUT(amplitude);
    CHECK_INPUT(beta);
    CHECK_INPUT(tile_offsets);
    CHECK_INPUT(flatten_ids);

    uint32_t N = means2d.size(0); // number of gaussians
    uint32_t tile_height = tile_offsets.size(0);
    uint32_t tile_width = tile_offsets.size(1);
    uint32_t n_isects = flatten_ids.size(0);

    // Each block covers a tile on the image. In total there are
    // 1 * tile_height * tile_width blocks.
    dim3 threads = {tile_size, tile_size, 1};
    dim3 blocks = {1, tile_height, tile_width};

    torch::Tensor renders = torch::empty({image_height, image_width},
                                         means2d.options().dtype(torch::kFloat32));

    at::cuda::CUDAStream stream = at::cuda::getCurrentCUDAStream();
    const uint32_t shared_mem =
        tile_size * tile_size * (sizeof(int32_t) + sizeof(float2) + sizeof(float3) + sizeof(float) + sizeof(float));

    if (cudaFuncSetAttribute(rasterize_to_pixels_fwd_kernel,
                             cudaFuncAttributeMaxDynamicSharedMemorySize,
                             shared_mem) != cudaSuccess) {
        AT_ERROR("Failed to set maximum shared memory size (requested ", shared_mem,
                 " bytes), try lowering tile_size.");
    }
    rasterize_to_pixels_fwd_kernel<<<blocks, threads, shared_mem, stream>>>(
        N, n_isects, (float2 *)means2d.data_ptr<float>(),
        (float3 *)conics.data_ptr<float>(), (float *)amplitude.data_ptr<float>(),
        (float *)beta.data_ptr<float>(),
        image_width, image_height, tile_size, tile_width, tile_height,
        tile_offsets.data_ptr<int32_t>(), flatten_ids.data_ptr<int32_t>(),
        renders.data_ptr<float>());

    return renders;
}

__global__ void rasterize_to_pixels_bwd_kernel(
    const uint32_t N, const uint32_t n_isects,
    const float2 *__restrict__ means2d,    // [N, 2]
    const float3 *__restrict__ conics,     // [N, 3]
    const float *__restrict__ amplitude,   // [N]
    const float *__restrict__ beta,   // [N]
    const uint32_t image_width, const uint32_t image_height, const uint32_t tile_size,
    const uint32_t tile_width, const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    const float *__restrict__ render_colors_grad, // [image_height, image_width]
    float2 *__restrict__ means2d_grad,    // [N, 2]
    float3 *__restrict__ conics_grad,     // [N, 3]
    float *__restrict__ amplitude_grad,    // [N]
    float *__restrict__ beta_grad    // [N]
) {
    // each thread processes one pixel, but also timeshares caching gaussians in a
    // shared tile

    auto block = cg::this_thread_block();
    int32_t tile_id = block.group_index().y * tile_width + block.group_index().z;
    uint32_t i = block.group_index().y * tile_size + block.thread_index().y;
    uint32_t j = block.group_index().z * tile_size + block.thread_index().x;

    float px = (float)j + 0.5f;
    float py = (float)i + 0.5f;
    int32_t pix_id = i * image_width + j;

    // return if out of bounds
    bool inside = (i < image_height && j < image_width);
    bool done = !inside;

    // have all threads in tile process the same gaussians in batches
    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    const uint32_t block_size = block.size();
    uint32_t num_batches = (range_end - range_start + block_size - 1) / block_size;

    extern __shared__ int s[];
    int32_t *id_batch = (int32_t *)s;                              // [block_size]
    float2 *xy_batch = (float2 *)&id_batch[block_size];    // [block_size]
    float3 *conic_batch = (float3 *)&xy_batch[block_size]; // [block_size]
    float *amplitude_batch = (float *)&conic_batch[block_size]; // [block_size]
    float *beta_batch = (float *)&amplitude_batch[block_size]; // [block_size]

    // index of most recent gaussian to write to this thread's pixel
    uint32_t cur_idx = 0;

    // collect and process batches of gaussians
    uint32_t tr = block.thread_rank();

    float grad_out = inside ? render_colors_grad[pix_id] : 0.f;
    for (uint32_t b = 0; b < num_batches; ++b) {
        if (__syncthreads_count(done) >= block_size) {
            break;
        }

        uint32_t batch_start = range_start + block_size * b;
        uint32_t idx = batch_start + tr;
        if (idx < range_end) {
            int32_t g = flatten_ids[idx];
            id_batch[tr] = g;
            xy_batch[tr] = means2d[g];
            conic_batch[tr] = conics[g];
            amplitude_batch[tr] = amplitude[g];
            beta_batch[tr] = beta[g];
        }

        block.sync();

        uint32_t batch_size = min(block_size, range_end - batch_start);
        for (uint32_t t = 0; (t < batch_size) && !done; ++t) {
            const float3 conic = conic_batch[t];
            const float2 xy = xy_batch[t];
            const float2 delta = {xy.x - px, xy.y - py};
            const float sigma =
                0.5f * (conic.x * delta.x * delta.x + conic.z * delta.y * delta.y) +
                conic.y * delta.x * delta.y;
            float pow = __powf(sigma, beta_batch[t]);
            float exp = __expf(-pow);
            float I = exp * amplitude_batch[t];

            if (sigma < 0.f) {
                continue;
            }

            int32_t g = id_batch[t];
            float vis_grad = grad_out;

            float I_grad = vis_grad;
            float p_grad = I_grad * -exp * amplitude_batch[t];
            float sigma_grad = p_grad * beta_batch[t] * __powf(sigma, beta_batch[t]-1);
            float2 delta_grad;
            delta_grad.x = sigma_grad * (conic.x * delta.x + conic.y * delta.y);
            delta_grad.y = sigma_grad * (conic.z * delta.y + conic.y * delta.x);

            float2 xy_grad;
            xy_grad.x = delta_grad.x;
            xy_grad.y = delta_grad.y;

            atomicAdd(&means2d_grad[g].x, xy_grad.x);
            atomicAdd(&means2d_grad[g].y, xy_grad.y);
            atomicAdd(&conics_grad[g].x, sigma_grad * 0.5f * delta.x * delta.x);
            atomicAdd(&conics_grad[g].y, sigma_grad * delta.x * delta.y);
            atomicAdd(&conics_grad[g].z, sigma_grad * 0.5f * delta.y * delta.y);
            atomicAdd(&amplitude_grad[g], I_grad * exp);
            atomicAdd(&beta_grad[g], I_grad * p_grad * pow * __logf(sigma));

            cur_idx = batch_start + t;
        }
    }
}

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor>
rasterize_to_pixels_bwd_tensor(
    // Gaussian parameters
    const torch::Tensor &means2d,                   // [N, 2]
    const torch::Tensor &conics,                    // [N, 3]
    const torch::Tensor &amplitude,                 // [N]
    const torch::Tensor &beta,                 // [N]
    // image size
    const uint32_t image_width, const uint32_t image_height, const uint32_t tile_size,
    // intersections
    const torch::Tensor &tile_offsets, // [tile_height, tile_width]
    const torch::Tensor &flatten_ids,  // [n_isects]
    // gradients of outputs
    const torch::Tensor &grad_render // [image_height, image_width]
    ) {
    DEVICE_GUARD(means2d);
    CHECK_INPUT(means2d);
    CHECK_INPUT(conics);
    CHECK_INPUT(beta);
    CHECK_INPUT(tile_offsets);
    CHECK_INPUT(flatten_ids);
    CHECK_INPUT(grad_render);

    uint32_t N = means2d.size(0); // number of gaussians
    uint32_t n_isects = flatten_ids.size(0);
    uint32_t tile_height = tile_offsets.size(0);
    uint32_t tile_width = tile_offsets.size(1);

    // Each block covers a tile on the image. In total there are
    // tile_height * tile_width blocks.
    dim3 threads = {tile_size, tile_size, 1};
    dim3 blocks = {1, tile_height, tile_width};

    torch::Tensor grad_means2d = torch::zeros_like(means2d);
    torch::Tensor grad_conics = torch::zeros_like(conics);
    torch::Tensor grad_amplitude = torch::zeros_like(amplitude);
    torch::Tensor grad_beta = torch::zeros_like(beta);

    if (n_isects) {
        const uint32_t shared_mem = tile_size * tile_size *
                                    (sizeof(int32_t) + sizeof(float2) + sizeof(float3) +
                                     sizeof(float) + sizeof(float));
        at::cuda::CUDAStream stream = at::cuda::getCurrentCUDAStream();
        if (cudaFuncSetAttribute(rasterize_to_pixels_bwd_kernel,
                                 cudaFuncAttributeMaxDynamicSharedMemorySize,
                                 shared_mem) != cudaSuccess) {
            AT_ERROR("Failed to set maximum shared memory size (requested ",
                     shared_mem, " bytes), try lowering tile_size.");
        }
        rasterize_to_pixels_bwd_kernel<<<blocks, threads, shared_mem, stream>>>(
            N, n_isects, (float2 *)means2d.data_ptr<float>(),
            (float3 *)conics.data_ptr<float>(), amplitude.data_ptr<float>(),
            beta.data_ptr<float>(),
            image_width, image_height, tile_size, tile_width, tile_height,
            tile_offsets.data_ptr<int32_t>(), flatten_ids.data_ptr<int32_t>(),
            grad_render.data_ptr<float>(),
            (float2 *)grad_means2d.data_ptr<float>(),
            (float3 *)grad_conics.data_ptr<float>(),
            grad_amplitude.data_ptr<float>(),
            grad_beta.data_ptr<float>()
            );
    }

    return std::make_tuple(grad_means2d, grad_conics, grad_amplitude, grad_beta);
}
