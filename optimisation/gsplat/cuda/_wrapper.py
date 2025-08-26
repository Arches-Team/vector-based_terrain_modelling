from typing import Callable, Optional, Tuple

import torch
from torch import Tensor


def _make_lazy_cuda_func(name: str) -> Callable:
    def call_cuda(*args, **kwargs):
        # pylint: disable=import-outside-toplevel
        from ._backend import _C

        return getattr(_C, name)(*args, **kwargs)

    return call_cuda


@torch.no_grad()
def isect_tiles(
    means2d: Tensor,  # [N, 2]
    radii: Tensor,  # [N]
    tile_size: int,
    tile_width: int,
    tile_height: int,
    sort: bool = True,
) -> Tuple[Tensor, Tensor, Tensor]:
    """Maps projected Gaussians to intersecting tiles.

    Args:
        means2d: Projected Gaussian means. [N, 2]
        radii: Maximum radii of the projected Gaussians. [N]
        tile_size: Tile size.
        tile_width: Tile width.
        tile_height: Tile height.
        sort: If True, the returned intersections will be sorted by the intersection ids. Default: True.

    Returns:
        A tuple:

        - **Tiles per Gaussian**. The number of tiles intersected by each Gaussian.
          Int32 [N]
        - **Intersection ids**. Each id is an 64-bit integer with the following
          information: tile_id.
        - **Flatten ids**. The global flatten indices in [N]
    """

    N, _ = means2d.shape
    assert means2d.shape == (N, 2), means2d.size()
    assert radii.shape == (N,), radii.size()

    tiles_per_gauss, isect_ids, flatten_ids = _make_lazy_cuda_func("isect_tiles")(
        means2d.contiguous(),
        radii.contiguous(),
        tile_size,
        tile_width,
        tile_height,
        sort,
        True,  # DoubleBuffer: memory efficient radixsort
    )
    return tiles_per_gauss, isect_ids, flatten_ids


@torch.no_grad()
def isect_offset_encode(
    isect_ids: Tensor, tile_width: int, tile_height: int
) -> Tensor:
    """Encodes intersection ids to offsets.

    Args:
        isect_ids: Intersection ids. [n_isects]
        tile_width: Tile width.
        tile_height: Tile height.

    Returns:
        Offsets. [tile_height, tile_width]
    """
    return _make_lazy_cuda_func("isect_offset_encode")(
        isect_ids.contiguous(), tile_width, tile_height
    )


def rasterize_to_pixels(
    means2d: Tensor,  # [N, 2]
    conics: Tensor,  # [N, 3]
    amplitude: Tensor,  # [N]
    beta: Tensor,  # [N]
    image_width: int,
    image_height: int,
    tile_size: int,
    isect_offsets: Tensor,  # [tile_height, tile_width]
    flatten_ids: Tensor,  # [n_isects]
) -> Tensor:
    """Rasterizes Gaussians to pixels.

    Args:
        means2d: Projected Gaussian means. [N, 2]
        conics: Inverse of the projected covariances with only upper triangle values. [N, 3]
        amplitude: Gaussian amplitude [N]
        beta: Gaussian amplitude [N]
        image_width: Image width.
        image_height: Image height.
        tile_size: Tile size.
        isect_offsets: Intersection offsets outputs from `isect_offset_encode()`. [tile_height, tile_width]
        flatten_ids: The global flatten indices in [N]

    Returns:
        A tensor:

        - **Rendered colors**. [image_height, image_width]
    """

    N = means2d.size(0)
    assert means2d.shape == (N, 2), means2d.shape
    assert conics.shape == (N, 3), conics.shape
    assert amplitude.shape == (N,), amplitude.shape

    tile_height, tile_width = isect_offsets.shape
    assert (
        tile_height * tile_size >= image_height
    ), f"Assert Failed: {tile_height} * {tile_size} >= {image_height}"
    assert (
        tile_width * tile_size >= image_width
    ), f"Assert Failed: {tile_width} * {tile_size} >= {image_width}"

    render_colors = _RasterizeToPixels.apply(
        means2d.contiguous(),
        conics.contiguous(),
        amplitude.contiguous(),
        beta.contiguous(),
        image_width,
        image_height,
        tile_size,
        isect_offsets.contiguous(),
        flatten_ids.contiguous(),
    )

    return render_colors


class _RasterizeToPixels(torch.autograd.Function):
    """Rasterize gaussians"""

    @staticmethod
    def forward(
        ctx,
        means2d: Tensor,  # [N, 2]
        conics: Tensor,  # [N, 3]
        amplitude: Tensor,  # [N]
        beta: Tensor,  # [N]
        width: int,
        height: int,
        tile_size: int,
        isect_offsets: Tensor,  # [tile_height, tile_width]
        flatten_ids: Tensor,  # [n_isects]
    ) -> Tensor:
        render_colors = _make_lazy_cuda_func(
            "rasterize_to_pixels_fwd"
        )(
            means2d,
            conics,
            amplitude,
            beta,
            width,
            height,
            tile_size,
            isect_offsets,
            flatten_ids,
        )

        ctx.save_for_backward(
            means2d,
            conics,
            amplitude,
            beta,
            isect_offsets,
            flatten_ids,
        )
        ctx.width = width
        ctx.height = height
        ctx.tile_size = tile_size

        return render_colors

    @staticmethod
    def backward(
        ctx,
        v_render_colors: Tensor,  # [H, W, 1]
    ):
        (
            means2d,
            conics,
            amplitude,
            beta,
            isect_offsets,
            flatten_ids,
        ) = ctx.saved_tensors
        width = ctx.width
        height = ctx.height
        tile_size = ctx.tile_size

        (
            v_means2d,
            v_conics,
            v_amplitude,
            v_beta,
        ) = _make_lazy_cuda_func("rasterize_to_pixels_bwd")(
            means2d,
            conics,
            amplitude,
            beta,
            width,
            height,
            tile_size,
            isect_offsets,
            flatten_ids,
            v_render_colors.contiguous(),
        )

        return (
            v_means2d,
            v_conics,
            v_amplitude,
            v_beta,
            None,
            None,
            None,
            None,
            None,
            None,
        )
