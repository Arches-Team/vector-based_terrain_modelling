import math
from typing import Dict, Optional, Tuple

import torch
from torch import Tensor
from typing_extensions import Literal

from .cuda._wrapper import (
    isect_offset_encode,
    isect_tiles,
    rasterize_to_pixels,
)


def rasterization(
    means2d: Tensor,  # [N, 2]
    theta: Tensor,  # [N]
    scales: Tensor,  # [N, 2]
    conics: Tensor,  # [N, 3]
    amplitude: Tensor,  # [N]
    beta: Tensor,  # [N]
    width: int,
    height: int,
    tile_size: int = 16,
) -> Tuple[Tensor, Dict]:
    """Rasterize a set of 2D Gaussians (N).

    This function provides a handful features for 2D Gaussian rasterization, which
    we detail in the following notes. A complete profiling of the these features
    can be found in the :ref:`profiling` page.

    .. note::
        **Batch Rasterization**: This function allows for rasterizing a set of 2D Gaussians

    Args:
        means2d: The 3D centers of the Gaussians. [N, 2]
        theta: The rotation in radian of the Gaussians. [N]
        scales: The scales of the Gaussians in pixels. [N, 2]
        conics: Inverse of the projected covariances with only upper triangle values. [N, 3]
        amplitude: The amplitude of the Gaussians. [N]
        beta: The beta of the Gaussians (GES). [N]
        width: The width of the image.
        height: The height of the image.
        tile_size: The size of the tiles for rasterization. Default is 16.
            (Note: other values are not tested)

    Returns:
        A tuple:

        **render_colors**: The rendered colors. [height, width]

        **meta**: A dictionary of intermediate results of the rasterization.

    """

    N = means2d.shape[0]
    assert means2d.shape == (N, 2), means2d.shape
    assert theta.shape == (N,), theta.shape
    assert scales.shape == (N, 2), scales.shape
    assert amplitude.shape == (N,), scales.shape

    radii = torch.max(scales, dim=1).values.int() * 3

    # Identify intersecting tiles
    tile_width = math.ceil(width / float(tile_size))
    tile_height = math.ceil(height / float(tile_size))
    tiles_per_gauss, isect_ids, flatten_ids = isect_tiles(
        means2d,
        radii,
        tile_size,
        tile_width,
        tile_height,
        sort=True,
    )
    isect_offsets = isect_offset_encode(isect_ids, tile_width, tile_height)

    render_colors = rasterize_to_pixels(
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

    meta = {
        "radii": radii,
        "means2d": means2d,
        "conics": conics,
        "amplitude": amplitude,
        "beta": beta,
        "tile_width": tile_width,
        "tile_height": tile_height,
        "tiles_per_gauss": tiles_per_gauss,
        "isect_ids": isect_ids,
        "flatten_ids": flatten_ids,
        "isect_offsets": isect_offsets,
        "width": width,
        "height": height,
        "tile_size": tile_size,
    }
    return render_colors, meta
