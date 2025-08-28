import torch
import utils
from gsplat import rasterization


def rasterize_cuda(gaussians_params, image_size, return_zero_det_gaussians: bool = False):

    sigma_x, sigma_y, mod_scale, theta, amplitude, x, y, beta = gaussians_params.unbind(dim=-1)

    sigma_x = sigma_x * mod_scale * image_size
    sigma_y = sigma_y * mod_scale * image_size

    # Cuda device is compulsory
    zero_tensor = torch.zeros(theta.shape, device='cuda')
    S = torch.stack(
        [torch.stack([sigma_x, zero_tensor], dim=-1),
         torch.stack([zero_tensor, sigma_y], dim=-1)],
        dim=-2
    )
    R = utils.build_rotation(theta)
    M = S @ R
    covariance = torch.transpose(M, 1, 2) @ M
    covariance = torch.stack([covariance[:, 0, 0], covariance[:, 0, 1], covariance[:, 1, 1]], dim=-1)

    # Check for positive semi-definiteness
    det = covariance[:, 0] * covariance[:, 2] - covariance[:, 1] * covariance[:, 1]
    keep_gaussians = det > 0.

    det = det[keep_gaussians]
    inv_det = 1. / det
    covariance = covariance[keep_gaussians]
    inv_covariance = torch.stack([covariance[:, 2] * inv_det, covariance[:, 1] * inv_det, covariance[:, 0] * inv_det],
                                 dim=-1)

    # change range
    x = ((x+1.)/2.)*image_size
    y = ((y+1.)/2.)*image_size

    img, _ = rasterization(torch.stack((x, y), dim=1).float()[keep_gaussians], theta.float()[keep_gaussians],
                           torch.stack((sigma_x, sigma_y), dim=1).float()[keep_gaussians],
                           inv_covariance.float(), amplitude.float()[keep_gaussians], beta.float()[keep_gaussians],
                           width=image_size, height=image_size, tile_size=16)

    if return_zero_det_gaussians:
        return img, ~keep_gaussians

    return img
