import os
import time

import pyrallis
import numpy as np
import torch.nn as nn
import torch
import cv2
import lpips

from gaussian_model import GaussianModel
from losses import d_ssim_loss
from config import Config
from optimise_utils import optimisation_loop, backup_code
from rasterizer import rasterize_cuda


def main():
    cfg = pyrallis.parse(Config)
    print('The generated config is:')
    print(pyrallis.dump(cfg))

    # Copy some file to the result folder for backup
    file_to_copy = [
        'optimise.py',
        'optimise_utils.py',
        'inference.py',
        'rasterizer.py',
        'gaussian_model.py',
        'utils.py',
        'losses.py',
        'config.py',
    ]

    os.makedirs(cfg.log.results_folder, exist_ok=True)
    backup_code(cfg.log.results_folder, file_to_copy, cfg)

    if cfg.optim.image_file_name.endswith('.png'):
        original_image = cv2.imread(cfg.optim.image_file_name, cv2.IMREAD_UNCHANGED)
        original_image = original_image / np.iinfo(original_image.dtype).max
    elif cfg.optim.image_file_name.endswith('.npy'):
        original_image = np.load(cfg.optim.image_file_name)
    else:
        print(f'{cfg.optim.image_file_name} format not recognized')
        exit(1)

    if cfg.optim.invert_gt:
        original_image = 1. - original_image

    if cfg.optim.image_gaussian_blur > 0:
        if cfg.optim.image_gaussian_blur % 2 == 0:
            raise "Gaussian blur must be an odd number"
        original_image = cv2.GaussianBlur(original_image, (cfg.optim.image_gaussian_blur, cfg.optim.image_gaussian_blur), 0)

    lpips_fn = lpips.LPIPS(net='vgg').to(device=torch.device(cfg.device))
    losses_fn = [
        ('L1', nn.L1Loss(), 0.4, False),
        ('SSIM', d_ssim_loss, 0.2, False),
        ('LPIPS', lpips_fn, 0.2, True),
    ]

    gaussians_models = []
    rasterized_img = None

    for i, res in enumerate(cfg.multiscale.resolutions):
        original_image_resized = cv2.resize(original_image, (res, res))
        cv2.imwrite(os.path.join(cfg.log.results_folder, f'gt_{res}.png'), (original_image_resized * 255).astype(np.uint8))

        init_img = original_image_resized

        if rasterized_img is not None:
            r = cv2.resize(rasterized_img.cpu().numpy(), (res, res))
            r = (r - np.min(r)) / (np.max(r) - np.min(r))
            init_img = init_img - r

        cv2.imwrite(os.path.join(cfg.log.results_folder, f'init_{res}.png'),
                    (((init_img-np.min(init_img))/(np.max(init_img)-np.min(init_img))) * 256).astype(np.uint8))

        init_img = torch.tensor(init_img, dtype=torch.float32, device=cfg.device)

        primary_samples = cfg.multiscale.samples_per_res[i] if cfg.multiscale.samples_per_res else cfg.gaussians.default_samples
        mod_scale_init = cfg.multiscale.mod_scale_per_res[i] if cfg.multiscale.mod_scale_per_res else 20/res

        gaussians = GaussianModel(primary_samples, init_img,
                                  device=cfg.device, lr=cfg.optim.learning_rate,
                                  mod_scale_init=mod_scale_init)
        gaussians_models.append(gaussians)

        with torch.no_grad():
            r = rasterize_cuda(gaussians.get_concat_params(apply_activation=True),
                               cfg.optim.image_size, return_zero_det_gaussians=False)
        if rasterized_img is not None:
            rasterized_img = rasterized_img + r
        else:
            rasterized_img = r
        raster_init = rasterized_img.cpu().numpy()
        cv2.imwrite(os.path.join(cfg.log.results_folder, f'gaussians_init_{res}.png'),
                    (((raster_init-np.min(raster_init))/(np.max(raster_init)-np.min(raster_init))) * 256).astype(
                        np.uint8))

    gaussians = torch.concat([g.get_concat_params(apply_activation=False) for g in gaussians_models], dim=0)
    gaussians = GaussianModel(init_from_array=gaussians.clone().detach(), device=cfg.device,
                              lr=cfg.optim.learning_rate)

    gaussians.add_densification_stats()
    gaussians.prune(min_amplitude=cfg.optim.amplitude_threshold,
                    max_sigma_ratio=cfg.optim.max_ratio_threshold,
                    image_size=cfg.optim.image_size)

    target_image = torch.tensor(cv2.resize(original_image, (cfg.optim.image_size, cfg.optim.image_size)),
                                dtype=torch.float32, device=cfg.device)

    start = time.time()
    gaussians, best_model = optimisation_loop(cfg, gaussians, losses_fn, target_image,
                                              cfg.log.results_folder, get_best_model=True)
    with open(os.path.join(cfg.log.results_folder, 'timing.txt'), 'w') as f:
        f.write(str(time.time() - start))
    gaussians.save_npy(cfg.log.results_folder, 'gaussians')
    best_model.save_npy(cfg.log.results_folder, 'best_model')


if __name__ == "__main__":
    main()
