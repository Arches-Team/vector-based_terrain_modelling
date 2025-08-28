import os
from dataclasses import dataclass

import numpy as np
import cv2
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import lpips

import losses

lpips_fn = lpips.LPIPS(net='vgg')

@dataclass
class Config:
    image_size: int = 1024
    folder: str = './results/xxx/'
    dem: str = 'final.png'
    gt: str = None
    gaussians: torch.tensor = None
    recursive: bool = True


def main():
    """Generate metrics for one optimisation result.

    Produces losses (if GT is given), heatmap, sigmas means.
    Change Config class accordingly to your data.
    """
    list_metrics = [gaussians_pos_and_heatmap, statistics]

    config = Config()

    config.gt = r'./data/xxx.png'

    if not config.recursive:
        launch_metrics(config, list_metrics)
    else:
        main_folder = config.folder
        for folder, dirs, files in os.walk(main_folder):
            if any([f == 'gaussians.npy' for f in files]):
                config.folder = folder
                print(config.folder)
                launch_metrics(config, list_metrics)


def launch_metrics(config: Config, metrics):
    gaussians = np.load(f'{config.folder}/gaussians.npy')
    config.gaussians = gaussians
    for metric in metrics:
        metric(config)


def gaussians_pos_and_heatmap(config: Config):
    image_size = config.image_size
    dem = np.zeros((config.image_size, config.image_size), dtype=np.float32)
    if os.path.isfile(os.path.join(config.folder, config.dem)):
        dem = ((cv2.imread(os.path.join(config.folder, config.dem), cv2.IMREAD_UNCHANGED) / 65535) * 255).astype(
            np.uint8)
        image_size = dem.shape[0]

    x = (((np.clip(config.gaussians[:, 5], -1, 1) + 1) / 2) * (image_size-1)).astype(np.uint16)
    y = (((np.clip(config.gaussians[:, 6], -1, 1) + 1) / 2) * (image_size-1)).astype(np.uint16)
    amplitude = config.gaussians[:, 4]
    colors = np.repeat(np.asarray(amplitude)[..., None], 3, axis=-1)
    colors[amplitude >= 0.] = (0, 0, 255)
    colors[amplitude < 0.] = (255, 0, 0)

    img = np.stack((dem,)*3, axis=-1)
    heatmap = np.zeros((image_size, image_size), dtype=np.uint8)
    img[y, x] = colors
    heatmap[y, x] = 1.
    kernel_size = 21
    kernel = np.ones((kernel_size, kernel_size))
    heatmap = cv2.filter2D(heatmap, ddepth=-1, kernel=kernel)
    heatmap = (((heatmap - np.min(heatmap)) / (np.max(heatmap) - np.min(heatmap))) * 255).astype(np.uint8)
    heatmap = cv2.applyColorMap(heatmap, cv2.COLORMAP_JET)

    cv2.imwrite(f'{config.folder}/gaussian_{config.gaussians.shape[0]}_pos.png', img)
    cv2.imwrite(f'{config.folder}/gaussian_distribution_heatmap.png', heatmap)


def statistics(config: Config):
    txt_file = os.path.join(config.folder, 'statistics.txt')

    gaussians = config.gaussians
    sigma_x = gaussians[:, 0]
    sigma_y = gaussians[:, 1]

    H, xedges, yedges = np.histogram2d(sigma_x, sigma_y)
    H = H.T
    plt.imshow(H, interpolation='nearest', origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    plt.savefig(os.path.join(config.folder, 'histogram_sigma.png'))

    mean_sigma_x = np.mean(sigma_x)
    mean_sigma_y = np.mean(sigma_y)

    with open(txt_file, 'w') as f:
        f.write(f'Name: {config.folder}\n')
        f.write(f'Nb primitives: {gaussians.shape[0]}\n')
        f.write(f'Mean sigma_x: {mean_sigma_x}\n')
        f.write(f'Mean sigma_y: {mean_sigma_y}\n')
        with open(os.path.join(config.folder, 'timing.txt')) as timing:
            t = timing.readlines()[0]
            f.write(f'Time (s): {t}\n')

    # Computing errors
    is_dem_exist = os.path.isfile(os.path.join(config.folder, config.dem))

    if config.gt is not None and is_dem_exist:
        dem = torch.tensor(cv2.imread(os.path.join(config.folder, config.dem), cv2.IMREAD_UNCHANGED) / 65535.)[
            ..., None]
        image_size = dem.shape[0]

        gt = cv2.imread(config.gt, cv2.IMREAD_UNCHANGED) / 65535.
        gt = cv2.resize(gt, (image_size, image_size))
        gt = torch.tensor(gt)[..., None]

        gt_rgb = gt.permute(2, 0, 1)[None, ...].repeat(1, 3, 1, 1).float()
        dem_rgb = dem.permute(2, 0, 1)[None, ...].repeat(1, 3, 1, 1).float()

        l1 = torch.nn.L1Loss()(gt, dem)
        l2 = torch.nn.MSELoss()(gt, dem)
        ssim = losses.d_ssim_loss(gt, dem)
        lpips_loss = lpips_fn(gt_rgb, dem_rgb)
        with open(txt_file, 'a') as f:
            f.write(f'L1: {float(l1):0.8f}\n')
            f.write(f'L2: {float(l2):0.8f}\n')
            f.write(f'SSIM: {float(ssim):0.8f}\n')
            f.write(f'LPIPS: {float(lpips_loss):0.8f}\n')


if __name__ == "__main__":
    main()
