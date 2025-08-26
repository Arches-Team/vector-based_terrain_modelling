import cv2
import numpy as np
import torch
import os

from gaussian_model import GaussianModel
from rasterizer import rasterize_cuda
import utils


def infer(npy_file, image_size, device):
    gaussians = GaussianModel(device=device)
    gaussians.load_npy(npy_file)

    with torch.no_grad():
        img = rasterize_cuda(gaussians.get_concat_params(apply_activation=False), image_size)

    img = (torch.clamp(img, 0., 1)/1.).cpu().detach().numpy()
    return img


def main():
    image_size = 1024
    folder = './results/gsplat/dem/SRTM/test_init/N38W003_0002_0000_start10k_display1_lr1e-3_noblur_initrandom/'
    device = utils.get_device()

    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.npy'):
                file_path = os.path.join(root, file)
                img = infer(file_path, image_size, device)
                cv2.imwrite(f'{root}/dem_final_{file[:-4].zfill(4)}.png', (img * 65535).astype(np.uint16))
                print(f'Saved {root}/final.png')


if __name__ == "__main__":
    main()
