import math
import os
from typing import List
import shutil
from collections import defaultdict

import pyrallis
import torch
from tqdm import tqdm
import matplotlib

import utils

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cv2
import numpy as np

from losses import combined_loss
from rasterizer import rasterize_cuda
from config import Config
from gaussian_model import GaussianModel


def optimisation_loop(cfg: Config, gaussians: GaussianModel, losses_fn: List,
                      target_image: torch.Tensor, out_folder: str,
                      img_to_add: torch.Tensor = None, get_best_model: bool = False):
    loss_dict_history = defaultdict(list)
    nb_gaussians_history = []
    densify_until_iter = int(cfg.optim.densify_before * cfg.optim.num_epochs)

    epoch_loss_not_improve = 0
    loss_min = 0
    last_epoch = cfg.optim.num_epochs
    early_stopping_patience = False
    best_model = None
    best_image = None
    best_model_epoch = None

    for epoch in tqdm(range(cfg.optim.num_epochs)):
        rasterized_img, zero_det_gaussians = rasterize_cuda(gaussians.get_concat_params(apply_activation=True),
                                                            cfg.optim.image_size, return_zero_det_gaussians=True)

        if img_to_add is not None:
            rasterized_img = rasterized_img + img_to_add

        loss, loss_dict = combined_loss(rasterized_img, target_image, losses_fn, return_dict=True)
        gaussians.optimizer.zero_grad()

        loss.backward()
        gaussians.optimizer.step()

        if epoch > 0:
            if loss.item() > loss_min + cfg.optim.early_stopping:
                epoch_loss_not_improve += 1
            else:
                epoch_loss_not_improve = 0
                best_model = gaussians.clone()
                best_image = rasterized_img
                best_model_epoch = epoch
        else:
            loss_min = loss.item()

        loss_min = min(loss_min, loss.item())
        loss_dict_history['main'].append(loss.item())
        for l in loss_dict:
            loss_dict_history[l].append(loss_dict[l])
        nb_gaussians_history.append(len(gaussians))

        with torch.no_grad():
            gaussians.add_densification_stats()
            gaussians.prune_points(zero_det_gaussians)

            if not early_stopping_patience:
                if epoch % cfg.optim.densification_interval == 0:
                    gaussians.prune(min_amplitude=cfg.optim.amplitude_threshold,
                                    max_sigma_ratio=cfg.optim.max_ratio_threshold,
                                    image_size=cfg.optim.image_size)

                if epoch % cfg.optim.densification_interval == 0 and 0 < epoch < densify_until_iter:
                    gaussians.densify(max_grad=cfg.optim.gradient_threshold)

            if epoch % cfg.log.display_interval == 0:
                save_plot(cfg, out_folder, rasterized_img, target_image, len(gaussians), loss_dict, loss_dict_history,
                          nb_gaussians_history, epoch, cfg.optim.num_epochs, densify_until_iter)
                if cfg.log.save_gaussians_at_display:
                    gaussians.save_npy(cfg.log.results_folder, f'{epoch}')

        if early_stopping_patience:
            cfg.optim.early_stopping_patience -= 1
            if cfg.optim.early_stopping_patience < 1:
                print(f'Stopped after patience.')
                last_epoch = epoch
                break

        if cfg.optim.early_stopping_epoch_not_improve is not None and epoch_loss_not_improve > cfg.optim.early_stopping_epoch_not_improve:
            print(f'Early stopping! {epoch_loss_not_improve} epochs with no improvement.')
            if cfg.optim.early_stopping_patience > 0:
                early_stopping_patience = True
            else:
                last_epoch = epoch
                break

    save_plot(cfg, out_folder, rasterized_img, target_image, len(gaussians), loss_dict, loss_dict_history,
              nb_gaussians_history, last_epoch - 1, cfg.optim.num_epochs, densify_until_iter)
    if get_best_model:
        generated_array = np.clip(best_image.cpu().detach().numpy(), 0., 1.)
        img = (generated_array * 255).astype(np.uint8)
        file_path_png = os.path.join(out_folder, f'best_model_epoch{best_model_epoch}.png')
        cv2.imwrite(file_path_png, img)
        return gaussians, best_model
    return gaussians


def save_plot(cfg: Config, directory, img, target, nb_gaussians, loss_dict, loss_dict_history, nb_gaussians_history,
              epoch, num_epochs, densify_until_iter):
    num_subplots = 2
    fig_size_width = 12

    fig, ax = plt.subplots(2, num_subplots, figsize=(fig_size_width, fig_size_width))

    generated_array = img.cpu().detach().numpy()

    ax[0, 0].imshow(img.cpu().detach().numpy())
    ax[0, 0].set_title('2D Gaussian Splatting')
    ax[0, 0].axis('off')

    ax[0, 1].imshow(target.cpu().detach().numpy())
    ax[0, 1].set_title('Ground Truth')
    ax[0, 1].axis('off')

    ax[1, 0].plot(range(epoch + 1), loss_dict_history['main'][:epoch + 1])
    ax[1, 0].set_title('Loss vs. Epochs')
    ax[1, 0].set_xlabel('Epoch')
    ax[1, 0].set_ylabel('Loss')
    ax[1, 0].set_xlim(0, num_epochs)

    ax[1, 1].plot(range(epoch + 1), nb_gaussians_history[:epoch + 1])
    ax[1, 1].set_title('Number of gaussians vs. Epochs')
    ax[1, 1].set_xlabel('Epoch')
    ax[1, 1].set_ylabel('Nb gaussians')
    ax[1, 1].set_xlim(0, num_epochs)
    ax[1, 1].axvline(
        math.floor(densify_until_iter / cfg.optim.densification_interval) * cfg.optim.densification_interval,
        linestyle="dashed", color='red')

    # Display the image
    plt.subplots_adjust(wspace=0.1)

    with open(os.path.join(directory, 'losses.txt'), 'w') as f:
        for l in loss_dict_history:
            f.write(f'{l}:{str(loss_dict_history[l])}\n')

    with open(os.path.join(directory, 'nb_gaussians_history.txt'), 'w') as f:
        f.write(str(nb_gaussians_history))

    generated_array = np.clip(generated_array, 0., 1.)
    img = (generated_array * 255).astype(np.uint8)

    filename_jpg = f"{epoch}.jpg"
    filename_png = f"dem_{epoch}.png"
    file_path_png = os.path.join(directory, filename_png)
    file_path_jpg = os.path.join(directory, filename_jpg)

    cv2.imwrite(file_path_png, img)
    fig.savefig(file_path_jpg, bbox_inches='tight')

    plt.clf()
    plt.close()

    print(
        f"Epoch {epoch + 1}/{num_epochs}, Loss: {loss_dict_history['main'][-1]}, on {nb_gaussians} points.\nLosses: {loss_dict}")


def backup_code(out_folder: str, files: List[str], cfg: Config):
    current_folder = os.path.dirname(__file__)
    backup_folder = os.path.join(out_folder, 'code')
    os.makedirs(backup_folder, exist_ok=True)
    for file in files:
        backup_file = os.path.join(backup_folder, file)
        recursive_folder = os.path.dirname(backup_file)
        os.makedirs(recursive_folder, exist_ok=True)
        shutil.copy(os.path.join(current_folder, file), backup_file)
    with open(os.path.join(backup_folder, 'config.yaml'), 'w') as f:
        f.write(pyrallis.dump(cfg))
