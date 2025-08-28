import os
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

import utils


class GaussianModel:
    def __init__(self, nb_gaussians: int = 0, image: torch.tensor = None, device: str = "cpu", lr: float = 0.01,
                 init_from_array: np.ndarray or torch.tensor = None, mod_scale_init: float = None, optimizable: bool = True):
        self._sigmas_activation = torch.sigmoid
        self._mod_scale_activation = F.softplus
        self._theta_activation = lambda x: torch.tanh(x) * utils.PI
        self._amplitude_activation = nn.Identity()
        self._xy_activation = torch.tanh
        self._beta_activation = lambda x: torch.sigmoid(x) + 1.

        self._inverse_xy_activation = torch.atanh
        self._inverse_amplitude_activation = nn.Identity()
        self._inverse_mod_scale_activation = lambda x: x + torch.log(-torch.expm1(-x))

        if image is not None:
            nb_gaussians = min(nb_gaussians, image.shape[0] * image.shape[1])
        self._sigmas = nn.Parameter(torch.rand(nb_gaussians, 2, device=device).requires_grad_(True))

        if mod_scale_init is None:
            self._mod_scale = nn.Parameter((torch.ones(nb_gaussians, device=device) * 0.2).requires_grad_(True))
        else:
            self._mod_scale = nn.Parameter((self._inverse_mod_scale_activation(torch.zeros(nb_gaussians, device=device)
                                            + mod_scale_init)).requires_grad_(True))
        self._theta = nn.Parameter(
            (2 * utils.PI * torch.rand(nb_gaussians, device=device) - utils.PI).requires_grad_(True))
        self._amplitude = nn.Parameter(torch.rand(nb_gaussians, device=device).requires_grad_(True) * 0.1)
        self._xy = nn.Parameter(self._inverse_xy_activation((2 * torch.rand(nb_gaussians, 2, device=device) - 1.)).requires_grad_(True))
        self._beta = nn.Parameter((torch.rand(nb_gaussians, device=device)).requires_grad_(True))

        self._xy_gradient_accum = torch.zeros((self._xy.shape[0]), device=device)
        self._denom = 0

        self._device = device

        # Overwrite all previous init
        if init_from_array is not None:
            self._init_from_array(init_from_array)

        self._init_param_accumulation()
        self._xy_gradient_accum = torch.zeros((self._xy.shape[0]), device=self._device)
        self._denom = 0

        params = [
            {'params': [self._sigmas], "name": "sigmas"},
            {'params': [self._mod_scale], "name": "mod_scale"},
            {'params': [self._theta], "name": "theta"},
            {'params': [self._amplitude], "name": "amplitude"},
            {'params': [self._xy], "name": "xy"},
            {'params': [self._beta], "name": "beta"}
        ]

        self.optimizer = torch.optim.Adam(params, lr=lr)

        if not optimizable:
            self.stop_optim(apply_activation=False)

    def _init_with_image(self, image, nb_gaussians):
        width, height = image.shape
        coords_int = np.random.randint(0, [width, height], size=(nb_gaussians, 2))
        self._xy = nn.Parameter(
            self._inverse_xy_activation(2. * torch.tensor(coords_int / [width, height],
                                                          device=self._device).float() - 1).requires_grad_(True))

        colors = [image[coord[1], coord[0]] * 0.2 for coord in coords_int]
        self._amplitude = nn.Parameter(
            self._inverse_amplitude_activation(torch.tensor(colors, device=self._device).float()).requires_grad_(True))

    def _init_from_array(self, array: np.ndarray or torch.tensor):
        self._sigmas = nn.Parameter(array[:, 0:2].clone().detach().to(self._device).requires_grad_(True).requires_grad_(True))
        self._mod_scale = nn.Parameter(array[:, 2].clone().detach().to(self._device).requires_grad_(True).requires_grad_(True))
        self._theta = nn.Parameter(array[:, 3].clone().detach().to(self._device).requires_grad_(True).requires_grad_(True))
        self._amplitude = nn.Parameter(array[:, 4].clone().detach().to(self._device).requires_grad_(True).requires_grad_(True))
        self._xy = nn.Parameter(array[:, 5:7].clone().detach().to(self._device).requires_grad_(True).requires_grad_(True))
        self._beta = nn.Parameter(array[:, 7].clone().detach().to(self._device).requires_grad_(True).requires_grad_(True))

    def _init_param_accumulation(self):
        self._amplitude_accum = torch.zeros(self._amplitude.shape, device=self._device)
        self._sigmas_accum = torch.zeros(self._sigmas.shape, device=self._device)
        self._mod_scale_accum = torch.zeros(self._mod_scale.shape, device=self._device)
        self._param_denom = 0

    @property
    def amplitude(self) -> torch.Tensor:
        return self._amplitude_activation(self._amplitude)

    @property
    def mod_scale(self) -> torch.Tensor:
        return self._mod_scale_activation(self._mod_scale)

    @property
    def sigmas(self) -> torch.Tensor:
        return self._sigmas_activation(self._sigmas)

    @property
    def xy(self) -> torch.Tensor:
        return self._xy_activation(self._xy)

    @property
    def theta(self) -> torch.Tensor:
        return self._theta_activation(self._theta)

    @property
    def beta(self) -> torch.Tensor:
        return self._beta_activation(self._beta)

    def clone(self):
        return GaussianModel(init_from_array=self.get_concat_params(apply_activation=False).clone().detach(),
                             device=self._device)

    def load_npy(self, npy_filepath, require_grad=False):
        npy_gaussians = np.load(npy_filepath)
        self._sigmas = nn.Parameter(
            torch.tensor(npy_gaussians[:, :2], device=self._device).requires_grad_(require_grad))
        self._mod_scale = nn.Parameter(
            torch.tensor(npy_gaussians[:, 2], device=self._device).requires_grad_(require_grad))
        self._theta = nn.Parameter(torch.tensor(npy_gaussians[:, 3], device=self._device).requires_grad_(require_grad))
        self._amplitude = nn.Parameter(
            torch.tensor(npy_gaussians[:, 4], device=self._device).requires_grad_(require_grad))
        self._xy = nn.Parameter(torch.tensor(npy_gaussians[:, 5:7], device=self._device).requires_grad_(require_grad))
        self._beta = nn.Parameter(torch.tensor(npy_gaussians[:, 7], device=self._device).requires_grad_(require_grad))

    def save_npy(self, folder, file):
        g = self.get_concat_params(apply_activation=True).detach().cpu().numpy()
        np.save(os.path.join(folder, f'{file}.npy'), g)

    def _get_params(self, apply_activation: bool = True):
        params = []
        if apply_activation:
            params.append(self.sigmas)
            params.append(self._mod_scale_activation(self._mod_scale))
            params.append(self._theta_activation(self._theta))
            params.append(self._amplitude_activation(self.amplitude))
            params.append(self._xy_activation(self._xy))
            params.append(self._beta_activation(self._beta))
        else:
            params.append(self._sigmas)
            params.append(self._mod_scale)
            params.append(self._theta)
            params.append(self._amplitude)
            params.append(self._xy)
            params.append(self._beta)

        return params

    def get_concat_params(self, apply_activation: bool = True):
        sigmas, mod_scale, theta, amplitude, xy, beta = self._get_params(apply_activation)
        mod_scale = mod_scale[..., None]
        theta = theta[..., None]
        amplitude = amplitude[..., None]
        beta = beta[..., None]

        return torch.cat((sigmas, mod_scale, theta, amplitude, xy, beta), dim=-1)

    def __len__(self):
        return self._sigmas.shape[0]

    def cat(self, gm: 'GaussianModel'):
        self.densification_postfix(gm._sigmas, gm._mod_scale, gm._theta, gm._amplitude, gm._xy, gm._beta)

    def translate(self, x: float, y: float, apply_activation: bool = True):
        p = torch.tensor([x, y], device=self._device)[None, ...]
        if apply_activation:
            self._xy = self._inverse_xy_activation(self.xy+p)
            self.prune_points(self._xy.isnan().any(dim=1))
        else:
            self._xy = self._xy+p

    def scale(self, scale: float, x: float, y: float, apply_activation: bool = True):
        self.translate(-x, -y, apply_activation)

        if apply_activation:
            self._xy = self._inverse_xy_activation(self.xy * scale)
            self._mod_scale = self._inverse_mod_scale_activation(self.mod_scale*scale)
        else:
            self._xy = self._xy * scale
            self._mod_scale = self._mod_scale * scale

        self.translate(x, y, apply_activation)

    # Set up the model after optimization: delete all activations and optimizer.
    # Scale and translate does not work without this due to propagation into the optimizer
    def stop_optim(self, apply_activation: bool = True):
        self.optimizer = None
        if apply_activation:
            self._sigmas = self.sigmas
            self._mod_scale = self.mod_scale
            self._theta = self.theta
            self._amplitude = self.amplitude
            self._xy = self.xy
            self._beta = self.beta

        self._sigmas_activation = nn.Identity()
        self._mod_scale_activation = nn.Identity()
        self._theta_activation = nn.Identity()
        self._amplitude_activation = nn.Identity()
        self._xy_activation = nn.Identity()
        self._beta_activation = nn.Identity()

        self._inverse_xy_activation = nn.Identity()
        self._inverse_amplitude_activation = nn.Identity()
        self._inverse_mod_scale_activation = nn.Identity()

    def get_inside_mask(self, x1: float, y1: float, x2: float, y2: float):
        x, y = self.xy.unbind(dim=1)
        mask_inside = torch.logical_and(torch.logical_and(x1 <= x, x <= x2),
                                        torch.logical_and(y1 <= y, y <= y2))
        return mask_inside

    # If the position is lower than step1, amplitude * 1., greater than step2, amplitude * 0., smooth step between
    def amplitude_blending(self, step1: float, step2: float, horizontal: bool):
        pos = self.xy[:, 1]
        if horizontal:
            pos = self.xy[:, 0]

        factor = utils.smoothstep(step2, step1, pos)
        self._amplitude = self._inverse_amplitude_activation(self.amplitude * factor)

    def prune_zero_amplitude(self):
        mask = torch.abs(self.amplitude) - utils.EPS < 0
        self.prune_points(mask)

    # Inspired by : Kerbl, B., Kopanas, G., Leimkuhler, T., & Drettakis, G. (2023).
    # 3D Gaussian Splatting for Real-Time Radiance Field Rendering. ACM Transactions on Graphics
    # -- START --
    def _prune_optimizer(self, mask):
        optimizable_tensors = {}
        for group in self.optimizer.param_groups:
            stored_state = self.optimizer.state.get(group['params'][0], None)
            if stored_state is not None:
                stored_state["exp_avg"] = stored_state["exp_avg"][mask]
                stored_state["exp_avg_sq"] = stored_state["exp_avg_sq"][mask]

                del self.optimizer.state[group['params'][0]]
                group["params"][0] = nn.Parameter((group["params"][0][mask].requires_grad_(True)))
                self.optimizer.state[group['params'][0]] = stored_state

                optimizable_tensors[group["name"]] = group["params"][0]
            else:
                group["params"][0] = nn.Parameter(group["params"][0][mask].requires_grad_(True))
                optimizable_tensors[group["name"]] = group["params"][0]
        return optimizable_tensors

    def prune_points(self, mask: torch.Tensor):
        valid_points_mask = ~mask
        if self.optimizer is not None:
            optimizable_tensors = self._prune_optimizer(valid_points_mask)

            self._sigmas = optimizable_tensors["sigmas"]
            self._mod_scale = optimizable_tensors["mod_scale"]
            self._theta = optimizable_tensors["theta"]
            self._amplitude = optimizable_tensors["amplitude"]
            self._xy = optimizable_tensors["xy"]
            self._beta = optimizable_tensors["beta"]

            self._xy_gradient_accum = self._xy_gradient_accum[valid_points_mask]
            self._amplitude_accum = self._amplitude_accum[valid_points_mask]
            self._mod_scale_accum = self._mod_scale_accum[valid_points_mask]
            self._sigmas_accum = self._sigmas_accum[valid_points_mask]
        else:
            self._sigmas = self._sigmas[valid_points_mask]
            self._mod_scale = self._mod_scale[valid_points_mask]
            self._theta = self._theta[valid_points_mask]
            self._amplitude = self._amplitude[valid_points_mask]
            self._xy = self._xy[valid_points_mask]
            self._beta = self._beta[valid_points_mask]

    def cat_tensors_to_optimizer(self, tensors_dict):
        optimizable_tensors = {}
        for group in self.optimizer.param_groups:
            assert len(group["params"]) == 1
            extension_tensor = tensors_dict[group["name"]]
            stored_state = self.optimizer.state.get(group['params'][0], None)
            if stored_state is not None:

                stored_state["exp_avg"] = torch.cat((stored_state["exp_avg"], torch.zeros_like(extension_tensor)),
                                                    dim=0)
                stored_state["exp_avg_sq"] = torch.cat((stored_state["exp_avg_sq"], torch.zeros_like(extension_tensor)),
                                                       dim=0)

                del self.optimizer.state[group['params'][0]]
                group["params"][0] = nn.Parameter(
                    torch.cat((group["params"][0], extension_tensor), dim=0).requires_grad_(True))
                self.optimizer.state[group['params'][0]] = stored_state

                optimizable_tensors[group["name"]] = group["params"][0]
            else:
                group["params"][0] = nn.Parameter(
                    torch.cat((group["params"][0], extension_tensor), dim=0).requires_grad_(True))
                optimizable_tensors[group["name"]] = group["params"][0]

        return optimizable_tensors

    def densification_postfix(self, new_sigmas, new_mod_scale, new_theta, new_amplitude, new_xy, new_beta):
        if self.optimizer is not None:
            d = {"sigmas": new_sigmas,
                 "mod_scale": new_mod_scale,
                 "theta": new_theta,
                 "amplitude": new_amplitude,
                 "xy": new_xy,
                 "beta": new_beta,
                 }

            optimizable_tensors = self.cat_tensors_to_optimizer(d)
            self._sigmas = optimizable_tensors["sigmas"]
            self._mod_scale = optimizable_tensors["mod_scale"]
            self._theta = optimizable_tensors["theta"]
            self._amplitude = optimizable_tensors["amplitude"]
            self._xy = optimizable_tensors["xy"]
            self._beta = optimizable_tensors["beta"]

            self._xy_gradient_accum = torch.zeros((self._xy.shape[0]), device=self._device)
            self._denom = 0
            self._init_param_accumulation()
        else:
            self._sigmas = torch.cat((self._sigmas, new_sigmas), dim=0)
            self._mod_scale = torch.cat((self._mod_scale, new_mod_scale), dim=0)
            self._theta = torch.cat((self._theta, new_theta), dim=0)
            self._amplitude = torch.cat((self._amplitude, new_amplitude), dim=0)
            self._xy = torch.cat((self._xy, new_xy), dim=0)
            self._beta = torch.cat((self._beta, new_beta), dim=0)

    def densify_and_split(self, grads, grad_threshold: float, N: int = 2):
        n_init_points = self.xy.shape[0]

        # Extract points that satisfy the gradient condition
        padded_grad = torch.zeros(n_init_points, device=self._device)
        padded_grad[:grads.shape[0]] = grads.squeeze()
        selected_pts_mask = padded_grad >= grad_threshold

        # Sample from 0 means gaussians distribution to apply rotation
        stds = self.sigmas[selected_pts_mask].repeat(N, 1)
        means = torch.zeros((stds.shape[0], 2), device=self._device)
        samples = torch.normal(mean=means, std=stds)
        rots = utils.build_rotation(self.theta[selected_pts_mask]).repeat(N, 1, 1)
        new_xy = self._inverse_xy_activation(torch.clamp(torch.bmm(rots, samples.unsqueeze(-1)).squeeze(-1)
                                         + self.xy[selected_pts_mask].repeat(N, 1),
                             -1 + 1e-8, 1 - 1e-8))
        new_sigmas = self._sigmas[selected_pts_mask].repeat(N, 1) / (0.8 * N)
        new_theta = self._theta[selected_pts_mask].repeat(N)
        new_mod_scale = self._mod_scale[selected_pts_mask].repeat(N)
        new_amplitude = self._amplitude[selected_pts_mask].repeat(N)
        new_beta = self._beta[selected_pts_mask].repeat(N)

        self.densification_postfix(new_sigmas, new_mod_scale, new_theta, new_amplitude, new_xy, new_beta)

        prune_filter = torch.cat(
            (selected_pts_mask, torch.zeros(N * selected_pts_mask.sum(), device=self._device, dtype=bool)))
        self.prune_points(prune_filter)

    def densify_and_clone(self, grads, grad_threshold: float):
        # Extract points that satisfy the gradient condition
        selected_pts_mask = grads >= grad_threshold

        new_sigmas = self._sigmas[selected_pts_mask]
        new_mod_scale = self._mod_scale[selected_pts_mask]
        new_theta = self._theta[selected_pts_mask]
        new_amplitude = self._amplitude[selected_pts_mask]
        new_xy = self._xy[selected_pts_mask]
        new_beta = self._beta[selected_pts_mask]

        self.densification_postfix(new_sigmas, new_mod_scale, new_theta, new_amplitude, new_xy, new_beta)

    def add_densification_stats(self):
        if self._xy.grad is not None:
            self._xy_gradient_accum += torch.linalg.vector_norm(self._xy.grad, dim=-1, keepdim=True).squeeze()
            self._denom += 1

        self._amplitude_accum += self._amplitude
        self._sigmas_accum += self._sigmas
        self._mod_scale_accum += self._mod_scale
        self._param_denom += 1
    # -- END -- Kerbl, B., Kopanas, G., Leimkuhler, T., & Drettakis, G. (2023).
    # 3D Gaussian Splatting for Real-Time Radiance Field Rendering. ACM Transactions on Graphics

    def densify(self, max_grad: float):
        grads = self._xy_gradient_accum / self._denom
        grads[grads.isnan()] = 0.0

        # TODO: add a maximum number ?
        self.densify_and_clone(grads, max_grad)
        self.densify_and_split(grads, max_grad)

        torch.cuda.empty_cache()

    def prune(self, min_amplitude: float, max_sigma_ratio: float, image_size: float):
        if self._param_denom == 0:
            return

        # TODO: sigmas norm instead of sigma ratio ?
        min_sigma_ratio = 1 / max_sigma_ratio
        sigmas = self._sigmas_activation(self._sigmas_accum / self._param_denom)
        mod_scale = self._mod_scale_activation(self._mod_scale_accum / self._param_denom)
        amplitude = self._amplitude_activation(self._amplitude_accum / self._param_denom)

        ratio = sigmas[:, 0] / sigmas[:, 1]

        # TODO: 2 because between -1 and 1. todo: remove hard coded value
        pixel_size = 2 / image_size
        size_sigmas = torch.max(3. * sigmas * mod_scale[..., None], dim=1).values

        prune_mask = torch.logical_or(torch.logical_or(torch.abs(amplitude) < min_amplitude, size_sigmas < pixel_size),
                                      torch.logical_or(ratio > max_sigma_ratio, ratio < min_sigma_ratio))
        self.prune_points(prune_mask)
        self._init_param_accumulation()

        torch.cuda.empty_cache()
