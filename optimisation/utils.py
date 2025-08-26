import time
import os

import numpy as np
import torch

PI = 3.1415927410125732
EPS = 1e-8


def get_device() -> str:
    if not torch.cuda.is_available():
        raise RuntimeError('CUDA not available.')
    return 'cuda'


def build_rotation(theta: torch.Tensor):
    R = torch.stack(
        [torch.stack([torch.cos(theta), -torch.sin(theta)], dim=-1),
         torch.stack([torch.sin(theta), torch.cos(theta)], dim=-1)],
        dim=-2
    )
    return R


# https://stackoverflow.com/questions/28889210/smoothstep-function
def smoothstep(edge0: float, edge1: float, x: torch.tensor):
    smooth = torch.clamp((x - edge0) / (edge1 - edge0), 0., 1.)
    return smooth * smooth * (3 - 2 * smooth)
