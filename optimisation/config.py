from dataclasses import dataclass, field
from typing import Optional, List

import utils


@dataclass
class GaussiansConfig:
    """ Config for gaussians model arguments """
    default_samples: int = 1000


@dataclass
class LogConfig:
    """ Config for logging arguments """
    results_folder: str = 'out'
    display_interval: int = 50
    save_gaussians_at_display: bool = False


@dataclass
class OptimizationConfig:
    """ Config for optimization arguments """
    num_epochs: int = 101
    densification_interval: int = 50
    densify_before: float = 0.8
    learning_rate: float = 0.01
    image_file_name: str = "../input.png"
    sketch_file_name: Optional[str] = None
    image_size: int = 256
    image_gaussian_blur: int = 0
    invert_gt: bool = False

    # Threshold of early stopping
    early_stopping: float = 0
    early_stopping_epoch_not_improve: int = None
    # Run for x epoch after early stopping without cloning/splitting
    early_stopping_patience: int = 0

    # Regularization threshold
    gradient_threshold: float = 0.
    amplitude_threshold: float = 0.
    max_ratio_threshold: float = 1e6


@dataclass
class MultiscaleConfig:
    """ Config for joined multiscale arguments """
    resolutions: List[int] = field(default_factory=list)
    samples_per_res: Optional[List[int]] = None
    mod_scale_per_res: Optional[List[float]] = None


@dataclass
class Config:
    log: LogConfig = field(default_factory=LogConfig)
    gaussians: GaussiansConfig = field(default_factory=GaussiansConfig)
    optim: OptimizationConfig = field(default_factory=OptimizationConfig)
    device: str = utils.get_device()
    multiscale: MultiscaleConfig = field(default_factory=MultiscaleConfig)
