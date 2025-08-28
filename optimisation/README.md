# Vector-based terrain modelling - Optimisation

## Installation

Install the required packages from the requirements : 
```
pip install -r requirements.txt
```
Please note that due to CUDA code, you will need a NVIDIA GPU with `nvcc` installed.

## Usage

The configuration of an optimization is handled by a `.yml` file. A default one is given as an example in the `configs/` folder.

After setting the parameters in the `.yml` file, launch the optimisation :
```
python optimise.py --config_path=./configs/default.yml
```

Results will be saved in the folder defined by `log.results_folder`.

## Notes

The `gsplat/` folder is adapted from ![3D Gaussian Splatting for Real-Time Radiance Field Rendering](https://github.com/graphdeco-inria/gaussian-splatting) official repository.