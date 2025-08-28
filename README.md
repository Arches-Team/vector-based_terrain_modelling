# Vector-based terrain modelling

![](docs/teaser.jpg)

[Website](https://simonperche.github.io/vector-based_terrain_modelling/) | [Paper](https://hal.science/hal-05104719v1/document) | [Video](https://www.youtube.com/watch?v=EOYhb4Z3j_E) 

This repository contains the official code and resources for the paper *Vector-based terrain modelling*.

## Abstract

Vector-based graphics offer numerous advantages over grid-based models, including resolution independence and ease of manipulation.
Despite these benefits, their use in landscape modeling remains uncommon because of a lack of direct editing and interactive feedback, essential for matching the artist's vision. 
We introduce a new vector-based model for creating digital terrains based on computationally efficient primitives. We propose a method to convert grid-based digital elevation maps to this representation with a user-defined level of accuracy. 
Once vectorized, the terrain can be authored using interactive high-level skeleton-based tools adapted to the primitive representation, allowing local deformations that automatically adapt to underlying geomorphological structures and landforms of the terrain.

## Structure

This work is divided onto two subfolders : 
- `/optimisation` containing the vectorization code (Python/PyTorch) allowing to convert a heightfield into a vector format as described in our paper
- `/edition` containing the edition app (C++/Qt) implementing the authoring of the vector terrain including geomorphology-based tool and landforms tools 

Each of these have their own README for the installation and usage.

## Citation
If you find our work useful, please consider citing:

```
@article{perche2025vectorterrains,
author  = {Perche, Simon and Guérin, Eric and Galin, Eric and Peytavie, Adrien},
title   = {Vector-Based Terrain Modelling},
journal = {Computer Graphics Forum},
year    = {2025},
}
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
