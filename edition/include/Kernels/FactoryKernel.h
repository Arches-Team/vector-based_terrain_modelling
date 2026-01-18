#pragma once

#include "Kernel.h"
#include "GaussianKernel.h"
#include "DetailsKernel.h"


namespace kernel
{
    inline std::unique_ptr<Kernel> create(KernelType type, std::vector<float> params) {
        auto scaleX = params[0];
        auto scaleY = params[1];
        auto theta = params[2];
        auto amplitude = params[3];
        auto posX = params[4];
        auto posY = params[5];
        switch (type) {
        case KernelType::GAUSSIAN:
            return std::make_unique<GaussianKernel>(scaleX, scaleY, 1.f, theta, amplitude, posX, posY, params[6]);
        case KernelType::DETAILS:
            return std::make_unique<DetailsKernel>(scaleX, scaleY, theta, amplitude, posX, posY, params[6], params[7], params[8], params[9], params[10], params[11]);
        default:
            return nullptr;
        }
    }
}
