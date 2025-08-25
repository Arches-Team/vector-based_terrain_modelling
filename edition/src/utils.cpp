#pragma once
#include "utils.h"

Vector2 utils::toVec2(QPointF p)
{
    return Vector2(p.x(), p.y());
}
