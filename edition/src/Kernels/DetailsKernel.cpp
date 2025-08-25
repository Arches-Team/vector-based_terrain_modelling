#include "Kernels/DetailsKernel.h"


void DetailsKernel::print(std::ostream& os) const
{
	os << "DetailsKernel(Scale: (" << m_scaleX << ", " << m_scaleY << "), Pos: (" << m_posX << ", " << m_posY <<
		"), ";
	os << ", Amplitude: " << m_amplitude << ", Theta: " << m_theta << ", Texture: (" << m_textureX << ", " << m_textureY << "), ";
	os << " Original(Scale: (" << m_scaleXOrig << ", " << m_scaleYOrig << "), Theta: " << m_thetaOrig << ")), W/O deformation theta: (" << m_thetaWODeform << ")";
}

