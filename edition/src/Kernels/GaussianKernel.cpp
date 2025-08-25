#include "Kernels/GaussianKernel.h"


void GaussianKernel::print(std::ostream& os) const
{
	os << "GaussianKernel(Scale: (" << m_scaleX << ", " << m_scaleY << "), Pos: (" << m_posX << ", " << m_posY <<
		"), ";
	os << ", Amplitude: " << m_amplitude << ", Theta: " << m_theta << ", Beta: " << m_beta << ")";
}
