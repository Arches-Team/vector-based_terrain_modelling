#include "Kernels/GaussianKernel.h"
#include "smooth.h"

//Kernel GaussianKernel::operator+(const Kernel& k) const
//{
//    return {
//        m_scaleX + k.m_scaleX, m_scaleY + k.m_scaleY, m_modScale + k.m_modScale, m_theta + k.m_theta,
//        m_amplitude * m_modAmplitude + k.m_amplitude * k.m_modAmplitude, m_posX + k.m_posX, m_posY + k.m_posY,
//        m_beta + k.m_beta
//    };
//}
//
//Kernel GaussianKernel::operator/(const float val) const
//{
//    return {
//        m_scaleX / val, m_scaleY / val, m_modScale / val, m_theta / val, (m_amplitude * m_modAmplitude) / val,
//        m_posX / val, m_posY / val, m_beta / val
//    };
//}
//
//Kernel GaussianKernel::operator*(const float val) const
//{
//    return {
//        m_scaleX * val, m_scaleY * val, m_modScale * val, m_theta * val, (m_amplitude * m_modAmplitude) * val,
//        m_posX * val, m_posY * val, m_beta * val
//    };
//}

void GaussianKernel::print(std::ostream& os) const
{
	os << "GaussianKernel(Scale: (" << m_scaleX << ", " << m_scaleY << "), Pos: (" << m_posX << ", " << m_posY <<
		"), ";
	os << ", Amplitude: " << m_amplitude << ", Theta: " << m_theta << ", Beta: " << m_beta << ")";
}
