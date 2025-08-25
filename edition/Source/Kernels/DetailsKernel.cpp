#include "Kernels/DetailsKernel.h"
#include "smooth.h"


//Kernel Kernel::operator+(const Kernel& k) const
//{
//    return {
//        m_scaleX + k.m_scaleX, m_scaleY + k.m_scaleY, m_modScale + k.m_modScale, m_theta + k.m_theta,
//        m_amplitude * m_modAmplitude + k.m_amplitude * k.m_modAmplitude, m_posX + k.m_posX, m_posY + k.m_posY,
//        m_beta + k.m_beta
//    };
//}
//
//Kernel Kernel::operator/(const float val) const
//{
//    return {
//        m_scaleX / val, m_scaleY / val, m_modScale / val, m_theta / val, (m_amplitude * m_modAmplitude) / val,
//        m_posX / val, m_posY / val, m_beta / val
//    };
//}
//
//Kernel Kernel::operator*(const float val) const
//{
//    return {
//        m_scaleX * val, m_scaleY * val, m_modScale * val, m_theta * val, (m_amplitude * m_modAmplitude) * val,
//        m_posX * val, m_posY * val, m_beta * val
//    };
//}

void DetailsKernel::print(std::ostream& os) const
{
	os << "DetailsKernel(Scale: (" << m_scaleX << ", " << m_scaleY << "), Pos: (" << m_posX << ", " << m_posY <<
		"), ";
	os << ", Amplitude: " << m_amplitude << ", Theta: " << m_theta << ", Texture: (" << m_textureX << ", " << m_textureY << "), ";
	os << " Original(Scale: (" << m_scaleXOrig << ", " << m_scaleYOrig << "), Theta: " << m_thetaOrig << ")), W/O deformation theta: (" << m_thetaWODeform << ")";
}

