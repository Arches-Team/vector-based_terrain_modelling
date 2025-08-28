#include "Kernels/Kernel.h"

#include <math.h>
#include <tuple>

#include "libs/smooth.h"
#include "Eigen/Dense"


double Kernel::distCenter(const Vector2& p) const
{
    return Norm(Vector2(m_posX, m_posY) - p);
}

Ellipse2 Kernel::getEllipse(const float sigma) const
{
    Vector2 majorAxis(1., 0.);

    const Vector2 center(m_posX, m_posY);
    const double angle = 3*M_PI / 2 - m_theta;
    const double cs = cos(angle);
    const double sn = sin(angle);
    Vector2 majorAxisRotate;
    majorAxisRotate[0] = majorAxis[0] * cs - majorAxis[1] * sn;
    majorAxisRotate[1] = majorAxis[0] * sn + majorAxis[1] * cs;
    return Ellipse2(center, m_scaleX * sigma, m_scaleY * sigma, majorAxisRotate);
}

bool Kernel::isInside(const Vector2& p) const
{
    const SmoothEllipse2 ellipse(getEllipse());
    return ellipse.Value(p) > 0.;
}

/*!
\brief Scale the kernel with o origin
*/
void Kernel::scale(const float scale, const Vector2& o)
{
    m_scaleX *= scale;
    m_scaleY *= scale;
    translate(-o);
    m_posX *= scale;
    m_posY *= scale;
    translate(o);
    m_amplitude *= scale;
}

void Kernel::scale(Vector2 axis, float factor)
{
    // Double has not enough precision
    float multiplier = 100.f;

    m_scaleX *= multiplier;
    m_scaleY *= multiplier;

    // Scaling needs an angle between [0 ; pi] and 'a' defined major axis
    axis = Normalized(axis);

    auto angleVector = [](Vector2 x, Eigen::Vector2f y) -> double
        {
            return atan2(y[1] * x[0] - y[0] * x[1], y[0] * x[0] + y[1] * x[1]);
        };

    auto x = [](double a, double b, double theta, double t) -> double
        {
            return a * cos(t) * cos(theta) - b * sin(t) * sin(theta);
        };
    auto y = [](double a, double b, double theta, double t) -> double
        {
            return a * cos(t) * sin(theta) + b * sin(t) * cos(theta);
        };

    auto ellipse_parameters = [](double A, double B, double C, double D, double E, double F) -> std::tuple<double, double, double>
        {
            double theta = (utils::pi / 2)+ 0.5 * atan2(-B, C - A);
            double det = B * B - 4 * A * C;

            if (det >= 0)
                std::cout << "Oops, det is not < 0, these points do not define an ellipse." << std::endl;

            double sqrt_part = sqrt((A - C) * (A - C) + B * B);

            double a = -sqrt(2 * ((det * F) * ((A + C) - sqrt_part))) / det;
            double b = -sqrt(2 * ((det * F) * ((A + C) + sqrt_part))) / det;

            return std::make_tuple(a, b, theta);
        };


    Eigen::Vector2f xAxis(0, 1);
    double angle = angleVector(axis, xAxis);
    m_theta -= static_cast<float>(angle);

    Eigen::MatrixXf A(5, 6);

    constexpr int n = 5;
    for (int i = 0; i < n; i++)
    {
        double t = i * sqrt(2) * utils::pi;
        double current_x = x(m_scaleX, m_scaleY, m_theta, t);
        double current_y = y(m_scaleX, m_scaleY, m_theta, t);

        current_x *= factor;
        A(i, 0) = current_x * current_x;
        A(i, 1) = current_x * current_y;
        A(i, 2) = current_y * current_y;
        A(i, 3) = current_x;
        A(i, 4) = current_y;
        A(i, 5) = 1.0;
    }

    Eigen::MatrixXf ATA = A.transpose() * A;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(ATA);

    // Extract the eigenvector associated with the minimum eigenvalue
    auto eigenvalues = solver.eigenvalues();
    double minEigenValue = std::numeric_limits<double>::max();
    int minEigenValueId = -1;
    for (int i = 0; i < eigenvalues.size(); ++i)
    {
        double absval = abs(eigenvalues[i]);
        if (absval < minEigenValue)
        {
            minEigenValue = absval;
            minEigenValueId = i;
        }
    }

    if (minEigenValueId == -1)
    {
        throw std::logic_error("No eigenvalues found. Please check the scaling of primitives.");
    }
        
    Eigen::VectorXf z = solver.eigenvectors().col(minEigenValueId);

    // Calculate ellipse parameters
    auto [a, b, theta] = ellipse_parameters(z(0), z(1), z(2), z(3), z(4), z(5));

    a = roundf(a * 10000) / 10000;
    b = roundf(b * 10000) / 10000;
    theta = roundf(theta * 10000) / 10000;

    // Check for NaN
    if (a == a && b == b && theta == theta)
    {
        m_scaleX = a;
        m_scaleY = b;
        m_theta = theta;
    }

    m_theta += angle;

    m_scaleX /= multiplier;
    m_scaleY /= multiplier;
}

void Kernel::translate(const Vector2& offset)
{
    m_posX += static_cast<float>(offset[0]);
    m_posY += static_cast<float>(offset[1]);
}

void Kernel::setPos(const Vector2& pos)
{
    m_posX = static_cast<float>(pos[0]);
    m_posY = static_cast<float>(pos[1]);
}

void Kernel::rotate(const float theta, const Vector2& o)
{
    const double cs = cos(theta);
    const double sn = sin(theta);

    translate(-o);
    Vector2 p(m_posX, m_posY);
    m_posX = p[0] * cs - p[1] * sn;
    m_posY = p[0] * sn + p[1] * cs;
    translate(o);
    m_theta -= theta;
}

void Kernel::normalize()
{
    if (m_scaleY > m_scaleX)
    {
        float tmp = m_scaleX;
        m_scaleX = m_scaleY;
        m_scaleY = tmp;
        m_theta += utils::pi / 2.f;
    }

    if (m_theta > M_PI || m_theta < 0.f)
        m_theta = std::fmod(m_theta, utils::pi);
}

void Kernel::print(std::ostream& os) const
{
    os << "Kernel(Scale: (" << m_scaleX << ", " << m_scaleY << "), Pos: (" << m_posX << ", " << m_posY <<
        "), ";
    os << ", Amplitude: " << m_amplitude << ", Theta: " << m_theta << ")";
}