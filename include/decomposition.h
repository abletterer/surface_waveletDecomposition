#ifndef _DECOMPOSITION_H_
#define _DECOMPOSITION_H_

#include <QImage>

#include <vector>

namespace CGoGN
{
namespace SCHNApps
{

enum Axis
{
    HORIZONTAL=0,
    VERTICAL
};

class NQRgb
{
public:
    NQRgb(int red=0, int green=0, int blue=0)
        : m_r(red), m_g(green), m_b(blue)
    {}

    NQRgb(QRgb color)
    {
        m_r = qRed(color);
        m_g = qGreen(color);
        m_b = qBlue(color);
    }

    ~NQRgb()
    {}

    int getRed() { return m_r; }
    void setRed(int red) { m_r = red; }

    int getGreen() { return m_g; }
    void setGreen(int green) { m_g = green; }

    int getBlue() { return m_b; }
    void setBlue(int blue) { m_b = blue; }

    void setMean(NQRgb a, NQRgb b)
    {
        m_r = round((a.getRed() + b.getRed())/2.f);
        m_g = round((a.getGreen() + b.getGreen())/2.f);
        m_b = round((a.getBlue() + b.getBlue())/2.f);
    }

    NQRgb operator+(NQRgb rgb)
    {
        return NQRgb(m_r+rgb.getRed(), m_g+rgb.getGreen(), m_b+rgb.getBlue());
    }

    NQRgb operator+=(NQRgb rgb)
    {
        m_r += rgb.getRed();
        m_g += rgb.getGreen();
        m_b += rgb.getBlue();
        return *this;
    }

    NQRgb operator-(NQRgb rgb)
    {
        return NQRgb(m_r-rgb.getRed(), m_g-rgb.getGreen(), m_b-rgb.getBlue());
    }

    NQRgb operator-=(NQRgb rgb)
    {
        m_r -= rgb.getRed();
        m_g -= rgb.getGreen();
        m_b -= rgb.getBlue();
        return *this;
    }

    NQRgb operator/(const float value)
    {
        return NQRgb(round(m_r/value), round(m_g/value), round(m_b/value));
    }

    NQRgb operator/=(const float value)
    {
        m_r = round(m_r/value);
        m_g = round(m_g/value);
        m_b = round(m_b/value);
        return *this;
    }

private:
    int m_r, m_g, m_b;
};

class Decomposition
{
public:
    Decomposition(const int width, const int height, const int level = 0)
        : m_width(width),
          m_height(height),
          m_matrix_decomposition(width*height),
          m_level(level),
          m_max_level(level)
    {
    }

    ~Decomposition()
    {}

    void setCoefficient(const int x, const int y, const int value)
    {
        if(x >= 0 && y >= 0 && x < m_width && y < m_height)
        {
            m_matrix_decomposition[x+m_width*y] = value;
        }
        else
        {
            CGoGNerr << "setCoefficient : Indices {" << x << ", " << y << "} not in the range [0; {" << m_width-1 << ", " << m_height-1 << "}]" << CGoGNendl;
        }
    }

    int getCoefficient(const int x, const int y)
    {
        if(x >= 0 && y >= 0 && x < m_width && y < m_height)
        {
            return m_matrix_decomposition[x+m_width*y];
        }
        else
        {
            CGoGNerr << "getCoefficient : Indices {" << x << ", " << y << "} not in the range [0; {" << m_width-1 << ", " << m_height-1 << "}]" << CGoGNendl;
            return 0;
        }
    }

    int getWidth() { return m_width; }

    int getHeight() { return m_height; }

    const int getLevel() { return m_level; }
    void setLevel(int level) { m_level = level; }

    const int getMaxLevel() { return m_max_level; }
    void setMaxLevel(int level) { m_max_level = level; }

    void moveUpDecomposition() { --m_level; }
    void moveDownDecomposition() { ++m_level; }

    void setCoefficientMatrix(const std::vector<int>& matrix)
    {
        m_matrix_decomposition = std::vector<int>(matrix);
    }

    std::vector<int>* getCoefficientMatrix() { return &m_matrix_decomposition; }

private:
    int m_width, m_height;
    std::vector<int> m_matrix_decomposition;
    int m_level;
    int m_max_level;
};

}
}

#endif
