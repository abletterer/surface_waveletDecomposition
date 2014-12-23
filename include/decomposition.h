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
        m_r = (a.getRed() + b.getRed())/2;
        m_g = (a.getGreen() + b.getGreen())/2;
        m_b = (a.getBlue() + b.getBlue())/2;
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
          m_level(level)
    {
    }

    ~Decomposition()
    {}

    void setValue(const int x, const int y, const NQRgb& value)
    {
        m_matrix_decomposition[x+m_width*y] = value;
    }

    NQRgb getValue(const int x, const int y)
    {
        if(x >= 0 && y >= 0 && x < m_width && y < m_height)
        {
            return m_matrix_decomposition[x+m_width*y];
        }
        else
        {
            CGoGNerr << "Get : Indices {" << x << ", " << y << "} not in the range [0; {" << m_width-1 << ", " << m_height-1 << "}]" << CGoGNendl;
            return NQRgb();
        }
    }

    int getWidth() { return m_width; }

    int getHeight() { return m_height; }

    const int getLevel() { return m_level; }
    void setLevel(int level) { m_level = level; }
    void getUpDecomposition() { --m_level; }
    void getDownDecomposition() { ++m_level; }

    void setMatrix(const std::vector<NQRgb>& matrix)
    {
        m_matrix_decomposition = std::vector<NQRgb>(matrix);
    }

    std::vector<NQRgb>* getMatrix() { return &m_matrix_decomposition; }

private:
    int m_width, m_height;
    std::vector<NQRgb> m_matrix_decomposition;
    int m_level;
};

}
}

#endif
