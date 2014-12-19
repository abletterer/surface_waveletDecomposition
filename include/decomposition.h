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

    ~NQRgb()
    {}

    int getRed() { return m_r; }
    void setRed(int red) { m_r = red; }

    int getGreen() { return m_g; }
    void setGreen(int green) { m_g = green; }

    int getBlue() { return m_b; }
    void setBlue(int blue) { m_b = blue; }
private:
    int m_r, m_g, m_b;
};

class Decomposition
{
public:
    Decomposition(const int x=0, const int y=0, const short int level=0, Decomposition* parent = NULL)
        : m_image(),
          m_correction_x(x), m_correction_y(y),
          m_horizontal_correction(x*y), m_vertical_correction(x*y),
          m_child(NULL),
          m_parent(parent),
          m_level(level)
    {
    }

    ~Decomposition()
    {
        if(m_child)
        {
            delete m_child;
        }
    }

    QImage& getImage() { return m_image; }
    void setImage(const QImage& image) { m_image = image; }

    int getCorrectionX() { return m_correction_x; }
    int getCorrectionY() { return m_correction_y; }

    NQRgb& getHorizontalCorrection(const int x, const int y)
    {
        if(x >= 0 && y >= 0 && x < m_correction_x && y < m_correction_y)
        {
            return m_horizontal_correction[x+y*m_correction_x];
        }
        else
        {
            CGoGNerr << "Get : Indices {" << x << ", " << y << "} not in the range [0; {" << m_correction_x << ", " << m_correction_y << "}]" << CGoGNendl;
        }
    }
    void setHorizontalCorrection(const int x, const int y, const NQRgb& correction)
    {
        if(x >= 0 && y >= 0 && x < m_correction_x && y < m_correction_y)
        {
            m_horizontal_correction[x+y*m_correction_x] = correction;
        }
        else
        {
            CGoGNerr << "Set : Indices {" << x << ", " << y << "} not in the range [0; {" << m_correction_x << ", " << m_correction_y << "}]" << CGoGNendl;
        }
    }

    NQRgb& getVerticalCorrection(const int x, const int y)
    {
        if(x >= 0 && y >= 0 && x < m_correction_x && y < m_correction_y)
        {
            return m_vertical_correction[x+y*m_correction_x];
        }
        else
        {
            CGoGNerr << "Get : Indices {" << x << ", " << y << "} not in the range [0; {" << m_correction_x << ", " << m_correction_y << "}]" << CGoGNendl;
        }
    }
    void setVerticalCorrection(const int x, const int y, const NQRgb& correction)
    {
        if(x >= 0 && y >= 0 && x < m_correction_x && y < m_correction_y)
        {
            m_vertical_correction[x+y*m_correction_x] = correction;
        }
        else
        {
            CGoGNerr << "Set : Indices {" << x << ", " << y << "} not in the range [0; {" << m_correction_x << ", " << m_correction_y << "}]" << CGoGNendl;
        }
    }

    Decomposition* getChild() { return m_child; }
    Decomposition* addChild()
    {
        if(!m_child)
        {
            m_child = new Decomposition(m_image.width()/2, m_image.height()/2, m_level+1, this);
            return m_child;
        }
        CGoGNerr << "Analysis already done for level " << m_level << CGoGNendl;
        return this;
    }

    Decomposition* getParent() { return m_parent; }

    short int getLevel() { return m_level; }

private:
    QImage m_image;
    int m_correction_x, m_correction_y;
    std::vector<NQRgb> m_horizontal_correction;
    std::vector<NQRgb> m_vertical_correction;
    Decomposition* m_child;
    Decomposition* m_parent;
    short int m_level;
};

}
}

#endif
