#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

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

class Decomposition
{
public:
    Decomposition(int x, int y, short int level, Axis transformation_type, Decomposition* parent = NULL)
        : m_image(),
          m_correction_x(x), m_correction_y(y),
          m_correction(m_correction_x*m_correction_y),
          m_child(NULL),
          m_parent(parent),
          m_level(level),
          m_transformation_type(transformation_type)
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

    QRgb& getCorrection(int x, int y)
    {
        if(x < m_correction_x && y < m_correction_y)
        {
            return m_correction[x+y*m_correction_x];
        }
    }
    void setCorrection(const QRgb& correction, int x, int y)
    {
        if(x < m_correction_x && y < m_correction_y)
        {
            m_correction[x+y*m_correction_x] = correction;
        }
    }

    Decomposition* getChild() { return m_child; }
    Decomposition* addChild()
    {
        if(!m_child)
        {
            if(getTransformationType()==VERTICAL)
            {
                m_child = new Decomposition(m_correction_x/2+m_correction_x%2, m_correction_y/2+m_correction_y%2, m_level+1, HORIZONTAL, this);
            }
            else
            {
                m_child = new Decomposition(m_correction_x/2+m_correction_x%2, m_correction_y/2+m_correction_y%2, m_level+1, VERTICAL, this);
            }
            return m_child;
        }
        CGoGNerr << "Analysis already done for level " << m_level << CGoGNendl;
        return this;
    }

    Decomposition* getParent() { return m_parent; }

    short int getLevel() { return m_level; }

    Axis getTransformationType() { return m_transformation_type; }
    void setTransformationType(Axis transformation_type) { m_transformation_type = transformation_type; }

private:
    QImage m_image;
    int m_correction_x, m_correction_y;
    std::vector<QRgb> m_correction;
    Decomposition* m_child;
    Decomposition* m_parent;
    short int m_level;
    Axis m_transformation_type;
};

}
}

#endif // DECOMPOSITION_H
