#ifndef _IMAGECOORDINATES_H_
#define _IMAGECOORDINATES_H_

namespace CGoGN
{
namespace SCHNApps
{

/*
* Classe définissant l'attribut de sommet Coordinates, définissant les coordonnées images associées à un point
*/
class ImageCoordinates
{
public:
    ImageCoordinates()
        : m_x(0.f), m_y(0.f)
    {}

    ~ImageCoordinates()
    {
    }

    void setCoordinates(float x, float y)
    {
        m_x = x;
        m_y = y;
    }

    float getXCoordinate() { return m_x; }
    void setXCoordinate(float x)
    {
        m_x = x;
    }

    float getYCoordinate() { return m_y; }
    void setYCoordinate(float y)
    {
        m_y = y;
    }

    friend std::ostream& operator<< (std::ostream& stream, const ImageCoordinates &coordinates)
    {
        return stream;
    }

    static std::string CGoGNnameOfType()
    {
        return "ImageCoordinates" ;
    }

private:
    float m_x, m_y;
};

} //namespace SCHNApps
} //namespace CGoGN

#endif
