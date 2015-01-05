#include "surface_waveletDecomposition.h"

#include "mapHandler.h"

#define NORMALE 5

namespace CGoGN
{

namespace SCHNApps
{

bool Surface_WaveletDecomposition_Plugin::enable()
{
    m_waveletDecompositionDialog = new Dialog_Surface_WaveletDecomposition(m_schnapps);

    m_waveletDecompositionAction = new QAction("Wavelet Decomposition", this);

    m_schnapps->addMenuAction(this, "Surface;Wavelet Decomposition", m_waveletDecompositionAction);

    connect(m_waveletDecompositionAction, SIGNAL(triggered()), this, SLOT(openWaveletDecompositionDialog()));

    connect(m_waveletDecompositionDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeWaveletDecompositionDialog()));
    connect(m_waveletDecompositionDialog->button_decompose, SIGNAL(clicked()), this, SLOT(decomposeFromDialog()));
    connect(m_waveletDecompositionDialog->button_saveImages, SIGNAL(clicked()), this, SLOT(saveImagesFromDialog()));

    m_decomposition = NULL;
    m_camera = NULL;

    return true;
}

void Surface_WaveletDecomposition_Plugin::disable()
{
    disconnect(m_waveletDecompositionAction, SIGNAL(triggered()), this, SLOT(openWaveletDecompositionDialog()));

    disconnect(m_waveletDecompositionDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeWaveletDecompositionDialog()));
    disconnect(m_waveletDecompositionDialog->button_decompose, SIGNAL(clicked()), this, SLOT(decomposeFromDialog()));
    disconnect(m_waveletDecompositionDialog->button_saveImages, SIGNAL(clicked()), this, SLOT(saveImagesFromDialog()));

    delete m_decomposition;
}

void Surface_WaveletDecomposition_Plugin::openWaveletDecompositionDialog()
{
    m_waveletDecompositionDialog->show();
}

void Surface_WaveletDecomposition_Plugin::closeWaveletDecompositionDialog()
{
    m_waveletDecompositionDialog->close();
}

void Surface_WaveletDecomposition_Plugin::decomposeFromDialog()
{
    decompose();
}

void Surface_WaveletDecomposition_Plugin::saveImagesFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty())
    {
        saveImages(currentItems[0]->text());
    }
}

const QString Surface_WaveletDecomposition_Plugin::initializeObject(const QString& view, QString& filename, const bool multiple)
{
    if(!filename.isEmpty() && !m_decomposition)
    {
        QString file = filename.left(filename.lastIndexOf('.'));
        file = file.mid(file.lastIndexOf('/')).remove(0, 1);

        QString extension = filename.mid(filename.lastIndexOf('.'));
        extension.toUpper();

        QImage image;
        if(!image.load(filename, extension.toUtf8().constData()))
        {
            CGoGNerr << "Image has not been loaded correctly" << CGoGNendl;
            return NULL;
        }
        image = image.convertToFormat(QImage::Format_RGB32);

        int width = image.width(), height = image.height();

        m_decomposition = new Decomposition(width, height);

        for(int i = 0 ; i < width; ++i)
        {
            for(int j = 0; j < height; ++j)
            {
                m_decomposition->setValue(i, height-j-1, NQRgb(image.pixel(i, height-j-1)));
            }
        }

        return file;
    }
    return NULL;
}

void Surface_WaveletDecomposition_Plugin::decompose()
{
    if(m_decomposition)
    {
        bool stop = false;
        int width = m_decomposition->getWidth(), height = m_decomposition->getHeight();
        int img_width = width, img_height = height;
        int img_width2, img_height2;
        while(!stop)
        {
            img_width2  = round(img_width/2.f);
            img_height2 = round(img_height/2.f);

            //Horizontal decomposition
            std::vector<NQRgb> tmp_matrix(*m_decomposition->getMatrix());
            for(int i = 0; i < img_width; ++i)
            {
                for(int j = 0; j < img_height; ++j)
                {
                    if(i%2==0)
                    {
                        tmp_matrix[i/2+width*j] = m_decomposition->getValue(i, j);
                    }
                    else
                    {
                        NQRgb left = m_decomposition->getValue(i-1, j);
                        NQRgb result;
                        if(i != img_width-1)
                        {
                            NQRgb right = m_decomposition->getValue(i+1, j);
                            result = m_decomposition->getValue(i, j);
                            result -= (left+right)/2.f;
                            tmp_matrix[img_width2+round(i/2.f)-1+width*j] = result;
                        }
                        else
                        {
                            result = m_decomposition->getValue(i, j);
                            result -= left;
                            tmp_matrix[img_width2+round(i/2.f)-1+width*j] = result;
                        }
                    }
                }
            }

            m_decomposition->setMatrix(tmp_matrix);

            //Vertical decomposition
            for(int i = 0; i < img_width; ++i)
            {
                for(int j = 0; j < img_height; ++j)
                {
                    if(j%2 == 0)
                    {
                        tmp_matrix[i+width*(j/2)] = m_decomposition->getValue(i, j);
                    }
                    else
                    {
                        NQRgb up = m_decomposition->getValue(i, j-1);
                        NQRgb result;
                        if(j != img_height-1)
                        {
                            NQRgb down = m_decomposition->getValue(i, j+1);
                            result = m_decomposition->getValue(i, j);
                            result -= (up+down)/2.f;
                            tmp_matrix[i+width*(img_height2+round(j/2.f)-1)] = result;
                        }
                        else
                        {
                            result = m_decomposition->getValue(i, j);
                            result -= up;
                            tmp_matrix[i+width*(img_height2+round(j/2.f)-1)] = result;
                        }
                    }
                }
            }

            m_decomposition->setMatrix(tmp_matrix);

            img_width  = img_width2;
            img_height = img_height2;
            stop = true;
            if(img_width < 2 || img_height < 2)
            {
                stop = true;
            }
            else
            {
                m_decomposition->getDownDecomposition();
            }
        }
    }
}

void Surface_WaveletDecomposition_Plugin::saveImages(const QString& name, const QString& directory)
{
    if(m_decomposition && !name.isEmpty() && !directory.isEmpty())
    {
        int width = m_decomposition->getWidth();
        int height = m_decomposition->getHeight();
        QImage image(width, height, QImage::Format_RGB32);
        QString filename(directory);
        filename.append(name);
        filename.append("-");
        filename.append(QString::number(m_decomposition->getLevel()));
        filename.append(".png");

        for(int i = 0; i < width; ++i)
        {
            for(int j = 0; j < height; ++j)
            {
                NQRgb color = m_decomposition->getValue(i, j);
                image.setPixel(i, j, qRgb(qAbs(color.getRed()), qAbs(color.getGreen()), qAbs(color.getBlue())));
            }
        }

        if(!image.save(filename))
        {
            CGoGNerr << "Image '" << filename.toStdString() << "' has not been saved" << CGoGNendl;
        }
    }
}

MapHandlerGen* Surface_WaveletDecomposition_Plugin::drawCoarseImage(const QString& mapName)
{
    if(m_decomposition)
    {
        CGoGNout << "Drawing image .." << CGoGNflush;
        Utils::Chrono chrono;
        chrono.start();

        MapHandlerGen* mhg_map = m_schnapps->addMap(mapName, 2);
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);
        PFP2::MAP* map = mh_map->getMap();

        VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinates = mh_map->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
        if(!planeCoordinates.isValid())
        {
            planeCoordinates = mh_map->addAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
        }

        VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_map->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
        if(!imageCoordinates.isValid())
        {
            imageCoordinates = mh_map->addAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
        }

        int image_width = m_decomposition->getWidth(), image_height = m_decomposition->getHeight();
        int image_min_width = image_width/pow(2, m_decomposition->getLevel()), image_min_height = image_height/pow(2, m_decomposition->getLevel());

        Algo::Surface::Tilings::Square::Grid<PFP2> grid(*map, image_min_width-1, image_min_height-1);
        grid.embedIntoGrid(planeCoordinates, image_width-1, image_height-1);

        std::vector<Dart> vDarts = grid.getVertexDarts();

        for(int i = 0; i < image_min_width; ++i)
        {
            for(int j = 0; j < image_min_height; ++j)
            {
                imageCoordinates[vDarts[j*image_min_width+i]].setCoordinates(i, image_min_height-j-1);
            }
        }

        mh_map->notifyConnectivityModification();
        mh_map->updateBB(planeCoordinates);
        mh_map->notifyAttributeModification(planeCoordinates);
        mh_map->notifyAttributeModification(imageCoordinates);

        CGoGNout << ".. finished in " << chrono.elapsed() << " ms." << CGoGNendl;

//        project2DImageTo3DSpace(mapName);

        return mhg_map;
    }
    return NULL;
}

void Surface_WaveletDecomposition_Plugin::project2DImageTo3DSpace(const QString& mapName)
{
    if(m_decomposition)
    {
        MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

        if(mh_map)
        {
            PFP2::MAP* map = mh_map->getMap();

            VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
            if(!position.isValid())
            {
                position = mh_map->addAttribute<PFP2::VEC3, VERTEX>("position");
            }

            VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinates = mh_map->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
            if(!planeCoordinates.isValid())
            {
                planeCoordinates = mh_map->addAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
            }

            VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_map->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
            if(!imageCoordinates.isValid())
            {
                imageCoordinates = mh_map->addAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
            }

            if(!m_camera)
            {
                m_camera = new qglviewer::Camera(*(m_schnapps->getSelectedView()->camera()));
                m_camera->setZNearCoefficient(1.5f);
            }

            qglviewer::Vec camera_p = m_camera->position();

            PFP2::VEC3 camera_position = PFP2::VEC3(camera_p.x, camera_p.y, camera_p.z);

            PFP2::VEC3 plane_center_position = camera_position;
            plane_center_position[2] -= m_camera->zNear()*200;

            float d_C_Xc = (plane_center_position-camera_position).norm2();

            //Reprojection des points => Théorème de Thalès
            TraversorV<PFP2::MAP> trav_vert_map(*map);
            for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
            {
                PFP2::VEC3 plane_point_position = plane_center_position;
                plane_point_position[0] += planeCoordinates[d][0];
                plane_point_position[1] += planeCoordinates[d][1];

                NQRgb color = m_decomposition->getValue(imageCoordinates[d].getXCoordinate(), imageCoordinates[d].getYCoordinate());

                float z_coordinate = (m_camera->zFar()-m_camera->zNear())*(1.f-(qAbs(color.getRed())/255.f))*200;

                PFP2::VEC3 projected_point_position = plane_center_position;
                projected_point_position[2] = z_coordinate;

                float d_C_Xp = (plane_point_position-camera_position).norm();
                float d_C_Xh = (projected_point_position-camera_position).norm2();

                float d_C_Xe = sqrt(d_C_Xh/d_C_Xc)*d_C_Xp;

                position[d] = (plane_point_position-camera_position)/d_C_Xp*(d_C_Xe+d_C_Xp);
//                position[d] = plane_point_position;
//                position[d][2] += (m_camera->zFar()-m_camera->zNear())*200;
            }

            mh_map->notifyAttributeModification(position);
            mh_map->updateBB(position);
            m_schnapps->getSelectedView()->updateGL();
        }
    }
}

void Surface_WaveletDecomposition_Plugin::projectNewPointsTo3DSpace(MapHandler<PFP2>* mh_map, const std::vector<Dart>& vertices, const std::vector<NQRgb>& matrix)
{
    if(mh_map)
    {
        VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
        if(!position.isValid())
        {
            position = mh_map->addAttribute<PFP2::VEC3, VERTEX>("position");
        }

        VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinates = mh_map->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
        if(!planeCoordinates.isValid())
        {
            planeCoordinates = mh_map->addAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
        }

        VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_map->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
        if(!imageCoordinates.isValid())
        {
            imageCoordinates = mh_map->addAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
        }

        qglviewer::Vec camera_p = m_camera->position();

        PFP2::VEC3 camera_position = PFP2::VEC3(camera_p.x, camera_p.y, camera_p.z);

        PFP2::VEC3 plane_center_position = camera_position;
        plane_center_position[2] -= m_camera->zNear()*200;

        float d_C_Xc = (plane_center_position-camera_position).norm2();

        for(std::vector<Dart>::const_iterator d = vertices.begin(); d != vertices.end(); ++d)
        {
            PFP2::VEC3 plane_point_position = plane_center_position;
            plane_point_position[0] += planeCoordinates[*d][0];
            plane_point_position[1] += planeCoordinates[*d][1];

            NQRgb color = matrix[imageCoordinates[*d].getXCoordinate()+m_decomposition->getWidth()*imageCoordinates[*d].getYCoordinate()];

            float z_coordinate = (m_camera->zFar()-m_camera->zNear())*(1.f-(qAbs(color.getRed())/255.f))*200;

            PFP2::VEC3 projected_point_position = plane_center_position;
            projected_point_position[2] = z_coordinate;

            float d_C_Xp = (plane_point_position-camera_position).norm();
            float d_C_Xh = (projected_point_position-camera_position).norm2();

            float d_C_Xe = sqrt(d_C_Xh/d_C_Xc)*d_C_Xp;

            position[*d] = (plane_point_position-camera_position)/d_C_Xp*(d_C_Xe+d_C_Xp);
        }

        mh_map->notifyAttributeModification(position);
        mh_map->updateBB(position);
        m_schnapps->getSelectedView()->updateGL();
    }
}

void Surface_WaveletDecomposition_Plugin::triangulateMap(const QString& mapName)
{
    if(m_decomposition)
    {
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

        if(mh_map)
        {
            PFP2::MAP* map = mh_map->getMap();

            DartMarker<PFP2::MAP> marker(*map);

            TraversorF<PFP2::MAP> trav_face_map(*map);
            for(Dart d = trav_face_map.begin(); d != trav_face_map.end(); d = trav_face_map.next())
            {
                if(!marker.isMarked(d))
                {
                    Dart phi_11_d;
                    if(map->vertexDegree(d)<map->vertexDegree(map->phi1(d)))
                    {
                        phi_11_d = map->phi_1(map->phi_1(d));
                    }
                    else
                    {
                        phi_11_d = map->phi<11>(d);
                    }
                    map->splitFace(d, phi_11_d);
                    marker.markOrbit<FACE>(d);
                    marker.markOrbit<FACE>(phi_11_d);
                }
            }

            mh_map->notifyConnectivityModification();

            m_schnapps->getSelectedView()->updateGL();
        }
    }
}

void Surface_WaveletDecomposition_Plugin::moveUpDecomposition(const QString& mapName)
{
    if(m_decomposition && m_decomposition->getLevel()>0)
    {
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));
        if(mh_map)
        {
            PFP2::MAP* map = mh_map->getMap();

            VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
            if(!position.isValid())
            {
                CGoGNerr << "position attribute is not valid" << CGoGNendl;
            }

            VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinates = mh_map->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
            if(!planeCoordinates.isValid())
            {
                CGoGNerr << "PlaneCoordinates attribute is not valid" << CGoGNendl;
            }

            VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_map->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
            if(!imageCoordinates.isValid())
            {
                CGoGNerr << "ImageCoordinates attribute is not valid" << CGoGNendl;
            }

            int width = m_decomposition->getWidth()/pow(2, m_decomposition->getLevel());
            int height = m_decomposition->getHeight()/pow(2, m_decomposition->getLevel());

            DartMarker<PFP2::MAP> marker(*map);
            DartMarker<PFP2::MAP> marker_face(*map);
            DartMarker<PFP2::MAP> marker_horizontal(*map), marker_vertical(*map), marker_diagonal(*map);

            std::vector<Dart> verticesAdded;
            verticesAdded.reserve(width*height);

            std::vector<Dart> verticesModified;
            verticesModified.reserve(width*height);

            TraversorF<PFP2::MAP> trav_face_map(*map);
            for(Dart d = trav_face_map.begin(); d != trav_face_map.end(); d = trav_face_map.next())
            {
                if(!marker_face.isMarked(d))
                {
                    Traversor2FE<PFP2::MAP> trav_edge_face_map(*map, d);
                    for(Dart dd = trav_edge_face_map.begin(); dd != trav_edge_face_map.end(); dd = trav_edge_face_map.next())
                    {
                        if(!marker.isMarked(dd))
                        {
                            Dart dd1 = map->phi1(dd);
                            Dart ddd = map->cutEdge(dd);
                            marker.markOrbit<EDGE>(dd);
                            marker.markOrbit<EDGE>(ddd);

                            if(imageCoordinates[dd].getXCoordinate()==imageCoordinates[dd1].getXCoordinate())
                            {
                                marker_vertical.markOrbit<VERTEX>(ddd);
                            }
                            else if(imageCoordinates[dd].getYCoordinate()==imageCoordinates[dd1].getYCoordinate())
                            {
                                marker_horizontal.markOrbit<VERTEX>(ddd);
                            }
                            else
                            {
                                marker_diagonal.markOrbit<VERTEX>(ddd);
                            }
                            verticesAdded.push_back(ddd);
                            verticesModified.push_back(dd);
                            planeCoordinates[ddd] = (planeCoordinates[dd]+planeCoordinates[dd1])/2.f;
                        }
                    }
                }
            }

//            TraversorV<PFP2::MAP> trav_vert_map(*map);
            for(std::vector<Dart>::const_iterator d = verticesAdded.begin(); d != verticesAdded.end(); ++d)
            {
                if(marker_horizontal.isMarked(*d))
                {
                    Dart d_1 = map->phi_1(*d), d1 = map->phi1(*d);
                    int min = imageCoordinates[d_1].getXCoordinate()>imageCoordinates[d1].getXCoordinate()?imageCoordinates[d1].getXCoordinate():imageCoordinates[d_1].getXCoordinate();
                    imageCoordinates[*d].setCoordinates((min*2+1), imageCoordinates[d_1].getYCoordinate()*2);
                }
                else if(marker_vertical.isMarked(*d))
                {
                    Dart d_1 = map->phi_1(*d), d1 = map->phi1(*d);
                    int min = imageCoordinates[d_1].getYCoordinate()>imageCoordinates[d1].getYCoordinate()?imageCoordinates[d1].getYCoordinate():imageCoordinates[d_1].getYCoordinate();
                    imageCoordinates[*d].setCoordinates(imageCoordinates[d_1].getXCoordinate()*2, (min*2+1));
                }
                else if(marker_diagonal.isMarked(*d))
                {
                    Dart d_1 = map->phi_1(map->phi2(*d)), d1 = map->phi2(*d);
                    int min_x = imageCoordinates[d_1].getXCoordinate()>imageCoordinates[d1].getXCoordinate()?imageCoordinates[d1].getXCoordinate():imageCoordinates[d_1].getXCoordinate();
                    int min_y = imageCoordinates[d_1].getYCoordinate()>imageCoordinates[d1].getYCoordinate()?imageCoordinates[d1].getYCoordinate():imageCoordinates[d_1].getYCoordinate();
                    imageCoordinates[*d].setCoordinates((min_x*2+1), (min_y*2+1));
                }
            }

            for(Dart d = trav_face_map.begin(); d != trav_face_map.end(); d = trav_face_map.next())
            {
                if(!marker_face.isMarked(d))
                {
                    Dart d1 = map->phi<11>(d);
                    Dart d11 = map->phi<11>(d1);

                    map->splitFace(map->phi1(d), map->phi_1(d));
                    map->splitFace(map->phi1(d1), map->phi_1(d1));
                    map->splitFace(map->phi1(d11), map->phi_1(d11));
                    marker_face.markOrbit<FACE>(d);
                    marker_face.markOrbit<FACE>(d1);
                    marker_face.markOrbit<FACE>(d11);
                    marker_face.markOrbit<FACE>(map->phi<12>(d11));
                }
            }

            for(std::vector<Dart>::const_iterator d = verticesModified.begin(); d != verticesModified.end(); ++d)
            {
                imageCoordinates[*d].setCoordinates(imageCoordinates[*d].getXCoordinate()*2, imageCoordinates[*d].getYCoordinate()*2);
            }

            std::vector<NQRgb> matrix = std::vector<NQRgb>(*m_decomposition->getMatrix());

            for(int i = 0; i < width*2; ++i)
            {
                for(int j = 0; j < height*2; ++j)
                {
                    int index = i+m_decomposition->getWidth()*j;
                    if(j%2 == 1)
                    {
                        NQRgb up = m_decomposition->getValue(i, round(j/2.f)-1);
                        NQRgb result = m_decomposition->getValue(i, height+round(j/2.f)-1);
                        if(j != height*2-1)
                        {
                            NQRgb down = m_decomposition->getValue(i, round(j/2.f)+1);
                            result += (up+down)/2.f;
                        }
                        else
                        {
                            result += up;
                        }
                        matrix[index] = result;
                    }
                    else
                    {
                        matrix[index] = m_decomposition->getValue(i, j/2);
                    }
                }
            }

            std::vector<NQRgb> matrix2 = std::vector<NQRgb>(matrix);

            for(int i = 0; i < width*2; ++i)
            {
                for(int j = 0; j < height*2 ; ++j)
                {
                    int index = i+m_decomposition->getWidth()*j;
                    if(i%2 == 1)
                    {
                        NQRgb left = matrix2[round(i/2.f)-1+m_decomposition->getWidth()*j];
                        NQRgb result = matrix2[width+round(i/2.f)-1+m_decomposition->getWidth()*j];
                        if(i != width*2-1)
                        {
                            NQRgb right = matrix2[round(i/2.f)+1+m_decomposition->getWidth()*j];
                            result += (left+right)/2.f;
                        }
                        else
                        {
                            result += left;
                        }
                        matrix[index] = result;
                    }
                    else
                    {
                        matrix[index] = matrix2[i/2+m_decomposition->getWidth()*j];
                    }
                }
            }

//            marker_face.unmarkAll();

//            for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
//            {
//                if(position[d]==PFP2::VEC3())
//                {
//                    //CGoGNout << imageCoordinates[d].getXCoordinate() << " ; " << imageCoordinates[d].getYCoordinate() << CGoGNendl;
//                }
//            }

//            QImage image(width*2, height*2, QImage::Format_RGB32);

//            for(int i = 0; i < width*2; ++i)
//            {
//                for(int j = 0; j < height*2 ; ++j)
//                {
//                    NQRgb value = matrix[i+m_decomposition->getWidth()*j];
//                    image.setPixel(i, j, qRgb(qAbs(value.getRed()), qAbs(value.getBlue()), qAbs(value.getGreen())));
//                }
//            }

//            QString filename("/home/blettere/Projets/Models/Test/");
//            filename.append(mapName);
//            filename.append(".png");

//            if(!image.save(filename))
//            {
//                CGoGNerr << "Image '" << filename.toStdString() << "' has not been saved" << CGoGNendl;
//            }

            mh_map->notifyAttributeModification(planeCoordinates);
            mh_map->notifyAttributeModification(imageCoordinates);
            mh_map->notifyConnectivityModification();

            projectNewPointsTo3DSpace(mh_map, verticesAdded, matrix);
        }
    }
}

void Surface_WaveletDecomposition_Plugin::moveDownDecomposition(const QString& mapName)
{
    if(m_decomposition && (m_decomposition->getWidth()>=2 && m_decomposition->getHeight()>=2))
    {

    }
}

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_Plugin, Surface_WaveletDecomposition_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_PluginD, Surface_WaveletDecomposition_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
