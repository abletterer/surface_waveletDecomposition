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

        int imageX = image.width(), imageY = image.height();

        m_decomposition = new Decomposition();
        m_decomposition->setImage(image);

        return file;
    }
    return NULL;
}

void Surface_WaveletDecomposition_Plugin::decompose()
{
    if(m_decomposition)
    {
        //We search the last level of decomposition already computed
        Decomposition* decomposition = m_decomposition;
        while(decomposition->getChild())
        {
            decomposition = decomposition->getChild();
        }

        while(decomposition->getImage().width()>2 && decomposition->getImage().height()>2)
        {
            //If a decomposition is still feasible

            QImage& image_parent = decomposition->getImage();

            decomposition = decomposition->addChild();

            int image_width = image_parent.width()/2 + image_parent.width()%2;
            int image_height = image_parent.height()/2 + image_parent.height()%2;

            //Creation of a new image composed of the even of parent image
            QImage image(image_width, image_height, image_parent.format());

            for(int i = 0; i < image_width; ++i)
            {
                for(int j = 0; j < image_height; ++j)
                {
                    image.setPixel(i, j, image_parent.pixel(i*2, j*2));
                }
            }

            decomposition->setImage(image);

            //Creation of a matrix composed of the difference between odd of parent image and prediction value (linear interpolation)
            QRgb cur_pixel;
            NQRgb hori_pixel, vert_pixel;
            for(int i = 0; i < image_parent.width(); ++i)
            {
                for(int j = 0; j < image_parent.height(); ++j)
                {
                    if(i%2!=0 || j%2!=0)
                    {
                        cur_pixel = image_parent.pixel(i, j);

                        hori_pixel = NQRgb();
                        vert_pixel = NQRgb();

                        if(i%2==1)
                        {
                            //Odd column number
                            if(i==image_parent.width()-1)
                            {
                                //Mirror effect at the border of the image
                                hori_pixel.setRed(qRed(cur_pixel) - qRed(image_parent.pixel(i-1,j)));
                                hori_pixel.setGreen(qGreen(cur_pixel) - qGreen(image_parent.pixel(i-1,j)));
                                hori_pixel.setBlue(qBlue(cur_pixel) - qBlue(image_parent.pixel(i-1,j)));
                            }
                            else
                            {
                                hori_pixel.setRed(qRed(cur_pixel) - (qRed(image_parent.pixel(i-1,j))+qRed(image_parent.pixel(i+1,j)))/2);
                                hori_pixel.setGreen(qGreen(cur_pixel) - (qGreen(image_parent.pixel(i-1,j))+qGreen(image_parent.pixel(i+1,j)))/2);
                                hori_pixel.setBlue(qBlue(cur_pixel) - (qBlue(image_parent.pixel(i-1,j))+qBlue(image_parent.pixel(i+1,j)))/2);
                            }
                        }
                        if(j%2==1)
                        {
                            //Odd line number
                            if(j==image_parent.height()-1)
                            {
                                //Mirror effect at the border of the image
                                vert_pixel.setRed(qRed(cur_pixel) - qRed(image_parent.pixel(i,j-1)));
                                vert_pixel.setGreen(qGreen(cur_pixel) - qGreen(image_parent.pixel(i,j-1)));
                                vert_pixel.setBlue(qBlue(cur_pixel) - qBlue(image_parent.pixel(i,j-1)));
                            }
                            else
                            {
                                vert_pixel.setRed(qRed(cur_pixel) - (qRed(image_parent.pixel(i,j-1))+qRed(image_parent.pixel(i,j+1)))/2);
                                vert_pixel.setGreen(qGreen(cur_pixel) - (qGreen(image_parent.pixel(i,j-1))+qGreen(image_parent.pixel(i,j+1)))/2);
                                vert_pixel.setBlue(qBlue(cur_pixel) - (qBlue(image_parent.pixel(i,j-1))+qBlue(image_parent.pixel(i,j+1)))/2);
                            }
                        }
                        decomposition->setHorizontalCorrection(i/2, j/2, hori_pixel);
                        decomposition->setVerticalCorrection(i/2, j/2, vert_pixel);
                    }
                }
            }
            CGoGNout << "New level of decomposition created" << CGoGNendl;
        }
    }
}

void Surface_WaveletDecomposition_Plugin::saveImages(const QString& name)
{
    if(m_decomposition && !name.isEmpty())
    {
        Decomposition* decomposition = m_decomposition;
        do
        {
            QImage image = decomposition->getImage();
            QString filename("/home/blettere/Projets/Models/Decomposition/");
            filename.append(name);
            filename.append("-");
            QString filename2(filename);
            filename.append(QString::number(decomposition->getLevel()));
            filename.append(".png");

            if(!image.save(filename))
            {
                CGoGNerr << "Image '" << filename.toStdString() << "' has not been saved" << CGoGNendl;
            }

            filename2.append(QString::number(decomposition->getLevel()+1));
            filename = filename2;
            filename.append("-horizontalCoherence.png");
            filename2.append("-verticalCoherence.png");

            decomposition = decomposition->getChild();

            if(decomposition)
            {
                QImage image2(decomposition->getCorrectionX(), decomposition->getCorrectionY(), image.format());
                image = QImage(decomposition->getCorrectionX(), decomposition->getCorrectionY(), image2.format());
                NQRgb correction;
                for(int i = 0; i < decomposition->getCorrectionX(); ++i)
                {
                    for(int j = 0; j < decomposition->getCorrectionY(); ++j)
                    {
                        correction = decomposition->getHorizontalCorrection(i, j);
                        image.setPixel(i, j, qRgb(qAbs(correction.getRed()), qAbs(correction.getGreen()), qAbs(correction.getBlue())));
                        correction = decomposition->getVerticalCorrection(i, j);
                        image2.setPixel(i, j, qRgb(qAbs(correction.getRed()), qAbs(correction.getGreen()), qAbs(correction.getBlue())));
                    }
                }

                if(!image.save(filename))
                {
                    CGoGNerr << "Image '" << filename2.toStdString() << "' has not been saved" << CGoGNendl;
                }

                if(!image2.save(filename2))
                {
                    CGoGNerr << "Image '" << filename2.toStdString() << "' has not been saved" << CGoGNendl;
                }
            }
        } while(decomposition);
    }
}

MapHandlerGen* Surface_WaveletDecomposition_Plugin::drawCoarseImage(const QString& mapName)
{
    if(m_decomposition)
    {
        MapHandlerGen* mhg_map = m_schnapps->addMap(mapName, 2);
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);
        PFP2::MAP* map = mh_map->getMap();

        VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
        if(!position.isValid())
        {
            position = mh_map->addAttribute<PFP2::VEC3, VERTEX>("position");
        }

        VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_map->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
        if(!imageCoordinates.isValid())
        {
            imageCoordinates = mh_map->addAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
        }

        Decomposition* decomposition = m_decomposition;
//        while(decomposition->getChild())
//        {
//            decomposition = decomposition->getChild();
//        }

        QImage image = decomposition->getImage();
        int imageX = image.width(), imageY = image.height();

        Algo::Surface::Tilings::Square::Grid<PFP2> grid(*map, imageX-1, imageY-1);
        grid.embedIntoGrid(position, imageX-1, imageY-1);

        std::vector<Dart> vDarts = grid.getVertexDarts();

        mh_map->updateBB(position);

        qglviewer::Vec bb_max = mh_map->getBBmax();
        qglviewer::Vec bb_min = mh_map->getBBmin();

        for(int i = 0; i < imageX; ++i)
        {
            for(int j = 0; j < imageY; ++j)
            {
                imageCoordinates[vDarts[j*imageX+i]].setCoordinates(i,imageY-j-1);
            }
        }

        QImage image_parent = m_decomposition->getImage();
        int image_parentX = image_parent.width(), image_parentY = image_parent.height();

        float width_step_parent = (bb_max.x-bb_min.x)/(image_parentX-1);
        float height_step_parent = (bb_max.y-bb_min.y)/(image_parentY-1);

        int shift_counterX = 0, shift_counterY = 0;

        while(image_parentX != imageX && image_parentY != imageY)
        {
            if(image_parentX%2==0)
            {
                ++shift_counterX;
            }
            if(image_parentY%2==0)
            {
                ++shift_counterY;
            }

            image_parentX = image_parentX/2+image_parentX%2;
            image_parentY = image_parentY/2+image_parentY%2;
        }

        shift_counterX = pow(2, shift_counterX);
        shift_counterY = pow(2, shift_counterY);

        PFP2::MATRIX44 transform_matrix;
        transform_matrix.identity();
        transform_matrix.setSubVectorV(0, 3, PFP2::VEC4(-width_step_parent*shift_counterX/2., height_step_parent*shift_counterY/2., 0., 1.));
        grid.transform(position, transform_matrix);

        mh_map->notifyAttributeModification(position);
        mh_map->notifyAttributeModification(imageCoordinates);
        mh_map->notifyConnectivityModification();
        mh_map->updateBB(position);

        m_schnapps->getSelectedView()->updateGL();

        m_decomposition = decomposition;

        return mhg_map;
    }
    return NULL;
}

void Surface_WaveletDecomposition_Plugin::project2DImageTo3DSpace(const QString& mapName)
{
    if(m_decomposition)
    {
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

        if(mh_map)
        {
            PFP2::MAP* map = mh_map->getMap();

            if(!m_camera)
            {
                m_camera = new qglviewer::Camera(*(m_schnapps->getSelectedView()->camera()));
                m_camera->setZNearCoefficient(1.f);
            }

            VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
            if(!position.isValid())
            {
                CGoGNerr << "Position attribute is not valid" << CGoGNendl;
            }

            VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_map->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
            if(!imageCoordinates.isValid())
            {
                CGoGNerr << "ImageCoordinates attribute is not valid" << CGoGNendl;
            }

            QImage image = m_decomposition->getImage();

            qglviewer::Vec camera_p = m_camera->position();

            PFP2::VEC3 camera_position = PFP2::VEC3(camera_p.x, camera_p.y, camera_p.z);

            PFP2::VEC3 plane_center_position = camera_position;
            plane_center_position[2] -= m_camera->zNear();

            float d_C_Xc = (plane_center_position-camera_position).norm2();

            //Reprojection des points => Théorème de Thalès
            TraversorV<PFP2::MAP> trav_vert_map(*map);
            for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
            {
                PFP2::VEC3 plane_point_position = plane_center_position;
                plane_point_position[0] += position[d][0];
                plane_point_position[1] += position[d][1];

                QRgb pixel = image.pixel(imageCoordinates[d].getXCoordinate(), imageCoordinates[d].getYCoordinate());

                float z_coordinate = (m_camera->zFar()-m_camera->zNear())*(1.f-(qRed(pixel)/255.f));

                PFP2::VEC3 projected_point_position = plane_center_position;
                projected_point_position[2] -= z_coordinate;

                float d_C_Xp = (plane_point_position-camera_position).norm();
                float d_C_Xh = (projected_point_position-camera_position).norm2();

                float d_C_Xe = sqrt(d_C_Xh/d_C_Xc)*d_C_Xp;

                position[d] = (plane_point_position-camera_position)/d_C_Xp*(d_C_Xe+d_C_Xp);
            }

            mh_map->notifyAttributeModification(position);
            mh_map->notifyConnectivityModification();
            mh_map->updateBB(position);

            //camera->setSceneBoundingBox(mh_map->getBBmin(), mh_map->getBBmax());

            m_schnapps->getSelectedView()->camera()->setSceneBoundingBox(mh_map->getBBmin(), mh_map->getBBmax());

            m_schnapps->getSelectedView()->updateGL();
        }
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
            VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
            if(!position.isValid())
            {
                CGoGNerr << "Position attribute is not valid" << CGoGNendl;
            }

            DartMarker<PFP2::MAP> marker(*map);

            TraversorF<PFP2::MAP> trav_face_map(*map);
            for(Dart d = trav_face_map.begin(); d != trav_face_map.end(); d = trav_face_map.next())
            {
                if(!marker.isMarked(d))
                {
                    Dart phi_11_d = map->phi<11>(d);
                    map->splitFace(d, phi_11_d);
                    marker.markOrbit<FACE>(d);
                    marker.markOrbit<FACE>(phi_11_d);
                }
            }

            mh_map->notifyAttributeModification(position);
            mh_map->notifyConnectivityModification();

            m_schnapps->getSelectedView()->updateGL();
        }
    }
}

void Surface_WaveletDecomposition_Plugin::moveUpDecomposition(const QString& mapName)
{
    if(m_decomposition->getParent())
    {
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

        if(mh_map)
        {
            PFP2::MAP* map = mh_map->getMap();

            VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
            VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_map->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");

            DartMarker<PFP2::MAP> marker(*map);
            DartMarker<PFP2::MAP> marker_face(*map);
            DartMarker<PFP2::MAP> marker_horizontal(*map);
            DartMarker<PFP2::MAP> marker_vertical(*map);
            DartMarker<PFP2::MAP> marker_diagonal(*map);
            DartMarker<PFP2::MAP> marker_update(*map);

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

                            if(!marker_update.isMarked(dd))
                            {
                                imageCoordinates[dd].setXCoordinate(imageCoordinates[dd].getXCoordinate()*2);
                                imageCoordinates[dd].setYCoordinate(imageCoordinates[dd].getYCoordinate()*2);
                                marker_update.markOrbit<VERTEX>(dd);
                            }

                            if(!marker_update.isMarked(dd1))
                            {
                                imageCoordinates[dd1].setXCoordinate(imageCoordinates[dd1].getXCoordinate()*2);
                                imageCoordinates[dd1].setYCoordinate(imageCoordinates[dd1].getYCoordinate()*2);
                                marker_update.markOrbit<VERTEX>(dd1);
                            }

                            //Position = projection du point du plan vers l'espace avec la caméra

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
                            position[ddd] = (position[dd]+position[dd1])/2.f;
                        }
                    }
                }
            }

            TraversorV<PFP2::MAP> trav_vert_map(*map);
            for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
            {
                if(marker_vertical.isMarked(d))
                {
                    Dart d_1 = map->phi_1(d);
                    Dart d1 = map->phi1(d);
                    imageCoordinates[d].setXCoordinate(imageCoordinates[d_1].getXCoordinate());
                    imageCoordinates[d].setYCoordinate((imageCoordinates[d_1].getYCoordinate()+imageCoordinates[d1].getYCoordinate())/2);
                }
                else if(marker_horizontal.isMarked(d))
                {
                    Dart d_1 = map->phi_1(d);
                    Dart d1 = map->phi1(d);
                    imageCoordinates[d].setXCoordinate((imageCoordinates[d_1].getXCoordinate()+imageCoordinates[d1].getYCoordinate())/2);
                    imageCoordinates[d].setYCoordinate(imageCoordinates[d_1].getYCoordinate());
                }
            }

            for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
            {
                if(marker_diagonal.isMarked(d))
                {
                    Dart d_1 = map->phi_1(map->phi2(d));
                    Dart d1 = map->phi_1(map->phi_1(d));
                    imageCoordinates[d].setXCoordinate((imageCoordinates[d_1].getXCoordinate()+imageCoordinates[d1].getXCoordinate())/2);
                    imageCoordinates[d].setYCoordinate(imageCoordinates[d_1].getYCoordinate());
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

            QImage image = m_decomposition->getImage();

            qglviewer::Vec camera_p = m_camera->position();
            PFP2::VEC3 camera_position = PFP2::VEC3(camera_p.x, camera_p.y, camera_p.z);

            PFP2::VEC3 plane_center_position = camera_position;
            plane_center_position[2] -= m_camera->zNear();

            float d_C_Xc = (plane_center_position-camera_position).norm2();

            qglviewer::Vec bb_min = mh_map->getBBmin();
            qglviewer::Vec bb_max = mh_map->getBBmax();

            std::swap(bb_min.y, bb_max.y);

            for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
            {
                if(marker_vertical.isMarked(d) || marker_horizontal.isMarked(d) || marker_diagonal.isMarked(d))
                {
                    int x = imageCoordinates[d].getXCoordinate();
                    int y = imageCoordinates[d].getYCoordinate();

                    int prediction;
                    int correction;

                    if(marker_vertical.isMarked(d))
                    {
                        prediction = (qRed(image.pixel(x/2, (y-1)/2)) + qRed(image.pixel(x/2, (y+1)/2)))/2;
                        correction = m_decomposition->getVerticalCorrection(x/2, y/2).getRed();
                    }
                    else
                    {
                        prediction = (qRed(image.pixel((x-1)/2, y/2)) + qRed(image.pixel((x+1)/2, y/2)))/2;
                        correction = m_decomposition->getHorizontalCorrection(x/2, y/2).getRed();
                    }

                    float z_coordinate = (m_camera->zFar()-m_camera->zNear())*(1.f-((prediction)/255.f));

                    PFP2::VEC3 plane_point_position = plane_center_position;
                    plane_point_position[0] += (x-(bb_min.x+bb_max.x)/2)/(bb_max.x-bb_min.x);
                    plane_point_position[1] += (y-(bb_min.y+bb_max.y)/2)/(bb_min.y-bb_max.y);

                    PFP2::VEC3 projected_point_position = plane_center_position;
                    projected_point_position[2] -= z_coordinate;

                    float d_C_Xp = (plane_point_position-camera_position).norm();
                    float d_C_Xh = (projected_point_position-camera_position).norm2();

                    float d_C_Xe = sqrt(d_C_Xh/d_C_Xc)*d_C_Xp;

                    position[d] = (plane_point_position-camera_position)/d_C_Xp*(d_C_Xe+d_C_Xp);
                }
            }

            mh_map->notifyAttributeModification(position);
            mh_map->notifyAttributeModification(imageCoordinates);
            mh_map->notifyConnectivityModification();
            mh_map->updateBB(position);

            m_decomposition = m_decomposition->getParent();

            m_schnapps->getSelectedView()->updateGL();
        }
    }
}

void Surface_WaveletDecomposition_Plugin::moveDownDecomposition(const QString& mapName)
{
    if(m_decomposition->getChild())
    {
        m_decomposition = m_decomposition->getChild();

        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

        if(mh_map)
        {
            PFP2::MAP* map = mh_map->getMap();
            TraversorF<PFP2::MAP> trav_face_map(*map);

            DartMarker<PFP2::MAP> marker(*map);

            marker.markAll();
        }
    }
}

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_Plugin, Surface_WaveletDecomposition_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_PluginD, Surface_WaveletDecomposition_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
