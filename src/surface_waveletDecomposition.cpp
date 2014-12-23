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
        int width = m_decomposition->getWidth();
        int height = m_decomposition->getHeight();
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
                            result.setMean(left, right);
                            result -= m_decomposition->getValue(i, j);
                            tmp_matrix[img_width2-1+round(i/2.f)+width*j] = result;
                        }
                        else
                        {
                            result = left;
                            result -= m_decomposition->getValue(i, j);
                            tmp_matrix[img_width2-1+round(i/2.f)+width*j] = result;
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
                            result.setMean(up, down);
                            result -= m_decomposition->getValue(i, j);
                            tmp_matrix[i+width*(img_height2-1+round(j/2.f))] = result;
                        }
                        else
                        {
                            result = up;
                            result -= m_decomposition->getValue(i, j);
                            tmp_matrix[i+width*(img_height2-1+round(j/2.f))] = result;
                        }
                    }
                }
            }

            m_decomposition->setMatrix(tmp_matrix);

            img_width  = img_width2;
            img_height = img_height2;
            m_decomposition->getDownDecomposition();
            if(img_width < 2 || img_height < 2)
            {
                stop = true;
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
//    if(m_decomposition)
//    {
//        MapHandlerGen* mhg_map = m_schnapps->addMap(mapName, 2);
//        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);
//        PFP2::MAP* map = mh_map->getMap();

//        VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
//        if(!position.isValid())
//        {
//            position = mh_map->addAttribute<PFP2::VEC3, VERTEX>("position");
//        }

//        VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_map->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
//        if(!imageCoordinates.isValid())
//        {
//            imageCoordinates = mh_map->addAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
//        }

//        Decomposition* decomposition = m_decomposition;
////        while(decomposition->getChild())
////        {
////            decomposition = decomposition->getChild();
////        }

//        QImage image = decomposition->getImage();
//        int imageX = image.width(), imageY = image.height();

//        Algo::Surface::Tilings::Square::Grid<PFP2> grid(*map, imageX-1, imageY-1);
//        grid.embedIntoGrid(position, imageX-1, imageY-1);

//        std::vector<Dart> vDarts = grid.getVertexDarts();

//        mh_map->updateBB(position);

//        qglviewer::Vec bb_max = mh_map->getBBmax();
//        qglviewer::Vec bb_min = mh_map->getBBmin();

//        for(int i = 0; i < imageX; ++i)
//        {
//            for(int j = 0; j < imageY; ++j)
//            {
//                imageCoordinates[vDarts[j*imageX+i]].setCoordinates(i,imageY-j-1);
//            }
//        }

//        QImage image_parent = m_decomposition->getImage();
//        int image_parentX = image_parent.width(), image_parentY = image_parent.height();

//        float width_step_parent = (bb_max.x-bb_min.x)/(image_parentX-1);
//        float height_step_parent = (bb_max.y-bb_min.y)/(image_parentY-1);

//        int shift_counterX = 0, shift_counterY = 0;

//        while(image_parentX != imageX && image_parentY != imageY)
//        {
//            if(image_parentX%2==0)
//            {
//                ++shift_counterX;
//            }
//            if(image_parentY%2==0)
//            {
//                ++shift_counterY;
//            }

//            image_parentX = image_parentX/2+image_parentX%2;
//            image_parentY = image_parentY/2+image_parentY%2;
//        }

//        shift_counterX = pow(2, shift_counterX);
//        shift_counterY = pow(2, shift_counterY);

//        PFP2::MATRIX44 transform_matrix;
//        transform_matrix.identity();
//        transform_matrix.setSubVectorV(0, 3, PFP2::VEC4(-width_step_parent*shift_counterX/2., height_step_parent*shift_counterY/2., 0., 1.));
//        grid.transform(position, transform_matrix);

//        mh_map->notifyAttributeModification(position);
//        mh_map->notifyAttributeModification(imageCoordinates);
//        mh_map->notifyConnectivityModification();
//        mh_map->updateBB(position);

//        m_schnapps->getSelectedView()->updateGL();

//        m_decomposition = decomposition;

//        return mhg_map;
//    }
    return NULL;
}

void Surface_WaveletDecomposition_Plugin::project2DImageTo3DSpace(const QString& mapName)
{
//    if(m_decomposition)
//    {
//        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

//        if(mh_map)
//        {
//            PFP2::MAP* map = mh_map->getMap();

//            if(!m_camera)
//            {
//                m_camera = new qglviewer::Camera(*(m_schnapps->getSelectedView()->camera()));
//                m_camera->setZNearCoefficient(1.f);
//            }

//            VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
//            if(!position.isValid())
//            {
//                CGoGNerr << "Position attribute is not valid" << CGoGNendl;
//            }

//            VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_map->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
//            if(!imageCoordinates.isValid())
//            {
//                CGoGNerr << "ImageCoordinates attribute is not valid" << CGoGNendl;
//            }

//            QImage image = m_decomposition->getImage();

//            qglviewer::Vec camera_p = m_camera->position();

//            PFP2::VEC3 camera_position = PFP2::VEC3(camera_p.x, camera_p.y, camera_p.z);

//            PFP2::VEC3 plane_center_position = camera_position;
//            plane_center_position[2] -= m_camera->zNear();

//            float d_C_Xc = (plane_center_position-camera_position).norm2();

//            //Reprojection des points => Théorème de Thalès
//            TraversorV<PFP2::MAP> trav_vert_map(*map);
//            for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
//            {
//                PFP2::VEC3 plane_point_position = plane_center_position;
//                plane_point_position[0] += position[d][0];
//                plane_point_position[1] += position[d][1];

//                QRgb pixel = image.pixel(imageCoordinates[d].getXCoordinate(), imageCoordinates[d].getYCoordinate());

//                float z_coordinate = (m_camera->zFar()-m_camera->zNear())*(1.f-(qRed(pixel)/255.f));

//                PFP2::VEC3 projected_point_position = plane_center_position;
//                projected_point_position[2] -= z_coordinate;

//                float d_C_Xp = (plane_point_position-camera_position).norm();
//                float d_C_Xh = (projected_point_position-camera_position).norm2();

//                float d_C_Xe = sqrt(d_C_Xh/d_C_Xc)*d_C_Xp;

//                position[d] = (plane_point_position-camera_position)/d_C_Xp*(d_C_Xe+d_C_Xp);
//            }

//            mh_map->notifyAttributeModification(position);
//            mh_map->notifyConnectivityModification();
//            mh_map->updateBB(position);

//            //camera->setSceneBoundingBox(mh_map->getBBmin(), mh_map->getBBmax());

//            m_schnapps->getSelectedView()->camera()->setSceneBoundingBox(mh_map->getBBmin(), mh_map->getBBmax());

//            m_schnapps->getSelectedView()->updateGL();
//        }
//    }
}

void Surface_WaveletDecomposition_Plugin::triangulateMap(const QString& mapName)
{
//    if(m_decomposition)
//    {
//        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

//        if(mh_map)
//        {
//            PFP2::MAP* map = mh_map->getMap();
//            VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
//            if(!position.isValid())
//            {
//                CGoGNerr << "Position attribute is not valid" << CGoGNendl;
//            }

//            DartMarker<PFP2::MAP> marker(*map);

//            TraversorF<PFP2::MAP> trav_face_map(*map);
//            for(Dart d = trav_face_map.begin(); d != trav_face_map.end(); d = trav_face_map.next())
//            {
//                if(!marker.isMarked(d))
//                {
//                    Dart phi_11_d = map->phi<11>(d);
//                    map->splitFace(d, phi_11_d);
//                    marker.markOrbit<FACE>(d);
//                    marker.markOrbit<FACE>(phi_11_d);
//                }
//            }

//            mh_map->notifyAttributeModification(position);
//            mh_map->notifyConnectivityModification();

//            m_schnapps->getSelectedView()->updateGL();
//        }
//    }
}

void Surface_WaveletDecomposition_Plugin::moveUpDecomposition(const QString& mapName)
{
    if(m_decomposition)
    {
    }
}

void Surface_WaveletDecomposition_Plugin::moveDownDecomposition(const QString& mapName)
{
    if(m_decomposition)
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
