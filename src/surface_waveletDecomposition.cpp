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

    connect(m_waveletDecompositionDialog->button_delete_background, SIGNAL(clicked()), this, SLOT(deleteBackgroundFromDialog()));
    connect(m_waveletDecompositionDialog->button_move_up, SIGNAL(clicked()), this, SLOT(moveUpFromDialog()));
    connect(m_waveletDecompositionDialog->button_move_down, SIGNAL(clicked()), this, SLOT(moveDownFromDialog()));

    m_decomposition = NULL;
    m_camera = NULL;
    m_matrix_coef = std::vector<int>();
    m_drawer = new Utils::Drawer();

    return true;
}

void Surface_WaveletDecomposition_Plugin::disable()
{
    disconnect(m_waveletDecompositionAction, SIGNAL(triggered()), this, SLOT(openWaveletDecompositionDialog()));

    disconnect(m_waveletDecompositionDialog->button_delete_background, SIGNAL(clicked()), this, SLOT(deleteBackgroundFromDialog()));
    disconnect(m_waveletDecompositionDialog->button_move_up, SIGNAL(clicked()), this, SLOT(moveUpFromDialog()));
    disconnect(m_waveletDecompositionDialog->button_move_down, SIGNAL(clicked()), this, SLOT(moveDownFromDialog()));

    if(m_decomposition)
    {
        delete m_decomposition;
    }
    if(m_camera)
    {
        delete m_camera;
    }
    if(m_drawer)
    {
        delete m_drawer;
    }
}

void Surface_WaveletDecomposition_Plugin::draw(View* view)
{
    m_drawer->callList();
}

void Surface_WaveletDecomposition_Plugin::openWaveletDecompositionDialog()
{
    m_waveletDecompositionDialog->show();
}

void Surface_WaveletDecomposition_Plugin::closeWaveletDecompositionDialog()
{
    m_waveletDecompositionDialog->close();
}

void Surface_WaveletDecomposition_Plugin::deleteBackgroundFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty())
    {
        deleteBackground(currentItems[0]->text());
    }
}

void Surface_WaveletDecomposition_Plugin::moveUpFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty())
    {
        moveUpDecomposition(currentItems[0]->text());
    }
}

void Surface_WaveletDecomposition_Plugin::moveDownFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty())
    {
        moveDownDecomposition(currentItems[0]->text());
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
                m_decomposition->setCoefficient(i, height-j-1, qRed(image.pixel(i, height-j-1)));
            }
        }

        return file;
    }
    return NULL;
}

void Surface_WaveletDecomposition_Plugin::decompose(const int max_counter)
{
    if(m_decomposition)
    {
        bool stop = false;
        int width = m_decomposition->getWidth(), height = m_decomposition->getHeight();
        int img_width = width, img_height = height;
        int img_width2, img_height2;
        int counter = 0;
        while(!stop && (max_counter == -1 || counter < max_counter))
        {
            img_width2  = round(img_width/2.f);
            img_height2 = round(img_height/2.f);

            //Horizontal decomposition
            std::vector<int> tmp_coef_matrix(*m_decomposition->getCoefficientMatrix());
            std::vector<int> tmp_diff_matrix(*m_decomposition->getDifferenceMatrix());
            for(int i = 0; i < img_width; ++i)
            {
                for(int j = 0; j < img_height; ++j)
                {
                    if(i%2==0)
                    {
                        tmp_coef_matrix[i/2+width*j] = m_decomposition->getCoefficient(i, j);
                        tmp_diff_matrix[i/2+width*j] = m_decomposition->getDifference(i, j);
                    }
                    else
                    {
                        int left = m_decomposition->getCoefficient(i-1, j);
                        int result;
                        if(i != img_width-1)
                        {
                            int right = m_decomposition->getCoefficient(i+1, j);
                            result = m_decomposition->getCoefficient(i, j);
                            result -= floor((left+right)/2.f + 1/2.f);
                            tmp_coef_matrix[img_width2+i/2+width*j] = result;
                            tmp_diff_matrix[img_width2+i/2+width*j] = left-right;
                        }
                        else
                        {
                            result = m_decomposition->getCoefficient(i, j);
                            result -= left;
                            tmp_coef_matrix[img_width2+i/2+width*j] = result;
                            tmp_diff_matrix[img_width2+i/2+width*j] = 0;
                        }
                    }
                }
            }

            m_decomposition->setCoefficientMatrix(tmp_coef_matrix);
            m_decomposition->setDifferenceMatrix(tmp_diff_matrix);

            //Vertical decomposition
            for(int i = 0; i < img_width; ++i)
            {
                for(int j = 0; j < img_height; ++j)
                {
                    if(j%2 == 0)
                    {
                        tmp_coef_matrix[i+width*(j/2)] = m_decomposition->getCoefficient(i, j);
                        tmp_diff_matrix[i+width*(j/2)] = m_decomposition->getDifference(i, j);
                    }
                    else
                    {
                        int up = m_decomposition->getCoefficient(i, j-1);
                        int result;
                        if(j != img_height-1)
                        {
                            int down = m_decomposition->getCoefficient(i, j+1);
                            result = m_decomposition->getCoefficient(i, j);
                            result -= floor((up+down)/2.f + 1/2.f);
                            tmp_coef_matrix[i+width*(img_height2+j/2)] = result;
                            tmp_diff_matrix[i+width*(img_height2+j/2)] = up-down;
                        }
                        else
                        {
                            result = m_decomposition->getCoefficient(i, j);
                            result -= up;
                            tmp_coef_matrix[i+width*(img_height2+j/2)] = result;
                            tmp_diff_matrix[i+width*(img_height2+j/2)] = 0;
                        }
                    }
                }
            }

            m_decomposition->setCoefficientMatrix(tmp_coef_matrix);
            m_decomposition->setDifferenceMatrix(tmp_diff_matrix);

            img_width  = img_width2;
            img_height = img_height2;
            m_decomposition->moveDownDecomposition();
            if(img_width < 3 || img_height < 3)
            {
                stop = true;
                m_decomposition->setMaxLevel(m_decomposition->getLevel());
                CGoGNout << m_decomposition->getLevel()+1 << " niveau(x) de décomposition" << CGoGNendl;
            }
            ++counter;
        }

        if(m_matrix_coef.empty())
        {
            m_matrix_coef = std::vector<int>(*m_decomposition->getCoefficientMatrix());
        }

//        updateDrawer();
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
                int color = m_decomposition->getCoefficient(i, j);
                image.setPixel(i, j, qRgb(qAbs(color), qAbs(color), qAbs(color)));
            }
        }

        if(!image.save(filename))
        {
            CGoGNerr << "Image '" << filename.toStdString() << "' has not been saved" << CGoGNendl;
        }
    }
}

void Surface_WaveletDecomposition_Plugin::saveAllImages(const QString& name, const QString& directory)
{
    if(m_decomposition && !name.isEmpty() && !directory.isEmpty())
    {
        int img_width = m_decomposition->getWidth();
        int img_height = m_decomposition->getHeight();

        std::vector<int> matrix(m_matrix_coef);
        std::vector<int> matrix2(matrix);

        int level = m_decomposition->getLevel();
        int counter = level;

        QString filename(directory);
        filename.append(name);
        filename.append("-Reconstructions-");

        while(counter >= 0)
        {
            matrix = std::vector<int>(m_matrix_coef);
            int width = img_width/pow(2, level);
            int height = img_height/pow(2, level);

            bool use_coef = true;

            while(width != img_width && height != img_height)
            {
                if(width == img_width/pow(2, counter) || height == img_height/pow(2, counter))
                {
                    use_coef = false;
                }

                matrix2 = std::vector<int>(matrix);

                for(int i = 0; i < width*2; ++i)
                {
                    for(int j = 0; j < height*2; ++j)
                    {
                        int index = i+img_width*j;
                        if(j%2 == 1)
                        {
                            int up = matrix2[i+img_width*(j/2)];
                            int result = 0;
                            if(use_coef)
                            {
                                result= matrix2[i+img_width*(height+j/2)];
                            }
                            if(j != height*2-1)
                            {
                                int down = matrix2[i+img_width*(j/2+1)];
                                result += floor((up+down)/2.f + 1/2.f);
                            }
                            else
                            {
                                result += up;
                            }
                            matrix[index] = result;
                        }
                        else
                        {
                            matrix[index] = matrix2[i+img_width*(j/2)];
                        }
                    }
                }

                matrix2 = std::vector<int>(matrix);

                for(int i = 0; i < width*2; ++i)
                {
                    for(int j = 0; j < height*2 ; ++j)
                    {
                        int index = i+img_width*j;
                        if(i%2 == 1)
                        {
                            int left = matrix2[i/2+img_width*j];
                            int result = 0;
                            if(use_coef)
                            {
                                //Ute of the wavelet coefficients to correct the prediction
                                result = matrix2[width+i/2+img_width*j];
                            }
                            if(i != width*2-1)
                            {
                                int right = matrix2[i/2+1+img_width*j];
                                result += floor((left+right)/2.f + 1/2.f);
                            }
                            else
                            {
                                result += left;
                            }
                            matrix[index] = result;
                        }
                        else
                        {
                            matrix[index] = matrix2[i/2+img_width*j];
                        }
                    }
                }
                width *= 2;
                height *= 2;
            }

            QString filename2(filename);

            filename2.append(QString::number(counter));
            filename2.append(".png");

            QImage image(width, height, QImage::Format_RGB32);

            for(int i = 0; i < width; ++i)
            {
                for(int j = 0; j < height; ++j)
                {
                    int color = matrix[i + img_width*j];
                    image.setPixel(i, j, qRgb(qAbs(color), qAbs(color), qAbs(color)));
                }
            }

            if(!image.save(filename2))
            {
                CGoGNerr << "Image '" << filename2.toStdString() << "' has not been saved" << CGoGNendl;
            }

            --counter;
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

                int color = m_decomposition->getCoefficient(imageCoordinates[d].getXCoordinate(), imageCoordinates[d].getYCoordinate());

                float z_coordinate = (m_camera->zFar()-m_camera->zNear())*(1.f-(qAbs(color)/255.f))*200;

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

void Surface_WaveletDecomposition_Plugin::projectNewPointsTo3DSpace(MapHandler<PFP2>* mh_map, const std::vector<Dart>& vertices)
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

            int color = m_matrix_coef[imageCoordinates[*d].getXCoordinate()+m_decomposition->getWidth()*imageCoordinates[*d].getYCoordinate()];

            float z_coordinate = (m_camera->zFar()-m_camera->zNear())*(1.f-(qAbs(color)/255.f))*200;

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
                    Dart phi_11_d = map->phi<11>(d);
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

void Surface_WaveletDecomposition_Plugin::deleteBackground(const QString& mapName)
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

            TraversorF<PFP2::MAP> trav_face_map(*map);
            Dart next;
            bool stop = false;
            for(Dart d = trav_face_map.begin(); d != trav_face_map.end(); d = next)
            {
                next = trav_face_map.next();
                stop = false;
                Traversor2FV<PFP2::MAP> trav_vert_face_map(*map, d);
                for(Dart dd = trav_vert_face_map.begin(); !stop && dd != trav_vert_face_map.end(); dd = trav_vert_face_map.next())
                {
                    int color = m_matrix_coef[imageCoordinates[dd].getXCoordinate()+m_decomposition->getWidth()*imageCoordinates[dd].getYCoordinate()];
                    if(color==0)
                    {
                        map->deleteFace(d);
                        stop = true;
                    }
                }
            }

            mh_map->notifyConnectivityModification();
            mh_map->notifyAttributeModification(position);
            mh_map->notifyAttributeModification(planeCoordinates);
            mh_map->notifyAttributeModification(imageCoordinates);

            mh_map->updateBB(position);

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

            int img_width = m_decomposition->getWidth();
            int img_height = m_decomposition->getHeight();

            int width = img_width/pow(2, m_decomposition->getLevel());
            int height = img_height/pow(2, m_decomposition->getLevel());

            DartMarker<PFP2::MAP> marker(*map);
            DartMarker<PFP2::MAP> marker_face(*map);
            DartMarker<PFP2::MAP> marker_horizontal(*map), marker_vertical(*map), marker_diagonal(*map);
            DartMarker<PFP2::MAP> marker_modified(*map);

            std::vector<Dart> vertices_added;
            vertices_added.reserve(width*height*4);

            std::vector<Dart> vertices_modified;
            vertices_modified.reserve(width*height);

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
                            vertices_added.push_back(ddd);
                            if(!marker_modified.isMarked(dd))
                            {
                                marker_modified.markOrbit<VERTEX>(dd);
                                vertices_modified.push_back(dd);
                            }
                            planeCoordinates[ddd] = (planeCoordinates[dd]+planeCoordinates[dd1])/2.f;
                        }
                    }
                }
            }

            for(std::vector<Dart>::const_iterator d = vertices_added.begin(); d != vertices_added.end(); ++d)
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
                    Dart d_1 = map->phi_1(*d), d1 = map->phi1(*d);
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

//            Dart lower_left_vertex;

            for(std::vector<Dart>::const_iterator d = vertices_modified.begin(); d != vertices_modified.end(); ++d)
            {
//                if(imageCoordinates[*d].getXCoordinate() == 0 && imageCoordinates[*d].getYCoordinate() == height-1)
//                {
//                    lower_left_vertex = *d;
//                }
                imageCoordinates[*d].setCoordinates(imageCoordinates[*d].getXCoordinate()*2, imageCoordinates[*d].getYCoordinate()*2);
            }

            //Adding exterior points
//            Dart d = lower_left_vertex, d1, d2;

//            std::vector<Dart> faces_added;
//            faces_added.reserve(width+height);

//            while(!map->isBoundaryEdge(d))
//            {
//                d = map->phi<21>(d);
//            }

//            d1 = map->newFace(4);
//            map->sewFaces(d, d1);
//            faces_added.push_back(d1);
//            d1 = map->phi_1(d1);

//            planeCoordinates[d1] = planeCoordinates[map->phi1(d)] + planeCoordinates[map->phi1(d)]-planeCoordinates[map->phi_1(d)];
//            imageCoordinates[d1].setCoordinates(imageCoordinates[map->phi1(d)].getXCoordinate(), imageCoordinates[map->phi1(d)].getYCoordinate()+1);
//            vertices_added.push_back(d1);

//            planeCoordinates[map->phi_1(d1)] = planeCoordinates[d] + planeCoordinates[d]-planeCoordinates[map->phi_1(map->phi2(map->phi_1(d)))];
//            imageCoordinates[map->phi_1(d1)].setCoordinates(imageCoordinates[map->phi1(d)].getXCoordinate(), imageCoordinates[map->phi1(d)].getYCoordinate()+1);
//            vertices_added.push_back(map->phi_1(d1));

//            bool stop = false, up = false;
//            while(!stop)
//            {
//                d = map->phi<12121>(d);

//                d2 = map->newFace(4);
//                map->sewFaces(d, d2);
//                map->sewFaces(d1, map->phi1(d2));
//                faces_added.push_back(d2);
//                d1 = map->phi_1(d2);

//                planeCoordinates[d1] = planeCoordinates[map->phi1(d)] + planeCoordinates[map->phi1(d)]-planeCoordinates[map->phi_1(d)];
//                if(up)
//                {
//                    imageCoordinates[d1].setCoordinates(imageCoordinates[map->phi1(d)].getXCoordinate()+1, imageCoordinates[map->phi1(d)].getYCoordinate());
//                }
//                else
//                {
//                    imageCoordinates[d1].setCoordinates(imageCoordinates[map->phi1(d)].getXCoordinate(), imageCoordinates[map->phi1(d)].getYCoordinate()+1);
//                }
//                vertices_added.push_back(d1);

//                if(imageCoordinates[map->phi1(d)].getYCoordinate()==0)
//                {
//                    stop = true;
//                }
//                else if(!up && imageCoordinates[map->phi1(d)].getXCoordinate() == width*2-2 && imageCoordinates[map->phi1(d)].getYCoordinate() == height*2-2)
//                {
//                    d = map->phi1(d);

//                    d2 = map->newFace(4);
//                    map->sewFaces(d1, d2);
//                    faces_added.push_back(d2);
//                    d1 = map->phi_1(d2);

//                    d2 = map->newFace(4);
//                    map->sewFaces(d1, d2);
//                    map->sewFaces(d, map->phi_1(d2));
//                    faces_added.push_back(d2);

//                    planeCoordinates[d1] = planeCoordinates[d] + planeCoordinates[d]-planeCoordinates[map->phi_1(d)];
//                    imageCoordinates[d1].setCoordinates(imageCoordinates[d].getXCoordinate()+1, imageCoordinates[d].getYCoordinate());
//                    vertices_added.push_back(d1);

//                    planeCoordinates[map->phi_1(d1)] = planeCoordinates[d] + planeCoordinates[d]-planeCoordinates[map->phi_1(map->phi2(map->phi1(d)))];
//                    imageCoordinates[map->phi_1(d1)].setCoordinates(imageCoordinates[d].getXCoordinate()+1, imageCoordinates[d].getYCoordinate()+1);
//                    vertices_added.push_back(map->phi_1(d1));

//                    d1 = map->phi<11>(d2);

//                    planeCoordinates[d1] = planeCoordinates[map->phi1(d1)] + planeCoordinates[map->phi1(d1)]-planeCoordinates[map->phi_1(map->phi2(map->phi1(d)))];
//                    imageCoordinates[d1].setCoordinates(imageCoordinates[map->phi1(d1)].getXCoordinate()+1, imageCoordinates[d].getYCoordinate());
//                    vertices_added.push_back(d1);
//                    up = true;
//                }
//            }

//            for(std::vector<Dart>::const_iterator d = faces_added.begin(); d != faces_added.end(); ++d)
//            {
//                Dart phi_11_d = map->phi<11>(*d);
//                map->splitFace(*d, phi_11_d);
//            }

            if(m_matrix_coef.empty())
            {
                m_matrix_coef = std::vector<int>(*m_decomposition->getCoefficientMatrix());
            }

            std::vector<int> coef_matrix2 = std::vector<int>(m_matrix_coef);

            for(int i = 0; i < width*2; ++i)
            {
                for(int j = 0; j < height*2; ++j)
                {
                    int index = i+img_width*j;
                    if(j%2 == 1)
                    {
                        int up = coef_matrix2[i+img_width*(j/2)];
                        int result = coef_matrix2[i+img_width*(height+j/2)];
                        if(j != height*2-1)
                        {
                            int down = coef_matrix2[i+img_width*(j/2+1)];
                            result += floor((up+down)/2.f + 1/2.f);
                        }
                        else
                        {
                            result += up;
                        }
                        m_matrix_coef[index] = result;
                    }
                    else
                    {
                        m_matrix_coef[index] = coef_matrix2[i+img_width*(j/2)];
                    }
                }
            }

            coef_matrix2 = std::vector<int>(m_matrix_coef);

            for(int i = 0; i < width*2; ++i)
            {
                for(int j = 0; j < height*2 ; ++j)
                {
                    int index = i+img_width*j;
                    if(i%2 == 1)
                    {
                        int left = coef_matrix2[i/2+img_width*j];
                        int result = coef_matrix2[width+i/2+img_width*j];
                        if(i != width*2-1)
                        {
                            int right = coef_matrix2[i/2+1+img_width*j];
                            result += floor((left+right)/2.f + 1/2.f);
                        }
                        else
                        {
                            result += left;
                        }
                        m_matrix_coef[index] = result;
                    }
                    else
                    {
                        m_matrix_coef[index] = coef_matrix2[i/2+img_width*j];
                    }
                }
            }

            QImage image(width*2, height*2, QImage::Format_RGB32);

            for(int i = 0; i < width*2; ++i)
            {
                for(int j = 0; j < height*2 ; ++j)
                {
                    int value = m_matrix_coef[i+m_decomposition->getWidth()*j];
                    image.setPixel(i, j, qRgb(qAbs(value), qAbs(value), qAbs(value)));
                }
            }

            QString filename("/home/blettere/Projets/Models/Test/");
            filename.append(mapName);
            filename.append(".png");

            if(!image.save(filename))
            {
                CGoGNerr << "Image '" << filename.toStdString() << "' has not been saved" << CGoGNendl;
            }

            m_decomposition->moveUpDecomposition();

            mh_map->notifyAttributeModification(planeCoordinates);
            mh_map->notifyAttributeModification(imageCoordinates);
            mh_map->notifyConnectivityModification();

            projectNewPointsTo3DSpace(mh_map, vertices_added);
        }
    }
}

void Surface_WaveletDecomposition_Plugin::moveDownDecomposition(const QString& mapName)
{
    if(m_decomposition && m_decomposition->getLevel() < m_decomposition->getMaxLevel())
    {
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));
        if(mh_map)
        {
            PFP2::MAP* map = mh_map->getMap();

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

            int img_width = m_decomposition->getWidth();
            int img_height = m_decomposition->getHeight();

            int width = img_width/pow(2, m_decomposition->getLevel());
            int height = img_height/pow(2, m_decomposition->getLevel());

            Dart starting_vertex;
            bool stop = false;

            TraversorV<PFP2::MAP> trav_vert_map(*map);
            for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end() && !stop; d = trav_vert_map.next())
            {
                if(imageCoordinates[d].getXCoordinate()==0 && imageCoordinates[d].getYCoordinate()==0)
                {
                    starting_vertex = d;
                    stop = true;
                }
            }

            Dart d = starting_vertex;

            while(imageCoordinates[d].getYCoordinate() != height-2)
            {
                while(!map->isBoundaryMarked<2>(d))
                {
                    Dart dd = map->deleteVertex(map->phi_1(map->phi<12>(d)));
                    map->mergeFaces(map->phi<111>(dd));
                    map->mergeFaces(dd);
                    d = map->phi<111112>(d);
                }
                d = map->phi2(d);
                while(!map->isBoundaryMarked<2>(d))
                {
                    d = map->phi<1112>(d);
                }
                d = map->phi<21121>(d);
            }

            DartMarker<PFP2::MAP> marker_edge(*map);

            for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
            {
                if(!marker_edge.isMarked(d))
                {
                    if(!map->uncutEdge(d))
                    {
                        d = map->phi_1(d);
                        map->uncutEdge(d);
                    }
                    marker_edge.markOrbit<EDGE>(d);
                }
            }

            std::vector<int> coef_matrix2 = std::vector<int>(m_matrix_coef);

            for(int i = 0; i < width; ++i)
            {
                for(int j = 0; j < height; ++j)
                {
                    if(i%2==0)
                    {
                        m_matrix_coef[i/2+img_width*j] = coef_matrix2[i+img_width*j];
                    }
                    else
                    {
                        int left = coef_matrix2[i-1+img_width*j];
                        int result;
                        if(i != width-1)
                        {
                            int right = coef_matrix2[i+1+img_width*j];
                            result = coef_matrix2[i+img_width*j];
                            result -= floor((left+right)/2.f + 1/2.f);
                            m_matrix_coef[width/2+i/2+img_width*j] = result;
                        }
                        else
                        {
                            result = coef_matrix2[i+img_width*j];
                            result -= left;
                            m_matrix_coef[width/2+i/2+img_width*j] = result;
                        }
                    }
                }
            }

            for(int i = 0; i < width; ++i)
            {
                for(int j = 0; j < height; ++j)
                {
                    if(j%2==0)
                    {
                        m_matrix_coef[i+img_width*(j/2)] = coef_matrix2[i+img_width*j];
                    }
                    else
                    {
                        int left = coef_matrix2[i+img_width*(j-1)];
                        int result;
                        if(j != height-1)
                        {
                            int right = coef_matrix2[i+img_width*(j+1)];
                            result = coef_matrix2[i+img_width*j];
                            result -= floor((left+right)/2.f + 1/2.f);
                            m_matrix_coef[i+img_width*(height/2+j/2)] = result;
                        }
                        else
                        {
                            result = coef_matrix2[i+img_width*j];
                            result -= left;
                            m_matrix_coef[i+img_width*(height/2+j/2)] = result;
                        }
                    }
                }
            }

            //triangulateMap(mapName);

            m_decomposition->moveDownDecomposition();

            mh_map->notifyConnectivityModification();
            mh_map->notifyAttributeModification(planeCoordinates);
            mh_map->notifyAttributeModification(imageCoordinates);
            m_schnapps->getSelectedView()->updateGL();

        }
    }
}

void Surface_WaveletDecomposition_Plugin::updateDrawer()
{
    if(m_decomposition)
    {
        int img_width = m_decomposition->getWidth();
        int img_height = m_decomposition->getHeight();

        int diff_level = pow(2, 0);

        int width = img_width/diff_level;
        int height = img_height/diff_level;

        m_drawer->newList(GL_COMPILE_AND_EXECUTE);
            m_drawer->begin(GL_POINTS);
                m_drawer->color3f(1.f, 1.f, 0.f);
                m_drawer->pointSize(10.f);

                for(int i = 0; i < width; ++i)
                {
                    for(int j = 0; j < height; ++j)
                    {
                        m_drawer->vertex3f(i*diff_level, j*diff_level, 0.f);
                    }
                }

            m_drawer->end();
        m_drawer->endList();

        m_schnapps->getSelectedView()->getCurrentCamera()->setSceneCenter(qglviewer::Vec(0.f, 0.f, 0.f));
        m_schnapps->getSelectedView()->updateGL();
    }
}

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_Plugin, Surface_WaveletDecomposition_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_PluginD, Surface_WaveletDecomposition_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
