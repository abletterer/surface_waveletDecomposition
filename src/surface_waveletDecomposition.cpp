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
    connect(m_waveletDecompositionDialog->button_triangulate_map, SIGNAL(clicked()), this, SLOT(triangulateMapFromDialog()));

    connect(m_waveletDecompositionDialog->button_move_up, SIGNAL(clicked()), this, SLOT(moveUpFromDialog()));
    connect(m_waveletDecompositionDialog->button_move_down, SIGNAL(clicked()), this, SLOT(moveDownFromDialog()));

    m_decomposition = NULL;
    m_camera = NULL;
    m_matrix_coef = std::vector<int>();

    return true;
}

void Surface_WaveletDecomposition_Plugin::disable()
{
    disconnect(m_waveletDecompositionAction, SIGNAL(triggered()), this, SLOT(openWaveletDecompositionDialog()));

    disconnect(m_waveletDecompositionDialog->button_delete_background, SIGNAL(clicked()), this, SLOT(deleteBackgroundFromDialog()));
    disconnect(m_waveletDecompositionDialog->button_triangulate_map, SIGNAL(clicked()), this, SLOT(triangulateMapFromDialog()));

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

void Surface_WaveletDecomposition_Plugin::triangulateMapFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty())
    {
        triangulateMap(currentItems[0]->text());
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
            for(int i = 0; i < img_width; ++i)
            {
                for(int j = 0; j < img_height; ++j)
                {
                    if(i%2==0)
                    {
                        tmp_coef_matrix[i/2+width*j] = m_decomposition->getCoefficient(i, j);
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
                        }
                        else
                        {
                            result = m_decomposition->getCoefficient(i, j);
                            result -= left;
                            tmp_coef_matrix[img_width2+i/2+width*j] = result;
                        }
                    }
                }
            }

            m_decomposition->setCoefficientMatrix(tmp_coef_matrix);

            //Vertical decomposition
            for(int i = 0; i < img_width; ++i)
            {
                for(int j = 0; j < img_height; ++j)
                {
                    if(j%2 == 0)
                    {
                        tmp_coef_matrix[i+width*(j/2)] = m_decomposition->getCoefficient(i, j);
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
                        }
                        else
                        {
                            result = m_decomposition->getCoefficient(i, j);
                            result -= up;
                            tmp_coef_matrix[i+width*(img_height2+j/2)] = result;
                        }
                    }
                }
            }

            m_decomposition->setCoefficientMatrix(tmp_coef_matrix);

            img_width  = img_width2;
            img_height = img_height2;
            m_decomposition->moveDownDecomposition();
            if(img_width < 32 || img_height < 32)
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
                                //Use of the wavelet coefficients to correct the prediction
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
                                //Use of the wavelet coefficients to correct the prediction
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

void Surface_WaveletDecomposition_Plugin::saveDecompositions(const QString& name, const QString& directory)
{
    if(m_decomposition && !name.isEmpty() && !directory.isEmpty())
    {

        int img_width = m_decomposition->getWidth();
        int img_height = m_decomposition->getHeight();

        int level = m_decomposition->getLevel();

        int width = img_width/pow(2, level);
        int height = img_height/pow(2, level);

        CGoGNout << "Enregistrement du niveau " << level << CGoGNendl;

        QString filename(directory);
        filename.append("/");
        filename.append(name);
        filename.append("/");

        mkdir(filename.toStdString().c_str(), 0777);

        filename.append(name);

        if(level>0)
        {
            filename.append("-");
            filename.append(QString::number(width));
            filename.append("x");
            filename.append(QString::number(height));
            filename.append("/");

            mkdir(filename.toStdString().c_str(), 0777);

            filename.append(name);
            filename.append("-");
            filename.append(QString::number(width));
            filename.append("x");
            filename.append(QString::number(height));

            std::vector<int> values_ll;
            values_ll.resize(256, 0);

            //Print LL values
            for(int i = 0; i < width; ++i)
            {
                for(int j = 0; j < height; ++j)
                {
                    ++values_ll[m_matrix_coef[i+img_width*j]];
                }
            }

            std::vector<int> values_hl;
            values_hl.resize(511, 0);

            //Print HL coefficients
            for(int i = width; i < width*2; ++i)
            {
                for(int j = 0; j < height; ++j)
                {
                    ++values_hl[255+m_matrix_coef[i+img_width*j]];
                }
            }

            std::vector<int> values_lh;
            values_lh.resize(511, 0);

            //Print LH coefficients
            for(int i = 0; i < width; ++i)
            {
                for(int j = height; j < height*2; ++j)
                {
                    ++values_lh[255+m_matrix_coef[i+img_width*j]];
                }
            }

            std::vector<int> values_hh;
            values_hh.resize(511, 0);

            //Print HH coefficients
            for(int i = width; i < width*2; ++i)
            {
                for(int j = height; j < height*2; ++j)
                {
                    ++values_hh[255+m_matrix_coef[i+img_width*j]];
                }
            }

            CGoGNStream::Out file_ll, file_hl, file_lh, file_hh;
            QString filename2(filename);
            filename2.append("-LL.dat");
            file_ll.toFile(filename2.toStdString());
            filename2 = QString(filename);
            filename2.append("-HL.dat");
            file_hl.toFile(filename2.toStdString());
            filename2 = QString(filename);
            filename2.append("-LH.dat");
            file_lh.toFile(filename2.toStdString());
            filename2 = QString(filename);
            filename2.append("-HH.dat");
            file_hh.toFile(filename2.toStdString());

            file_ll.toStd(false);
            file_hl.toStd(false);
            file_lh.toStd(false);
            file_hh.toStd(false);

            for(int i = 0; i < 256; ++i)
            {
                file_ll << i << " " << values_ll[i]/(float)(width*height) << CGoGNendl;
            }

            file_ll.close();

            for(int i = 0; i < 511; ++i)
            {
                file_hl << i-255 << " " << values_hl[i]/(float)(width*height) << CGoGNendl;
            }

            file_hl.close();

            for(int i = 0; i < 511; ++i)
            {
                file_lh << i-255 << " " << values_lh[i]/(float)(width*height) << CGoGNendl;
            }

            file_lh.close();

            for(int i = 0; i < 511; ++i)
            {
                file_hh << i-255 << " " << values_hh[i]/(float)(width*height) << CGoGNendl;
            }

            file_hh.close();
        }
        else
        {
            std::vector<int> values;
            values.resize(256, 0);

            //Print HL coefficients
            for(int i = 0; i < width; ++i)
            {
                for(int j = 0; j < height; ++j)
                {
                    ++values[m_matrix_coef[i+img_width*j]];
                }
            }

            CGoGNStream::Out file;
            filename.append(".dat");
            file.toFile(filename.toStdString());
            file.toStd(false);

            for(int i = 0; i < 256; ++i)
            {
                file << i << " " << values[i]/(float)(width*height) << CGoGNendl;
            }

            file.close();
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

        int img_width = m_decomposition->getWidth(), img_height = m_decomposition->getHeight();
        int width = img_width/pow(2, m_decomposition->getLevel()), height = img_height/pow(2, m_decomposition->getLevel());

        Algo::Surface::Tilings::Square::Grid<PFP2> grid(*map, width-1, height-1);
//        grid.embedIntoGrid(planeCoordinates, img_width-1, img_height-1, 0, img_width-1, img_height-1);
        grid.embedIntoGrid(planeCoordinates,
                           img_width-1-pow(2, m_decomposition->getLevel()),
                           img_height-1-pow(2, m_decomposition->getLevel()),
                           0.f, img_width-1, img_height-1);

        std::vector<Dart> vDarts = grid.getVertexDarts();

        for(int i = 0; i < width; ++i)
        {
            for(int j = 0; j < height; ++j)
            {
                imageCoordinates[vDarts[j*width+i]].setCoordinates(i, height-j-1);
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

                int color = m_matrix_coef[imageCoordinates[d].getXCoordinate()+m_decomposition->getWidth()*imageCoordinates[d].getYCoordinate()];

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
                    if(color<128)
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

bool Surface_WaveletDecomposition_Plugin::moveUpDecomposition(const QString& mapName)
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

            map->clear(false);

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

            m_decomposition->moveUpDecomposition();

            width *= 2;
            height *= 2;

            Algo::Surface::Tilings::Square::Grid<PFP2> grid(*map, width-1, height-1);
//            grid.embedIntoGrid(planeCoordinates, img_width-1, img_height-1, 0, img_width-1, img_height-1);
            grid.embedIntoGrid(planeCoordinates,
                               img_width-1-pow(2, m_decomposition->getLevel()),
                               img_height-1-pow(2, m_decomposition->getLevel()),
                               0.f, img_width-1, img_height-1);

            std::vector<Dart> vDarts = grid.getVertexDarts();

            for(int i = 0; i < width; ++i)
            {
                for(int j = 0; j < height; ++j)
                {
                    imageCoordinates[vDarts[i+width*j]].setCoordinates(i, height-j-1);
                }
            }

            mh_map->notifyAttributeModification(planeCoordinates);
            mh_map->notifyAttributeModification(imageCoordinates);
            mh_map->notifyConnectivityModification();

            //triangulateMap(mapName);

            project2DImageTo3DSpace(mapName);
        }
    }
    else
    {
        return false;
    }
    return true;
}

bool Surface_WaveletDecomposition_Plugin::moveDownDecomposition(const QString& mapName)
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

            map->clear(false);

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
                        int result = coef_matrix2[i+img_width*j];
                        if(i != width-1)
                        {
                            int right = coef_matrix2[i+1+img_width*j];
                            result -= floor((left+right)/2.f + 1/2.f);
                            m_matrix_coef[width/2+i/2+img_width*j] = result;
                        }
                        else
                        {
                            result -= left;
                            m_matrix_coef[width/2+i/2+img_width*j] = result;
                        }
                    }
                }
            }

            coef_matrix2 = std::vector<int>(m_matrix_coef);

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
                        int result = coef_matrix2[i+img_width*j];
                        if(j != height-1)
                        {
                            int right = coef_matrix2[i+img_width*(j+1)];
                            result -= floor((left+right)/2.f + 1/2.f);
                            m_matrix_coef[i+img_width*(height/2+j/2)] = result;
                        }
                        else
                        {
                            result -= left;
                            m_matrix_coef[i+img_width*(height/2+j/2)] = result;
                        }
                    }
                }
            }

            m_decomposition->moveDownDecomposition();

            width /= 2;
            height /= 2;

            Algo::Surface::Tilings::Square::Grid<PFP2> grid(*map, width-1, height-1);
//            grid.embedIntoGrid(planeCoordinates, img_width-1, img_height-1, 0, img_width-1, img_height-1);
            grid.embedIntoGrid(planeCoordinates,
                               img_width-1-pow(2, m_decomposition->getLevel()),
                               img_height-1-pow(2, m_decomposition->getLevel()),
                               0.f, img_width-1, img_height-1);

            std::vector<Dart> vDarts = grid.getVertexDarts();

            for(int i = 0; i < width; ++i)
            {
                for(int j = 0; j < height; ++j)
                {
                    imageCoordinates[vDarts[i+width*j]].setCoordinates(i, height-j-1);
                }
            }

            mh_map->notifyAttributeModification(planeCoordinates);
            mh_map->notifyAttributeModification(imageCoordinates);
            mh_map->notifyConnectivityModification();

            //triangulateMap(mapName);

            project2DImageTo3DSpace(mapName);
        }
    }
    else
    {
        return false;
    }
    return true;
}

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_Plugin, Surface_WaveletDecomposition_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_PluginD, Surface_WaveletDecomposition_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
