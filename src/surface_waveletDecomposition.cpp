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

    m_colorPerVertexShader = new CGoGN::Utils::ShaderColorPerVertex();
    registerShader(m_colorPerVertexShader);

    m_positionVBO = new Utils::VBO();
    m_colorVBO = new Utils::VBO();

    m_toDraw = false;

    connect(m_waveletDecompositionAction, SIGNAL(triggered()), this, SLOT(openWaveletDecompositionDialog()));

    connect(m_waveletDecompositionDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeWaveletDecompositionDialog()));
    connect(m_waveletDecompositionDialog->button_synthesize, SIGNAL(clicked()), this, SLOT(synthesizeFromDialog()));
    connect(m_waveletDecompositionDialog->button_analyze, SIGNAL(clicked()), this, SLOT(analyzeFromDialog()));
    connect(m_waveletDecompositionDialog->button_saveImages, SIGNAL(clicked()), this, SLOT(saveImagesFromDialog()));
    connect(m_waveletDecompositionDialog->button_completeAnalysis, SIGNAL(clicked()), this, SLOT(completeAnalysisFromDialog()));

    return true;
}

void Surface_WaveletDecomposition_Plugin::disable()
{
    delete m_colorPerVertexShader;
    delete m_positionVBO;
    delete m_colorVBO;

    disconnect(m_waveletDecompositionAction, SIGNAL(triggered()), this, SLOT(openWaveletDecompositionDialog()));

    disconnect(m_waveletDecompositionDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeWaveletDecompositionDialog()));
    disconnect(m_waveletDecompositionDialog->button_synthesize, SIGNAL(clicked()), this, SLOT(synthesizeFromDialog()));
    disconnect(m_waveletDecompositionDialog->button_analyze, SIGNAL(clicked()), this, SLOT(analyzeFromDialog()));
    disconnect(m_waveletDecompositionDialog->button_saveImages, SIGNAL(clicked()), this, SLOT(saveImagesFromDialog()));
    disconnect(m_waveletDecompositionDialog->button_completeAnalysis, SIGNAL(clicked()), this, SLOT(completeAnalysisFromDialog()));
}

void Surface_WaveletDecomposition_Plugin::drawMap(View* view, MapHandlerGen *map)
{
    if(m_toDraw && m_schnapps->getSelectedView() == view)
    {
        //If VBO are initialized
        glPolygonMode(GL_FRONT, GL_FILL);
        glEnable(GL_LIGHTING);
        glEnable(GL_POLYGON_OFFSET_FILL);
        m_colorPerVertexShader->setAttributePosition(m_positionVBO);
        m_colorPerVertexShader->setAttributeColor(m_colorVBO);
        m_colorPerVertexShader->setOpacity(1.0);
        map->draw(m_colorPerVertexShader, CGoGN::Algo::Render::GL2::TRIANGLES);
        glDisable(GL_POLYGON_OFFSET_FILL);
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

void Surface_WaveletDecomposition_Plugin::synthesizeFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty() && m_decomposition)
    {
        const QString& mapName = currentItems[0]->text();
        const QString& positionName = m_waveletDecompositionDialog->combo_positionAttribute->currentText();

        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

        if(mh_map && m_decomposition)
        {
            PFP2::MAP* map = mh_map->getMap();
            VertexAttribute<PFP2::VEC3, PFP2::MAP> positionMap = mh_map->getAttribute<PFP2::VEC3, VERTEX>(positionName);

            Decomposition* decomposition = m_decomposition;
        }
    }
}

void Surface_WaveletDecomposition_Plugin::analyzeFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty() && m_decomposition)
    {
        const QString& mapName = currentItems[0]->text();
        const QString& positionName = m_waveletDecompositionDialog->combo_positionAttribute->currentText();

        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapName));

        if(mh_map && m_decomposition)
        {
            PFP2::MAP* map = mh_map->getMap();
            VertexAttribute<PFP2::VEC3, PFP2::MAP> positionMap = mh_map->getAttribute<PFP2::VEC3, VERTEX>(positionName);

            //We search the last level of decomposition
            Decomposition* decomposition = m_decomposition;
            while(decomposition->getChild())
            {
                decomposition = decomposition->getChild();
            }

            if(decomposition->getCorrectionX()>2 && decomposition->getCorrectionY()>2)
            {
                //If a decomposition is still feasible
                decomposition = decomposition->addChild();

                Decomposition* parent = decomposition->getParent();

                //Création d'une nouvelle image composée uniquement des pairs
                QImage& image_parent = parent->getImage();

                int image_width = 0, image_height = 0;
                int width_step = 0, height_step = 0;
                switch(parent->getTransformationType())
                {
                case HORIZONTAL:
                    image_width = image_parent.width()/2 + image_parent.width()%2;
                    width_step = 2;
                    image_height = image_parent.height();
                    height_step = 1;
                    break;
                case VERTICAL:
                    image_width = image_parent.width();
                    width_step = 1;
                    image_height = image_parent.height()/2 + image_parent.height()%2;
                    height_step = 2;
                    break;
                default:
                    break;
                }

                //Creation of a new image composed of the even of parent image
                QImage image(image_width, image_height, image_parent.format());

                for(int i = 0; i < image_width; ++i)
                {
                    for(int j = 0; j < image_height; ++j)
                    {
                        image.setPixel(i, j, image_parent.pixel(i*width_step, j*height_step));
                    }
                }

                decomposition->setImage(image);

                //Creation of a matrix composed of the difference between odd of parent image and prediction value (linear interpolation)
                int diff_red, diff_green, diff_blue;

                QRgb cur_pixel;
                switch(decomposition->getTransformationType())
                {
                case HORIZONTAL:
                    for(int i = 1; i < image_parent.width(); i += width_step)
                    {
                        for(int j = 0; j < image_parent.height(); j += height_step)
                        {
                            cur_pixel = image_parent.pixel(i, j);
                            if(i==image_parent.width()-1)
                            {
                                //Mirror effect at the border of the image
                                diff_red = qRed(cur_pixel) - qRed(image_parent.pixel(i-1,j));
                                diff_green = qGreen(cur_pixel) - qGreen(image_parent.pixel(i-1,j));
                                diff_blue = qBlue(cur_pixel) - qBlue(image_parent.pixel(i-1,j));
                            }
                            else
                            {
                                diff_red = qRed(cur_pixel) - 0.5*(qRed(image_parent.pixel(i-1,j))+qRed(image_parent.pixel(i+1,j)));
                                diff_green = qGreen(cur_pixel) - 0.5*(qGreen(image_parent.pixel(i-1,j))+qGreen(image_parent.pixel(i+1,j)));
                                diff_blue = qBlue(cur_pixel) - 0.5*(qBlue(image_parent.pixel(i-1,j))+qBlue(image_parent.pixel(i+1,j)));
                            }
                            parent->setCorrection(i, j, NQRgb(diff_red, diff_green, diff_blue));
                        }
                    }
                    break;
                case VERTICAL:
                    for(int i = 0; i < image_parent.width(); i += width_step)
                    {
                        for(int j = 1; j < image_parent.height(); j += height_step)
                        {
                            cur_pixel = image_parent.pixel(i, j);
                            if(j==image_parent.height()-1)
                            {
                                //Mirror effect at the border of the image
                                diff_red = qRed(cur_pixel) - qRed(image_parent.pixel(i,j-1));
                                diff_green = qGreen(cur_pixel) - qGreen(image_parent.pixel(i,j-1));
                                diff_blue = qBlue(cur_pixel) - qBlue(image_parent.pixel(i,j-1));
                            }
                            else
                            {
                                diff_red = qRed(cur_pixel) - 0.5*(qRed(image_parent.pixel(i,j-1))+qRed(image_parent.pixel(i,j+1)));
                                diff_green = qGreen(cur_pixel) - 0.5*(qGreen(image_parent.pixel(i,j-1))+qGreen(image_parent.pixel(i,j+1)));
                                diff_blue = qBlue(cur_pixel) - 0.5*(qBlue(image_parent.pixel(i,j-1))+qBlue(image_parent.pixel(i,j+1)));
                            }
                            parent->setCorrection(i, j, NQRgb(diff_red, diff_green, diff_blue));
                        }
                    }
                    break;
                default:
                    break;
                }
                CGoGNout << "New level of decomposition created" << CGoGNendl;
            }
            else
            {
                CGoGNerr << "No more decomposition possible" << CGoGNendl;
            }
        }
    }
}

void Surface_WaveletDecomposition_Plugin::saveImagesFromDialog()
{
    QList<QListWidgetItem*> currentItems = m_waveletDecompositionDialog->list_maps->selectedItems();
    if(!currentItems.empty() && m_decomposition)
    {
        const QString& mapName = currentItems[0]->text();

        Decomposition* decomposition = m_decomposition;
        do
        {
            QImage image = decomposition->getImage();
            QString filename("/home/blettere/Projets/Models/Decomposition/");
            filename.append(mapName);
            filename.append("-");
            filename.append(QString::number(decomposition->getLevel()));
            filename.append(".png");
            if(!image.save(filename))
            {
                CGoGNerr << "Image '" << filename.toStdString() << "' has not been saved" << CGoGNendl;
            }

            decomposition = decomposition->getChild();
            if(decomposition && decomposition->getChild())
            {
                decomposition = decomposition->getChild();
            }
            else
            {
                break;
            }
        } while(decomposition);
    }
}

void Surface_WaveletDecomposition_Plugin::completeAnalysisFromDialog()
{
    if(m_decomposition)
    {
        Decomposition* decomposition = m_decomposition;
        do
        {
            analyzeFromDialog();
            decomposition = decomposition->getChild();
        } while(decomposition->getCorrectionX()>2 && decomposition->getCorrectionY()>2);
    }
}

MapHandlerGen* Surface_WaveletDecomposition_Plugin::initializeObject(const QString& view, const QString& filename)
{
    if(!filename.isEmpty())
    {
        QString file = filename.left(filename.lastIndexOf('.'));
        file = file.mid(file.lastIndexOf('/')).remove(0, 1);
        MapHandlerGen* mhg_map = m_schnapps->addMap(file, 2);
        MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);
        PFP2::MAP* map = mh_map->getMap();

        VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
        if(!position.isValid())
        {
            position = mh_map->addAttribute<PFP2::VEC3, VERTEX>("position");
        }

        VertexAttribute <PFP2::VEC3, PFP2::MAP> color = mh_map->getAttribute<PFP2::VEC3, VERTEX>("color");
        if(!color.isValid())
        {
            color = mh_map->addAttribute<PFP2::VEC3, VERTEX>("color");
        }

        QString extension = filename.mid(filename.lastIndexOf('.'));
        extension.toUpper();

        QImage image;
        if(!image.load(filename, extension.toUtf8().constData()))
        {
            CGoGNout << "Image has not been loaded correctly" << CGoGNendl;
            return NULL;
        }

        image = image.convertToFormat(QImage::Format_RGB32);

        int imageX = image.width(), imageY = image.height();

        Algo::Surface::Tilings::Square::Grid<PFP2> grid(*map, imageX-1, imageY-1);
        grid.embedIntoGrid(position, imageX-1, imageY-1);

        std::vector<Dart> vDarts = grid.getVertexDarts();

        QRgb pixel;

        for(int i = 0; i < imageX; ++i)
        {
            for(int j = 0; j < imageY; ++j)
            {
                pixel = image.pixel(i,(imageY-j)-1);
                color[vDarts[j*imageX+i]] = PFP2::VEC3(qRed(pixel)/255.f, qGreen(pixel)/255.f, qBlue(pixel)/255.f);
            }
        }

        m_decomposition = new Decomposition(imageX/2, imageY, 0, HORIZONTAL);

        m_decomposition->setImage(image);

        mh_map->updateBB(position);
        mh_map->notifyAttributeModification(position);
        mh_map->notifyAttributeModification(color);
        mh_map->notifyConnectivityModification();
        m_colorVBO->updateData(color);
        m_positionVBO->updateData(position);

        m_toDraw = true;

        map->enableQuickTraversal<PFP2::MAP, VERTEX>();

        if(view == m_schnapps->getSelectedView()->getName())
        {
            m_schnapps->getSelectedView()->updateGL();
        }

        return mhg_map;
    }
    return NULL;
}


#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_Plugin, Surface_WaveletDecomposition_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_WaveletDecomposition_PluginD, Surface_WaveletDecomposition_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
